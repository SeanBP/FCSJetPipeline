#ifndef GAUSSCORE_H
#define GAUSSCORE_H

// Iterative Gaussian-core fit of a weighted sample, for response
// distributions with long tails and skewed peaks (EM-only jets: a sharp upper
// edge near 1 at low energy, leakage below, hadron showers above). Seed: the
// peak, not the median -- the maximum of a weighted histogram over
// median +- 4 MAD widths (bins of a tenth of the MAD width, smoothed over three
// bins), with sigma_HM = (narrower half-width at half maximum) / 1.1774, so a
// long tail on one side does not widen the window. Each iteration fills a
// weighted histogram on centre +- 5 sigma and fits a Gaussian on
// [centre - n_lo sigma, centre + n_hi sigma], then re-centres on the fitted
// mean and width; it stops when the mean moves by less than tol x sigma or
// after max_iter fits. The fitted mean is held inside the first (peak-centred)
// window and the width within 0.3-3 x sigma_HM throughout, so the window
// cannot walk away from the peak. A degenerate sample (MAD width at or below sigma_floor,
// e.g. more than half the entries identical) has no Gaussian core: its "core
// mean" is then the weighted mean within centre +- n sigma_floor, flagged
// degenerate. The JES maps and the linearity test share it, so "Gaussian
// mean" means the same estimator everywhere.

#include <TH1D.h>
#include <TF1.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

struct GaussCoreConfig
{
    double n_lo = 1.0, n_hi = 1.0;  // fit window in sigmas below / above the centre
    int n_bins = 50;                // histogram bins over centre +- 5 sigma
    int max_iter = 3;
    double tol = 0.25;              // converged when |delta mean| < tol x sigma (the mean is held in the first window anyway)
    double sigma_floor = 0.002;     // absolute, in the variable's units
};

struct GaussCoreResult
{
    double mu = NAN, mu_err = NAN, sigma = NAN, sigma_err = NAN, chi2_ndf = NAN;
    double median = NAN;
    int n_iter = 0;
    bool converged = false;
    bool degenerate = false;
};

// Weighted median: first value whose cumulative weight reaches half the total
inline double GaussCoreMedian(std::vector<std::pair<double, double>> vw)
{
    if (vw.empty()) return NAN;
    std::sort(vw.begin(), vw.end());
    double tot = 0;
    for (auto &p : vw) tot += p.second;
    double run = 0;
    for (auto &p : vw) { run += p.second; if (run >= 0.5 * tot) return p.first; }
    return vw.back().first;
}

// keep_hist / keep_fit, if given, receive the last iteration's histogram and
// function (owned by the caller), for plotting
inline GaussCoreResult GaussCoreFit(const std::vector<double> &v, const std::vector<double> &w, const GaussCoreConfig &cfg = GaussCoreConfig(),
                                    TH1D **keep_hist = 0, TF1 **keep_fit = 0)
{
    GaussCoreResult r;
    if (v.size() < 5) return r;
    std::vector<std::pair<double, double>> vw(v.size());
    double sum_w = 0;
    for (size_t i = 0; i < v.size(); i++) { vw[i] = {v[i], w[i]}; sum_w += w[i]; }
    if (!(sum_w > 0)) return r;
    r.median = GaussCoreMedian(vw);
    std::vector<std::pair<double, double>> dev(v.size());
    for (size_t i = 0; i < v.size(); i++) dev[i] = {std::fabs(v[i] - r.median), w[i]};
    const double mad_sigma = 1.4826 * GaussCoreMedian(dev);
    double centre = r.median, sigma = std::max(mad_sigma, cfg.sigma_floor);
    const double w_norm = v.size() / sum_w; // histogram contents at the scale of entry counts
    if (!(mad_sigma > cfg.sigma_floor))
    {
        double s = 0, sw = 0, s2 = 0;
        for (size_t i = 0; i < v.size(); i++)
            if (v[i] >= centre - cfg.n_lo * cfg.sigma_floor && v[i] <= centre + cfg.n_hi * cfg.sigma_floor) { s += w[i] * v[i]; sw += w[i]; s2 += w[i] * w[i]; }
        r.mu = sw > 0 ? s / sw : centre;
        r.sigma = mad_sigma;
        r.mu_err = sw > 0 ? cfg.sigma_floor / std::sqrt(sw * sw / s2) : NAN;
        r.degenerate = true;
        r.converged = true;
        return r;
    }

    // Peak seed: maximum of a smoothed weighted histogram, width from the
    // narrower side's half width at half maximum
    {
        const double bw = 0.1 * sigma;
        const int nb = (int)std::ceil(8 * sigma / bw);
        std::vector<double> c(nb, 0.0), sm(nb, 0.0);
        const double x0 = r.median - 4 * sigma;
        for (size_t i = 0; i < v.size(); i++)
        {
            const int b = (int)std::floor((v[i] - x0) / bw);
            if (b >= 0 && b < nb) c[b] += w[i];
        }
        for (int b = 0; b < nb; b++)
        {
            double s3 = 0; int n3 = 0;
            for (int d = -1; d <= 1; d++) if (b + d >= 0 && b + d < nb) { s3 += c[b + d]; n3++; }
            sm[b] = s3 / n3;
        }
        int bmax = (int)(std::max_element(sm.begin(), sm.end()) - sm.begin());
        const double half = 0.5 * sm[bmax];
        int bl = bmax, br = bmax;
        while (bl > 0 && sm[bl] > half) bl--;
        while (br < nb - 1 && sm[br] > half) br++;
        const double hl = (bmax - bl) * bw, hr = (br - bmax) * bw;
        const double hw = std::min(hl > 0 ? hl : hr, hr > 0 ? hr : hl);
        if (hw > 0)
        {
            centre = x0 + (bmax + 0.5) * bw;
            sigma = std::max(hw / 1.1774, cfg.sigma_floor);
        }
    }
    const double sigma0 = sigma, mu_lo = centre - cfg.n_lo * sigma0, mu_hi = centre + cfg.n_hi * sigma0;
    for (int it = 0; it < cfg.max_iter; it++)
    {
        TH1D *h = new TH1D(Form("gc_h_%d_%p", it, (void *)&r), "", cfg.n_bins, centre - 5 * sigma, centre + 5 * sigma);
        h->Sumw2();
        h->SetDirectory(0);
        for (size_t i = 0; i < v.size(); i++) h->Fill(v[i], w[i] * w_norm);
        const double lo = centre - cfg.n_lo * sigma, hi = centre + cfg.n_hi * sigma;
        TF1 *f = new TF1(Form("gc_f_%d_%p", it, (void *)&r), "gaus", lo, hi);
        f->SetParameters(h->GetBinContent(h->FindBin(centre)), centre, sigma);
        f->SetParLimits(1, mu_lo, mu_hi);
        f->SetParLimits(2, 0.3 * sigma0, 3 * sigma0);
        const int st = h->Fit(f, "RQN0");
        r.n_iter = it + 1;
        if (st != 0 || !(f->GetParameter(2) > 0))
        {
            delete h; delete f;
            r.converged = false;
            return r;
        }
        const double mu = f->GetParameter(1), s = std::max(f->GetParameter(2), cfg.sigma_floor);
        r.mu = mu; r.mu_err = f->GetParError(1); r.sigma = f->GetParameter(2); r.sigma_err = f->GetParError(2);
        r.chi2_ndf = f->GetNDF() > 0 ? f->GetChisquare() / f->GetNDF() : NAN;
        // converged when the mean moved by less than tol x sigma, or by less
        // than twice its own statistical error (low statistics)
        const bool done = std::fabs(mu - centre) < std::max(cfg.tol * sigma, 2 * r.mu_err);
        centre = mu; sigma = s;
        if (done || it == cfg.max_iter - 1)
        {
            r.converged = done;
            if (keep_hist) *keep_hist = h; else delete h;
            if (keep_fit) *keep_fit = f; else delete f;
            return r;
        }
        delete h; delete f;
    }
    return r;
}

#endif
