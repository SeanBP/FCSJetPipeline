// Prints a side-by-side comparison of the diagnostics captured by this
// folder's JetEnergyScaleFineGrid/CalibrationMap forks, for two samples
// (e.g. EM-only production vs. the ECAL+HCAL reference). No plotting --
// console summary only, matching this pipeline's other temporary
// diagnostic tools.
//
// Usage: ./summarize_diagnostics <label_a> <bins_a.root> <fitdiag_a.root> \
//                                 <label_b> <bins_b.root> <fitdiag_b.root>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

double percentile_sorted(const std::vector<double> &sorted, double p)
{
    if (sorted.empty()) return NAN;
    if (sorted.size() == 1) return sorted[0];
    double rank = (sorted.size() - 1) * p / 100.0;
    size_t lo = (size_t)std::floor(rank), hi = (size_t)std::ceil(rank);
    if (lo == hi) return sorted[lo];
    double frac = rank - lo;
    return sorted[lo] + frac * (sorted[hi] - sorted[lo]);
}

struct Stats { double mean, median, p90; long n; };

Stats summarize(std::vector<double> v)
{
    Stats s;
    v.erase(std::remove_if(v.begin(), v.end(), [](double x){ return !std::isfinite(x); }), v.end());
    s.n = (long)v.size();
    if (v.empty()) { s.mean = s.median = s.p90 = NAN; return s; }
    double sum = 0; for (double x : v) sum += x;
    s.mean = sum / v.size();
    std::sort(v.begin(), v.end());
    s.median = percentile_sorted(v, 50.0);
    s.p90 = percentile_sorted(v, 90.0);
    return s;
}

void print_stats(const std::string &name, const Stats &s)
{
    std::cout << "    " << std::left << std::setw(28) << name
              << " n=" << std::setw(7) << s.n
              << " mean=" << std::setw(10) << s.mean
              << " median=" << std::setw(10) << s.median
              << " p90=" << s.p90 << std::endl;
}

void analyze_bins_tree(const std::string &label, const std::string &path)
{
    TFile f(path.c_str());
    if (f.IsZombie()) { std::cerr << "Cannot open " << path << std::endl; return; }
    TTree *t = (TTree *)f.Get("bins");
    if (!t) { std::cerr << "No 'bins' tree in " << path << std::endl; return; }

    double scale_mu, resp_mu, scale_sigma, resp_sigma;
    double scale_chi2_ndf, resp_chi2_ndf;
    double scale_raw_mean, scale_raw_sigma, resp_raw_mean, resp_raw_sigma;
    int n_entries;

    t->SetBranchAddress("scale_mu", &scale_mu);
    t->SetBranchAddress("response_mu", &resp_mu);
    t->SetBranchAddress("scale_sigma", &scale_sigma);
    t->SetBranchAddress("response_sigma", &resp_sigma);
    t->SetBranchAddress("scale_chi2_ndf", &scale_chi2_ndf);
    t->SetBranchAddress("response_chi2_ndf", &resp_chi2_ndf);
    t->SetBranchAddress("scale_raw_mean", &scale_raw_mean);
    t->SetBranchAddress("scale_raw_sigma", &scale_raw_sigma);
    t->SetBranchAddress("response_raw_mean", &resp_raw_mean);
    t->SetBranchAddress("response_raw_sigma", &resp_raw_sigma);
    t->SetBranchAddress("n_entries", &n_entries);

    Long64_t n = t->GetEntries();
    long n_candidate = 0, n_converged = 0;
    std::vector<double> scale_chi2, resp_chi2, scale_drift_sigma, resp_drift_sigma, scale_sigma_ratio, resp_sigma_ratio;

    for (Long64_t i = 0; i < n; i++)
    {
        t->GetEntry(i);
        if (n_entries < 5) continue;
        n_candidate++;

        if (scale_mu < 0 && resp_mu < 0) continue; // BinData() default sentinel = fit failed
        n_converged++;

        scale_chi2.push_back(scale_chi2_ndf);
        resp_chi2.push_back(resp_chi2_ndf);

        if (scale_raw_sigma > 0)
            scale_drift_sigma.push_back(std::fabs(scale_mu - scale_raw_mean) / scale_raw_sigma);
        if (resp_raw_sigma > 0)
            resp_drift_sigma.push_back(std::fabs(resp_mu - resp_raw_mean) / resp_raw_sigma);

        if (scale_raw_sigma > 0)
            scale_sigma_ratio.push_back(scale_sigma / scale_raw_sigma);
        if (resp_raw_sigma > 0)
            resp_sigma_ratio.push_back(resp_sigma / resp_raw_sigma);
    }

    std::cout << "=== [" << label << "] JetEnergyScaleFineGrid (per-(E,x,y)-bin Gaussian fit) ===" << std::endl;
    std::cout << "  Candidate bins (n_entries>=5): " << n_candidate
              << "   Converged fits: " << n_converged
              << "   (" << (100.0 * n_converged / std::max(1L, n_candidate)) << "%)" << std::endl;
    print_stats("scale chi2/ndf", summarize(scale_chi2));
    print_stats("response chi2/ndf", summarize(resp_chi2));
    print_stats("|scale mu - raw mean| / raw sigma", summarize(scale_drift_sigma));
    print_stats("|resp mu - raw mean| / raw sigma", summarize(resp_drift_sigma));
    print_stats("scale fit sigma / raw sigma", summarize(scale_sigma_ratio));
    print_stats("resp fit sigma / raw sigma", summarize(resp_sigma_ratio));
    std::cout << std::endl;

    f.Close();
}

void analyze_fit_diagnostics(const std::string &label, const std::string &path)
{
    TFile f(path.c_str());
    if (f.IsZombie()) { std::cerr << "Cannot open " << path << std::endl; return; }
    TTree *t = (TTree *)f.Get("fit_diagnostics");
    if (!t) { std::cerr << "No 'fit_diagnostics' tree in " << path << std::endl; return; }

    int reason, n_points, d_at_bound, e0_at_bound, sigma_at_bound;
    double chi2, bg_chi2, delta_chi2;
    int ndof, bg_ndof;

    t->SetBranchAddress("reason", &reason);
    t->SetBranchAddress("n_points", &n_points);
    t->SetBranchAddress("d_at_bound", &d_at_bound);
    t->SetBranchAddress("e0_at_bound", &e0_at_bound);
    t->SetBranchAddress("sigma_at_bound", &sigma_at_bound);
    t->SetBranchAddress("chi2", &chi2);
    t->SetBranchAddress("bg_chi2", &bg_chi2);
    t->SetBranchAddress("delta_chi2", &delta_chi2);
    t->SetBranchAddress("ndof", &ndof);
    t->SetBranchAddress("bg_ndof", &bg_ndof);

    Long64_t n = t->GetEntries();
    long n_bump_kept = 0, n_few_points = 0, n_bound_hit = 0, n_chi2_not_improved = 0;
    long n_d_bound = 0, n_e0_bound = 0, n_sigma_bound = 0;
    std::vector<double> bg_reduced_chi2_all, kept_reduced_chi2, delta_chi2_attempted;

    for (Long64_t i = 0; i < n; i++)
    {
        t->GetEntry(i);

        switch (reason)
        {
            case 0: n_bump_kept++; break;
            case 1: n_few_points++; break;
            case 2: n_bound_hit++; break;
            case 3: n_chi2_not_improved++; break;
        }

        if (reason == 2)
        {
            if (d_at_bound) n_d_bound++;
            if (e0_at_bound) n_e0_bound++;
            if (sigma_at_bound) n_sigma_bound++;
        }

        if (bg_ndof > 0)
            bg_reduced_chi2_all.push_back(bg_chi2 / bg_ndof);

        if (ndof > 0)
            kept_reduced_chi2.push_back(chi2 / ndof);

        if (reason != 1) // bump was attempted (not skipped for too-few-points)
            delta_chi2_attempted.push_back(delta_chi2);
    }

    std::cout << "=== [" << label << "] CalibrationMap (per-(x,y)-bin background+bump fit) ===" << std::endl;
    std::cout << "  Total position bins fit: " << n << std::endl;
    std::cout << "  Bump kept:            " << n_bump_kept << " (" << (100.0 * n_bump_kept / n) << "%)" << std::endl;
    std::cout << "  Background-only, why:" << std::endl;
    std::cout << "    too few E points:    " << n_few_points << " (" << (100.0 * n_few_points / n) << "%)" << std::endl;
    std::cout << "    bump hit a bound:    " << n_bound_hit << " (" << (100.0 * n_bound_hit / n) << "%)"
              << "  [D=" << n_d_bound << " E0=" << n_e0_bound << " sigma=" << n_sigma_bound << "]" << std::endl;
    std::cout << "    chi2 not improved:   " << n_chi2_not_improved << " (" << (100.0 * n_chi2_not_improved / n) << "%)" << std::endl;
    print_stats("background-only reduced chi2 (ALL bins)", summarize(bg_reduced_chi2_all));
    print_stats("kept-model reduced chi2", summarize(kept_reduced_chi2));
    print_stats("delta_chi2 (bump attempted bins)", summarize(delta_chi2_attempted));
    std::cout << std::endl;

    f.Close();
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc < 7)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <label_a> <bins_a.root> <fitdiag_a.root> <label_b> <bins_b.root> <fitdiag_b.root>" << std::endl;
        return 1;
    }

    analyze_bins_tree(argv[1], argv[2]);
    analyze_fit_diagnostics(argv[1], argv[3]);
    analyze_bins_tree(argv[4], argv[5]);
    analyze_fit_diagnostics(argv[4], argv[6]);

    return 0;
}
