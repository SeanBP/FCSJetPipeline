// Diagnostics-only, standalone plotting tool: picks a few representative
// CalibrationMap (x,y) position bins -- high total jet statistics, and
// where the Gaussian bump hit the OLD (production) fit bounds -- and
// plots each bin's scale-vs-response calibration points with the
// background-only and background+bump curves overlaid, refit under both
// the OLD (production) bounds and the loosened diagnostics bounds. Lets
// us see, per bin, whether loosening the bounds recovers a converged,
// sensible bump shape that actually matches the data, or whether the data
// just doesn't support one.
//
// Self-contained: duplicates the small set of CalibrationMap.cpp helpers
// it needs (BackgroundChi2/GaussChi2/fit machinery) rather than sharing a
// header, since this is throwaway diagnostic tooling.
//
// Usage:
//   ./plot_bump_bins <jet_calibration_diag.root> <old_fit_diagnostics.root> \
//                     <n_bins_to_plot> <output_dir>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TError.h>
#include <TLatex.h>
#include <TString.h>
#include <TColor.h>
#include <TAxis.h>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <algorithm>

// OLD = production CalibrationMap.cpp bounds. NEW = this folder's
// loosened diagnostics bounds (see CalibrationMap.cpp's own comment for
// why). Kept as separate named sets here so both curves can be drawn on
// the same plot regardless of what's currently hardcoded in the .cpp.
struct Bounds { double d_min, d_max, e0_min, e0_max, sigma_min, sigma_max; };
static const Bounds OLD_BOUNDS = {0.0, 0.2, 10.0, 35.0, 1e-3, 25.0};
static const Bounds NEW_BOUNDS = {0.0, 1.0, 5.0, 70.0, 1e-3, 50.0};

const int MIN_POINTS_FOR_BUMP = 7;
const double RESPONSE_SANITY_FACTOR = 1.5;
const double SCALE_SANITY_MAX = 1.5;
const double BUMP_SIGNIFICANCE_CL = 0.5;
const double BOUND_HIT_RTOL = 0.03;

struct BinRow
{
    double x_center, y_center;
    double E_low, E_high;
    double scale_mu, response_mu;
    double scale_mu_err, response_mu_err;
    int n_entries;
};

double log_quad_model(double E, double A, double B, double C)
{
    double lnE = std::log(E);
    return A + B * lnE + C * lnE * lnE;
}

double log_quad_gauss_model(double E, double A, double B, double C, double D, double E0, double sigma)
{
    return log_quad_model(E, A, B, C) + D * std::exp(-((E - E0) * (E - E0)) / (2 * sigma * sigma));
}

struct BackgroundChi2
{
    const std::vector<double> &response, &scale, &scale_err;
    double operator()(const double *p) const
    {
        double sum = 0.0;
        for (size_t i = 0; i < response.size(); i++)
        {
            double model = log_quad_model(response[i], p[0], p[1], p[2]);
            double r = (scale[i] - model) / scale_err[i];
            sum += r * r;
        }
        return sum;
    }
};

struct GaussChi2
{
    const std::vector<double> &response, &scale, &scale_err;
    double operator()(const double *p) const
    {
        double sum = 0.0;
        for (size_t i = 0; i < response.size(); i++)
        {
            double model = log_quad_gauss_model(response[i], p[0], p[1], p[2], p[3], p[4], p[5]);
            double r = (scale[i] - model) / scale_err[i];
            sum += r * r;
        }
        return sum;
    }
};

bool minimize_once(ROOT::Math::Functor &func, const std::string &algo,
                    const std::vector<double> &x0, const std::vector<double> &lower,
                    const std::vector<double> &upper, std::vector<double> &xOut, double &fOut)
{
    std::unique_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit2", algo.c_str()));
    min->SetMaxFunctionCalls(20000);
    min->SetMaxIterations(20000);
    min->SetTolerance(1e-6);
    min->SetPrintLevel(-1);
    min->SetFunction(func);

    for (size_t i = 0; i < x0.size(); i++)
    {
        std::string name = "p" + std::to_string(i);
        double step = std::max(1e-3, std::fabs(x0[i]) * 0.1);
        bool hasLo = std::isfinite(lower[i]), hasHi = std::isfinite(upper[i]);
        if (hasLo && hasHi) min->SetLimitedVariable(i, name.c_str(), x0[i], step, lower[i], upper[i]);
        else if (hasLo) min->SetLowerLimitedVariable(i, name.c_str(), x0[i], step, lower[i]);
        else if (hasHi) min->SetUpperLimitedVariable(i, name.c_str(), x0[i], step, upper[i]);
        else min->SetVariable(i, name.c_str(), x0[i], step);
    }

    bool ok = min->Minimize();
    xOut.assign(min->X(), min->X() + x0.size());
    fOut = min->MinValue();
    return ok;
}

void fit_with_restarts(ROOT::Math::Functor &func, const std::vector<std::vector<double>> &x0_list,
                        const std::vector<double> &lower, const std::vector<double> &upper,
                        std::vector<double> &bestX, double &bestF)
{
    bestF = std::numeric_limits<double>::infinity();
    for (const auto &x0 : x0_list)
    {
        std::vector<double> x; double val;
        bool ok = minimize_once(func, "Migrad", x0, lower, upper, x, val);
        if (!ok) minimize_once(func, "Simplex", x0, lower, upper, x, val);
        if (val < bestF) { bestF = val; bestX = x; }
    }
}

bool at_bound(double value, double lo, double hi, bool check_lower, double rtol = BOUND_HIT_RTOL)
{
    double span = hi - lo;
    double eps = std::max(rtol * span, 1e-8);
    if (check_lower && value <= lo + eps) return true;
    return value >= hi - eps;
}

std::vector<double> linspace(double lo, double hi, int n)
{
    std::vector<double> v(n);
    if (n == 1) { v[0] = lo; return v; }
    for (int i = 0; i < n; i++) v[i] = lo + (hi - lo) * i / (n - 1);
    return v;
}

struct FitOutcome
{
    double A, B, C, D, E0, sigma;
    double chi2;
    int ndof;
    bool bump_kept;
};

// Same decision logic as CalibrationMap.cpp's fit_log_gauss_function_constrained,
// parameterized on Bounds so both OLD and NEW can be run back to back here.
FitOutcome fit_one(const std::vector<double> &response, const std::vector<double> &scale,
                    const std::vector<double> &scale_err, double chi2_improvement_min,
                    const Bounds &b)
{
    FitOutcome out;

    BackgroundChi2 bgChi2{response, scale, scale_err};
    ROOT::Math::Functor bgFunc(bgChi2, 3);
    std::vector<double> lower_bg = {0.0, NAN, NAN}, upper_bg = {NAN, NAN, NAN};
    std::vector<double> bgX; double bgF;
    fit_with_restarts(bgFunc, {{0.2, 0.01, 0.01}}, lower_bg, upper_bg, bgX, bgF);

    auto background_only = [&]()
    {
        out.A = bgX[0]; out.B = bgX[1]; out.C = bgX[2];
        out.D = 0.0; out.E0 = NAN; out.sigma = NAN;
        out.chi2 = bgF; out.ndof = (int)response.size() - 3;
        out.bump_kept = false;
    };

    if ((int)response.size() < MIN_POINTS_FOR_BUMP) { background_only(); return out; }

    double A0 = bgX[0], B0 = bgX[1], C0 = bgX[2];
    std::vector<double> e0_candidates = linspace(b.e0_min, b.e0_max, 3);
    std::vector<double> sigma_candidates = {b.sigma_max * 0.3, b.sigma_max * 0.6};

    GaussChi2 gaussChi2{response, scale, scale_err};
    ROOT::Math::Functor gaussFunc(gaussChi2, 6);
    std::vector<double> lower_g = {0.0, NAN, NAN, b.d_min, b.e0_min, b.sigma_min};
    std::vector<double> upper_g = {NAN, NAN, NAN, b.d_max, b.e0_max, b.sigma_max};

    std::vector<std::vector<double>> x0_list;
    for (double e0 : e0_candidates)
        for (double sigma : sigma_candidates)
            x0_list.push_back({A0, B0, C0, 0.5 * b.d_max, e0, sigma});

    std::vector<double> gX; double gF;
    fit_with_restarts(gaussFunc, x0_list, lower_g, upper_g, gX, gF);

    double D = gX[3], E0v = gX[4], sigma = gX[5];
    double delta_chi2 = bgF - gF;

    bool bump_at_bound = at_bound(D, b.d_min, b.d_max, false) ||
                          at_bound(E0v, b.e0_min, b.e0_max, true) ||
                          at_bound(sigma, b.sigma_min, b.sigma_max, true);

    if (bump_at_bound || delta_chi2 < chi2_improvement_min) { background_only(); return out; }

    out.A = gX[0]; out.B = gX[1]; out.C = gX[2];
    out.D = gX[3]; out.E0 = gX[4]; out.sigma = gX[5];
    out.chi2 = gF; out.ndof = (int)response.size() - 6;
    out.bump_kept = true;
    return out;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;

    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <jet_calibration_diag.root> <old_fit_diagnostics.root> <n_bins> <output_dir>" << std::endl;
        return 1;
    }

    std::string bins_file = argv[1];
    std::string olddiag_file = argv[2];
    int n_to_plot = std::atoi(argv[3]);
    std::string out_dir = argv[4];

    // --- Load per-(E,x,y) calibration points, grouped by (x,y) ---
    TFile fbins(bins_file.c_str());
    TTree *tbins = (TTree *)fbins.Get("bins");
    double x_center, y_center, E_low, E_high, scale_mu, response_mu, scale_mu_err, response_mu_err;
    int n_entries;
    tbins->SetBranchAddress("x_center", &x_center);
    tbins->SetBranchAddress("y_center", &y_center);
    tbins->SetBranchAddress("E_low", &E_low);
    tbins->SetBranchAddress("E_high", &E_high);
    tbins->SetBranchAddress("scale_mu", &scale_mu);
    tbins->SetBranchAddress("response_mu", &response_mu);
    tbins->SetBranchAddress("scale_mu_err", &scale_mu_err);
    tbins->SetBranchAddress("response_mu_err", &response_mu_err);
    tbins->SetBranchAddress("n_entries", &n_entries);

    std::map<std::pair<double, double>, std::vector<BinRow>> groups;
    std::map<std::pair<double, double>, long> total_entries;

    Long64_t nrows = tbins->GetEntries();
    for (Long64_t i = 0; i < nrows; i++)
    {
        tbins->GetEntry(i);
        auto key = std::make_pair(x_center, y_center);
        total_entries[key] += n_entries;
        if (scale_mu > 0 && response_mu > 0) // only store fitted (non-default) rows
            groups[key].push_back({x_center, y_center, E_low, E_high, scale_mu, response_mu,
                                    scale_mu_err, response_mu_err, n_entries});
    }
    fbins.Close();

    // --- Load OLD-bounds fit_diagnostics to find bound-hit bins ---
    TFile folddiag(olddiag_file.c_str());
    TTree *tdiag = (TTree *)folddiag.Get("fit_diagnostics");
    double dx, dy;
    int d_at_bound;
    tdiag->SetBranchAddress("x_center", &dx);
    tdiag->SetBranchAddress("y_center", &dy);
    tdiag->SetBranchAddress("d_at_bound", &d_at_bound);

    std::set<std::pair<double, double>> bound_hit_keys;
    Long64_t ndiag = tdiag->GetEntries();
    for (Long64_t i = 0; i < ndiag; i++)
    {
        tdiag->GetEntry(i);
        if (d_at_bound) bound_hit_keys.insert({dx, dy});
    }
    folddiag.Close();

    std::cout << "Loaded " << groups.size() << " position bins, "
              << bound_hit_keys.size() << " were D-bound-hit under OLD bounds" << std::endl;

    // --- Rank D-bound-hit bins by total jet statistics ---
    std::vector<std::pair<long, std::pair<double, double>>> ranked;
    for (const auto &key : bound_hit_keys)
    {
        if (groups.count(key) == 0) continue;
        ranked.push_back({total_entries[key], key});
    }
    std::sort(ranked.begin(), ranked.end(),
              std::greater<std::pair<long, std::pair<double, double>>>());

    if ((int)ranked.size() < n_to_plot)
        std::cout << "WARNING: only " << ranked.size() << " D-bound-hit bins available, requested " << n_to_plot << std::endl;

    double chi2_improvement_min = TMath::ChisquareQuantile(BUMP_SIGNIFICANCE_CL, 3);

    {
        std::string cmd = "mkdir -p " + out_dir;
        system(cmd.c_str());
    }

    int n_plotted = 0;
    for (const auto &entry : ranked)
    {
        if (n_plotted >= n_to_plot) break;
        long n_stat = entry.first;
        double xc = entry.second.first, yc = entry.second.second;

        const std::vector<BinRow> &rows = groups[{xc, yc}];

        std::vector<double> response, scale, scale_err, response_err;
        for (const auto &r : rows)
        {
            bool valid = std::isfinite(r.scale_mu) && std::isfinite(r.response_mu) &&
                         std::isfinite(r.scale_mu_err) && r.scale_mu_err > 0 &&
                         r.scale_mu <= SCALE_SANITY_MAX && r.response_mu <= RESPONSE_SANITY_FACTOR * r.E_high;
            if (!valid) continue;
            response.push_back(r.response_mu);
            scale.push_back(r.scale_mu);
            scale_err.push_back(r.scale_mu_err);
            response_err.push_back(r.response_mu_err);
        }

        if (response.size() < 4) continue;

        FitOutcome old_fit = fit_one(response, scale, scale_err, chi2_improvement_min, OLD_BOUNDS);
        FitOutcome new_fit = fit_one(response, scale, scale_err, chi2_improvement_min, NEW_BOUNDS);

        std::cout << "Bin (x=" << xc << ", y=" << yc << "), n_entries=" << n_stat
                  << ", n_points=" << response.size()
                  << " | OLD bump_kept=" << old_fit.bump_kept << " chi2/ndof=" << (old_fit.chi2 / old_fit.ndof)
                  << " | NEW bump_kept=" << new_fit.bump_kept << " chi2/ndof=" << (new_fit.chi2 / new_fit.ndof)
                  << (new_fit.bump_kept ? "" : "") << std::endl;
        if (new_fit.bump_kept)
            std::cout << "    NEW bump params: D=" << new_fit.D << " E0=" << new_fit.E0 << " sigma=" << new_fit.sigma << std::endl;

        // --- Plot ---
        double rmin = *std::min_element(response.begin(), response.end());
        double rmax = *std::max_element(response.begin(), response.end());
        double pad = 0.1 * (rmax - rmin);

        TGraphErrors g((int)response.size(), response.data(), scale.data(), response_err.data(), scale_err.data());
        g.SetMarkerStyle(20);
        g.SetMarkerColor(kBlack);
        g.SetLineColor(kBlack);
        g.SetTitle(Form("Bin x=%.1f y=%.1f  (n_entries=%ld, n_points=%d);response_mu (GeV);scale_mu",
                         xc, yc, n_stat, (int)response.size()));

        const int NCURVE = 200;
        std::vector<double> curve_x(NCURVE);
        std::vector<double> curve_old(NCURVE), curve_new(NCURVE);
        for (int i = 0; i < NCURVE; i++)
        {
            double e = rmin - pad + (rmax - rmin + 2 * pad) * i / (NCURVE - 1);
            curve_x[i] = e;
            curve_old[i] = old_fit.bump_kept
                ? log_quad_gauss_model(e, old_fit.A, old_fit.B, old_fit.C, old_fit.D, old_fit.E0, old_fit.sigma)
                : log_quad_model(e, old_fit.A, old_fit.B, old_fit.C);
            curve_new[i] = new_fit.bump_kept
                ? log_quad_gauss_model(e, new_fit.A, new_fit.B, new_fit.C, new_fit.D, new_fit.E0, new_fit.sigma)
                : log_quad_model(e, new_fit.A, new_fit.B, new_fit.C);
        }

        TGraph g_old(NCURVE, curve_x.data(), curve_old.data());
        g_old.SetLineColor(kRed);
        g_old.SetLineStyle(2);
        g_old.SetLineWidth(2);

        TGraph g_new(NCURVE, curve_x.data(), curve_new.data());
        g_new.SetLineColor(kBlue);
        g_new.SetLineStyle(1);
        g_new.SetLineWidth(2);

        // Base axis ranges on the data VALUES, not the error-bar extents --
        // TGraphErrors otherwise auto-scales to include full error bars, so
        // one poorly-constrained point (huge response_mu_err/scale_mu_err)
        // can squash every other point into an unreadable clump.
        double smin = *std::min_element(scale.begin(), scale.end());
        double smax = *std::max_element(scale.begin(), scale.end());
        double spad = 0.15 * (smax - smin);
        if (spad <= 0) spad = 0.1 * std::max(1e-3, std::fabs(smax));

        TCanvas c("c", "", 900, 650);
        g.Draw("AP");
        g.GetXaxis()->SetLimits(rmin - pad, rmax + pad);
        g.GetYaxis()->SetRangeUser(smin - spad, smax + spad);
        g_old.Draw("L SAME");
        g_new.Draw("L SAME");

        TLegend leg(0.55, 0.15, 0.88, 0.35);
        leg.AddEntry(&g, "data (scale_mu vs response_mu)", "p");
        leg.AddEntry(&g_old, Form("OLD bounds (%s, #chi^{2}/ndof=%.2f)",
                                   old_fit.bump_kept ? "bump kept" : "background-only", old_fit.chi2 / old_fit.ndof), "l");
        leg.AddEntry(&g_new, Form("NEW loosened bounds (%s, #chi^{2}/ndof=%.2f)",
                                   new_fit.bump_kept ? "bump kept" : "background-only", new_fit.chi2 / new_fit.ndof), "l");
        leg.Draw();

        std::string outpath = out_dir + Form("/bump_bin_x%.0f_y%.0f.png", xc, yc);
        c.Print(outpath.c_str());
        std::cout << "  -> " << outpath << std::endl;

        n_plotted++;
    }

    std::cout << "Plotted " << n_plotted << " bins to " << out_dir << std::endl;
    return 0;
}
