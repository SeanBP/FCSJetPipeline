// Step 2 of the EM-only JES calibration: builds the JES lookup JSON from the
// "bins" tree of JetEnergyScaleFineGrid_em_only (per-bin weighted medians).
// Per (x,y) position bin, fits scale vs response with
//   A + B/E + C/E^2 + D exp(-(E-E0)^2 / 2 sigma^2),  sigma = sigmaL below E0, sigmaR above,
// then smooths the parameter maps over the grid and masks cells outside the
// buffered ECAL fiducial region. JSON bins carry {i, j, A, B, C, D, E0, sigmaL,
// sigmaR} (read by jet_calibration/ApplyCorrections_em_only.cpp).
//
// Usage: ./CalibrationMap_em_only <bins.root> <out.json>
//        ./CalibrationMap_em_only <bins.root> <out.json> x x --shard <i> <n> <shard.csv>
//        ./CalibrationMap_em_only <bins.root> <out.json> x x --merge <shard_list> <out.json>
// (argv[3] and argv[4] are ignored placeholders)
// --merge also writes <out>.chi2.csv and <out>.direct_fits.csv (pre-smoothing fits).

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMath.h>
#include <TError.h>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <string>
#include <cmath>
#include <limits>
#include <memory>
#include <algorithm>
#include <cstdio>
#include <cstdlib>


// CONFIGURATION

std::string ROOT_FILE = "jet_calibration_grid_test.root";
std::string OUTPUT_JSON = "JetEnergyScale_lookup_GridTest.json";

// Below this many energy points only the 3-parameter background is fit
const int MIN_POINTS_FOR_SKEW_GAUSSIAN = 8;

// Sanity cuts dropping the noise tail of low-statistics bins
const double RESPONSE_SANITY_FACTOR = 1.5;
const double SCALE_SANITY_MAX = 1.5;

// Skew-Gaussian bounds. E0_MAX/SIGMA_MAX are set in main() for the EM-only
// response domain (Emax = 100). The bounds sit well above the optima found in
// the hottest bins, so a fit at a bound signals a problem, not a clipped peak.
const double D_MIN = 0.0, D_MAX = 0.5;
const double E0_MIN = 1.0;
double E0_MAX = 35.0;
const double E0_MAX_ecal_hcal = 35.0, E0_MAX_em_only = 40.0;
// About the finest E-bin spacing: narrower widths fit single-point steps
const double SIGMA_MIN = 1.0;
double SIGMA_MAX = 25.0;
const double SIGMA_MAX_ecal_hcal = 25.0, SIGMA_MAX_em_only = 35.0;

// The skew Gaussian is kept only if its chi2 improvement beats this quantile
// of a chi-square with 4 dof (D, E0, sigmaL, sigmaR)
const double SKEW_GAUSSIAN_SIGNIFICANCE_CL = 0.5;

// Copied from shared/JetParameters.h (keep in sync; that header needs
// fastjet). Upstream, JetEnergyScaleFineGrid only fills grid points inside
// this region, so cells outside it have no statistics; the saved map is
// masked to it so smoothing never extrapolates past the acceptance.
struct FiducialRect
{
    float x_inner, x_outer, y_min, y_max;
};

const float cut_x_inner = 24.15591212f, cut_x_outer = 131.82035143f;
const float cut_y_min = -83.19687629f, cut_y_max = 86.41376856f;
const float cut_x_inner_ecal = 20.53604584f, cut_x_outer_ecal = 137.59500000f;
const float cut_y_min_ecal = -96.72439178f, cut_y_max_ecal = 86.25565167f;

const FiducialRect kFiducialRect = {cut_x_inner, cut_x_outer, cut_y_min, cut_y_max};
const FiducialRect kFiducialRectEcal = {cut_x_inner_ecal, cut_x_outer_ecal, cut_y_min_ecal, cut_y_max_ecal};

const float reco_fiducial_buffer = 10.0f;
const float reco_fiducial_buffer_ecal = 20.0f;

// Set in main(): the EM-only rectangle and buffer
const FiducialRect *ACTIVE_FIDUCIAL_RECT = &kFiducialRect;
float ACTIVE_FIDUCIAL_BUFFER = reco_fiducial_buffer;

inline bool pass_fiducial_cut(float jetXE, float jetYE, float buffer, const FiducialRect &rect)
{
    const float ax = std::fabs(jetXE);
    const bool in_x = ax >= (rect.x_inner + buffer) && ax <= (rect.x_outer - buffer);
    const bool in_y = jetYE >= (rect.y_min + buffer) && jetYE <= (rect.y_max - buffer);
    return in_x && in_y;
}

// A fit within this fraction of a parameter range from a bound counts as at the bound
const double BOUND_HIT_RTOL = 0.03;

// Below this, D counts as "no skew Gaussian" for reporting purposes.
const double SKEW_GAUSSIAN_PRESENT_THRESHOLD = 1e-4;

// Map smoothing parameters
const int SMOOTH_N_ITER = 6;
const int SMOOTH_RADIUS = 2;
// -1 = unfit bins take the diffusion average (cells far outside the data stay NaN)
const int SMOOTH_MAX_FILL_DISTANCE = -1;

// D is only meaningful where a direct fit supports it, but diffusion would
// let unfit cells inherit a raised D from each other. D decays exponentially
// with distance (in 2.5 cm grid cells) from the nearest fitted cell.
const double D_DECAY_LENGTH_BINS = 3.0;


struct BinRow
{
    double x_center, y_center;
    double E_low, E_high;
    double scale_mu, response_mu;
    double scale_mu_err, response_mu_err;
};

struct FitRow
{
    double x_center, y_center;
    double A, B, C, D, E0, sigmaL, sigmaR;
    double chi2;
    int n_params;
    int ndof;
};

typedef std::vector<std::vector<double>> Grid2D;

std::vector<BinRow> load_bins(const std::string &root_file)
{
    TFile f(root_file.c_str());
    if (f.IsZombie())
    {
        std::cerr << "ERROR: cannot open " << root_file << std::endl;
        exit(1);
    }

    TTree *tree = (TTree *)f.Get("bins");
    if (!tree)
    {
        std::cerr << "ERROR: no 'bins' tree in " << root_file << std::endl;
        exit(1);
    }

    double x_center, y_center, E_low, E_high;
    double scale_mu, response_mu, scale_mu_err, response_mu_err;

    tree->SetBranchAddress("x_center", &x_center);
    tree->SetBranchAddress("y_center", &y_center);
    tree->SetBranchAddress("E_low", &E_low);
    tree->SetBranchAddress("E_high", &E_high);
    tree->SetBranchAddress("scale_mu", &scale_mu);
    tree->SetBranchAddress("response_mu", &response_mu);
    tree->SetBranchAddress("scale_mu_err", &scale_mu_err);
    tree->SetBranchAddress("response_mu_err", &response_mu_err);

    std::vector<BinRow> rows;
    rows.reserve(tree->GetEntries());

    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; i++)
    {
        tree->GetEntry(i);
        rows.push_back({x_center, y_center, E_low, E_high,
                         scale_mu, response_mu, scale_mu_err, response_mu_err});
    }

    f.Close();
    return rows;
}

// Inverse-power background: tends to the high-E plateau A as E -> infinity
double bg_model(double E, double A, double B, double C)
{
    return A + B / E + C / (E * E);
}

// Asymmetric Gaussian: independent rise (E < E0) and fall widths, since many
// bins show a fast rise and slow decline
double bg_gauss_model(double E, double A, double B, double C, double D, double E0, double sigmaL, double sigmaR)
{
    double sigma = (E < E0) ? sigmaL : sigmaR;
    return bg_model(E, A, B, C) + D * std::exp(-((E - E0) * (E - E0)) / (2 * sigma * sigma));
}

struct BackgroundChi2
{
    const std::vector<double> &response, &scale, &scale_err;

    double operator()(const double *p) const
    {
        double sum = 0.0;
        for (size_t i = 0; i < response.size(); i++)
        {
            double model = bg_model(response[i], p[0], p[1], p[2]);
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
            double model = bg_gauss_model(response[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
            double r = (scale[i] - model) / scale_err[i];
            sum += r * r;
        }
        return sum;
    }
};

// Bounded chi2 minimization: Migrad (primary), Simplex fallback if it doesn't converge.
bool minimize_once(
    ROOT::Math::Functor &func,
    const std::string &algo,
    const std::vector<double> &x0,
    const std::vector<double> &lower,
    const std::vector<double> &upper,
    std::vector<double> &xOut,
    double &fOut)
{
    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", algo.c_str()));

    min->SetMaxFunctionCalls(20000);
    min->SetMaxIterations(20000);
    min->SetTolerance(1e-6);
    min->SetPrintLevel(-1);
    min->SetFunction(func);

    for (size_t i = 0; i < x0.size(); i++)
    {
        std::string name = "p" + std::to_string(i);
        double step = std::max(1e-3, std::fabs(x0[i]) * 0.1);

        bool hasLo = std::isfinite(lower[i]);
        bool hasHi = std::isfinite(upper[i]);

        if (hasLo && hasHi)
            min->SetLimitedVariable(i, name.c_str(), x0[i], step, lower[i], upper[i]);
        else if (hasLo)
            min->SetLowerLimitedVariable(i, name.c_str(), x0[i], step, lower[i]);
        else if (hasHi)
            min->SetUpperLimitedVariable(i, name.c_str(), x0[i], step, upper[i]);
        else
            min->SetVariable(i, name.c_str(), x0[i], step);
    }

    bool ok = min->Minimize();

    xOut.assign(min->X(), min->X() + x0.size());
    fOut = min->MinValue();

    return ok;
}

void fit_with_restarts(
    ROOT::Math::Functor &func,
    const std::vector<std::vector<double>> &x0_list,
    const std::vector<double> &lower,
    const std::vector<double> &upper,
    std::vector<double> &bestX,
    double &bestF)
{
    bestF = std::numeric_limits<double>::infinity();

    for (const auto &x0 : x0_list)
    {
        std::vector<double> x;
        double val;

        bool ok = minimize_once(func, "Migrad", x0, lower, upper, x, val);
        if (!ok)
            minimize_once(func, "Simplex", x0, lower, upper, x, val);

        if (val < bestF)
        {
            bestF = val;
            bestX = x;
        }
    }
}

// True if value sits at (within rtol of) lo or hi. check_lower=false skips the
// lower check -- for D, lower bound 0 is a legitimate "no skew Gaussian" outcome, not degenerate.
bool at_bound(double value, double lo, double hi, bool check_lower, double rtol = BOUND_HIT_RTOL)
{
    double span = hi - lo;
    double eps = std::max(rtol * span, 1e-8);

    if (check_lower && value <= lo + eps)
        return true;

    return value >= hi - eps;
}

std::vector<double> linspace(double lo, double hi, int n)
{
    std::vector<double> v(n);
    if (n == 1)
    {
        v[0] = lo;
        return v;
    }
    for (int i = 0; i < n; i++)
        v[i] = lo + (hi - lo) * i / (n - 1);
    return v;
}

// Fits background, then background+skew-Gaussian if enough points (keeping the
// skew Gaussian only if significant and not bound-pinned). False if too few points to fit at all.
bool fit_log_gauss_function_constrained(
    const std::vector<double> &response_in,
    const std::vector<double> &scale_in,
    const std::vector<double> &scale_err_in,
    double chi2_improvement_min,
    FitRow &out)
{
    std::vector<double> response, scale, scale_err;

    for (size_t i = 0; i < response_in.size(); i++)
    {
        if (std::isfinite(response_in[i]) && std::isfinite(scale_in[i]) &&
            std::isfinite(scale_err_in[i]) && response_in[i] > 0 && scale_err_in[i] > 0)
        {
            response.push_back(response_in[i]);
            scale.push_back(scale_in[i]);
            scale_err.push_back(scale_err_in[i]);
        }
    }

    if (response.size() < 3)
        return false;

    BackgroundChi2 bgChi2{response, scale, scale_err};
    ROOT::Math::Functor bgFunc(bgChi2, 3);

    std::vector<double> lower_bg = {0.0, NAN, NAN};
    std::vector<double> upper_bg = {NAN, NAN, NAN};

    std::vector<double> bgX;
    double bgF;
    // A is the high-E scale, near 1
    fit_with_restarts(bgFunc, {{1.0, 1.0, -1.0}}, lower_bg, upper_bg, bgX, bgF);

    auto background_only = [&](const std::string &reason)
    {
        (void)reason;
        out.A = bgX[0];
        out.B = bgX[1];
        out.C = bgX[2];
        out.D = 0.0;
        out.E0 = std::numeric_limits<double>::quiet_NaN();
        out.sigmaL = std::numeric_limits<double>::quiet_NaN();
        out.sigmaR = std::numeric_limits<double>::quiet_NaN();
        out.chi2 = bgF;
        out.n_params = 3;
        out.ndof = (int)response.size() - 3;
    };

    if ((int)response.size() < MIN_POINTS_FOR_SKEW_GAUSSIAN)
    {
        background_only("few_points");
        return true;
    }

    double A0 = bgX[0], B0 = bgX[1], C0 = bgX[2];

    // E0 seeds: the union of a 3- and a 6-point grid. Some bins converge to a
    // poor local optimum from a sparse grid; the union keeps every seed of the
    // smaller grid, so no bin can do worse than with it.
    std::vector<double> e0_candidates = linspace(E0_MIN, E0_MAX, 3);
    {
        std::vector<double> e0_dense = linspace(E0_MIN, E0_MAX, 6);
        e0_candidates.insert(e0_candidates.end(), e0_dense.begin(), e0_dense.end());
        // drop the duplicated endpoints
        std::sort(e0_candidates.begin(), e0_candidates.end());
        e0_candidates.erase(std::unique(e0_candidates.begin(), e0_candidates.end()), e0_candidates.end());
    }
    // Every (sigmaL, sigmaR) combination, including asymmetric ones, so
    // Migrad can start near a fast-rise/slow-fall shape
    std::vector<double> sigma_candidates = {SIGMA_MAX * 0.1, SIGMA_MAX * 0.15, SIGMA_MAX * 0.2, SIGMA_MAX * 0.3,
                                             SIGMA_MAX * 0.35, SIGMA_MAX * 0.55, SIGMA_MAX * 0.6};

    GaussChi2 gaussChi2{response, scale, scale_err};
    ROOT::Math::Functor gaussFunc(gaussChi2, 7);

    std::vector<double> lower_g = {0.0, NAN, NAN, D_MIN, E0_MIN, SIGMA_MIN, SIGMA_MIN};
    std::vector<double> upper_g = {NAN, NAN, NAN, D_MAX, E0_MAX, SIGMA_MAX, SIGMA_MAX};

    std::vector<std::vector<double>> x0_list;
    for (double e0 : e0_candidates)
        for (double sigmaL : sigma_candidates)
            for (double sigmaR : sigma_candidates)
                x0_list.push_back({A0, B0, C0, 0.5 * D_MAX, e0, sigmaL, sigmaR});

    // Data-driven seed: the largest positive residual above the background
    // fit, with widths from its half-maximum crossings on either side
    {
        int peak_idx = -1;
        double peak_resid = 0.0;
        std::vector<double> resid(response.size());
        for (size_t i = 0; i < response.size(); i++)
        {
            resid[i] = scale[i] - bg_model(response[i], A0, B0, C0);
            if (response[i] >= E0_MIN && response[i] <= E0_MAX && resid[i] > peak_resid)
            {
                peak_resid = resid[i];
                peak_idx = (int)i;
            }
        }
        if (peak_idx >= 0)
        {
            double half = peak_resid / 2.0;
            double e_lo = response[0], e_hi = response.back();
            for (int i = peak_idx; i >= 0; i--)
                if (resid[i] < half) { e_lo = response[i]; break; }
            for (size_t i = peak_idx; i < response.size(); i++)
                if (resid[i] < half) { e_hi = response[i]; break; }
            double e0_est = std::min(E0_MAX, std::max(E0_MIN, response[peak_idx]));
            double sigmaL_est = std::max(SIGMA_MIN * 2, (e0_est - e_lo) / 1.1775);
            double sigmaR_est = std::max(SIGMA_MIN * 2, (e_hi - e0_est) / 1.1775);
            sigmaL_est = std::min(sigmaL_est, SIGMA_MAX * 0.9);
            sigmaR_est = std::min(sigmaR_est, SIGMA_MAX * 0.9);
            double D_est = std::min(D_MAX * 0.9, std::max(D_MIN + 1e-3, peak_resid));
            x0_list.push_back({A0, B0, C0, D_est, e0_est, sigmaL_est, sigmaR_est});
        }
    }

    // Keep the best fit overall and the best fit not at a bound: the global
    // best is sometimes a very broad, bound-pinned Gaussian mimicking
    // background curvature, while a well-localized interior fit exists.
    std::vector<double> gX, bestValidX;
    double gF = std::numeric_limits<double>::infinity();
    double bestValidF = std::numeric_limits<double>::infinity();
    bool haveValid = false;

    for (const auto &x0 : x0_list)
    {
        std::vector<double> x;
        double val;
        bool ok = minimize_once(gaussFunc, "Migrad", x0, lower_g, upper_g, x, val);
        if (!ok)
            minimize_once(gaussFunc, "Simplex", x0, lower_g, upper_g, x, val);

        if (val < gF)
        {
            gF = val;
            gX = x;
        }

        // Only the upper width bound is degenerate; a width at the lower
        // bound is a sharp rise or fall
        bool pinned =
            at_bound(x[3], D_MIN, D_MAX, false) ||
            at_bound(x[4], E0_MIN, E0_MAX, true) ||
            at_bound(x[5], SIGMA_MIN, SIGMA_MAX, false) ||
            at_bound(x[6], SIGMA_MIN, SIGMA_MAX, false);

        if (!pinned && val < bestValidF)
        {
            bestValidF = val;
            bestValidX = x;
            haveValid = true;
        }
    }

    double D = gX[3], E0v = gX[4], sigmaL = gX[5], sigmaR = gX[6];
    double delta_chi2 = bgF - gF;

    bool skew_gaussian_at_bound =
        at_bound(D, D_MIN, D_MAX, false) ||
        at_bound(E0v, E0_MIN, E0_MAX, true) ||
        at_bound(sigmaL, SIGMA_MIN, SIGMA_MAX, false) ||
        at_bound(sigmaR, SIGMA_MIN, SIGMA_MAX, false);

    // Fall back to the best fit not at a bound
    if (skew_gaussian_at_bound && haveValid)
    {
        gX = bestValidX;
        gF = bestValidF;
        D = gX[3];
        E0v = gX[4];
        sigmaL = gX[5];
        sigmaR = gX[6];
        delta_chi2 = bgF - gF;
        skew_gaussian_at_bound = false;
    }

    if (skew_gaussian_at_bound)
    {
        background_only("bound_hit");
        return true;
    }

    if (delta_chi2 < chi2_improvement_min)
    {
        background_only("chi2_not_improved");
        return true;
    }

    out.A = gX[0];
    out.B = gX[1];
    out.C = gX[2];
    out.D = gX[3];
    out.E0 = gX[4];
    out.sigmaL = gX[5];
    out.sigmaR = gX[6];
    out.chi2 = gF;
    out.n_params = 7;
    out.ndof = (int)response.size() - 7;

    return true;
}

// Position bins are fit independently (bins with < 4 valid entries are
// dropped). To parallelize, run separate processes with --shard / --merge;
// ROOT 5.34's Minuit2 is safe neither with threads nor after fork().
struct GroupEntry
{
    double x_center, y_center;
    std::vector<int> idxs;
};

// Groups rows by (x_center, y_center), in deterministic map order, so the
// partition k % n_shards == shard is the same in every process
std::vector<GroupEntry> build_position_groups(const std::vector<BinRow> &rows)
{
    std::map<std::pair<double, double>, std::vector<int>> groups;
    for (size_t i = 0; i < rows.size(); i++)
        groups[{rows[i].x_center, rows[i].y_center}].push_back((int)i);

    std::vector<GroupEntry> entries;
    entries.reserve(groups.size());
    for (const auto &kv : groups)
        entries.push_back({kv.first.first, kv.first.second, kv.second});
    return entries;
}

// Fits every group k with k % n_shards == shard_index
std::vector<FitRow> compute_shard(const std::vector<BinRow> &rows, double chi2_improvement_min,
                                   unsigned shard_index, unsigned n_shards)
{
    std::vector<GroupEntry> entries = build_position_groups(rows);

    std::vector<FitRow> local;
    for (size_t k = shard_index; k < entries.size(); k += n_shards)
    {
        const GroupEntry &g = entries[k];

        std::vector<double> response, scale, response_err, scale_err;

        for (int idx : g.idxs)
        {
            const BinRow &r = rows[idx];

            bool valid =
                std::isfinite(r.scale_mu) && std::isfinite(r.response_mu) &&
                std::isfinite(r.scale_mu_err) && std::isfinite(r.response_mu_err) &&
                r.scale_mu > 0 && r.response_mu > 0 &&
                r.scale_mu <= SCALE_SANITY_MAX &&
                r.response_mu <= RESPONSE_SANITY_FACTOR * r.E_high;

            if (!valid)
                continue;

            response.push_back(r.response_mu);
            scale.push_back(r.scale_mu);
            response_err.push_back(r.response_mu_err);
            scale_err.push_back(r.scale_mu_err);
        }

        if (response.size() < 4)
            continue;

        std::vector<size_t> order(response.size());
        for (size_t i = 0; i < order.size(); i++)
            order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](size_t a, size_t b)
                  { return response[a] < response[b]; });

        std::vector<double> response_sorted(order.size()), scale_sorted(order.size()),
            response_err_sorted(order.size()), scale_err_sorted(order.size());

        for (size_t i = 0; i < order.size(); i++)
        {
            response_sorted[i] = response[order[i]];
            scale_sorted[i] = scale[order[i]];
            response_err_sorted[i] = response_err[order[i]];
            scale_err_sorted[i] = scale_err[order[i]];
        }

        FitRow fit;
        bool ok = fit_log_gauss_function_constrained(
            response_sorted, scale_sorted, scale_err_sorted, chi2_improvement_min, fit);

        if (!ok)
            continue;

        fit.x_center = g.x_center;
        fit.y_center = g.y_center;
        local.push_back(fit);
    }
    return local;
}

// Fits every position bin in this process
std::vector<FitRow> compute_all_regions(const std::vector<BinRow> &rows, double chi2_improvement_min)
{
    return compute_shard(rows, chi2_improvement_min, 0, 1);
}

Grid2D make_nan_grid(int nx, int ny)
{
    return Grid2D(nx, std::vector<double>(ny, std::numeric_limits<double>::quiet_NaN()));
}

// The saved map must span every position bin defined in the ROOT tree, not
// just the ones a region fit succeeded for.
std::map<std::string, Grid2D> build_coarse_maps(
    const std::vector<FitRow> &fit_rows,
    const std::vector<double> &xs,
    const std::vector<double> &ys,
    const std::vector<std::string> &params)
{
    std::map<std::pair<double, double>, const FitRow *> lookup;
    for (const auto &r : fit_rows)
        lookup[{r.x_center, r.y_center}] = &r;

    std::map<std::string, Grid2D> maps;
    for (const auto &param : params)
        maps[param] = make_nan_grid((int)xs.size(), (int)ys.size());

    for (size_t i = 0; i < xs.size(); i++)
    {
        for (size_t j = 0; j < ys.size(); j++)
        {
            auto it = lookup.find({xs[i], ys[j]});
            if (it == lookup.end())
                continue;

            const FitRow &r = *it->second;

            for (const auto &param : params)
            {
                double val;
                if (param == "A")
                    val = r.A;
                else if (param == "B")
                    val = r.B;
                else if (param == "C")
                    val = r.C;
                else if (param == "D")
                    val = r.D;
                else if (param == "E0")
                    val = r.E0;
                else if (param == "sigmaL")
                    val = r.sigmaL;
                else if (param == "sigmaR")
                    val = r.sigmaR;
                else // reduced_chi2
                    val = r.chi2 / r.ndof;

                maps[param][i][j] = val;
            }
        }
    }

    return maps;
}

// Fill every NaN cell inside footprint_mask with the single global mean of
// all real (non-NaN) values in Z. Cells outside footprint_mask are left NaN.
Grid2D fill_nan_global_mean(const Grid2D &Z, const std::vector<std::vector<bool>> &footprint_mask)
{
    double sum = 0.0;
    long count = 0;

    for (const auto &row : Z)
        for (double v : row)
            if (std::isfinite(v))
            {
                sum += v;
                count++;
            }

    double global_mean = count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();

    Grid2D out = Z;
    for (size_t i = 0; i < Z.size(); i++)
        for (size_t j = 0; j < Z[i].size(); j++)
            if (!std::isfinite(Z[i][j]) && footprint_mask[i][j])
                out[i][j] = global_mean;

    return out;
}

std::vector<std::vector<bool>> dilate(const std::vector<std::vector<bool>> &mask)
{
    int nx = (int)mask.size();
    int ny = nx > 0 ? (int)mask[0].size() : 0;

    std::vector<std::vector<bool>> out(nx, std::vector<bool>(ny, false));

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        {
            bool any = false;
            for (int di = -1; di <= 1 && !any; di++)
                for (int dj = -1; dj <= 1 && !any; dj++)
                {
                    int ni = i + di, nj = j + dj;
                    if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && mask[ni][nj])
                        any = true;
                }
            out[i][j] = any;
        }

    return out;
}

// Chebyshev (8-connected) grid distance from each cell to the nearest cell
// where Z is finite -- 0 at a finite cell itself, growing by 1 per ring of
// dilate() outward. Used to pull the skew-Gaussian amplitude D back toward 0
// away from real data (see D_DECAY_LENGTH_BINS).
Grid2D chebyshev_distance_to_finite(const Grid2D &Z)
{
    int nx = (int)Z.size();
    int ny = nx > 0 ? (int)Z[0].size() : 0;

    std::vector<std::vector<bool>> mask(nx, std::vector<bool>(ny, false));
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            mask[i][j] = std::isfinite(Z[i][j]);

    Grid2D dist(nx, std::vector<double>(ny, -1.0));
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            if (mask[i][j])
                dist[i][j] = 0.0;

    std::vector<std::vector<bool>> frontier = mask;
    int max_steps = nx + ny; // grid diagonal, an upper bound on any reachable distance
    for (int step = 1; step <= max_steps; step++)
    {
        std::vector<std::vector<bool>> dilated = dilate(frontier);
        bool any_new = false;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                if (dilated[i][j] && dist[i][j] < 0.0)
                {
                    dist[i][j] = (double)step;
                    any_new = true;
                }
        frontier = dilated;
        if (!any_new)
            break;
    }

    return dist;
}

// Neighbor-averaging smoother with MAD-based outlier rejection
Grid2D smooth_map(
    const Grid2D &Z,
    int n_iter = SMOOTH_N_ITER,
    int radius = SMOOTH_RADIUS,
    double outlier_sigma = 3.0,
    int min_neighbors = 0,
    double neighbor_weight = 0.5,
    int max_fill_distance = SMOOTH_MAX_FILL_DISTANCE)
{
    int nx = (int)Z.size();
    int ny = nx > 0 ? (int)Z[0].size() : 0;

    std::vector<std::vector<bool>> original_finite(nx, std::vector<bool>(ny));
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            original_finite[i][j] = std::isfinite(Z[i][j]);

    std::vector<std::pair<int, int>> offsets;
    for (int di = -radius; di <= radius; di++)
        for (int dj = -radius; dj <= radius; dj++)
            if (!(di == 0 && dj == 0))
                offsets.push_back({di, dj});

    Grid2D Z_smooth = Z;

    for (int iter = 0; iter < n_iter; iter++)
    {
        std::vector<std::vector<bool>> outlier(nx, std::vector<bool>(ny, false));

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                double center = Z_smooth[i][j];
                if (!std::isfinite(center))
                    continue;

                std::vector<double> neighbor_only;
                for (const auto &off : offsets)
                {
                    int ni = i + off.first, nj = j + off.second;
                    if (ni >= 0 && ni < nx && nj >= 0 && nj < ny)
                    {
                        double val = Z_smooth[ni][nj];
                        if (std::isfinite(val))
                            neighbor_only.push_back(val);
                    }
                }

                if ((int)neighbor_only.size() < min_neighbors)
                    continue;

                std::vector<double> sorted_vals = neighbor_only;
                std::sort(sorted_vals.begin(), sorted_vals.end());
                size_t n = sorted_vals.size();
                double median = n == 0 ? 0.0 : (n % 2 == 1 ? sorted_vals[n / 2]
                                                            : 0.5 * (sorted_vals[n / 2 - 1] + sorted_vals[n / 2]));

                std::vector<double> abs_dev(neighbor_only.size());
                for (size_t k = 0; k < neighbor_only.size(); k++)
                    abs_dev[k] = std::fabs(neighbor_only[k] - median);
                std::sort(abs_dev.begin(), abs_dev.end());
                double mad = abs_dev.empty() ? 0.0 : (abs_dev.size() % 2 == 1
                                                           ? abs_dev[abs_dev.size() / 2]
                                                           : 0.5 * (abs_dev[abs_dev.size() / 2 - 1] + abs_dev[abs_dev.size() / 2]));
                double sigma = 1.4826 * mad;

                if (sigma < 1e-8)
                {
                    if (std::fabs(center - median) > 1e-6)
                        outlier[i][j] = true;
                }
                else if (std::fabs(center - median) > outlier_sigma * sigma)
                {
                    outlier[i][j] = true;
                }
            }
        }

        Grid2D Z_new = Z_smooth;

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                std::vector<double> values;

                for (const auto &off : offsets)
                {
                    int ni = i + off.first, nj = j + off.second;
                    if (ni >= 0 && ni < nx && nj >= 0 && nj < ny)
                    {
                        if (outlier[ni][nj])
                            continue;

                        double val = Z_smooth[ni][nj];
                        if (std::isfinite(val))
                            values.push_back(val);
                    }
                }

                if (values.empty())
                    continue;

                double neighbor_average = 0.0;
                for (double v : values)
                    neighbor_average += v;
                neighbor_average /= values.size();

                if (!std::isfinite(Z_smooth[i][j]))
                    Z_new[i][j] = neighbor_average;
                else if (outlier[i][j])
                    Z_new[i][j] = neighbor_average;
                else
                    // Every finite, non-outlier cell is blended toward its
                    // neighborhood too: the map is interpolated downstream and
                    // must vary smoothly
                    Z_new[i][j] = (1 - neighbor_weight) * Z_smooth[i][j] + neighbor_weight * neighbor_average;
            }
        }

        Z_smooth = Z_new;
    }

    // Force the finite footprint to exactly match the original (or a small
    // margin around it, if max_fill_distance > 0).
    if (max_fill_distance >= 0)
    {
        std::vector<std::vector<bool>> fill_mask = original_finite;

        for (int k = 0; k < max_fill_distance; k++)
            fill_mask = dilate(fill_mask);

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                if (!fill_mask[i][j])
                    Z_smooth[i][j] = std::numeric_limits<double>::quiet_NaN();
    }

    return Z_smooth;
}

void write_double(std::ostream &os, double v)
{
    char buf[64];
    snprintf(buf, sizeof(buf), "%.17g", v);
    os << buf;
}

void write_calibration_json(
    const std::vector<FitRow> &fit_rows_in,
    const std::vector<double> &xs,
    const std::vector<double> &ys,
    const std::string &save_file)
{
    std::vector<FitRow> fit_rows = fit_rows_in;

    std::vector<std::string> all_params = {"A", "B", "C", "D", "E0", "sigmaL", "sigmaR", "reduced_chi2"};

    std::map<std::string, Grid2D> coarse_maps = build_coarse_maps(fit_rows, xs, ys, all_params);

    long n_skew_gaussian_cells = 0, n_fit_cells = 0;
    for (const auto &row : coarse_maps["D"])
        for (double v : row)
            if (std::isfinite(v))
            {
                n_fit_cells++;
                if (v > SKEW_GAUSSIAN_PRESENT_THRESHOLD)
                    n_skew_gaussian_cells++;
            }

    std::cout << "A real skew Gaussian was fit for " << n_skew_gaussian_cells << "/" << n_fit_cells
              << " calibrated grid cells (rest are background-only)" << std::endl;

    // Unfit cells inherit the neighborhood average via smooth_map. A/B/C/D
    // share one NaN footprint (a FitRow sets all four or none); E0/sigma
    // have a larger one (also NaN for background-only fits), so may not
    // always find a skew-Gaussian-bearing neighbor to average from.
    std::map<std::string, Grid2D> smoothed_maps;
    for (const auto &param : all_params)
        smoothed_maps[param] = smooth_map(coarse_maps[param]);

    // Taper D with distance from the nearest directly fitted cell (a
    // background-only fit counts, with D = 0); fitted cells are unchanged
    {
        Grid2D d_dist = chebyshev_distance_to_finite(coarse_maps["D"]);
        for (size_t i = 0; i < xs.size(); i++)
            for (size_t j = 0; j < ys.size(); j++)
            {
                if (d_dist[i][j] <= 0.0)
                    continue;
                double decay = std::exp(-d_dist[i][j] / D_DECAY_LENGTH_BINS);
                smoothed_maps["D"][i][j] *= decay;
            }
    }

    long n_abcd_nan_before = 0, n_abcd_filled = 0;
    for (size_t i = 0; i < xs.size(); i++)
        for (size_t j = 0; j < ys.size(); j++)
            if (!std::isfinite(coarse_maps["A"][i][j]))
            {
                n_abcd_nan_before++;
                if (std::isfinite(smoothed_maps["A"][i][j]))
                    n_abcd_filled++;
            }

    std::cout << "Filled " << n_abcd_filled << "/" << n_abcd_nan_before
              << " previously-empty position bins with the neighborhood average of nearby fits" << std::endl;

    // Where A/B/C/D are finite but E0/sigma aren't (no skew-Gaussian-bearing neighbor
    // reachable), fall back to the global mean skew-Gaussian E0/sigma rather than
    // leaving NaN (which would zero out a real, filled D skew-Gaussian term at save time).
    std::vector<std::vector<bool>> has_abcd(xs.size(), std::vector<bool>(ys.size()));
    long n_e0_nan_before = 0;
    for (size_t i = 0; i < xs.size(); i++)
        for (size_t j = 0; j < ys.size(); j++)
        {
            has_abcd[i][j] = std::isfinite(smoothed_maps["A"][i][j]) &&
                              std::isfinite(smoothed_maps["B"][i][j]) &&
                              std::isfinite(smoothed_maps["C"][i][j]) &&
                              std::isfinite(smoothed_maps["D"][i][j]);
            if (has_abcd[i][j] && !std::isfinite(smoothed_maps["E0"][i][j]))
                n_e0_nan_before++;
        }

    smoothed_maps["E0"] = fill_nan_global_mean(smoothed_maps["E0"], has_abcd);
    smoothed_maps["sigmaL"] = fill_nan_global_mean(smoothed_maps["sigmaL"], has_abcd);
    smoothed_maps["sigmaR"] = fill_nan_global_mean(smoothed_maps["sigmaR"], has_abcd);

    std::cout << "Filled E0/sigmaL/sigmaR for " << n_e0_nan_before
              << " remaining cells with the global mean of real skew-Gaussian-active values" << std::endl;

    const std::map<std::string, Grid2D> &maps_to_save = smoothed_maps;

    std::vector<std::tuple<int, int, double, double, double, double, double, double, double>> bins_out;
    long n_skew_gaussian_saved = 0, n_unexpected_nonfinite = 0, n_outside_fiducial = 0;

    for (size_t i = 0; i < xs.size(); i++)
    {
        for (size_t j = 0; j < ys.size(); j++)
        {
            // Outside the buffered fiducial region the map is undefined
            if (!pass_fiducial_cut((float)xs[i], (float)ys[j], ACTIVE_FIDUCIAL_BUFFER, *ACTIVE_FIDUCIAL_RECT))
            {
                n_outside_fiducial++;
                continue;
            }

            double A_val = maps_to_save.at("A")[i][j];
            double B_val = maps_to_save.at("B")[i][j];
            double C_val = maps_to_save.at("C")[i][j];

            if (!(std::isfinite(A_val) && std::isfinite(B_val) && std::isfinite(C_val)))
                continue;

            double D_val = maps_to_save.at("D")[i][j];
            double E0_val = maps_to_save.at("E0")[i][j];
            double sigmaL_val = maps_to_save.at("sigmaL")[i][j];
            double sigmaR_val = maps_to_save.at("sigmaR")[i][j];

            if (!(std::isfinite(D_val) && std::isfinite(E0_val) && std::isfinite(sigmaL_val) && std::isfinite(sigmaR_val)))
            {
                n_unexpected_nonfinite++;
                D_val = 0.0;
                E0_val = 0.0;
                sigmaL_val = 0.0;
                sigmaR_val = 0.0;
            }

            if (D_val > SKEW_GAUSSIAN_PRESENT_THRESHOLD)
                n_skew_gaussian_saved++;

            bins_out.emplace_back((int)i, (int)j, A_val, B_val, C_val, D_val, E0_val, sigmaL_val, sigmaR_val);
        }
    }

    if (n_unexpected_nonfinite > 0)
        std::cout << "WARNING: " << n_unexpected_nonfinite
                  << " bins had non-finite D/E0/sigmaL/sigmaR despite finite A/B/C -- forced to 0" << std::endl;

    std::cout << "Excluded " << n_outside_fiducial
              << " grid cells outside the buffer-restricted fiducial acceptance (never defined, regardless of smoothing)" << std::endl;

    long n_total = (long)xs.size() * (long)ys.size();
    std::cout << "Saving " << bins_out.size() << "/" << n_total
              << " bins with finite A,B,C (dropping " << (n_total - (long)bins_out.size())
              << " with no calibration)" << std::endl;
    std::cout << n_skew_gaussian_saved << "/" << bins_out.size()
              << " saved bins include an actual skew Gaussian term" << std::endl;

    std::ofstream out(save_file.c_str());
    if (!out.is_open())
    {
        std::cerr << "ERROR: cannot open " << save_file << " for writing" << std::endl;
        exit(1);
    }

    out << "{\n";
    out << "  \"nx\": " << xs.size() << ",\n";
    out << "  \"ny\": " << ys.size() << ",\n";

    out << "  \"x_grid\": [";
    for (size_t i = 0; i < xs.size(); i++)
    {
        if (i > 0)
            out << ", ";
        write_double(out, xs[i]);
    }
    out << "],\n";

    out << "  \"y_grid\": [";
    for (size_t j = 0; j < ys.size(); j++)
    {
        if (j > 0)
            out << ", ";
        write_double(out, ys[j]);
    }
    out << "],\n";

    out << "  \"bins\": [\n";
    for (size_t k = 0; k < bins_out.size(); k++)
    {
        const auto &b = bins_out[k];
        out << "    {\"i\": " << std::get<0>(b) << ", \"j\": " << std::get<1>(b)
            << ", \"A\": ";
        write_double(out, std::get<2>(b));
        out << ", \"B\": ";
        write_double(out, std::get<3>(b));
        out << ", \"C\": ";
        write_double(out, std::get<4>(b));
        out << ", \"D\": ";
        write_double(out, std::get<5>(b));
        out << ", \"E0\": ";
        write_double(out, std::get<6>(b));
        out << ", \"sigmaL\": ";
        write_double(out, std::get<7>(b));
        out << ", \"sigmaR\": ";
        write_double(out, std::get<8>(b));
        out << "}" << (k + 1 < bins_out.size() ? "," : "") << "\n";
    }
    out << "  ]\n";
    out << "}\n";

    out.close();

    std::cout << "Lookup table saved to " << save_file << std::endl;
}

// Full-precision FitRow serialization for --shard / --merge
void write_fitrow_csv(const std::vector<FitRow> &fit_rows, const std::string &save_file)
{
    std::ofstream out(save_file.c_str());
    if (!out.is_open())
    {
        std::cerr << "ERROR: cannot open " << save_file << " for writing" << std::endl;
        exit(1);
    }
    out << "x_center,y_center,A,B,C,D,E0,sigmaL,sigmaR,chi2,n_params,ndof\n";
    for (const auto &r : fit_rows)
    {
        write_double(out, r.x_center);
        out << ",";
        write_double(out, r.y_center);
        out << ",";
        write_double(out, r.A);
        out << ",";
        write_double(out, r.B);
        out << ",";
        write_double(out, r.C);
        out << ",";
        write_double(out, r.D);
        out << ",";
        write_double(out, r.E0);
        out << ",";
        write_double(out, r.sigmaL);
        out << ",";
        write_double(out, r.sigmaR);
        out << ",";
        write_double(out, r.chi2);
        out << "," << r.n_params << "," << r.ndof << "\n";
    }
    out.close();
}

std::vector<FitRow> read_fitrow_csv(const std::string &path)
{
    std::ifstream in(path.c_str());
    if (!in.is_open())
    {
        std::cerr << "ERROR: cannot open shard file " << path << " for reading" << std::endl;
        exit(1);
    }
    std::vector<FitRow> rows;
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string tok;
        std::vector<double> v;
        while (std::getline(ss, tok, ','))
            v.push_back(std::stod(tok));
        if (v.size() != 12)
        {
            std::cerr << "ERROR: malformed row in " << path << " (expected 12 fields, got " << v.size() << ")" << std::endl;
            exit(1);
        }
        FitRow r;
        r.x_center = v[0];
        r.y_center = v[1];
        r.A = v[2];
        r.B = v[3];
        r.C = v[4];
        r.D = v[5];
        r.E0 = v[6];
        r.sigmaL = v[7];
        r.sigmaR = v[8];
        r.chi2 = v[9];
        r.n_params = (int)v[10];
        r.ndof = (int)v[11];
        rows.push_back(r);
    }
    return rows;
}

// Side CSV of per-bin chi2/ndof of the final model (diagnostics only)
void write_chi2_csv(const std::vector<FitRow> &fit_rows, const std::string &save_file)
{
    std::ofstream out(save_file.c_str());
    out << "x,y,chi2,ndof,n_params,reduced_chi2\n";
    for (const auto &r : fit_rows)
    {
        if (r.ndof <= 0)
            continue;
        out << r.x_center << "," << r.y_center << "," << r.chi2 << "," << r.ndof << ","
            << r.n_params << "," << (r.chi2 / r.ndof) << "\n";
    }
    out.close();
    std::cout << "chi2/ndof CSV saved to " << save_file << std::endl;
}

// Side CSV of the direct (pre-smoothing) per-position fit parameters
// (diagnostics only)
void write_direct_fits_csv(const std::vector<FitRow> &fit_rows, const std::string &save_file)
{
    std::ofstream out(save_file.c_str());
    out << "x,y,A,B,C,D,E0,sigmaL,sigmaR\n";
    for (const auto &r : fit_rows)
    {
        out << r.x_center << "," << r.y_center << "," << r.A << "," << r.B << "," << r.C << ","
            << r.D << "," << r.E0 << "," << r.sigmaL << "," << r.sigmaR << "\n";
    }
    out.close();
    std::cout << "Direct (unsmoothed) fit CSV saved to " << save_file << std::endl;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning; // silence Minuit2's per-fit convergence chatter

    std::string input_file = ROOT_FILE;
    std::string output_file = OUTPUT_JSON;

    if (argc > 1)
        input_file = argv[1];
    if (argc > 2)
        output_file = argv[2];

    // argv[3] and argv[4] are ignored
    const bool em_only = true;
    E0_MAX = em_only ? E0_MAX_em_only : E0_MAX_ecal_hcal;
    SIGMA_MAX = em_only ? SIGMA_MAX_em_only : SIGMA_MAX_ecal_hcal;
    ACTIVE_FIDUCIAL_RECT = em_only ? &kFiducialRectEcal : &kFiducialRect;
    ACTIVE_FIDUCIAL_BUFFER = em_only ? reco_fiducial_buffer_ecal : reco_fiducial_buffer;

    // argv[5]: --shard <i> <n> <shard.csv> fits only shard i of n and writes
    // it; --merge <shard_list> <out.json> concatenates the shard CSVs and runs
    // the smoothing and save steps. Omitted = fit everything in one process.
    if (argc > 5 && std::string(argv[5]) == "--shard")
    {
        if (argc < 9)
        {
            std::cerr << "ERROR: --shard requires <shard_index> <n_shards> <shard_output_csv>" << std::endl;
            return 1;
        }
        unsigned shard_index = (unsigned)std::atoi(argv[6]);
        unsigned n_shards = (unsigned)std::atoi(argv[7]);
        std::string shard_output_csv = argv[8];

        if (n_shards == 0 || shard_index >= n_shards)
        {
            std::cerr << "ERROR: need 0 <= shard_index < n_shards (got shard_index=" << shard_index
                      << " n_shards=" << n_shards << ")" << std::endl;
            return 1;
        }

        std::cout << "Loading " << input_file << " ..." << std::endl;
        std::vector<BinRow> rows = load_bins(input_file);
        std::cout << "Loaded " << rows.size() << " (E,x,y) bin entries" << std::endl;

        double chi2_improvement_min = TMath::ChisquareQuantile(SKEW_GAUSSIAN_SIGNIFICANCE_CL, 4);

        std::cout << "Fitting shard " << shard_index << "/" << n_shards << " ..." << std::endl;
        std::vector<FitRow> shard_fit_rows = compute_shard(rows, chi2_improvement_min, shard_index, n_shards);
        std::cout << "Fit " << shard_fit_rows.size() << " position bins in this shard" << std::endl;

        write_fitrow_csv(shard_fit_rows, shard_output_csv);
        std::cout << "Shard results saved to " << shard_output_csv << std::endl;

        return 0;
    }

    if (argc > 5 && std::string(argv[5]) == "--merge")
    {
        if (argc < 8)
        {
            std::cerr << "ERROR: --merge requires <shard_list_file> <output_file>" << std::endl;
            return 1;
        }
        std::string shard_list_file = argv[6];
        std::string merge_output_file = argv[7];

        std::vector<FitRow> fit_rows;
        std::ifstream list_in(shard_list_file.c_str());
        if (!list_in.is_open())
        {
            std::cerr << "ERROR: cannot open shard list file " << shard_list_file << std::endl;
            return 1;
        }
        std::string shard_path;
        int n_shard_files = 0;
        while (std::getline(list_in, shard_path))
        {
            if (shard_path.empty())
                continue;
            std::vector<FitRow> shard_rows = read_fitrow_csv(shard_path);
            fit_rows.insert(fit_rows.end(), shard_rows.begin(), shard_rows.end());
            n_shard_files++;
        }
        std::cout << "Merged " << fit_rows.size() << " position-bin fits from " << n_shard_files << " shard files" << std::endl;

        std::cout << "Loading " << input_file << " for the full x/y grid ..." << std::endl;
        std::vector<BinRow> rows = load_bins(input_file);

        std::set<double> xs_set, ys_set;
        for (const auto &r : rows)
        {
            xs_set.insert(r.x_center);
            ys_set.insert(r.y_center);
        }
        std::vector<double> xs(xs_set.begin(), xs_set.end());
        std::vector<double> ys(ys_set.begin(), ys_set.end());

        write_calibration_json(fit_rows, xs, ys, merge_output_file);
        write_chi2_csv(fit_rows, merge_output_file + ".chi2.csv");
        write_direct_fits_csv(fit_rows, merge_output_file + ".direct_fits.csv");

        return 0;
    }

    std::cout << "Loading " << input_file << " ..." << std::endl;
    std::vector<BinRow> rows = load_bins(input_file);
    std::cout << "Loaded " << rows.size() << " (E,x,y) bin entries" << std::endl;

    // 4 dof: D, E0, sigmaL, sigmaR over the background
    double chi2_improvement_min = TMath::ChisquareQuantile(SKEW_GAUSSIAN_SIGNIFICANCE_CL, 4);

    std::cout << "Fitting position bins ..." << std::endl;
    std::vector<FitRow> fit_rows = compute_all_regions(rows, chi2_improvement_min);
    std::cout << "Fit " << fit_rows.size() << " position bins" << std::endl;

    std::set<double> xs_set, ys_set;
    for (const auto &r : rows)
    {
        xs_set.insert(r.x_center);
        ys_set.insert(r.y_center);
    }
    std::vector<double> xs(xs_set.begin(), xs_set.end());
    std::vector<double> ys(ys_set.begin(), ys_set.end());

    write_calibration_json(fit_rows, xs, ys, output_file);
    write_chi2_csv(fit_rows, output_file + ".chi2.csv");

    return 0;
}
