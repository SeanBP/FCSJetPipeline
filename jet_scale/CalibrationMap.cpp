// Builds the JES lookup JSON from a jet_calibration_grid*.root "bins" tree:
// per (x,y) bin, fit <scale> vs <response> (log-quadratic + optional
// Gaussian bump), then smooth the A/B/C/D/E0/sigma maps over the grid.
// C++ port of CalibrationMap.py -- same algorithm, no plotting.

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

// CONFIGURATION (mirrors CalibrationMap.py)

std::string ROOT_FILE = "jet_calibration_grid_test.root";
std::string OUTPUT_JSON = "JetEnergyScale_lookup_GridTest.json";

// Below this many valid energy points, a 6-parameter bump fit would have
// non-positive dof, so only the 3-parameter background is fit instead.
const int MIN_POINTS_FOR_BUMP = 7;

// Sanity cuts dropping the noise tail from bins with too few jets for a
// meaningful fitted mean (see CalibrationMap.py for derivation).
const double RESPONSE_SANITY_FACTOR = 1.5;
const double SCALE_SANITY_MAX = 1.5;

// Gaussian bump bounds; sigma's lower bound is pinned above zero (it's a denominator).
const double D_MIN = 0.0, D_MAX = 0.2;
const double E0_MIN = 10.0, E0_MAX = 35.0;
const double SIGMA_MIN = 1e-3, SIGMA_MAX = 25.0;

// Bump kept only if its chi2 improvement beats the chance level for a
// chi-square(3 dof) distribution (Wilks' theorem).
const double BUMP_SIGNIFICANCE_CL = 0.5;

// "a few percent" of each parameter's own range; a fit within this of
// either edge is treated as having hit the wall, not converged.
const double BOUND_HIT_RTOL = 0.03;

// Below this, D counts as "no bump" for reporting purposes.
const double BUMP_PRESENT_THRESHOLD = 1e-4;

// Map smoothing parameters (see plot_and_save_smoothed's call in
// CalibrationMap.py).
const int SMOOTH_N_ITER = 6;
const int SMOOTH_RADIUS = 2;
// -1 = no fill-distance clip: an unfit bin inherits smooth_map's diffusion
// average instead of staying NaN (bins outside acceptance still end up NaN,
// since diffusion has nothing finite nearby to average from).
const int SMOOTH_MAX_FILL_DISTANCE = -1;

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
    double A, B, C, D, E0, sigma;
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
// lower check -- for D, lower bound 0 is a legitimate "no bump" outcome, not degenerate.
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

// Fits background, then background+bump if enough points (keeping the bump
// only if significant and not bound-pinned). False if too few points to fit at all.
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
    fit_with_restarts(bgFunc, {{0.2, 0.01, 0.01}}, lower_bg, upper_bg, bgX, bgF);

    auto background_only = [&](const std::string &reason)
    {
        (void)reason;
        out.A = bgX[0];
        out.B = bgX[1];
        out.C = bgX[2];
        out.D = 0.0;
        out.E0 = std::numeric_limits<double>::quiet_NaN();
        out.sigma = std::numeric_limits<double>::quiet_NaN();
        out.chi2 = bgF;
        out.n_params = 3;
        out.ndof = (int)response.size() - 3;
    };

    if ((int)response.size() < MIN_POINTS_FOR_BUMP)
    {
        background_only("few_points");
        return true;
    }

    double A0 = bgX[0], B0 = bgX[1], C0 = bgX[2];

    std::vector<double> e0_candidates = linspace(E0_MIN, E0_MAX, 3);
    std::vector<double> sigma_candidates = {SIGMA_MAX * 0.3, SIGMA_MAX * 0.6};

    GaussChi2 gaussChi2{response, scale, scale_err};
    ROOT::Math::Functor gaussFunc(gaussChi2, 6);

    std::vector<double> lower_g = {0.0, NAN, NAN, D_MIN, E0_MIN, SIGMA_MIN};
    std::vector<double> upper_g = {NAN, NAN, NAN, D_MAX, E0_MAX, SIGMA_MAX};

    std::vector<std::vector<double>> x0_list;
    for (double e0 : e0_candidates)
        for (double sigma : sigma_candidates)
            x0_list.push_back({A0, B0, C0, 0.5 * D_MAX, e0, sigma});

    std::vector<double> gX;
    double gF;
    fit_with_restarts(gaussFunc, x0_list, lower_g, upper_g, gX, gF);

    double D = gX[3], E0v = gX[4], sigma = gX[5];
    double delta_chi2 = bgF - gF;

    bool bump_at_bound =
        at_bound(D, D_MIN, D_MAX, false) ||
        at_bound(E0v, E0_MIN, E0_MAX, true) ||
        at_bound(sigma, SIGMA_MIN, SIGMA_MAX, true);

    if (bump_at_bound)
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
    out.sigma = gX[5];
    out.chi2 = gF;
    out.n_params = 6;
    out.ndof = (int)response.size() - 6;

    return true;
}

// Per (x,y) bin, fit <scale> vs <response>; bins with <4 valid entries
// (or <3 finite after masking) are dropped entirely.
std::vector<FitRow> compute_all_regions(const std::vector<BinRow> &rows, double chi2_improvement_min)
{
    std::map<std::pair<double, double>, std::vector<int>> groups;

    for (size_t i = 0; i < rows.size(); i++)
        groups[{rows[i].x_center, rows[i].y_center}].push_back((int)i);

    std::vector<FitRow> fit_rows;
    fit_rows.reserve(groups.size());

    for (const auto &kv : groups)
    {
        double x_center = kv.first.first;
        double y_center = kv.first.second;

        std::vector<double> response, scale, response_err, scale_err;

        for (int idx : kv.second)
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

        fit.x_center = x_center;
        fit.y_center = y_center;
        fit_rows.push_back(fit);
    }

    return fit_rows;
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
                else if (param == "sigma")
                    val = r.sigma;
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

// Neighbor-averaging smoother with MAD-based outlier rejection (see
// CalibrationMap.py's smooth_map docstring for full parameter rationale).
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

    std::vector<std::string> all_params = {"A", "B", "C", "D", "E0", "sigma", "reduced_chi2"};

    std::map<std::string, Grid2D> coarse_maps = build_coarse_maps(fit_rows, xs, ys, all_params);

    long n_bump_cells = 0, n_fit_cells = 0;
    for (const auto &row : coarse_maps["D"])
        for (double v : row)
            if (std::isfinite(v))
            {
                n_fit_cells++;
                if (v > BUMP_PRESENT_THRESHOLD)
                    n_bump_cells++;
            }

    std::cout << "A real bump was fit for " << n_bump_cells << "/" << n_fit_cells
              << " calibrated grid cells (rest are background-only)" << std::endl;

    // Unfit cells inherit the neighborhood average via smooth_map. A/B/C/D
    // share one NaN footprint (a FitRow sets all four or none); E0/sigma
    // have a larger one (also NaN for background-only fits), so may not
    // always find a bump-bearing neighbor to average from.
    std::map<std::string, Grid2D> smoothed_maps;
    for (const auto &param : all_params)
        smoothed_maps[param] = smooth_map(coarse_maps[param]);

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

    // Where A/B/C/D are finite but E0/sigma aren't (no bump-bearing neighbor
    // reachable), fall back to the global mean bump E0/sigma rather than
    // leaving NaN (which would zero out a real, filled D bump term at save time).
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
    smoothed_maps["sigma"] = fill_nan_global_mean(smoothed_maps["sigma"], has_abcd);

    std::cout << "Filled E0/sigma for " << n_e0_nan_before
              << " remaining cells with the global mean of real bump-active values" << std::endl;

    const std::map<std::string, Grid2D> &maps_to_save = smoothed_maps;

    std::vector<std::tuple<int, int, double, double, double, double, double, double>> bins_out;
    long n_bump_saved = 0, n_unexpected_nonfinite = 0;

    for (size_t i = 0; i < xs.size(); i++)
    {
        for (size_t j = 0; j < ys.size(); j++)
        {
            double A_val = maps_to_save.at("A")[i][j];
            double B_val = maps_to_save.at("B")[i][j];
            double C_val = maps_to_save.at("C")[i][j];

            if (!(std::isfinite(A_val) && std::isfinite(B_val) && std::isfinite(C_val)))
                continue;

            double D_val = maps_to_save.at("D")[i][j];
            double E0_val = maps_to_save.at("E0")[i][j];
            double sigma_val = maps_to_save.at("sigma")[i][j];

            if (!(std::isfinite(D_val) && std::isfinite(E0_val) && std::isfinite(sigma_val)))
            {
                n_unexpected_nonfinite++;
                D_val = 0.0;
                E0_val = 0.0;
                sigma_val = 0.0;
            }

            if (D_val > BUMP_PRESENT_THRESHOLD)
                n_bump_saved++;

            bins_out.emplace_back((int)i, (int)j, A_val, B_val, C_val, D_val, E0_val, sigma_val);
        }
    }

    if (n_unexpected_nonfinite > 0)
        std::cout << "WARNING: " << n_unexpected_nonfinite
                  << " bins had non-finite D/E0/sigma despite finite A/B/C -- forced to 0" << std::endl;

    long n_total = (long)xs.size() * (long)ys.size();
    std::cout << "Saving " << bins_out.size() << "/" << n_total
              << " bins with finite A,B,C (dropping " << (n_total - (long)bins_out.size())
              << " with no calibration)" << std::endl;
    std::cout << n_bump_saved << "/" << bins_out.size()
              << " saved bins include an actual bump term" << std::endl;

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
        out << ", \"sigma\": ";
        write_double(out, std::get<7>(b));
        out << "}" << (k + 1 < bins_out.size() ? "," : "") << "\n";
    }
    out << "  ]\n";
    out << "}\n";

    out.close();

    std::cout << "Lookup table saved to " << save_file << std::endl;
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

    std::cout << "Loading " << input_file << " ..." << std::endl;
    std::vector<BinRow> rows = load_bins(input_file);
    std::cout << "Loaded " << rows.size() << " (E,x,y) bin entries" << std::endl;

    double chi2_improvement_min = TMath::ChisquareQuantile(BUMP_SIGNIFICANCE_CL, 3);

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

    return 0;
}
