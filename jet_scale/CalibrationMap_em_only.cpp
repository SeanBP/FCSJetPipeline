// EM-only version of CalibrationMap.cpp (same folder) -- step 2 of the
// EM-only JES calibration (see the README's jet_scale_em_only section). Kept as a
// separate file so the ECAL+HCAL production step is untouched.
//
// Builds the JES lookup JSON from the "bins" tree written by
// JetEnergyScaleFineGrid_em_only: per (x,y) position bin, fit <scale> vs
// <response> with an inverse-power background (A + B/E + C/E^2) plus an
// asymmetric skew Gaussian (D * exp(-(E-E0)^2 / 2 sigma^2), sigma = sigmaL
// below E0 and sigmaR above), then smooth the A/B/C/D/E0/sigmaL/sigmaR maps
// over the grid and mask cells outside the buffer-restricted ECAL
// fiducial acceptance. The JSON's bins carry {i, j, A, B, C, D, E0, sigmaL,
// sigmaR} -- NOT the production {..., sigma} log-quadratic schema, so this
// JSON cannot be read by jet_calibration/ApplyCorrections.cpp.
//
// Differs from the production CalibrationMap.cpp in: the model form above;
// EM-only skew-Gaussian position/width domain (E0_MAX/SIGMA_MAX, hardwired
// to the EM-only values in main()); the fiducial mask; and the position-bin
// fit being shardable across condor jobs (--shard/--merge in main(), driven
// by JetScalePlayground/submit_calibration_map_diag_condor.sh; the pipeline's
// jet_scale_em_only.xml just runs it serially, ~25 min -- in-process fork()
// parallelism is deliberately disabled, see compute_all_regions).
// --merge also writes <output>.chi2.csv and <output>.direct_fits.csv (the
// pre-smoothing per-position fit) next to the JSON, for diagnostics.
//
// Its input scale/response values are per-bin MEDIANS, not Gaussian-fit
// modes -- this file just fits whatever (response_mu, scale_mu) points it
// is given, so that choice lives entirely in
// JetEnergyScaleFineGrid_em_only.cpp (Cukierman & Nachman,
// arXiv:1609.05195).

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
#include <thread>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/wait.h>

// Forward declaration: defined below, needed here for compute_all_regions'
// lossless (round-trippable) shard-file serialization.
void write_double(std::ostream &os, double v);

// CONFIGURATION (mirrors CalibrationMap.py)

std::string ROOT_FILE = "jet_calibration_grid_test.root";
std::string OUTPUT_JSON = "JetEnergyScale_lookup_GridTest.json";

// Below this many valid energy points, a 7-parameter skew-Gaussian fit would have
// non-positive dof, so only the 3-parameter background is fit instead.
const int MIN_POINTS_FOR_SKEW_GAUSSIAN = 8;

// Sanity cuts dropping the noise tail from bins with too few jets for a
// meaningful fitted mean (see CalibrationMap.py for derivation).
const double RESPONSE_SANITY_FACTOR = 1.5;
const double SCALE_SANITY_MAX = 1.5;

// Skew Gaussian bounds; sigma's lower bound is pinned above zero (it's a
// denominator). D_MIN/D_MAX (skew-Gaussian amplitude, a multiplicative scale
// correction) and E0_MIN/SIGMA_MIN aren't tied to the absolute response
// energy scale, so they're unconditional; E0_MAX/SIGMA_MAX bound where in
// the response energy domain the skew Gaussian can sit, so they must fit within
// whatever Emax that domain actually has -- set in main() to the EM-only
// values (JetEnergyScaleFineGrid_em_only.cpp's Emax=100).
//
// All three of D_MAX/E0_MIN/SIGMA_MAX_em_only were quietly clipping a real,
// highly significant skew Gaussian: a direct chi2 scan (outside this pipeline, on
// the r10 production sample's hottest bin, x=45.5 y=0.8) walked the bounds
// out one at a time and each relaxation immediately pushed the converged
// fit further, meaning Migrad kept landing exactly on whichever wall was
// closest rather than at a genuine interior optimum -- and at_bound()
// correctly discarded every one of those as boundary-pinned, silently
// falling back to background-only. The true unconstrained optimum for that
// bin is D~0.206, E0~5.6-6.1 GeV, sigma~17.1-17.3 GeV (chi2 182.6->13.7 for
// +3 params, chi2/ndof~1.0 -- a genuine, well-converged minimum, not
// another wall). E0~6 sitting below our lowest E bin (~11.6 GeV) is a
// plausible low-energy threshold/turn-on effect extending past the edge of
// our data, not a fit degeneracy. Set with real headroom above that
// observed optimum (not just barely past it) so other bins' skew Gaussians, which
// may differ, also get a fair chance to converge to an interior value
// instead of hitting the same trap. If bins start pinning at these new
// walls too, that's the same signal to investigate/relax further, not
// evidence the fit is broken.
const double D_MIN = 0.0, D_MAX = 0.5;
const double E0_MIN = 1.0;
double E0_MAX = 35.0;
const double E0_MAX_ecal_hcal = 35.0, E0_MAX_em_only = 40.0;
// 1.0 GeV, not ~0: adjacent E-bins near the low-E turn-on are only ~1-2.5
// GeV apart, so an unconstrained sigma can shrink well below that spacing
// and fit a near-step transition pinned to a single data point -- smooth to
// the eye of a Gaussian only in the mathematical sense, not a genuine
// resolvable feature (seen directly on x=75.5,y=-39.2: sigmaL collapsed to
// 0.026, producing a visibly sharp kink at E0 despite a technically
// continuous formula). Flooring sigma at roughly the finest bin spacing
// keeps the skew Gaussian to features the data can actually constrain.
const double SIGMA_MIN = 1.0;
double SIGMA_MAX = 25.0;
const double SIGMA_MAX_ecal_hcal = 25.0, SIGMA_MAX_em_only = 35.0;

// Skew Gaussian kept only if its chi2 improvement beats the chance level for a
// chi-square(3 dof) distribution (Wilks' theorem).
const double SKEW_GAUSSIAN_SIGNIFICANCE_CL = 0.5;

// Duplicated from shared/JetParameters.h (source of truth -- keep these in
// sync if the geometry survey or buffer values there ever change) rather
// than #including that header directly, since it pulls in fastjet, a
// dependency this ROOT/Minuit2-only diagnostic has no other reason to link.
//
// JetEnergyScaleFineGrid.cpp (the aggregation step upstream of this file)
// already gates every grid point through this exact same
// rect+buffer-restricted pass_fiducial_cut before it will accumulate a
// single raw entry there (see its grid_valid array) -- so any grid point
// failing this same check is guaranteed, by construction, to have zero
// statistics, not just an unlucky zero this run. Without this, CalibrationMap
// had no way to tell "genuinely outside the reconstructable region" apart
// from "just empty by chance," and smooth_map's undecayed diffusion (see
// SMOOTH_MAX_FILL_DISTANCE) extrapolated confident-looking A/B/C values (and,
// before D_DECAY_LENGTH_BINS/clamp_isolated_D_islands, D too) tens of cm past
// the true buffer-restricted acceptance -- caught directly by comparing the
// r10 production constants map (extrapolated out to y=-96.7, the outer edge
// of the defined grid) against the chi2 heatmap (whose real per-bin fits
// stop at y=-76.7, matching cut_y_min_ecal+reco_fiducial_buffer_ecal
// almost exactly). Applied as the final filter on the saved calibration
// (see write_calibration_json), independent of and in addition to the
// existing D-specific decay/clamp logic, which still does useful work
// *within* this true acceptance region (filling statistically-empty-by-chance
// cells from real neighbors, excising the real-but-unstable patches) --
// this mask only stops the map from being defined at all past the boundary
// those mechanisms were never told about.
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

// Set in main() from the em_only flag, alongside E0_MAX/SIGMA_MAX -- points
// at kFiducialRect/reco_fiducial_buffer (ECAL+HCAL) or
// kFiducialRectEcal/reco_fiducial_buffer_ecal (EM-only).
const FiducialRect *ACTIVE_FIDUCIAL_RECT = &kFiducialRect;
float ACTIVE_FIDUCIAL_BUFFER = reco_fiducial_buffer;

inline bool pass_fiducial_cut(float jetXE, float jetYE, float buffer, const FiducialRect &rect)
{
    const float ax = std::fabs(jetXE);
    const bool in_x = ax >= (rect.x_inner + buffer) && ax <= (rect.x_outer - buffer);
    const bool in_y = jetYE >= (rect.y_min + buffer) && jetYE <= (rect.y_max - buffer);
    return in_x && in_y;
}

// "a few percent" of each parameter's own range; a fit within this of
// either edge is treated as having hit the wall, not converged.
const double BOUND_HIT_RTOL = 0.03;

// Below this, D counts as "no skew Gaussian" for reporting purposes.
const double SKEW_GAUSSIAN_PRESENT_THRESHOLD = 1e-4;

// Map smoothing parameters (see plot_and_save_smoothed's call in
// CalibrationMap.py).
const int SMOOTH_N_ITER = 6;
const int SMOOTH_RADIUS = 2;
// -1 = no fill-distance clip: an unfit bin inherits smooth_map's diffusion
// average instead of staying NaN (bins outside acceptance still end up NaN,
// since diffusion has nothing finite nearby to average from).
const int SMOOTH_MAX_FILL_DISTANCE = -1;

// The skew-Gaussian amplitude D is only meaningful where there's a genuine
// direct fit behind it; unlike A/B/C (a smoothly-varying background that's
// reasonable to extrapolate everywhere), a bump with no supporting data has
// no business persisting. smooth_map's diffusion has no zero prior though,
// so cells with zero raw statistics (never fit at all -- distinct from a
// real background-only fit, which already has finite D=0) can inherit an
// elevated D purely from repeated neighbor-averaging against other equally
// unfit cells, worst at grid corners where that self-reinforcement has
// nothing real nearby to pull it back down. Confirmed directly on the r10
// production map: its highest D values (~0.25, vs. a real interior median
// ~0.11) sit exactly at cells with n_entries=0. D_DECAY_LENGTH_BINS sets how
// many grid cells (2.5cm each) the exponential decay takes to fall off once
// away from the nearest cell that had an actual fit.
const double D_DECAY_LENGTH_BINS = 3.0;

// The r10 production map's two position-bin patches (upper-right
// ~x:100-118,y:55-70 and lower-right ~x:100-118,y:-56..-80) that once looked
// topologically separated from the main D blob, with elevated but still
// real chi2 (up to ~19 vs. ~2-3 in the main region), turned out to just be
// a real feature of the calibration surface, not a spurious rival region --
// confirmed once statistics were doubled (see git history around
// clamp_isolated_D_islands/fill_isolated_D_holes for the watershed-clamp +
// chi2-rescue + hole-fill machinery this used to require; all removed once
// the "confined to one region" premise it was built on turned out to be
// wrong). D is now smoothed continuously across the whole grid exactly like
// A/B/C, with no spatial exclusion.

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

// Inverse-power background (A + B/E + C/E^2), not the log-quadratic
// (A + B*ln(E) + C*ln(E)^2) used previously. A log-quadratic is a parabola
// in ln(E): once its curvature C is negative it MUST keep monotonically
// decreasing past its vertex, so it can never asymptote to the flat
// high-E plateau the data actually shows (confirmed by direct comparison
// on the x=75.5,y=-39.2 bin: log-quad predicted scale=0.895 at E=150 vs.
// the data's ~0.98-1.0 plateau; inverse-power predicted 0.982, matching,
// and fit chi2 dropped 57.2->32.6 on the same bin). A naturally
// asymptotes to the high-E plateau value as E->infinity, which is exactly
// the behavior needed here. See EM-only-truth-jet-fix session notes for
// the full derivation; this is diagnostics-fork-only for now, not yet
// carried into production CalibrationMap.cpp or ApplyCorrections.cpp.
double bg_model(double E, double A, double B, double C)
{
    return A + B / E + C / (E * E);
}

// Asymmetric skew Gaussian: independent rise (E<E0) and fall (E>=E0) widths. Survey
// of high-statistics/high-chi2 interior bins (>=20cm from the fiducial
// edge) found 33% dominated by one shape -- a fast rise then a much slower
// decline -- that a symmetric Gaussian cannot represent (confirmed via
// direct chi2 scan: for such bins the optimizer keeps pushing E0 toward its
// lower bound rather than settling on a genuine interior optimum, since a
// wide/low-E0 symmetric Gaussian is the closest symmetric approximation to
// an asymmetric rise/fall and that approximation only improves monotonically
// as E0->bound, with no separate competitive local minimum away from it --
// i.e. a real structural mismatch, not a seeding failure). Splitting sigma
// into sigmaL/sigmaR resolved it cleanly on every test bin (reduced chi2
// dropped by 2-3x, landing on non-degenerate, non-bound-pinned parameters).
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
    // A is now the high-E asymptotic scale (should sit near 1.0 for a
    // roughly-calibrated detector), not a log-quadratic constant term --
    // seeded accordingly, with B/C given O(1) steps matching the scale
    // found in the inverse-power vs. log-quadratic comparison test.
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

    // Densified after finding x=30.5,y=0.8's neighborhood (28.0/25.5/33.0,0.8)
    // converging to a narrow, early (E0~11-15) local optimum that fit none of
    // those bins' own data well (raw reduced chi2 up to ~28, one falling back
    // to background-only entirely), while the original sparse E0 grid --
    // {1, 20.5, 40} for em_only's [1,40] domain -- had no candidate anywhere
    // near the true, later-peaking optimum (E0~34, reduced chi2 3.3) that a
    // couple of neighboring bins did manage to find. Smoothing then diffused
    // that one good fit toward the bad neighborhood consensus, producing a
    // saved curve that missed x=30.5,y=0.8's own data by 58x in chi2 (see
    // session notes).
    //
    // A first attempt REPLACED the 3-point {1,20.5,40}/{0.15,0.3,0.6}*SIGMA_MAX
    // grid with a denser 6-point/4-point one -- that improved the mean raw
    // reduced chi2 across all 3443 bins (4.77->3.67), but because the new
    // discretization doesn't contain the old one, a couple of previously-good
    // bins (33.0,0.8: 8.85->14.0; 28.0,0.8: 17.9->20.5) lost the exact seed
    // that used to work for them and landed somewhere worse instead -- a
    // seed *swap* isn't guaranteed monotonic. Using the UNION of the old and
    // new candidate sets instead is: every bin can still start from any seed
    // the original grid offered, plus the new denser/wider ones, so no bin
    // can do worse than the original 3x3x3 grid already found.
    std::vector<double> e0_candidates = linspace(E0_MIN, E0_MAX, 3);
    {
        std::vector<double> e0_dense = linspace(E0_MIN, E0_MAX, 6);
        e0_candidates.insert(e0_candidates.end(), e0_dense.begin(), e0_dense.end());
        // Both linspace calls share the same E0_MIN/E0_MAX endpoints, so the
        // first/last entries duplicate exactly -- drop them rather than
        // wastefully refitting the same seed twice.
        std::sort(e0_candidates.begin(), e0_candidates.end());
        e0_candidates.erase(std::unique(e0_candidates.begin(), e0_candidates.end()), e0_candidates.end());
    }
    // Narrow/wide pair, not just a single scale, and every (sigmaL,sigmaR)
    // combination including the two asymmetric ones -- a plain sigma==sigma
    // grid never seeds a fast-rise/slow-fall (or the reverse) starting
    // point, so Migrad can converge to a spurious narrow "spike" near
    // closely-spaced low-E points instead of the genuine, broader skew Gaussian
    // (confirmed on x=75.5,y=-39.2: the spike's chi2/ndof was 6.2 while an
    // asymmetric-seeded refit of the same data found chi2/ndof 3.0 at the
    // visually correct peak -- a real, better optimum the old grid simply
    // never gave Migrad a starting point near).
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

    // Data-driven seed: find the point of maximum positive residual above
    // the background-only fit, and estimate the skew Gaussian's width from where
    // the residual falls to half that peak on either side (FWHM). In
    // high-statistics bins the chi2 landscape is sharp enough that Migrad
    // won't wander far from a mediocre starting point, so the generic grid
    // seeds above can converge to a spuriously small D even when a real,
    // well-localized skew Gaussian is clearly present in the residuals -- this seed
    // starts already close to the true optimum.
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
            // Half-max crossings on each side of the peak, estimated
            // independently -- the whole point of this seed is that these
            // two crossings need not be symmetric about the peak (that
            // asymmetry is exactly the shape the fixed-sigma model can't
            // capture; see bg_gauss_model's header comment).
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

    // Track the best fit overall AND the best fit that isn't boundary-
    // pinned, separately. The chi2-global-best restart is sometimes a
    // degenerate fit -- an extremely broad, wall-pinned skew Gaussian that
    // mimics extra background curvature rather than describing a genuine
    // localized skew Gaussian (confirmed by direct inspection: relaxing the wall
    // further just chases it lower with near-zero chi2 improvement,
    // instead of settling at an interior optimum). Picking only the
    // single global-best candidate and discarding everything on a bound
    // hit throws away a perfectly good, non-degenerate, still highly
    // significant local minimum (a well-localized skew Gaussian right at the
    // visible peak) whenever it happens to fit slightly worse than the
    // degenerate one -- so scan every restart and keep the best VALID
    // (non-bound-pinned) candidate as a fallback instead of giving up
    // entirely.
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

        // sigmaL/sigmaR: only the UPPER bound is degenerate (an extremely
        // wide skew Gaussian mimicking background curvature, the original failure
        // mode this check was built for). Hitting the LOWER bound just means
        // a genuinely sharp rise/fall -- a real, common feature now that the
        // two widths are independent, not a sign of a bad fit (same
        // rationale as D's own check_lower=false below).
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

    // The global-best candidate is boundary-pinned (likely degenerate) but
    // a non-pinned alternative exists -- use it instead of giving up.
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

// Per (x,y) bin, fit <scale> vs <response>; bins with <4 valid entries
// (or <3 finite after masking) are dropped entirely.
//
// *** n_procs_hint > 1 (fork-based parallel path below) IS CURRENTLY UNRELIABLE.
// *** DO NOT rely on it for a real calibration run. Both callers
// *** (run_calibration_map_diag.sh, submit_calibration_map_diag_condor.sh)
// *** have been reverted to force single-threaded (n_procs=1) until this is
// *** actually root-caused and re-validated. Two independent runs (one via
// *** condor with 8 shards, one interactive with 2) both silently lost
// *** every shard's contribution except the parent's own (shard 0) --
// *** merged bin counts landed suspiciously exactly at "total/n_procs" both
// *** times (446 with n_procs=8; 1721 with n_procs=2), with every child
// *** reporting a clean exit and, in the n_procs=2 case, the new
// *** read_count==0 sanity check below NOT firing -- meaning the shard file
// *** wasn't simply empty. A synthetic standalone reproduction of this exact
// *** fork/write/waitpid/read sequence (same code shape, but fabricated
// *** FitRow data instead of real Minuit2 fits, no ROOT involved at all)
// *** round-tripped correctly (3478/3478). That isolates the bug to
// *** something about forking specifically AFTER ROOT/Minuit2 has been
// *** initialized and used (fit_log_gauss_function_constrained calls
// *** ROOT::Math::Factory::CreateMinimizer) -- ROOT was never designed with
// *** fork() safety in mind, and this is consistent with some kind of
// *** post-fork corruption in its internal state that doesn't crash but
// *** produces wrong results. Not yet understood well enough to trust.
//
// Each position bin's fit is completely independent of every other (it only
// reads its own slice of `rows`), so this WOULD be an embarrassingly
// parallel workload -- worth exploiting once fixed, now that the
// restart-grid densification (see fit_log_gauss_function_constrained's
// seed-grid comment) made each bin's fit several times more expensive: with
// the seed grid alone, this was pushing full-production runs (3443 bins)
// toward 40 minutes single-threaded.
//
// Parallelized via fork(), not std::thread: an earlier std::thread version
// of this function (one thread per shard, each independently calling
// ROOT::Math::Factory::CreateMinimizer) segfaulted in every worker within
// ~10 seconds of starting, every time -- consistent with a race in ROOT
// 5.34/Minuit2's lazy plugin-loading machinery (TPluginManager and friends),
// which predates ROOT's later official thread-safety support and isn't
// documented as safe to call concurrently from multiple threads in one
// process. Forking separate OS processes was meant to sidestep that (each
// child gets its own independent copy of the process image, so every child
// does its own lazy initialization single-threaded, in total isolation --
// no shared mutable memory left to race on) -- and it does avoid the
// segfault, but has instead exposed the silent-corruption bug above.
struct GroupEntry
{
    double x_center, y_center;
    std::vector<int> idxs;
};

// Groups `rows` by (x_center, y_center) -- same deterministic std::map
// ordering every caller gets, which is what makes "entries[k] with
// k % n_shards == shard_index" a stable, reproducible partition across
// completely independent process invocations (condor shards included).
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

// Fits every position-bin group at flat index k where k % n_shards ==
// shard_index (see build_position_groups). Used two ways: (1) directly from
// main()'s --shard mode, one fresh OS process per condor job -- the
// currently-recommended parallel path, since each shard starts clean with
// no shared post-ROOT-init state to corrupt (see compute_all_regions'
// header comment for why in-process fork() is NOT safe here); (2) by
// compute_all_regions' own (disabled) fork()-based fit_shard, so that
// dead-but-documented path and this one share one source of truth instead
// of two copies of the same fit logic silently drifting apart.
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

std::vector<FitRow> compute_all_regions(const std::vector<BinRow> &rows, double chi2_improvement_min, unsigned n_procs_hint = 0)
{
    std::vector<GroupEntry> entries = build_position_groups(rows);

    // Fits every entries[k] with k % n_procs == shard, via compute_shard
    // (see its header comment). Called post-fork in each child (and
    // directly for shard 0 in the parent), so this always runs
    // single-threaded within whichever process calls it.
    auto fit_shard = [&](unsigned shard, unsigned n_procs) -> std::vector<FitRow>
    {
        return compute_shard(rows, chi2_improvement_min, shard, n_procs);
    };

    // n_procs_hint (from main()'s argv[4]) takes priority over auto-detection:
    // on a shared/multi-tenant machine, or a condor execute node whose
    // granted core count isn't reliably reflected by
    // hardware_concurrency(), auto-detecting the FULL local core count and
    // forking that many processes can over-fork relative to what's actually
    // available -- explicit is safer whenever the caller knows the real
    // number (e.g. a condor job passing exactly its own request_cpus).
    unsigned n_procs = n_procs_hint;
    if (n_procs == 0)
    {
        n_procs = std::thread::hardware_concurrency();
        if (n_procs == 0)
            n_procs = 1;
    }
    n_procs = std::min<unsigned>(n_procs, (unsigned)std::max<size_t>(entries.size(), 1));

    if (n_procs == 1)
        return fit_shard(0, 1);

    // One shard file per child (shard 0 is handled by the parent directly,
    // no file needed), written into the CURRENT WORKING DIRECTORY rather
    // than /tmp. A first version of this used mkstemp() under /tmp, and on
    // a real condor execute node (not this machine -- a local repro of the
    // exact same mkstemp/fork/ofstream/ifstream sequence round-tripped fine
    // here) it silently produced EMPTY shard files: every child ran and
    // burned real CPU (confirmed via the job's own resource usage log), and
    // every child exited 0, yet the parent merged in zero rows from any of
    // them -- only shard 0 (computed directly in the parent) survived,
    // silently discarding ~87% of all position bins' fits with no error of
    // any kind. That corrupted a saved calibration lookup table before
    // anyone noticed. Root cause not fully pinned down (the leading
    // suspect is /tmp on that specific execute node behaving differently
    // under the singularity bind-mount than a plain login-node /tmp does),
    // but rather than keep guessing at a remote filesystem's behavior, this
    // writes into the same GPFS directory every other read/write in this
    // program already uses successfully (the input ROOT file, the output
    // JSON/CSV) -- a path already proven reliable in this exact job,
    // instead of a new, unverified one. Never let this fail silently again:
    // every open is checked, and a shard whose child had entries to fit but
    // produced zero readable rows back is now a hard error, not something
    // that just quietly ships a wrong answer.
    char cwd_buf[4096];
    std::string cwd = getcwd(cwd_buf, sizeof(cwd_buf)) ? cwd_buf : ".";
    pid_t self_pid = getpid();

    std::vector<std::string> shard_files(n_procs);
    for (unsigned s = 1; s < n_procs; s++)
    {
        std::ostringstream path;
        path << cwd << "/.calibmap_shard_" << self_pid << "_" << s << ".tmp";
        shard_files[s] = path.str();
    }

    std::vector<pid_t> pids;
    for (unsigned s = 1; s < n_procs; s++)
    {
        pid_t pid = fork();
        if (pid < 0)
        {
            std::cerr << "ERROR: fork failed for shard " << s << std::endl;
            exit(1);
        }
        if (pid == 0)
        {
            // Child: fit this shard, serialize the results, then exit
            // immediately -- must not fall through and return into the rest
            // of main() a second time.
            std::vector<FitRow> local = fit_shard(s, n_procs);
            std::ofstream out(shard_files[s].c_str());
            if (!out.is_open())
            {
                std::cerr << "ERROR: shard " << s << " could not open " << shard_files[s]
                          << " for writing" << std::endl;
                _exit(1);
            }
            for (const auto &r : local)
            {
                write_double(out, r.x_center); out << " ";
                write_double(out, r.y_center); out << " ";
                write_double(out, r.A); out << " ";
                write_double(out, r.B); out << " ";
                write_double(out, r.C); out << " ";
                write_double(out, r.D); out << " ";
                write_double(out, r.E0); out << " ";
                write_double(out, r.sigmaL); out << " ";
                write_double(out, r.sigmaR); out << " ";
                write_double(out, r.chi2); out << " ";
                out << r.n_params << " " << r.ndof << "\n";
            }
            out.close();
            if (out.fail())
            {
                std::cerr << "ERROR: shard " << s << " failed writing " << shard_files[s] << std::endl;
                _exit(1);
            }
            _exit(0);
        }
        pids.push_back(pid);
    }

    // Parent takes shard 0 itself while the children run.
    std::vector<FitRow> fit_rows = fit_shard(0, n_procs);

    for (pid_t pid : pids)
    {
        int status = 0;
        waitpid(pid, &status, 0);
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        {
            std::cerr << "ERROR: a fit worker process (pid " << pid << ") failed" << std::endl;
            exit(1);
        }
    }

    for (unsigned s = 1; s < n_procs; s++)
    {
        // Cheap re-derivation (no fitting) of how many of this shard's
        // candidate bins existed, purely to sanity-check the read-back
        // below -- not every candidate produces a row (some have <4 valid
        // points), so this is an upper bound, not an exact expectation.
        size_t raw_assigned = 0;
        for (size_t k = s; k < entries.size(); k += n_procs)
            raw_assigned++;

        std::ifstream in(shard_files[s].c_str());
        if (!in.is_open())
        {
            std::cerr << "ERROR: could not read back shard " << s << " file " << shard_files[s] << std::endl;
            exit(1);
        }

        size_t read_count = 0;
        FitRow r;
        while (in >> r.x_center >> r.y_center >> r.A >> r.B >> r.C >> r.D >> r.E0 >>
               r.sigmaL >> r.sigmaR >> r.chi2 >> r.n_params >> r.ndof)
        {
            fit_rows.push_back(r);
            read_count++;
        }
        in.close();
        std::remove(shard_files[s].c_str());

        // The exact failure mode this whole rewrite exists to catch: a
        // shard with real candidate bins to fit that came back with
        // nothing readable at all.
        if (read_count == 0 && raw_assigned > 0)
        {
            std::cerr << "ERROR: shard " << s << " had " << raw_assigned
                      << " candidate bins but 0 were read back -- refusing to silently"
                      << " produce a corrupted calibration (see compute_all_regions'"
                      << " header comment)" << std::endl;
            exit(1);
        }
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
                    // Matches production CalibrationMap.cpp exactly: every
                    // finite, non-outlier cell still gets blended toward its
                    // neighborhood average. This is deliberate, not the
                    // dilution bug it was mistaken for earlier -- the
                    // downstream consumer linearly interpolates between grid
                    // points, so the saved surface needs to vary smoothly;
                    // full per-bin independence (no blending) would leave
                    // physically-unmotivated bin-to-bin discontinuities. A
                    // direct fit that disagrees sharply with every neighbor
                    // is more likely a statistical fluctuation in that one
                    // bin than a genuine sub-grid-spacing feature.
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

    // Pull D back toward 0 with distance from the nearest cell that actually
    // had a direct fit (real skew Gaussian OR a real background-only D=0
    // fit -- both count, since coarse_maps["D"] is finite=0.0 for
    // background-only bins, not NaN). Cells with a direct fit are left
    // exactly as smooth_map computed them (distance 0, decay factor 1); only
    // the pure-diffusion fill for cells with no fit at all is tapered. See
    // D_DECAY_LENGTH_BINS' comment for why this is needed.
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
            // The calibration must not be defined past the true
            // buffer-restricted fiducial acceptance -- see ACTIVE_FIDUCIAL_RECT's
            // comment. Checked first so nothing smooth_map extrapolated out
            // here (for ANY parameter, not just D) ever reaches the saved
            // table, regardless of how confident-looking that extrapolation is.
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

// Full-precision FitRow serialization for the condor shard/merge path (see
// main()'s --shard/--merge modes): one shard's worth of independent
// position-bin fits, written by a condor job that exits immediately after
// (no smoothing/save), later concatenated across every shard by --merge
// before running the existing smoothing/decay/clamp/fiducial-mask/save
// pipeline exactly as the single-process path does.
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

// Side CSV of per-bin chi2/ndof (final model kept -- background+skew Gaussian
// where significant, background-only otherwise), for chi2/dof heatmap
// plots. Not consumed by anything downstream of the JSON lookup.
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

// Side CSV of the raw per-position DIRECT fit (pre-smoothing/pre-decay/
// pre-clamp A,B,C,D,E0,sigmaL,sigmaR -- exactly what build_coarse_maps
// reads into coarse_maps before smooth_map() runs), for diagnostic plots
// that compare "direct fit" vs "final (smoothed)" per position. Not
// consumed by anything downstream of the JSON lookup.
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

    // argv[3] is ignored (kept only so the argv positions below stay the
    // same as the production tool's): this build is always EM-only, so the
    // skew-Gaussian position/width domain and fiducial mask are the
    // EM-only ones, matching JetEnergyScaleFineGrid_em_only's input.
    const bool em_only = true;
    E0_MAX = em_only ? E0_MAX_em_only : E0_MAX_ecal_hcal;
    SIGMA_MAX = em_only ? SIGMA_MAX_em_only : SIGMA_MAX_ecal_hcal;
    ACTIVE_FIDUCIAL_RECT = em_only ? &kFiducialRectEcal : &kFiducialRect;
    ACTIVE_FIDUCIAL_BUFFER = em_only ? reco_fiducial_buffer_ecal : reco_fiducial_buffer;

    // argv[5]: optional mode for true multi-PROCESS parallelism via separate
    // condor job invocations (see compute_shard's header comment for why
    // this is safe where in-process fork() is not). Omitted/"" = existing
    // single-process behavior below, unchanged.
    //
    //   --shard <shard_index> <n_shards> <shard_output_csv>
    //       Fits only this shard's position bins (see build_position_groups)
    //       and writes them to shard_output_csv, then exits immediately --
    //       no smoothing/decay/clamp/fiducial-mask/save. One condor job per
    //       shard, run with identical argv[1-3] (same input file, same
    //       em_only flag) across all n_shards jobs.
    //   --merge <shard_list_file> <output_file>
    //       shard_list_file is a plain text file, one shard CSV path per
    //       line (as written by every --shard job above). Concatenates all
    //       of them into one fit_rows vector, then runs the exact same
    //       smoothing/decay/clamp/fiducial-mask/save pipeline the
    //       single-process path uses -- this is the only place that
    //       pipeline runs, so results are identical in shape/logic to a
    //       single-process run over the same input, just assembled from
    //       shards instead of one big in-process fit. Fast (no fitting): run
    //       this locally right after condor_wait confirms every shard job
    //       finished, no need to route it through condor itself.
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

    // argv[4]: optional explicit worker-process count for compute_all_regions'
    // fork()-based parallel fit (see its header comment for why fork, not
    // std::thread). 0/omit = auto-detect via hardware_concurrency(), which is
    // right for an interactive box you have to yourself but WRONG on a
    // shared login node (would fork one process per physical core on the
    // whole machine, not just your fair share) or on a condor execute node
    // whose granted core count isn't guaranteed to be reflected by
    // hardware_concurrency(). Pass this explicitly -- matching whatever was
    // actually requested/allotted -- in either of those cases.
    //
    // HARD-DISABLED (>1 forced back to 1) until compute_all_regions' known
    // silent-corruption bug is actually root-caused: two separate runs
    // (n_procs=8 via condor, n_procs=2 interactive) each silently dropped
    // every shard's results except the parent's own, both times exiting 0
    // with no error. Do not remove this clamp without first fixing that --
    // see compute_all_regions' header comment for what's already been ruled
    // out. Auto-detect (argv[4] omitted -> hardware_concurrency()) is
    // clamped too: it hits the exact same bug, not just the login-node
    // fairness issue the comment above describes.
    unsigned n_procs_hint = (argc > 4) ? (unsigned)std::atoi(argv[4]) : 0;
    if (n_procs_hint != 1)
    {
        if (n_procs_hint > 1)
            std::cout << "WARNING: n_procs=" << n_procs_hint
                       << " requested but the fork-based parallel fit is currently"
                       << " known to silently corrupt results -- forcing single-threaded (1)."
                       << " See compute_all_regions' header comment." << std::endl;
        n_procs_hint = 1;
    }

    std::cout << "Loading " << input_file << " ..." << std::endl;
    std::vector<BinRow> rows = load_bins(input_file);
    std::cout << "Loaded " << rows.size() << " (E,x,y) bin entries" << std::endl;

    // 4 dof: the skew Gaussian adds D, E0, sigmaL, sigmaR over the 3-param background.
    double chi2_improvement_min = TMath::ChisquareQuantile(SKEW_GAUSSIAN_SIGNIFICANCE_CL, 4);

    std::cout << "Fitting position bins ..." << std::endl;
    std::vector<FitRow> fit_rows = compute_all_regions(rows, chi2_improvement_min, n_procs_hint);
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
