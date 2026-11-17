// Step 2 of the EM-only JES calibration: builds the JES lookup JSON from the
// "bins" tree of JetEnergyScaleFineGrid_em_only (Gaussian core or kernel peak).
// The scale vs response at each (x,y) position bin is described by
//   scale(E) = A + B lnE + C ln^2 E + D ln^3 E
// with A and B per position and C and D shared by the whole face: the shape
// in ln E is universal and only the level and slope vary across the ECal.
// C and D come from one least-squares fit over every position's points (A and
// B free per position, see fit_global_CD); A and B are then the per-position
// fits with C and D held there, which is the joint optimum since the model is
// linear. The A and B maps are smoothed over the grid (as the scale at two
// anchor energies, see ANCHOR_E) and masked to the buffered ECal fiducial
// region. The JSON records the range of all fitted responses [Emin, Emax].
// Since 2026-11-06 the evaluators (ApplyCorrections_em_only, jes_closure) do
// not use it: they hold E at 2 GeV below (JetParameters.h kJesEvalMinGeV, near
// the curve's minimum) and follow the curve above Emax up to sqrt(s)/2: on
// jesmap_pixring_v1 the log-cubic levels off past Emax = 96.8 GeV (median
// scale 0.907 there, 0.897 at 250 GeV), so holding it at Emax (as CMSSW's
// SimpleJetCorrector clamps pT to each record's range) was no longer needed.
// A map whose curve does not level off would need the Emax hold back. The
// range is the whole face's, not each position's: the shape is constrained
// over all of it, and past a position's own last point only its level and
// slope are carried along that shape (a per-position clamp held the outer
// positions, with no reliable points above about 40 GeV, flat from there,
// overcorrecting their high-E jets by up to 5%).
//
// A log-polynomial form is used, following ATLAS Collaboration, "Jet energy
// measurement with the ATLAS detector in proton-proton collisions at sqrt(s) =
// 7 TeV", Eur. Phys. J. C 73 (2013) 2304, arXiv:1112.6426, Eq. 9. The order is
// the lowest that describes the shape: on jesmap_pixring_v1 removing
// the cubic term worsens the global chi2 by 956 and shifts the scale by about
// 1% at the ends of the energy range (D is small only because ln^3 E is large:
// the D term is -0.295 at 100 GeV).
//
// The JSON carries "model": "logcubic", where it came from ("sample", "k", "tune",
// "commit": the pipeline commit it was built with), the global "C", "D", "Emin"
// and "Emax", and bins {i, j, A, B} (read by jet_calibration/ApplyCorrections_em_only.cpp
// and jet_scale/jes_closure.cpp). It also writes <out>.chi2.csv (chi2/ndof of
// the direct fits), <out>.chi2_smoothed.csv (the same points against the
// smoothed map) and <out>.direct_fits.csv (the unsmoothed per-position fits).
//
// Full (ECAL+HCAL) jets (2026-11-16): argv[7] "full" masks to the full-jet
// fiducial region (HCAL rectangle clipped to the ECAL, 20 cm buffer;
// JetParameters.h kFiducialRectFull) instead of the ECAL one; the model is
// the same log-cubic (Sean: as EM-only, no Gaussian term). The JSON records
// "jet_type".
//
// Which parameters can be face-wide: <out>.uniformity.txt compares nested
// models on the same points, by how many of A, B, C, D are free per position
// (the rest face-wide): A only; A, B (the map's choice); A, B, C; all four.
// Each is an exact linear least-squares fit (no start values). Per model:
// chi2, ndof, and the F-test of each step against the next (the extra
// per-position parameters' chi2 drop per extra parameter over the larger
// model's chi2/ndof, which absorbs the points' excess scatter).
// <out>.free_fits.csv: the all-free per-position fits (positions with at
// least 6 points), to map how C and D vary.
//
// Free parameters per position (2026-11-17, Sean: C, maybe D, may not be
// constant over the face; the map must be testable both ways): argv[8] = 2
// (default; A, B per position, C, D face-wide, the map above), 3 (A, B, C per
// position, D face-wide) or 4 (all four per position). The face-wide ones
// come from the nested fit (fit_nested, as the uniformity test); a position
// needs at least free + 2 points to be fit (others are filled by the
// smoothing). The free part is smoothed as the scale at `free` anchor
// energies (ANCHOR_E_FREE), each anchor map separately, and solved back per
// cell, as A and B always were. The JSON then carries per-bin "C" (3, 4) and
// "D" (4) next to the face-wide values, and "free_params"; the evaluators
// (ApplyCorrections_em_only, jes_closure) interpolate per-bin C and D
// bilinearly when present and use the face-wide ones otherwise.
//
// Usage: ./CalibrationMap_em_only <bins.root> <out.json> [sample] [k] [tune] [commit] [em|full] [2|3|4]
//        (filled in by submit_jet_scale_em_only.sh; "unknown" if left out)
// Build: g++ -m64 -O2 CalibrationMap_em_only.cpp -o CalibrationMap_em_only `root-config --cflags --libs` -std=c++11

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdio>
#include <cstdlib>


// CONFIGURATION

std::string ROOT_FILE = "jet_calibration_grid_test.root";
std::string OUTPUT_JSON = "JetEnergyScale_lookup_GridTest.json";

// Provenance written to the JSON, from the command line
const char *const MAP_INFO_KEY[4] = {"sample", "k", "tune", "commit"};
std::string map_info[4] = {"unknown", "unknown", "unknown", "unknown"};

// Sanity cuts dropping the noise tail of low-statistics bins
const double RESPONSE_SANITY_FACTOR = 1.5;
const double SCALE_SANITY_MAX = 1.5;

// Positions with fewer valid energy bins are not fit (they are filled by the smoothing)
const int MIN_POINTS = 4;

// Energy bins with a smaller Kish effective count are not used: an estimate
// from one or two effective jets is biased (2026-10 no-UE map: the 417
// bins below 3 sit 0.5-0.9% off the map, about 4 sigma; above 3 the pulls are
// about 1 at every n_eff). Bins down to n_eff 3 need the MAD floor of
// JetEnergyScaleFineGrid_em_only (27 Oct 2026); with the bare weighted MAD
// their errors collapse and n_eff >= 10 was needed instead.
const double MIN_NEFF = 3.0;

// A and B are strongly correlated, so they are smoothed as the scale at these
// two energies (inside the fitted range; the map to them is exact and
// invertible) rather than coefficient by coefficient, which could pair one
// position's A with its neighbours' B and put spots in the scale
const double ANCHOR_E[2] = {3.0, 30.0};
// With C (and D) free per position, as many anchors as free parameters,
// spread in ln E over the range most positions have points in (the outer
// ones have few reliable points above ~40 GeV)
const double ANCHOR_E3[3] = {3.0, 10.0, 35.0};
const double ANCHOR_E4[4] = {3.0, 7.0, 18.0, 45.0};
int g_free = 2; // argv[8]

// Copied from shared/JetParameters.h (keep in sync; that header needs
// fastjet). Upstream, JetEnergyScaleFineGrid only fills grid points inside
// this region, so cells outside it have no statistics; the saved map is
// masked to it so smoothing never extrapolates past the acceptance.
const float cut_x_inner_ecal = 20.53604584f, cut_x_outer_ecal = 137.59500000f;
const float cut_y_min_ecal = -96.72439178f, cut_y_max_ecal = 86.25565167f;
const float reco_fiducial_buffer_ecal = 20.0f;
// Full jets (kFiducialRectFull: the HCAL rectangle clipped to the ECAL one)
const float cut_x_inner_full = 24.15591212f, cut_x_outer_full = 131.82035143f;
const float cut_y_min_full = -83.19687629f, cut_y_max_full = 86.25565167f;
const float reco_fiducial_buffer_full = 20.0f;
bool g_full_jets = false; // argv[7] "full"

inline bool pass_fiducial_cut(float jetXE, float jetYE)
{
    const float buffer = g_full_jets ? reco_fiducial_buffer_full : reco_fiducial_buffer_ecal;
    const float xi = g_full_jets ? cut_x_inner_full : cut_x_inner_ecal, xo = g_full_jets ? cut_x_outer_full : cut_x_outer_ecal;
    const float y0 = g_full_jets ? cut_y_min_full : cut_y_min_ecal, y1 = g_full_jets ? cut_y_max_full : cut_y_max_ecal;
    const float ax = std::fabs(jetXE);
    const bool in_x = ax >= (xi + buffer) && ax <= (xo - buffer);
    const bool in_y = jetYE >= (y0 + buffer) && jetYE <= (y1 - buffer);
    return in_x && in_y;
}

// Map smoothing parameters
const int SMOOTH_N_ITER = 6;
const int SMOOTH_RADIUS = 2;
// -1 = unfit bins take the diffusion average (cells far outside the data stay NaN)
const int SMOOTH_MAX_FILL_DISTANCE = -1;


struct BinRow
{
    double x_center, y_center;
    double E_low, E_high;
    double scale_mu, response_mu;
    double scale_mu_err, response_mu_err;
    double n_eff;
};

// The points of one position: ln(response), scale, weight 1/scale_err^2
struct Cell
{
    double x_center, y_center;
    std::vector<double> l, y, w;
    double Emin, Emax; // range of the fitted responses (diagnostics; the clamp is the face's)
};

struct FitRow
{
    double x_center, y_center;
    double A, B, C, D; // C, D: the face-wide values unless free per position
    double Emin, Emax;
    double chi2;
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
    double scale_mu, response_mu, scale_mu_err, response_mu_err, n_eff;

    tree->SetBranchAddress("x_center", &x_center);
    tree->SetBranchAddress("y_center", &y_center);
    tree->SetBranchAddress("E_low", &E_low);
    tree->SetBranchAddress("E_high", &E_high);
    tree->SetBranchAddress("scale_mu", &scale_mu);
    tree->SetBranchAddress("response_mu", &response_mu);
    tree->SetBranchAddress("scale_mu_err", &scale_mu_err);
    tree->SetBranchAddress("response_mu_err", &response_mu_err);
    tree->SetBranchAddress("n_eff", &n_eff);

    std::vector<BinRow> rows;
    rows.reserve(tree->GetEntries());

    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; i++)
    {
        tree->GetEntry(i);
        rows.push_back({x_center, y_center, E_low, E_high,
                         scale_mu, response_mu, scale_mu_err, response_mu_err, n_eff});
    }

    f.Close();
    return rows;
}

// The (E,x,y) bins a position fit uses
bool valid_bin(const BinRow &r)
{
    return std::isfinite(r.scale_mu) && std::isfinite(r.response_mu) &&
           std::isfinite(r.scale_mu_err) && std::isfinite(r.response_mu_err) &&
           r.scale_mu > 0 && r.response_mu > 0 && r.scale_mu_err > 0 && r.n_eff >= MIN_NEFF &&
           r.scale_mu <= SCALE_SANITY_MAX &&
           r.response_mu <= RESPONSE_SANITY_FACTOR * r.E_high;
}

// Groups the valid bins by (x_center, y_center); positions with fewer than
// MIN_POINTS of them are dropped
std::vector<Cell> build_cells(const std::vector<BinRow> &rows)
{
    std::map<std::pair<double, double>, Cell> groups;
    for (const BinRow &r : rows)
    {
        if (!valid_bin(r))
            continue;
        Cell &c = groups[{r.x_center, r.y_center}];
        if (c.l.empty())
        {
            c.x_center = r.x_center;
            c.y_center = r.y_center;
            c.Emin = c.Emax = r.response_mu;
        }
        c.l.push_back(std::log(r.response_mu));
        c.y.push_back(r.scale_mu);
        c.w.push_back(1.0 / (r.scale_mu_err * r.scale_mu_err));
        c.Emin = std::min(c.Emin, r.response_mu);
        c.Emax = std::max(c.Emax, r.response_mu);
    }

    std::vector<Cell> cells;
    for (const auto &kv : groups)
        if ((int)kv.second.l.size() >= MIN_POINTS)
            cells.push_back(kv.second);
    return cells;
}

// Weighted sums G[a][b] = sum w col_a col_b over a cell's points, with columns
// 0 = 1, 1 = l, 2 = l^2, 3 = l^3, 4 = scale
void cell_sums(const Cell &c, double G[5][5])
{
    for (int a = 0; a < 5; a++)
        for (int b = 0; b < 5; b++)
            G[a][b] = 0.0;
    for (size_t i = 0; i < c.l.size(); i++)
    {
        const double l = c.l[i];
        const double col[5] = {1.0, l, l * l, l * l * l, c.y[i]};
        for (int a = 0; a < 5; a++)
            for (int b = 0; b < 5; b++)
                G[a][b] += c.w[i] * col[a] * col[b];
    }
}

// The face-wide C and D: minimizes the sum over positions of
// chi2 = sum ((scale - A - B l - C l^2 - D l^3) / scale_err)^2 with A and B
// free per position. A and B are profiled out exactly: with W the weights and
// X = [1, l] in a position, the residual projector
// M = W - W X (X^T W X)^-1 X^T W gives
//   (C, D) = (sum Z^T M Z)^-1 sum Z^T M y,  Z = [l^2, l^3].
void fit_global_CD(const std::vector<Cell> &cells, double &C, double &D)
{
    double S[2][2] = {{0, 0}, {0, 0}}, t[2] = {0, 0}, yMy = 0.0;
    long n_points = 0;
    for (const Cell &c : cells)
    {
        double G[5][5];
        cell_sums(c, G);
        const double det = G[0][0] * G[1][1] - G[0][1] * G[0][1];
        // u^T M v = G[u][v] - [G[u][0], G[u][1]] (X^T W X)^-1 [G[0][v], G[1][v]]^T
        auto uMv = [&](int u, int v)
        {
            const double a0 = G[u][0], a1 = G[u][1], b0 = G[0][v], b1 = G[1][v];
            return G[u][v] - (a0 * (G[1][1] * b0 - G[0][1] * b1) +
                              a1 * (G[0][0] * b1 - G[0][1] * b0)) / det;
        };
        for (int a = 0; a < 2; a++)
        {
            for (int b = 0; b < 2; b++)
                S[a][b] += uMv(2 + a, 2 + b);
            t[a] += uMv(2 + a, 4);
        }
        yMy += uMv(4, 4);
        n_points += (long)c.l.size();
    }

    const double det = S[0][0] * S[1][1] - S[0][1] * S[1][0];
    if (cells.empty() || !(det > 0))
    {
        std::cerr << "ERROR: global C, D fit is singular" << std::endl;
        exit(1);
    }
    const double cov[2][2] = {{S[1][1] / det, -S[0][1] / det}, {-S[1][0] / det, S[0][0] / det}};
    C = cov[0][0] * t[0] + cov[0][1] * t[1];
    D = cov[1][0] * t[0] + cov[1][1] * t[1];
    const double chi2 = yMy - (C * t[0] + D * t[1]);
    const long ndof = n_points - 2L * (long)cells.size() - 2;
    // errors scaled by sqrt(chi2/ndof): the points scatter beyond their errors
    const double inflate = ndof > 0 ? std::sqrt(std::max(1.0, chi2 / ndof)) : 1.0;
    std::printf("Global C, D fit over %zu positions, %ld points: chi2/ndof = %.1f/%ld\n"
                "  C = %.6f +- %.6f, D = %.6f +- %.6f (errors x %.2f), corr = %.3f\n",
                cells.size(), n_points, chi2, ndof,
                C, inflate * std::sqrt(cov[0][0]), D, inflate * std::sqrt(cov[1][1]),
                inflate, cov[0][1] / std::sqrt(cov[0][0] * cov[1][1]));
}

// A and B per position with C and D fixed: a weighted linear fit of
// scale - C l^2 - D l^3 on [1, l]
std::vector<FitRow> fit_cells(const std::vector<Cell> &cells, double C, double D)
{
    std::vector<FitRow> fits;
    for (const Cell &c : cells)
    {
        double s0 = 0, s1 = 0, s2 = 0, r0 = 0, r1 = 0;
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i], r = c.y[i] - l * l * (C + l * D);
            s0 += c.w[i];
            s1 += c.w[i] * l;
            s2 += c.w[i] * l * l;
            r0 += c.w[i] * r;
            r1 += c.w[i] * l * r;
        }
        const double det = s0 * s2 - s1 * s1;
        FitRow f;
        f.x_center = c.x_center;
        f.y_center = c.y_center;
        f.A = (s2 * r0 - s1 * r1) / det;
        f.B = (s0 * r1 - s1 * r0) / det;
        f.C = C;
        f.D = D;
        f.Emin = c.Emin;
        f.Emax = c.Emax;
        f.chi2 = 0.0;
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i];
            const double d = c.y[i] - (f.A + l * (f.B + l * (C + l * D)));
            f.chi2 += c.w[i] * d * d;
        }
        f.ndof = (int)c.l.size() - 2;
        fits.push_back(f);
    }
    return fits;
}

// Solves the n x n system M x = v (n <= 4) by Gaussian elimination with
// partial pivoting; false if singular
bool solve_small(int n, double M[4][4], const double v[4], double x[4])
{
    double a[4][5];
    for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) a[i][j] = M[i][j]; a[i][n] = v[i]; }
    for (int c = 0; c < n; c++)
    {
        int piv = c;
        for (int r = c + 1; r < n; r++) if (std::fabs(a[r][c]) > std::fabs(a[piv][c])) piv = r;
        if (!(std::fabs(a[piv][c]) > 1e-300)) return false;
        for (int j = 0; j <= n; j++) std::swap(a[c][j], a[piv][j]);
        for (int r = 0; r < n; r++)
        {
            if (r == c) continue;
            const double f = a[r][c] / a[c][c];
            for (int j = c; j <= n; j++) a[r][j] -= f * a[c][j];
        }
    }
    for (int i = 0; i < n; i++) x[i] = a[i][n] / a[i][i];
    return true;
}

// Nested model with the first k of (A, B, C, D) free per position and the
// rest face-wide: the face-wide ones from the profiled normal equations
// (per position, the residual projector of its free columns), then the free
// ones per position. Returns chi2 and ndof; glob gets the face-wide values.
bool fit_nested(const std::vector<Cell> &cells, int k, double &chi2, long &ndof, double glob[4])
{
    const int g = 4 - k;
    double S[4][4] = {{0}}, t[4] = {0};
    long n_points = 0;
    std::vector<std::vector<double>> Gs;
    for (const Cell &c : cells)
    {
        double G[5][5];
        cell_sums(c, G);
        n_points += (long)c.l.size();
        // u^T M v with M the projector off the first k columns
        double Gkk[4][4], inv[4][4];
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) Gkk[a][b] = G[a][b];
        for (int col = 0; col < k; col++)
        {
            double e[4] = {0, 0, 0, 0}, x[4];
            e[col] = 1;
            double tmp[4][4];
            for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) tmp[a][b] = Gkk[a][b];
            if (!solve_small(k, tmp, e, x)) return false;
            for (int a = 0; a < k; a++) inv[a][col] = x[a];
        }
        auto uMv = [&](int u, int v)
        {
            double r = G[u][v];
            for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) r -= G[u][a] * inv[a][b] * G[b][v];
            return r;
        };
        for (int a = 0; a < g; a++)
        {
            for (int b = 0; b < g; b++) S[a][b] += uMv(k + a, k + b);
            t[a] += uMv(k + a, 4);
        }
    }
    double gv[4] = {0, 0, 0, 0};
    if (g > 0 && !solve_small(g, S, t, gv)) return false;
    for (int a = 0; a < 4; a++) glob[a] = a < k ? NAN : gv[a - k];
    // per position: the free ones with the face-wide ones held, then chi2
    chi2 = 0;
    for (const Cell &c : cells)
    {
        double M[4][4] = {{0}}, v[4] = {0}, x[4] = {0};
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i], col[4] = {1.0, l, l * l, l * l * l};
            double r = c.y[i];
            for (int a = k; a < 4; a++) r -= gv[a - k] * col[a];
            for (int a = 0; a < k; a++) { v[a] += c.w[i] * col[a] * r; for (int b = 0; b < k; b++) M[a][b] += c.w[i] * col[a] * col[b]; }
        }
        if (k > 0 && !solve_small(k, M, v, x)) return false;
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i], col[4] = {1.0, l, l * l, l * l * l};
            double m = 0;
            for (int a = 0; a < 4; a++) m += (a < k ? x[a] : gv[a - k]) * col[a];
            chi2 += c.w[i] * (c.y[i] - m) * (c.y[i] - m);
        }
    }
    ndof = n_points - (long)k * (long)cells.size() - g;
    return true;
}

// The map's fits with the first k of (A, B, C, D) free per position (k = 3, 4)
// and the rest at the face-wide values glob (from fit_nested on the same
// cells): per position a weighted linear fit of scale minus the face-wide
// part on the free columns. Positions with fewer than k + 2 points are
// skipped (filled by the smoothing).
std::vector<FitRow> fit_cells_free(const std::vector<Cell> &cells, int k, const double glob[4])
{
    std::vector<FitRow> fits;
    for (const Cell &c : cells)
    {
        if ((int)c.l.size() < k + 2) continue;
        double M[4][4] = {{0}}, v[4] = {0}, x[4] = {0};
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i], col[4] = {1.0, l, l * l, l * l * l};
            double r = c.y[i];
            for (int a = k; a < 4; a++) r -= glob[a] * col[a];
            for (int a = 0; a < k; a++) { v[a] += c.w[i] * col[a] * r; for (int b = 0; b < k; b++) M[a][b] += c.w[i] * col[a] * col[b]; }
        }
        if (!solve_small(k, M, v, x)) continue;
        double p[4];
        for (int a = 0; a < 4; a++) p[a] = a < k ? x[a] : glob[a];
        FitRow f;
        f.x_center = c.x_center;
        f.y_center = c.y_center;
        f.A = p[0]; f.B = p[1]; f.C = p[2]; f.D = p[3];
        f.Emin = c.Emin;
        f.Emax = c.Emax;
        f.chi2 = 0.0;
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i];
            const double d = c.y[i] - (p[0] + l * (p[1] + l * (p[2] + l * p[3])));
            f.chi2 += c.w[i] * d * d;
        }
        f.ndof = (int)c.l.size() - k;
        fits.push_back(f);
    }
    return fits;
}

// The uniformity test (see the header): nested models on the positions with
// at least 6 points (so all four can be free everywhere), and the all-free fits
void write_uniformity(const std::vector<Cell> &all_cells, const std::string &txt, const std::string &csv)
{
    std::vector<Cell> cells;
    for (const Cell &c : all_cells) if ((int)c.l.size() >= 6) cells.push_back(c);
    std::ofstream out(txt);
    const char *nm[5] = {"", "A per position; B, C, D face-wide", "A, B per position; C, D face-wide (the map)",
                         "A, B, C per position; D face-wide", "A, B, C, D per position"};
    double chi[5]; long nd[5]; bool ok[5];
    out << "Positions with >= 6 points: " << cells.size() << " (of " << all_cells.size() << " fitted)\n";
    for (int k = 1; k <= 4; k++)
    {
        double glob[4];
        ok[k] = fit_nested(cells, k, chi[k], nd[k], glob);
        char line[400];
        if (!ok[k]) { std::snprintf(line, sizeof line, "%-48s singular\n", nm[k]); out << line; continue; }
        std::snprintf(line, sizeof line, "%-48s chi2/ndof = %.1f/%ld = %.3f", nm[k], chi[k], nd[k], chi[k] / nd[k]);
        out << line;
        for (int a = k; a < 4; a++) { std::snprintf(line, sizeof line, "  %c = %.6g", "ABCD"[a], glob[a]); out << line; }
        out << "\n";
    }
    out << "\nF-tests (each step's chi2 drop per extra parameter / the larger model's chi2/ndof):\n";
    for (int k = 1; k < 4; k++)
    {
        if (!ok[k] || !ok[k + 1]) continue;
        const long dn = nd[k] - nd[k + 1];
        const double F = ((chi[k] - chi[k + 1]) / dn) / (chi[k + 1] / nd[k + 1]);
        char line[300];
        std::snprintf(line, sizeof line, "  free %c per position: chi2 drop %.1f for %ld parameters, F = %.2f\n", "ABCD"[k], chi[k] - chi[k + 1], dn, F);
        out << line;
    }
    out << "(F near 1: the extra per-position freedom fits noise; well above 1: the parameter varies over the face)\n";
    // what it does to the scale: per position, the RMS over its points of
    // (all-free - shared-C,D) / shared, against the points' median relative error
    {
        double glob[4], ch; long nd2;
        if (ok[2] && fit_nested(cells, 2, ch, nd2, glob))
        {
            std::vector<FitRow> shared = fit_cells(cells, glob[2], glob[3]);
            std::vector<double> rms, relerr;
            for (size_t ci = 0; ci < cells.size(); ci++)
            {
                const Cell &c = cells[ci];
                double M[4][4] = {{0}}, v[4] = {0}, x[4];
                for (size_t i = 0; i < c.l.size(); i++)
                {
                    const double l = c.l[i], col[4] = {1.0, l, l * l, l * l * l};
                    for (int q = 0; q < 4; q++) { v[q] += c.w[i] * col[q] * c.y[i]; for (int r = 0; r < 4; r++) M[q][r] += c.w[i] * col[q] * col[r]; }
                }
                if (!solve_small(4, M, v, x)) continue;
                double s2 = 0;
                for (size_t i = 0; i < c.l.size(); i++)
                {
                    const double l = c.l[i];
                    const double pf = x[0] + l * (x[1] + l * (x[2] + l * x[3]));
                    const double ps = shared[ci].A + l * (shared[ci].B + l * (glob[2] + l * glob[3]));
                    s2 += (pf - ps) * (pf - ps) / (ps * ps);
                    relerr.push_back(1.0 / std::sqrt(c.w[i]) / c.y[i]);
                }
                rms.push_back(std::sqrt(s2 / c.l.size()));
            }
            auto q = [](std::vector<double> v, double f) { if (v.empty()) return 0.0; std::sort(v.begin(), v.end()); return v[(size_t)(f * (v.size() - 1))]; };
            char line[400];
            std::snprintf(line, sizeof line, "\nScale change from freeing C and D per position (RMS over each position's points, relative):\n"
                          "  median %.2f%%, 90%% of positions below %.2f%%, max %.2f%%; the points' median relative error is %.2f%%\n",
                          100 * q(rms, 0.5), 100 * q(rms, 0.9), 100 * q(rms, 1.0), 100 * q(relerr, 0.5));
            out << line;
        }
    }
    out.close();
    // all-free fits per position
    std::ofstream cf(csv);
    cf << "x_center,y_center,npoints,A,B,C,D,chi2,ndof\n";
    for (const Cell &c : cells)
    {
        double M[4][4] = {{0}}, v[4] = {0}, x[4];
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i], col[4] = {1.0, l, l * l, l * l * l};
            for (int a = 0; a < 4; a++) { v[a] += c.w[i] * col[a] * c.y[i]; for (int b = 0; b < 4; b++) M[a][b] += c.w[i] * col[a] * col[b]; }
        }
        if (!solve_small(4, M, v, x)) continue;
        double ch = 0;
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i];
            const double d = c.y[i] - (x[0] + l * (x[1] + l * (x[2] + l * x[3])));
            ch += c.w[i] * d * d;
        }
        char line[300];
        std::snprintf(line, sizeof line, "%.3f,%.3f,%zu,%.6g,%.6g,%.6g,%.6g,%.3f,%d\n", c.x_center, c.y_center, c.l.size(), x[0], x[1], x[2], x[3], ch, (int)c.l.size() - 4);
        cf << line;
    }
    std::cout << "Uniformity test: " << txt << std::endl;
    std::ifstream back(txt); std::string ln; while (std::getline(back, ln)) std::cout << "  " << ln << std::endl;
}

Grid2D make_nan_grid(int nx, int ny)
{
    return Grid2D(nx, std::vector<double>(ny, std::numeric_limits<double>::quiet_NaN()));
}

// The fitted positions' values of one quantity on the full x/y grid (NaN
// where no position was fit)
Grid2D coarse_map(const std::vector<FitRow> &fits, const std::vector<double> &xs,
                  const std::vector<double> &ys, double (*value)(const FitRow &))
{
    std::map<double, int> ix, iy;
    for (size_t i = 0; i < xs.size(); i++)
        ix[xs[i]] = (int)i;
    for (size_t j = 0; j < ys.size(); j++)
        iy[ys[j]] = (int)j;
    Grid2D Z = make_nan_grid((int)xs.size(), (int)ys.size());
    for (const FitRow &f : fits)
        Z[ix.at(f.x_center)][iy.at(f.y_center)] = value(f);
    return Z;
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
    const std::vector<FitRow> &fits,
    double C, double D,
    double Emin, double Emax,
    const std::vector<double> &xs,
    const std::vector<double> &ys,
    const std::string &save_file,
    Grid2D Ps_out[4])
{
    const size_t nx = xs.size(), ny = ys.size();
    const int k = g_free;
    const double *anchor = k == 2 ? ANCHOR_E : (k == 3 ? ANCHOR_E3 : ANCHOR_E4);
    double la[4];
    for (int m = 0; m < k; m++) la[m] = std::log(anchor[m]);

    // The free part sum_{p<k} par_p l^p at the k anchor energies (the face-wide
    // part is the same everywhere); with k = 2 this is A + B l at 3 and 30 GeV
    Grid2D P[4];
    P[0] = coarse_map(fits, xs, ys, [](const FitRow &f) { return f.A; });
    P[1] = coarse_map(fits, xs, ys, [](const FitRow &f) { return f.B; });
    P[2] = coarse_map(fits, xs, ys, [](const FitRow &f) { return f.C; });
    P[3] = coarse_map(fits, xs, ys, [](const FitRow &f) { return f.D; });
    Grid2D an[4];
    for (int m = 0; m < k; m++)
    {
        an[m] = make_nan_grid((int)nx, (int)ny);
        for (size_t i = 0; i < nx; i++)
            for (size_t j = 0; j < ny; j++)
            {
                double v = 0, lp = 1; // NaN where the position has no fit
                for (int q = 0; q < k; q++) { v += P[q][i][j] * lp; lp *= la[m]; }
                an[m][i][j] = v;
            }
        an[m] = smooth_map(an[m]);
    }
    // back from the anchors: solve sum_q par_q la_m^q = an_m (Vandermonde, exact)
    double V[4][4];
    for (int m = 0; m < k; m++) { double lp = 1; for (int q = 0; q < k; q++) { V[m][q] = lp; lp *= la[m]; } }

    long n_unfit = 0, n_filled = 0, n_outside_fiducial = 0;
    struct SavedBin { int i, j; double p[4]; };
    std::vector<SavedBin> bins_out;
    for (int q = 0; q < 4; q++) Ps_out[q] = make_nan_grid((int)nx, (int)ny);
    for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
        {
            double rhs[4], x[4], M[4][4];
            bool defined = true;
            for (int m = 0; m < k; m++) { rhs[m] = an[m][i][j]; defined = defined && std::isfinite(rhs[m]); }
            for (int m = 0; m < k; m++) for (int q = 0; q < k; q++) M[m][q] = V[m][q];
            if (defined) defined = solve_small(k, M, rhs, x);
            double ps[4];
            for (int q = 0; q < 4; q++)
                ps[q] = !defined ? std::numeric_limits<double>::quiet_NaN() : (q < k ? x[q] : (q == 2 ? C : D));
            for (int q = 0; q < 4; q++) Ps_out[q][i][j] = ps[q];
            if (!std::isfinite(P[0][i][j]))
            {
                n_unfit++;
                if (defined)
                    n_filled++;
            }
            // Outside the buffered fiducial region the map is undefined
            if (!pass_fiducial_cut((float)xs[i], (float)ys[j]))
            {
                n_outside_fiducial++;
                continue;
            }
            if (defined)
                bins_out.push_back({(int)i, (int)j, {ps[0], ps[1], ps[2], ps[3]}});
        }

    std::cout << "Filled " << n_filled << "/" << n_unfit
              << " position bins without a fit with the neighbourhood average of nearby fits" << std::endl;
    std::cout << "Excluded " << n_outside_fiducial
              << " grid cells outside the buffer-restricted fiducial acceptance" << std::endl;
    const long n_total = (long)nx * (long)ny;
    std::cout << "Saving " << bins_out.size() << "/" << n_total << " bins" << std::endl;

    std::ofstream out(save_file.c_str());
    if (!out.is_open())
    {
        std::cerr << "ERROR: cannot open " << save_file << " for writing" << std::endl;
        exit(1);
    }

    out << "{\n";
    // Evaluators refuse a map without this tag (earlier forms)
    out << "  \"model\": \"logcubic\",\n";
    out << "  \"jet_type\": \"" << (g_full_jets ? "full" : "em") << "\",\n";
    // how many of A, B, C, D are per position (2: C, D face-wide); with 3 or 4
    // the face-wide "C" / "D" below are the nested fit's (unused where per-bin)
    out << "  \"free_params\": " << k << ",\n";
    for (int i = 0; i < 4; i++)
        out << "  \"" << MAP_INFO_KEY[i] << "\": \"" << map_info[i] << "\",\n";
    out << "  \"C\": ";
    write_double(out, C);
    out << ",\n  \"D\": ";
    write_double(out, D);
    out << ",\n  \"Emin\": ";
    write_double(out, Emin);
    out << ",\n  \"Emax\": ";
    write_double(out, Emax);
    out << ",\n";
    out << "  \"nx\": " << nx << ",\n";
    out << "  \"ny\": " << ny << ",\n";

    out << "  \"x_grid\": [";
    for (size_t i = 0; i < nx; i++)
    {
        if (i > 0)
            out << ", ";
        write_double(out, xs[i]);
    }
    out << "],\n";

    out << "  \"y_grid\": [";
    for (size_t j = 0; j < ny; j++)
    {
        if (j > 0)
            out << ", ";
        write_double(out, ys[j]);
    }
    out << "],\n";

    out << "  \"bins\": [\n";
    for (size_t k = 0; k < bins_out.size(); k++)
    {
        const SavedBin &b = bins_out[k];
        out << "    {\"i\": " << b.i << ", \"j\": " << b.j << ", \"A\": ";
        write_double(out, b.p[0]);
        out << ", \"B\": ";
        write_double(out, b.p[1]);
        if (k >= 3) { out << ", \"C\": "; write_double(out, b.p[2]); }
        if (k >= 4) { out << ", \"D\": "; write_double(out, b.p[3]); }
        out << "}" << (k + 1 < bins_out.size() ? "," : "") << "\n";
    }
    out << "  ]\n";
    out << "}\n";

    out.close();

    std::cout << "Lookup table saved to " << save_file << std::endl;
}

// Side CSV of per-position chi2/ndof (diagnostics only)
void write_chi2_csv(const std::vector<FitRow> &fits, const std::string &save_file)
{
    std::ofstream out(save_file.c_str());
    out << "x,y,chi2,ndof,n_params,reduced_chi2\n";
    for (const auto &r : fits)
    {
        if (r.ndof <= 0)
            continue;
        out << r.x_center << "," << r.y_center << "," << r.chi2 << "," << r.ndof << "," << g_free << ","
            << (r.chi2 / r.ndof) << "\n";
    }
    out.close();
    std::cout << "chi2/ndof CSV saved to " << save_file << std::endl;
}

// Side CSV of each fitted position's chi2 against the smoothed map (the A, B
// of the JSON at that grid point), over the same points as the direct fit.
// No parameter is fitted to the position's points there, so ndof = points.
void write_chi2_smoothed_csv(const std::vector<Cell> &cells, const Grid2D Ps[4], const std::vector<double> &xs,
                             const std::vector<double> &ys, const std::string &save_file)
{
    std::map<double, int> ix, iy;
    for (size_t i = 0; i < xs.size(); i++)
        ix[xs[i]] = (int)i;
    for (size_t j = 0; j < ys.size(); j++)
        iy[ys[j]] = (int)j;
    std::ofstream out(save_file.c_str());
    out << "x,y,chi2,ndof,n_params,reduced_chi2\n";
    for (const Cell &c : cells)
    {
        const int ci = ix.at(c.x_center), cj = iy.at(c.y_center);
        const double A = Ps[0][ci][cj], B = Ps[1][ci][cj], C = Ps[2][ci][cj], D = Ps[3][ci][cj];
        if (!std::isfinite(A) || !std::isfinite(B) || !std::isfinite(C) || !std::isfinite(D))
            continue;
        double chi2 = 0.0;
        for (size_t i = 0; i < c.l.size(); i++)
        {
            const double l = c.l[i];
            const double d = c.y[i] - (A + l * (B + l * (C + l * D)));
            chi2 += c.w[i] * d * d;
        }
        const int ndof = (int)c.l.size();
        out << c.x_center << "," << c.y_center << "," << chi2 << "," << ndof << ",0," << (chi2 / ndof) << "\n";
    }
    out.close();
    std::cout << "Smoothed-map chi2/ndof CSV saved to " << save_file << std::endl;
}

// Side CSV of the direct (pre-smoothing) per-position fits (diagnostics only)
void write_direct_fits_csv(const std::vector<FitRow> &fits, const std::string &save_file)
{
    std::ofstream out(save_file.c_str());
    out << "x,y,A,B,C,D,Emin,Emax\n";
    for (const auto &r : fits)
        out << r.x_center << "," << r.y_center << "," << r.A << "," << r.B << "," << r.C << ","
            << r.D << "," << r.Emin << "," << r.Emax << "\n";
    out.close();
    std::cout << "Direct (unsmoothed) fit CSV saved to " << save_file << std::endl;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);

    std::string input_file = ROOT_FILE;
    std::string output_file = OUTPUT_JSON;

    if (argc > 1)
        input_file = argv[1];
    if (argc > 2)
        output_file = argv[2];
    for (int i = 0; i < 4 && i + 3 < argc; i++)
        map_info[i] = argv[i + 3];
    g_full_jets = argc > 7 && std::string(argv[7]) == "full";
    std::cout << "Jet type: " << (g_full_jets ? "full (ECAL+HCAL)" : "EM-only") << std::endl;
    if (argc > 8) g_free = std::atoi(argv[8]);
    if (g_free < 2 || g_free > 4) { std::cerr << "ERROR: free parameters per position must be 2, 3 or 4" << std::endl; return 1; }
    std::cout << "Free per position: " << (g_free == 2 ? "A, B (C, D face-wide)" : g_free == 3 ? "A, B, C (D face-wide)" : "A, B, C, D") << std::endl;

    std::cout << "Loading " << input_file << " ..." << std::endl;
    std::vector<BinRow> rows = load_bins(input_file);
    std::cout << "Loaded " << rows.size() << " (E,x,y) bin entries" << std::endl;

    std::vector<Cell> cells = build_cells(rows);
    double C, D;
    fit_global_CD(cells, C, D);
    std::vector<FitRow> fits;
    if (g_free == 2)
        fits = fit_cells(cells, C, D);
    else
    {
        // face-wide values of the non-free parameters from the nested fit on
        // the positions that can take k free (k + 2 points); C, D in the JSON
        // are then those (for k = 4: the k = 2 global fit, for reference)
        std::vector<Cell> fc;
        for (const Cell &c : cells) if ((int)c.l.size() >= g_free + 2) fc.push_back(c);
        double glob[4], ch; long nd;
        if (!fit_nested(fc, g_free, ch, nd, glob)) { std::cerr << "ERROR: nested fit failed" << std::endl; return 1; }
        if (g_free == 3) D = glob[3];
        const double gv[4] = {NAN, NAN, C, D};
        fits = fit_cells_free(fc, g_free, gv);
        std::printf("Nested fit (%d free per position) over %zu positions: chi2/ndof = %.1f/%ld%s\n", g_free, fc.size(), ch, nd,
                    g_free == 3 ? Form(", face-wide D = %.6f", D) : "");
    }
    std::cout << "Fit " << fits.size() << " position bins" << std::endl;
    double Emin = fits.front().Emin, Emax = fits.front().Emax;
    for (const FitRow &f : fits)
    {
        Emin = std::min(Emin, f.Emin);
        Emax = std::max(Emax, f.Emax);
    }
    std::printf("Fitted response range (the evaluators' clamp): %.2f - %.2f GeV\n", Emin, Emax);

    std::set<double> xs_set, ys_set;
    for (const auto &r : rows)
    {
        xs_set.insert(r.x_center);
        ys_set.insert(r.y_center);
    }
    std::vector<double> xs(xs_set.begin(), xs_set.end());
    std::vector<double> ys(ys_set.begin(), ys_set.end());

    Grid2D Ps[4];
    write_calibration_json(fits, C, D, Emin, Emax, xs, ys, output_file, Ps);
    write_chi2_csv(fits, output_file + ".chi2.csv");
    write_chi2_smoothed_csv(cells, Ps, xs, ys, output_file + ".chi2_smoothed.csv");
    write_direct_fits_csv(fits, output_file + ".direct_fits.csv");
    write_uniformity(cells, output_file + ".uniformity.txt", output_file + ".free_fits.csv");

    return 0;
}
