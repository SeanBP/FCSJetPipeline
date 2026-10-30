// Step 2 of the EM-only JES calibration: builds the JES lookup JSON from the
// "bins" tree of JetEnergyScaleFineGrid_em_only (per-bin weighted medians).
// The scale vs response at each (x,y) position bin is described by
//   scale(E) = A + B lnE + C ln^2 E + D ln^3 E
// with A and B per position and C and D shared by the whole face: the shape
// in ln E is universal and only the level and slope vary across the ECal.
// C and D come from one least-squares fit over every position's points (A and
// B free per position, see fit_global_CD); A and B are then the per-position
// fits with C and D held there, which is the joint optimum since the model is
// linear. The A and B maps are smoothed over the grid (as the scale at two
// anchor energies, see ANCHOR_E) and masked to the buffered ECal fiducial
// region. The polynomial has no plateau, so evaluators clamp E to the range
// of all fitted responses [Emin, Emax], holding the scale flat outside the
// data (as CMSSW's SimpleJetCorrector clamps pT to each record's range). The
// range is the whole face's, not each position's: the shape is constrained
// over all of it, and past a position's own last point only its level and
// slope are carried along that shape (a per-position clamp held the outer
// positions, with no reliable points above about 40 GeV, flat from there,
// overcorrecting their high-E jets by up to 5%).
//
// The JSON carries "model": "logcubic", the global "C", "D", "Emin" and
// "Emax", and bins {i, j, A, B} (read by jet_calibration/ApplyCorrections_em_only.cpp
// and jet_scale/jes_closure.cpp). It also writes <out>.chi2.csv and
// <out>.direct_fits.csv (the unsmoothed per-position fits).
//
// Usage: ./CalibrationMap_em_only <bins.root> <out.json>
// Build: g++ -m64 -O2 CalibrationMap_em_only.cpp -o CalibrationMap_em_only `root-config --cflags --libs` -std=c++11

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

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

// Sanity cuts dropping the noise tail of low-statistics bins
const double RESPONSE_SANITY_FACTOR = 1.5;
const double SCALE_SANITY_MAX = 1.5;

// Positions with fewer valid energy bins are not fit (they are filled by the smoothing)
const int MIN_POINTS = 4;

// Energy bins with a smaller Kish effective count are not used: a weighted
// median of one or two effective jets is biased (jesmap_nofilt_raw: the 417
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

// Copied from shared/JetParameters.h (keep in sync; that header needs
// fastjet). Upstream, JetEnergyScaleFineGrid only fills grid points inside
// this region, so cells outside it have no statistics; the saved map is
// masked to it so smoothing never extrapolates past the acceptance.
const float cut_x_inner_ecal = 20.53604584f, cut_x_outer_ecal = 137.59500000f;
const float cut_y_min_ecal = -96.72439178f, cut_y_max_ecal = 86.25565167f;
const float reco_fiducial_buffer_ecal = 20.0f;

inline bool pass_fiducial_cut(float jetXE, float jetYE)
{
    const float buffer = reco_fiducial_buffer_ecal;
    const float ax = std::fabs(jetXE);
    const bool in_x = ax >= (cut_x_inner_ecal + buffer) && ax <= (cut_x_outer_ecal - buffer);
    const bool in_y = jetYE >= (cut_y_min_ecal + buffer) && jetYE <= (cut_y_max_ecal - buffer);
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
    double A, B;
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
    const std::string &save_file)
{
    const size_t nx = xs.size(), ny = ys.size();
    const double l0 = std::log(ANCHOR_E[0]), l1 = std::log(ANCHOR_E[1]);

    // A + B l at the anchor energies (the C, D part is the same everywhere)
    Grid2D A = coarse_map(fits, xs, ys, [](const FitRow &f) { return f.A; });
    Grid2D B = coarse_map(fits, xs, ys, [](const FitRow &f) { return f.B; });
    Grid2D a0 = make_nan_grid((int)nx, (int)ny), a1 = make_nan_grid((int)nx, (int)ny);
    for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
        {
            a0[i][j] = A[i][j] + B[i][j] * l0; // NaN where the position has no fit
            a1[i][j] = A[i][j] + B[i][j] * l1;
        }
    a0 = smooth_map(a0);
    a1 = smooth_map(a1);

    long n_unfit = 0, n_filled = 0, n_outside_fiducial = 0;
    struct SavedBin { int i, j; double A, B; };
    std::vector<SavedBin> bins_out;
    for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
        {
            const double Bs = (a1[i][j] - a0[i][j]) / (l1 - l0);
            const double As = a0[i][j] - Bs * l0;
            const bool defined = std::isfinite(As) && std::isfinite(Bs);
            if (!std::isfinite(A[i][j]))
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
                bins_out.push_back({(int)i, (int)j, As, Bs});
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
        write_double(out, b.A);
        out << ", \"B\": ";
        write_double(out, b.B);
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
        out << r.x_center << "," << r.y_center << "," << r.chi2 << "," << r.ndof << ",2,"
            << (r.chi2 / r.ndof) << "\n";
    }
    out.close();
    std::cout << "chi2/ndof CSV saved to " << save_file << std::endl;
}

// Side CSV of the direct (pre-smoothing) per-position fits (diagnostics only)
void write_direct_fits_csv(const std::vector<FitRow> &fits, double C, double D,
                           const std::string &save_file)
{
    std::ofstream out(save_file.c_str());
    out << "x,y,A,B,C,D,Emin,Emax\n";
    for (const auto &r : fits)
        out << r.x_center << "," << r.y_center << "," << r.A << "," << r.B << "," << C << ","
            << D << "," << r.Emin << "," << r.Emax << "\n";
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

    std::cout << "Loading " << input_file << " ..." << std::endl;
    std::vector<BinRow> rows = load_bins(input_file);
    std::cout << "Loaded " << rows.size() << " (E,x,y) bin entries" << std::endl;

    std::vector<Cell> cells = build_cells(rows);
    double C, D;
    fit_global_CD(cells, C, D);
    std::vector<FitRow> fits = fit_cells(cells, C, D);
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

    write_calibration_json(fits, C, D, Emin, Emax, xs, ys, output_file);
    write_chi2_csv(fits, output_file + ".chi2.csv");
    write_direct_fits_csv(fits, C, D, output_file + ".direct_fits.csv");

    return 0;
}
