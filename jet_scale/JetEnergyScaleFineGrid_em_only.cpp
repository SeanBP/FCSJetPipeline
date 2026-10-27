// Step 1 of the EM-only JES calibration: per position on a fine x-y grid and
// per truth-energy bin, the weighted median scale and response of matched
// EM-only jets, written to a "bins" tree for CalibrationMap_em_only.
// Medians give exact closure under numerical inversion for any response
// shape (Cukierman & Nachman, arXiv:1609.05195, Eq. 35).
//
// The underlying event is subtracted at both levels by default:
// scale = (reco_E - reco_ue_E) / (truth_E - truth_ue_E), binned in subtracted
// truth energy, using pairs where both jets have an estimate. Input should be
// a calibration sample (sim_to_jet -C 1, UE rerun in minbias mode). A file
// without the UE branches is an error unless "noue" is given.
//
// Usage: JetEnergyScaleFineGrid_em_only <jettrees_dir/> [match_radius_cm]
//        [E0_GeV] [Emax_GeV] [noue]   (an empty string "" keeps that slot's default)

#include <TFile.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TROOT.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"

// CONFIGURATION

// Hardcoded default input; normally given as argv[1]
std::string input_dir = "/star/data01/pwg/seanp/tunes/pythia8_nnpdf23lo/JetTrees/";
std::string tree_name = "jetTree";

// reflection behavior
bool reflect_to_positive_half = true;

// thresholds
int min_total_entries = 5;

// Truth-energy binning, GeV; edges follow a power law in energy_curvature
// (finer at low E). EM-only truth jets are mostly below ~32 GeV; E0 = 3 GeV
// reaches below the response bump near 6 GeV. Jets above Emax are dropped.
double E0 = 3;
double Emax = 100;
int n_bins_E = 20;
double energy_curvature = 1.5;

// spatial grid -- set in main() from the ECAL-only fiducial rectangle
double x_min, x_max, y_min, y_max;

// fine grid of points across the x-y surface (adjustable spacing, same units as x/y);
// kept smaller than match_radius so neighboring points' circles overlap
double grid_spacing = 2.5;

// jets within this x-y (detector surface) distance of a grid point are
// included in its histogram
double match_radius = 5.0;

struct FitResult
{
    double mu;
    double mu_err;
    double sigma;
    double sigma_err;
    double chi2_ndf;

    FitResult()
        : mu(NAN),
          mu_err(NAN),
          sigma(NAN),
          sigma_err(NAN),
          chi2_ndf(NAN) {}
};

inline int FindBin(double v, const std::vector<double> &b)
{
    if (v < b.front() || v >= b.back())
        return -1;

    int low = 0;
    int high = (int)b.size() - 1;

    while (low < high)
    {
        int mid = (low + high) / 2;
        if (v >= b[mid + 1])
            low = mid + 1;
        else
            high = mid;
    }
    return low;
}

struct BinKey
{
    int q, ix, iy, ie;

    bool operator==(const BinKey &o) const
    {
        return q == o.q && ix == o.ix && iy == o.iy && ie == o.ie;
    }
};

struct BinKeyHash
{
    std::size_t operator()(const BinKey &k) const
    {
        return ((size_t)k.q << 24) ^
               ((size_t)k.ix << 16) ^
               ((size_t)k.iy << 8) ^
               (size_t)k.ie;
    }
};

struct BinData
{
    double scale_mu;
    double resp_mu;
    double scale_mu_err;
    double resp_mu_err;

    BinData()
        : scale_mu(-1),
          resp_mu(-1),
          scale_mu_err(NAN),
          resp_mu_err(NAN) {}

    BinData(double s, double r, double se, double re)
        : scale_mu(s),
          resp_mu(r),
          scale_mu_err(se),
          resp_mu_err(re) {}
};

// w_vals is each entry's MC cross-section weight (see JetXSecWeights.h),
// parallel to resp_vals/scale_vals.
struct BinAccum
{
    std::vector<float> resp_vals;
    std::vector<float> scale_vals;
    std::vector<float> w_vals;
};

// Weighted median (cumulative weight reaches half the total; with equal
// weights, the ordinary median) of each bin: f_med(x) = median[Y|X=x],
// R_med(x) = median[Y/x|X=x] (arXiv:1609.05195, Eqs. 29-30).
double median_of(const std::vector<double> &v, const std::vector<double> &w)
{
    std::vector<size_t> order(v.size());
    for (size_t i = 0; i < order.size(); i++)
        order[i] = i;
    std::sort(order.begin(), order.end(),
              [&v](size_t a, size_t b) { return v[a] < v[b]; });

    double half = 0.5 * std::accumulate(w.begin(), w.end(), 0.0);

    double cum = 0.0;
    for (size_t i = 0; i < order.size(); i++)
    {
        cum += w[order[i]];
        if (cum > half)
            return v[order[i]];
        if (cum == half && i + 1 < order.size())
            return 0.5 * (v[order[i]] + v[order[i + 1]]);
    }
    return v[order.back()];
}

FitResult median_fit(const std::vector<double> &v, const std::vector<double> &w)
{
    FitResult r;

    if (v.size() < (size_t)min_total_entries)
        return r;

    double med = median_of(v, w);

    // MAD-based sigma, only for an approximate standard error of the median
    // (1.2533 sigma / sqrt(n)), used downstream as a fit weight. A bin where
    // one entry holds over half the weight gets 0 and is rejected. The
    // weighted MAD also collapses when one entry holds nearly half the weight
    // and a neighbour is close to it (errors of 0.01-0.2% for a spread of a
    // few percent, at n_eff of 3-7), so it is floored by the unweighted MAD:
    // the resolution at a point does not depend on which ptHat range a jet
    // came from, and no single jet can hold the unweighted one.
    std::vector<double> absdev(v.size());
    for (size_t i = 0; i < v.size(); i++)
        absdev[i] = std::fabs(v[i] - med);
    double mad_sigma = 1.4826 * median_of(absdev, w);
    if (mad_sigma > 0.0)
        mad_sigma = std::max(mad_sigma, 1.4826 * median_of(absdev, std::vector<double>(v.size(), 1.0)));

    if (!std::isfinite(mad_sigma) || mad_sigma <= 0.0)
        return r;

    r.mu = med;
    // n is the Kish effective count
    r.mu_err = 1.2533 * mad_sigma / std::sqrt(EffectiveEntries(w));
    r.sigma = mad_sigma;
    // sigma_err/chi2_ndf stay NaN (no fit)

    return r;
}

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc > 1)
        input_dir = argv[1];

    // Optional overrides of the file-scope defaults above; an empty string
    // keeps the default, so a later slot can be set without an earlier one.
    //   argv[2]: match_radius (cm) -- x-y pooling radius around each grid
    //            point (see its declaration above)
    //   argv[3]: E0 (GeV) -- truth-energy binning floor
    //   argv[4]: Emax (GeV) -- truth-energy binning ceiling
    if (argc > 2 && argv[2][0] != '\0')
        match_radius = std::atof(argv[2]);
    if (argc > 3 && argv[3][0] != '\0')
        E0 = std::atof(argv[3]);
    if (argc > 4 && argv[4][0] != '\0')
        Emax = std::atof(argv[4]);
    //   argv[5]: "noue" -- no underlying-event subtraction (see header)
    const bool subtract_ue = !(argc > 5 && std::string(argv[5]) == "noue");
    std::cout << (subtract_ue ? "Underlying event subtracted from reco and truth jet energies"
                              : "No underlying-event subtraction (noue)") << std::endl;
    long n_pairs = 0, n_no_ue = 0;

    const FiducialRect &rect = kFiducialRectEcal;
    float fid_buffer = reco_fiducial_buffer_ecal;

    x_min = rect.x_inner + cut_buffer;
    x_max = rect.x_outer - cut_buffer;

    y_min = rect.y_min + cut_buffer;
    y_max = rect.y_max - cut_buffer;

    int qmax = reflect_to_positive_half ? 1 : 2;

    TSystemDirectory dir("input_dir", input_dir.c_str());
    TList *files = dir.GetListOfFiles();

    std::vector<std::string> root_files;

    TIter next(files);
    TSystemFile *file;

    while ((file = (TSystemFile *)next()))
    {
        std::string fname = file->GetName();
        if (!file->IsDirectory() && fname.find(".root") != std::string::npos)
            root_files.push_back(input_dir + fname);
    }

    // Precompute bin edges

    std::vector<double> e_bins(n_bins_E + 1);
    for (int i = 0; i <= n_bins_E; i++)
        e_bins[i] = E0 + (Emax - E0) * pow((double)i / n_bins_E, energy_curvature);

    // fine grid of points across the x-y surface
    int n_grid_x = (int)std::floor((x_max - x_min) / grid_spacing) + 1;
    int n_grid_y = (int)std::floor((y_max - y_min) / grid_spacing) + 1;

    std::vector<double> x_grid(n_grid_x);
    std::vector<double> y_grid(n_grid_y);

    for (int i = 0; i < n_grid_x; i++)
        x_grid[i] = x_min + i * grid_spacing;

    for (int i = 0; i < n_grid_y; i++)
        y_grid[i] = y_min + i * grid_spacing;

    // eta/phi of every grid point, and whether it's in the reconstructable
    // fiducial region (same buffer cut as JetMatcher's reco jet selection)
    std::vector<std::vector<double>> grid_eta(n_grid_x, std::vector<double>(n_grid_y));
    std::vector<std::vector<double>> grid_phi(n_grid_x, std::vector<double>(n_grid_y));
    std::vector<std::vector<bool>> grid_valid(n_grid_x, std::vector<bool>(n_grid_y));

    for (int ix = 0; ix < n_grid_x; ix++)
    for (int iy = 0; iy < n_grid_y; iy++)
    {
        computeEtaPhi(x_grid[ix], y_grid[iy], grid_eta[ix][iy], grid_phi[ix][iy]);
        grid_valid[ix][iy] = pass_fiducial_cut((float)x_grid[ix], (float)y_grid[iy], fid_buffer, rect);
    }

    // Per-file MC cross-section weight; every jet below is weighted by its
    // file's value instead of counting evenly (see JetXSecWeights.h).
    std::vector<float> file_weight;
    if (!ComputeXSecWeights(root_files, tree_name, file_weight))
        return 1;

    // Accumulate per-bin values directly (memory-efficient)

    std::unordered_map<BinKey, BinAccum, BinKeyHash> bin_accum;

    for (size_t fidx = 0; fidx < root_files.size(); fidx++)
    {
        if (file_weight[fidx] < 0)
            continue;

        TFile *f = TFile::Open(root_files[fidx].c_str());
        if (!f || f->IsZombie())
            continue;

        TTree *tree = (TTree *)f->Get(tree_name.c_str());
        if (!tree)
        {
            f->Close();
            delete f;
            continue;
        }

        std::vector<float> *reco_E = nullptr;
        std::vector<float> *reco_x = nullptr;
        std::vector<float> *reco_y = nullptr;
        std::vector<float> *truth_E = nullptr;
        std::vector<int> *match_truth_idx = nullptr;
        std::vector<int> *match_reco_idx = nullptr;

        tree->SetBranchAddress("reco_E", &reco_E);
        tree->SetBranchAddress("reco_x", &reco_x);
        tree->SetBranchAddress("reco_y", &reco_y);
        tree->SetBranchAddress("truth_E", &truth_E);
        tree->SetBranchAddress("match_truth_idx", &match_truth_idx);
        tree->SetBranchAddress("match_reco_idx", &match_reco_idx);

        std::vector<float> *reco_ue_E = nullptr, *truth_ue_E = nullptr;
        std::vector<int> *reco_ue_nproxy = nullptr, *truth_ue_ncone = nullptr;
        if (subtract_ue)
        {
            if (!tree->GetBranch("reco_ue_E") || !tree->GetBranch("truth_ue_E"))
            {
                std::cerr << "ERROR: " << root_files[fidx] << " has no reco_ue_E/truth_ue_E branches;"
                          << " run underlying_event/UnderlyingEvent on the sample first, or pass noue" << std::endl;
                return 1;
            }
            tree->SetBranchAddress("reco_ue_E", &reco_ue_E);
            tree->SetBranchAddress("truth_ue_E", &truth_ue_E);
            tree->SetBranchAddress("reco_ue_nproxy", &reco_ue_nproxy);
            tree->SetBranchAddress("truth_ue_ncone", &truth_ue_ncone);
        }

        Long64_t n_events = tree->GetEntries();

        for (Long64_t ev = 0; ev < n_events; ev++)
        {
            tree->GetEntry(ev);

            for (size_t m = 0; m < match_reco_idx->size(); m++)
            {
                int ri = match_reco_idx->at(m);
                int ti = match_truth_idx->at(m);

                double reco = reco_E->at(ri);
                double truth = truth_E->at(ti);
                n_pairs++;
                if (subtract_ue)
                {
                    if (reco_ue_nproxy->at(ri) <= 0 || truth_ue_ncone->at(ti) <= 0)
                    {
                        n_no_ue++;
                        continue;
                    }
                    reco -= reco_ue_E->at(ri);
                    truth -= truth_ue_E->at(ti);
                }

                double x = reflect_to_positive_half ? std::fabs(reco_x->at(ri)) : reco_x->at(ri);
                double y = reco_y->at(ri);

                int eb = FindBin(truth, e_bins);
                if (eb < 0)
                    continue;

                int q = (!reflect_to_positive_half && reco_x->at(ri) < 0) ? 1 : 0;

                double scale = reco / truth;

                // Include this jet at every grid point within match_radius (index
                // range computed directly since the grid is evenly spaced).
                int ix_lo = std::max(0, (int)std::floor((x - match_radius - x_min) / grid_spacing));
                int ix_hi = std::min(n_grid_x - 1, (int)std::ceil((x + match_radius - x_min) / grid_spacing));
                int iy_lo = std::max(0, (int)std::floor((y - match_radius - y_min) / grid_spacing));
                int iy_hi = std::min(n_grid_y - 1, (int)std::ceil((y + match_radius - y_min) / grid_spacing));

                for (int ix = ix_lo; ix <= ix_hi; ix++)
                for (int iy = iy_lo; iy <= iy_hi; iy++)
                {
                    if (!grid_valid[ix][iy])
                        continue;

                    double dx = x - x_grid[ix];
                    double dy = y - y_grid[iy];
                    double dr2 = dx * dx + dy * dy;

                    if (dr2 > match_radius * match_radius)
                        continue;

                    BinKey key = {q, ix, iy, eb};

                    BinAccum &acc = bin_accum[key];

                    acc.resp_vals.push_back((float)reco);
                    acc.scale_vals.push_back((float)scale);
                    acc.w_vals.push_back(file_weight[fidx]);
                }
            }
        }

        f->Close();
        delete f;
    }

    if (subtract_ue)
        std::cout << n_pairs << " matched pairs, " << n_no_ue << " left out for lack of a UE estimate" << std::endl;

    // Fit per-bin

    std::unordered_map<BinKey, BinData, BinKeyHash> bin_mu_map;

    for (auto it = bin_accum.begin(); it != bin_accum.end(); ++it)
    {
        const std::vector<float> &resp_f = it->second.resp_vals;
        const std::vector<float> &scale_f = it->second.scale_vals;
        const std::vector<float> &w_f = it->second.w_vals;

        if (resp_f.size() < (size_t)min_total_entries)
            continue;

        std::vector<double> resp_vals(resp_f.begin(), resp_f.end());
        std::vector<double> scale_vals(scale_f.begin(), scale_f.end());
        std::vector<double> w_vals(w_f.begin(), w_f.end());

        FitResult resp = median_fit(resp_vals, w_vals);
        FitResult sc = median_fit(scale_vals, w_vals);

        if (!std::isfinite(sc.mu) || !std::isfinite(resp.mu))
            continue;

        bin_mu_map[it->first] = BinData(sc.mu, resp.mu, sc.mu_err, resp.mu_err);
    }

    // Output

    TFile fout("jet_calibration.root", "RECREATE");
    TTree t("bins", "calibration bins");

    double E_low, E_high;
    double x_center, y_center;
    double eta_center, phi_center;

    double scale_mu, resp_mu;
    double scale_mu_err, resp_mu_err;
    int n_entries;
    double n_eff;

    t.Branch("E_low", &E_low);
    t.Branch("E_high", &E_high);
    t.Branch("x_center", &x_center);
    t.Branch("y_center", &y_center);
    t.Branch("eta_center", &eta_center);
    t.Branch("phi_center", &phi_center);

    t.Branch("scale_mu", &scale_mu);
    t.Branch("response_mu", &resp_mu);
    t.Branch("scale_mu_err", &scale_mu_err);
    t.Branch("response_mu_err", &resp_mu_err);

    // n_entries is the raw (unweighted) jet count; n_eff is the Kish
    // effective count under the cross-section weights, (sum w)^2 / sum w^2
    t.Branch("n_entries", &n_entries);
    t.Branch("n_eff", &n_eff);

    for (int ie = 0; ie < n_bins_E; ie++)
    for (int ix = 0; ix < n_grid_x; ix++)
    for (int iy = 0; iy < n_grid_y; iy++)
    for (int q = 0; q < (reflect_to_positive_half ? 1 : 2); q++)
    {
        BinKey key = {q, ix, iy, ie};

        auto it_acc = bin_accum.find(key);
        n_entries = (it_acc != bin_accum.end()) ? (int)it_acc->second.resp_vals.size() : 0;
        n_eff = 0.0;
        if (it_acc != bin_accum.end())
            n_eff = EffectiveEntries(std::vector<double>(it_acc->second.w_vals.begin(), it_acc->second.w_vals.end()));

        E_low = e_bins[ie];
        E_high = e_bins[ie + 1];

        x_center = x_grid[ix];
        y_center = y_grid[iy];

        eta_center = grid_eta[ix][iy];
        phi_center = grid_phi[ix][iy];

        const BinData &d =
            bin_mu_map.count(key) ? bin_mu_map[key] : BinData();

        scale_mu = d.scale_mu;
        resp_mu = d.resp_mu;
        scale_mu_err = d.scale_mu_err;
        resp_mu_err = d.resp_mu_err;

        t.Fill();
    }

    t.Write();
    fout.Close();

    std::cout << "Output written to jet_calibration.root" << std::endl;
}