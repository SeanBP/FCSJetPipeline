// DIAGNOSTICS FORK of ../jet_scale/JetEnergyScaleFineGrid.cpp (production) -- same
// algorithm, plus an EM-only-conditional Emax (see Emax_em_only below).
// ECAL+HCAL behavior (em_only=0) is byte-for-byte unchanged from
// production. Once validated against the new EM-only production, merge
// this Emax split back into the production file.

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

// CONFIGURATION

std::string input_dir = "/star/data01/pwg/seanp/tunes/pythia8_nnpdf23lo/JetTrees/";
std::string tree_name = "jetTree";
// input_dir may be overridden via argv[1] in main()

// reflection behavior
bool reflect_to_positive_half = true;

// thresholds
int min_total_entries = 5;

// energy binning -- Emax is set per em_only in main() (see there for why)
double E0 = 10;
double Emax = 200;
int n_bins_E = 20;
double energy_curvature = 1.5;

// EM-only truth jets (JetMatcher.cpp em_only mode clusters truth
// particles from photons only, produced within a tight vertex-
// displacement window of the primary vertex) carry only their jet's
// prompt-photon energy, so they're systematically much lower energy than
// a full hadronic truth jet at the same pT-hat -- the production Emax=200
// is mostly empty space for this sample. Data-driven from the complete
// 21-bin/420-job EM-only pT-hat scan (FCSJetAnalysis/
// compare_energy_distributions_temp.cpp): 98% of EM-only truth jets sit
// below 32.4 GeV vs. 96.4 GeV for the ECAL+HCAL reference (ratio 0.34).
// The cap itself is pushed out to the 99.9% cliff (71.0 GeV) rather than
// the 98% cliff, matching FCSJetAnalysis/JetStatistics.cpp's own
// Emax=70 -- an optimizer sweep there found job coverage barely moved
// (99.74% -> 99.62%) scanning Emax from 32.4 up to 70 GeV, so there was
// no reason to leave that tail out of the fit range. Keep the two files'
// values in sync by hand if either changes.
const double Emax_ecal_hcal = 200;
const double Emax_em_only = 70;

// spatial grid -- set in main() once em_only is known (default rectangle,
// or the ECAL-only one for EM-only JetTrees)
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

struct BinAccum
{
    std::vector<float> resp_vals;
    std::vector<float> scale_vals;
};

// f_med(x) = median[Y|X=x], R_med(x) = median[Y/x|X=x] (Cukierman & Nachman,
// "Mathematical Properties of Numerical Inversion for Jet Calibrations",
// arXiv:1609.05195, Eqs. 29-30) -- replaces the earlier Gaussian-fit-mode
// approach (f_mo/R_mo, Eqs. 3-4). That paper proves (Sec. 3.3.1, Eq. 35)
// that median-based numerical inversion achieves EXACT closure for any
// response function f, whereas mean- or mode-based calibration only closes
// if f is linear (Eqs. 10, 20-23) -- and the non-closure grows with both
// the response curve's curvature f''(x) and the resolution sigma(x) (Eq.
// 24), which is exactly why the mode-based map's downstream closure test
// showed its worst undershoot right where CalibrationMap.cpp's skew-
// Gaussian bump (large f'') coincides with poor low-energy resolution.
// Median is a strictly stronger fix than the median/MAD-SEEDED Gaussian
// fit this replaces: that only fixed the seed against tail bias, it didn't
// change which central-tendency definition was actually being reported.
double median_of(std::vector<double> v)
{
    size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    double median = v[mid];
    if (v.size() % 2 == 0)
    {
        std::nth_element(v.begin(), v.begin() + mid - 1, v.end());
        median = 0.5 * (median + v[mid - 1]);
    }
    return median;
}

FitResult median_fit(const std::vector<double> &v)
{
    FitResult r;

    if (v.size() < (size_t)min_total_entries)
        return r;

    double med = median_of(v);

    // MAD-based robust sigma, used only to get an approximate standard
    // error on the median below (Kenney & Keeping's asymptotic result for
    // an approximately-Gaussian-core distribution: SE_median ~
    // sqrt(pi/2)*sigma/sqrt(n) = 1.2533*sigma/sqrt(n)) -- this feeds
    // CalibrationMap.cpp's downstream fit only as a relative point weight,
    // not as a rigorous confidence interval, so the Gaussian-core
    // approximation is good enough here.
    std::vector<double> absdev(v.size());
    for (size_t i = 0; i < v.size(); i++)
        absdev[i] = std::fabs(v[i] - med);
    double mad_sigma = 1.4826 * median_of(absdev);

    if (!std::isfinite(mad_sigma) || mad_sigma <= 0.0)
        return r;

    r.mu = med;
    r.mu_err = 1.2533 * mad_sigma / std::sqrt((double)v.size());
    r.sigma = mad_sigma;
    // sigma_err/chi2_ndf are meaningless for a plain median (no fit was
    // performed) and aren't consumed downstream -- left at their NaN
    // defaults from FitResult's constructor.

    return r;
}

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc > 1)
        input_dir = argv[1];

    // argv[2]: 1 = these JetTrees are EM-only (ECAL hits only, no HCAL) --
    // use the ECAL-only fiducial rectangle/buffer instead of the default
    // ECAL+HCAL one. 0/omit = previous default behavior.
    bool em_only = (argc > 2) && (std::string(argv[2]) == "1");

    // argv[3]: optional match_radius override, for testing statistics-vs-
    // smoothing tradeoffs without editing this file each time (see
    // match_radius's declaration above for what it controls).
    if (argc > 3)
        match_radius = std::atof(argv[3]);

    // argv[4]: optional E0 (truth-energy binning floor) override. Default
    // 10 GeV left the lowest response points around ~11.6 GeV, which
    // turned out to be above where CalibrationMap.cpp's Gaussian bump
    // fits actually want their peak (E0~6 GeV, found via direct chi2 scan
    // of r10 production data) -- the fit was extrapolating into a region
    // with no real data. Lowering this gives actual low-E statistics
    // there instead of relying on extrapolation.
    if (argc > 4)
        E0 = std::atof(argv[4]);

    const FiducialRect &rect = em_only ? kFiducialRectEcal : kFiducialRect;
    float fid_buffer = em_only ? reco_fiducial_buffer_ecal : reco_fiducial_buffer;

    Emax = em_only ? Emax_em_only : Emax_ecal_hcal;

    // argv[5]: optional Emax (truth-energy binning ceiling) override. Truth
    // jets above Emax are dropped outright (FindBin returns -1), not merged
    // into the top bin -- so raising this recovers real statistics rather
    // than just extending an empty axis. min_total_entries still gates each
    // resulting (position, E) bin, so any new high-E bins that end up too
    // sparse are simply skipped rather than polluting the fit.
    if (argc > 5)
        Emax = std::atof(argv[5]);

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

    // Accumulate per-bin values directly (memory-efficient)

    std::unordered_map<BinKey, BinAccum, BinKeyHash> bin_accum;

    for (size_t fidx = 0; fidx < root_files.size(); fidx++)
    {
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

                    bin_accum[key].resp_vals.push_back((float)reco);
                    bin_accum[key].scale_vals.push_back((float)scale);
                }
            }
        }

        f->Close();
        delete f;
    }

    // Fit per-bin

    std::unordered_map<BinKey, BinData, BinKeyHash> bin_mu_map;

    for (auto it = bin_accum.begin(); it != bin_accum.end(); ++it)
    {
        const std::vector<float> &resp_f = it->second.resp_vals;
        const std::vector<float> &scale_f = it->second.scale_vals;

        if (resp_f.size() < (size_t)min_total_entries)
            continue;

        std::vector<double> resp_vals(resp_f.begin(), resp_f.end());
        std::vector<double> scale_vals(scale_f.begin(), scale_f.end());

        FitResult resp = median_fit(resp_vals);
        FitResult sc = median_fit(scale_vals);

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

    t.Branch("n_entries", &n_entries);

    for (int ie = 0; ie < n_bins_E; ie++)
    for (int ix = 0; ix < n_grid_x; ix++)
    for (int iy = 0; iy < n_grid_y; iy++)
    for (int q = 0; q < (reflect_to_positive_half ? 1 : 2); q++)
    {
        BinKey key = {q, ix, iy, ie};

        auto it_acc = bin_accum.find(key);
        n_entries = (it_acc != bin_accum.end()) ? (int)it_acc->second.resp_vals.size() : 0;

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