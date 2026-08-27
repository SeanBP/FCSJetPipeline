#include <TFile.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TF1.h>

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
int n_hist_bins = 50;
double asym = 1.0;
double n_sigma = 1.5;

// energy binning
double E0 = 10;
double Emax = 200;
int n_bins_E = 20;
double energy_curvature = 1.5;

// spatial grid
double x_min = cut_x_inner + cut_buffer;
double x_max = cut_x_outer - cut_buffer;

double y_min = cut_y_min + cut_buffer;
double y_max = cut_y_max - cut_buffer;

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

FitResult gaussian_fit(const std::vector<double> &v)
{
    FitResult r;

    if (v.size() < (size_t)min_total_entries)
        return r;

    double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();

    double var = 0.0;
    for (size_t i = 0; i < v.size(); i++)
        var += (v[i] - mean) * (v[i] - mean);

    double sigma = std::sqrt(var / std::max<size_t>(1, v.size() - 1));
    if (!std::isfinite(sigma) || sigma <= 0.0)
        return r;

    TH1D h("h", "", n_hist_bins, mean - 5 * sigma, mean + 5 * sigma);

    for (size_t i = 0; i < v.size(); i++)
        h.Fill(v[i]);

    double xmin = mean - n_sigma * sigma;
    double xmax = mean + n_sigma * sigma * asym;

    TF1 f("f", "gaus", xmin, xmax);
    f.SetParameters(h.GetMaximum(), mean, sigma);

    h.Fit(&f, "RQ");

    r.mu = f.GetParameter(1);
    r.mu_err = f.GetParError(1);
    r.sigma = f.GetParameter(2);
    r.sigma_err = f.GetParError(2);

    if (f.GetNDF() > 0)
        r.chi2_ndf = f.GetChisquare() / f.GetNDF();

    return r;
}

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc > 1)
        input_dir = argv[1];

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
        grid_valid[ix][iy] = pass_fiducial_cut((float)x_grid[ix], (float)y_grid[iy], reco_fiducial_buffer);
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

        FitResult resp = gaussian_fit(resp_vals);
        FitResult sc = gaussian_fit(scale_vals);

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