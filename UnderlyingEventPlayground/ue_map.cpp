// Underlying-event (UE) energy map for EM-only reco jets, with a truth-level
// jet-contamination diagnostic. Workshop version -- see README.md.
//
// Input: a folder of sim_to_jet output made with -E 1 -F 1 -s 1, i.e.
// jet_output_<suffix>.root (JetTree) and SimpleTree_<suffix>.root side by
// side. The JetTree gives reco jets and truth-jet constituents; the
// SimpleTree gives every FCS point in the event (the JetTree only keeps jet
// constituents). The two trees are read entry-by-entry in step, which is
// valid when sim_to_jet ran without a trigger filter (one jetTree entry per
// SimpleTree event); a pair with differing entry counts is skipped.
//
// The map covers half the FCS ECAL (x > 0), on the same fine x-y grid as the
// JES calibration (JetEnergyScaleFineGrid_em_only.cpp: same origin, spacing
// and fiducial validity). Every event is used twice, once as is and once
// left-right mirrored, which is done here by also evaluating the mirrored
// grid point (-x, y) in the unmirrored event -- the same thing.
//
// For every grid point P and every event (and its mirror):
//   1. veto:  no reco jet within 2R of P
//   2. tag:   at least one reco jet whose distance, min(dR(J,P), dR(J,P')),
//             is within ue_tag_tolerance_in_R * R of the best achievable one
//             for P (UESelection.h; P' = 180 deg in azimuth from P)
//   3. UE:    E_UE = summed energy of the FCS points within R of P, using the
//             same point selection as the jets (ECAL, pass_point_energy_cut)
//   4. contamination (diagnostic): E_contam = summed energy of TRUTH-JET
//             constituents within R of P
// Every entry is weighted by its event's cross-section weight
// (jet_scale/JetXSecWeights.h, the same weights the JES fits use).
//
// Per point, the E_UE distribution gives a mean and sigma, and the E_contam
// distribution gives the fraction of events with any contamination plus the
// mean and sigma of the nonzero part. Each mean/sigma comes from a Gaussian
// fit when one can be trusted, otherwise from the weighted numerical mean
// and standard deviation; the *_method grids say which (1 = fit, 0 =
// numerical) and the summary printed at the end counts them.
//
// Build and run from this folder:
//   g++ -m64 ue_map.cpp ../shared/VectorDict.cxx -o ue_map -I. -I../shared -I../jet_scale \
//       `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//   ./ue_map <input_dir/> [output.root] [grid_spacing_cm]
//
// Writes <output.root> (default ue_map.root: grids as TH2D, a "points" tree
// with everything per grid point, and the per-point distributions under
// hists/) and a matching .pdf of the grids.

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TDirectory.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"
#include "UESelection.h"

// CONFIGURATION

// x-y grid spacing (cm); same default as the JES fine grid
double grid_spacing = 2.5;

// sampling of the candidate jet-axis positions used to find each point's
// best achievable jet distance (cm)
double candidate_spacing = 1.0;

// cone radius around P for the UE and contamination sums
const double cone_radius = R;

// distributions: energy axis (GeV). Sums beyond the range still enter the
// numerical mean/std, which are computed from running sums, not the histogram.
int n_hist_bins = 300;
double hist_max = 150.0;

// A Gaussian fit is used only if the point has at least this many raw
// entries and effective entries, the fit converges to finite, in-range
// parameters, and its chi2/ndf is no worse than max_fit_chi2_ndf. Otherwise
// the weighted numerical mean/std are used and the point is flagged.
int min_fit_entries = 50;
double min_fit_n_eff = 30.0;
double max_fit_chi2_ndf = 5.0;

struct Deposit
{
    double eta, phi, E;
};

// Weighted distribution of one quantity at one grid point
struct Accum
{
    TH1D *hist;
    int n;              // raw entries
    double sw, sw2;     // sum of weights, of squared weights
    double swx, swxx;   // weighted sums of value, of value squared

    Accum() : hist(0), n(0), sw(0), sw2(0), swx(0), swxx(0) {}

    void Fill(double x, double w)
    {
        hist->Fill(x, w);
        n++;
        sw += w;
        sw2 += w * w;
        swx += w * x;
        swxx += w * x * x;
    }

    double Mean() const { return sw > 0 ? swx / sw : NAN; }
    double Std() const
    {
        if (!(sw > 0)) return NAN;
        double var = swxx / sw - Mean() * Mean();
        return var > 0 ? std::sqrt(var) : 0.0;
    }
    double NEff() const { return sw2 > 0 ? sw * sw / sw2 : 0.0; }
};

struct Summary
{
    double mean, sigma;       // adopted values
    int method;               // 1 = Gaussian fit, 0 = numerical, -1 = no entries
    double num_mean, num_std; // weighted numerical values
    double fit_mean, fit_sigma, fit_chi2_ndf;

    Summary() : mean(NAN), sigma(NAN), method(-1), num_mean(NAN), num_std(NAN),
                fit_mean(NAN), fit_sigma(NAN), fit_chi2_ndf(NAN) {}
};

// Gaussian fit of the accumulated distribution, falling back to the
// numerical mean/std when the fit cannot be trusted (see the thresholds
// above).
Summary Summarize(Accum &a)
{
    Summary s;
    if (a.n == 0)
        return s;

    s.num_mean = a.Mean();
    s.num_std = a.Std();
    s.mean = s.num_mean;
    s.sigma = s.num_std;
    s.method = 0;

    if (a.n < min_fit_entries || a.NEff() < min_fit_n_eff || !(s.num_std > 0))
        return s;

    // Fit a copy rescaled to mean weight 1, so the contents (and the fit's
    // amplitude) are at the scale of entry counts; the stored histogram
    // keeps its cross-section weights.
    TH1D h(*a.hist);
    h.Scale(a.n / a.sw);

    TF1 f("f_gaus", "gaus", 0.0, hist_max);
    f.SetParameters(h.GetMaximum(), s.num_mean, s.num_std);

    int status = h.Fit(&f, "RQN");

    s.fit_mean = f.GetParameter(1);
    s.fit_sigma = std::fabs(f.GetParameter(2));
    if (f.GetNDF() > 0)
        s.fit_chi2_ndf = f.GetChisquare() / f.GetNDF();

    bool ok = status == 0 &&
              std::isfinite(s.fit_mean) && std::isfinite(s.fit_sigma) &&
              s.fit_sigma > 0 && s.fit_mean >= 0.0 && s.fit_mean <= hist_max &&
              std::isfinite(s.fit_chi2_ndf) && s.fit_chi2_ndf <= max_fit_chi2_ndf;

    if (ok)
    {
        s.mean = s.fit_mean;
        s.sigma = s.fit_sigma;
        s.method = 1;
    }
    return s;
}

struct GridPoint
{
    int ix, iy;
    double x, y, eta, phi;
    double best_distance;   // best achievable jet distance (UESelection.h)
    double min_tag_distance; // best_distance - tolerance

    int n_events;           // raw entries passing veto + tag (both mirrors)
    Accum ue;               // E_UE, all selected events
    double sw_ue_zero;      // weight of selected events with E_UE = 0
    Accum contam;           // E_contam, selected events with E_contam > 0 only
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);

    if (argc < 2)
    {
        std::cerr << "Usage: ue_map <input_dir/> [output.root] [grid_spacing_cm]" << std::endl;
        return 1;
    }

    std::string input_dir = argv[1];
    if (input_dir.empty() || input_dir[input_dir.size() - 1] != '/')
        input_dir += "/";

    std::string out_name = "ue_map.root";
    if (argc > 2 && argv[2][0] != '\0')
        out_name = argv[2];
    if (argc > 3 && argv[3][0] != '\0')
        grid_spacing = std::atof(argv[3]);

    // ---------------- input files ----------------

    const std::string jet_prefix = "jet_output_";
    const std::string simple_prefix = "SimpleTree_";

    std::vector<std::string> jet_files, simple_files;
    {
        TSystemDirectory dir("input_dir", input_dir.c_str());
        TList *files = dir.GetListOfFiles();
        TIter next(files);
        TSystemFile *file;

        std::vector<std::string> names;
        while ((file = (TSystemFile *)next()))
        {
            std::string fname = file->GetName();
            if (!file->IsDirectory() && fname.compare(0, jet_prefix.size(), jet_prefix) == 0 &&
                fname.size() > 5 && fname.compare(fname.size() - 5, 5, ".root") == 0)
                names.push_back(fname);
        }
        std::sort(names.begin(), names.end());

        for (size_t i = 0; i < names.size(); i++)
        {
            jet_files.push_back(input_dir + names[i]);
            simple_files.push_back(input_dir + simple_prefix + names[i].substr(jet_prefix.size()));
        }
    }

    std::cout << "Found " << jet_files.size() << " JetTree files in " << input_dir << std::endl;

    // Per-file cross-section weight, same as the JES fits (-1 = skip file)
    std::vector<float> file_weight;
    if (!ComputeXSecWeights(jet_files, "jetTree", file_weight))
        return 1;

    // ---------------- grid ----------------

    const FiducialRect &rect = kFiducialRectEcal;
    const float fid_buffer = reco_fiducial_buffer_ecal;

    // Same grid definition as JetEnergyScaleFineGrid_em_only.cpp
    const double x_min = rect.x_inner + cut_buffer;
    const double x_max = rect.x_outer - cut_buffer;
    const double y_min = rect.y_min + cut_buffer;
    const double y_max = rect.y_max - cut_buffer;

    const int n_grid_x = (int)std::floor((x_max - x_min) / grid_spacing) + 1;
    const int n_grid_y = (int)std::floor((y_max - y_min) / grid_spacing) + 1;

    std::vector<SurfacePoint> candidates = FiducialJetAxisCandidates(candidate_spacing);

    std::vector<GridPoint> points;
    for (int ix = 0; ix < n_grid_x; ix++)
    for (int iy = 0; iy < n_grid_y; iy++)
    {
        GridPoint p;
        p.ix = ix;
        p.iy = iy;
        p.x = x_min + ix * grid_spacing;
        p.y = y_min + iy * grid_spacing;

        if (!pass_fiducial_cut((float)p.x, (float)p.y, fid_buffer, rect))
            continue;

        computeEtaPhi(p.x, p.y, p.eta, p.phi);
        p.best_distance = BestJetDistance(p.eta, p.phi, candidates);
        if (p.best_distance < 0)
            continue;
        p.min_tag_distance = p.best_distance - ue_tag_tolerance_in_R * R;

        p.n_events = 0;
        p.sw_ue_zero = 0.0;
        p.ue.hist = new TH1D(Form("h_ue_%d_%d", ix, iy),
                             Form("E_{UE} at (%.1f, %.1f) cm;E_{UE} (GeV);weighted events", p.x, p.y),
                             n_hist_bins, 0.0, hist_max);
        p.contam.hist = new TH1D(Form("h_contam_%d_%d", ix, iy),
                                 Form("Nonzero E_{contam} at (%.1f, %.1f) cm;E_{contam} (GeV);weighted events", p.x, p.y),
                                 n_hist_bins, 0.0, hist_max);
        p.ue.hist->Sumw2();
        p.contam.hist->Sumw2();

        points.push_back(p);
    }

    std::cout << "Grid: " << n_grid_x << " x " << n_grid_y << " points at " << grid_spacing
              << " cm spacing, " << points.size() << " inside the fiducial region" << std::endl;
    std::cout << "Selection: no reco jet within " << ue_min_dr_to_point / R << " R of the point; a reco jet within "
              << ue_tag_tolerance_in_R << " R of the best achievable distance; cone radius " << cone_radius << std::endl;

    // ---------------- event loop ----------------

    const int MAX_PTS = 200; // matches StSimpleReaderMaker.h's Cal_pt_* array size

    Long64_t n_events_read = 0, n_events_with_jet = 0;
    int n_pairs_used = 0, n_pairs_skipped = 0;

    for (size_t fidx = 0; fidx < jet_files.size(); fidx++)
    {
        if (file_weight[fidx] < 0)
            continue;
        const double w = file_weight[fidx];

        TFile *fj = TFile::Open(jet_files[fidx].c_str());
        TFile *fs = TFile::Open(simple_files[fidx].c_str());
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;

        if (!jt || !st || jt->GetEntries() != st->GetEntries() || !st->GetBranch("Cal_npts"))
        {
            std::cerr << "WARNING: skipping " << jet_files[fidx]
                      << " (missing/unreadable SimpleTree partner, no point branches, or entry counts differ: "
                      << (jt ? jt->GetEntries() : -1) << " vs " << (st ? st->GetEntries() : -1) << ")" << std::endl;
            n_pairs_skipped++;
            if (fj) { fj->Close(); delete fj; }
            if (fs) { fs->Close(); delete fs; }
            continue;
        }
        n_pairs_used++;

        std::vector<float> *reco_eta = 0, *reco_phi = 0;
        std::vector<float> *truth_const_px = 0, *truth_const_py = 0, *truth_const_pz = 0, *truth_const_E = 0;

        jt->SetBranchStatus("*", 0);
        const char *jet_branches[] = {"reco_eta", "reco_phi", "truth_const_px", "truth_const_py",
                                      "truth_const_pz", "truth_const_E"};
        for (int i = 0; i < 6; i++)
            jt->SetBranchStatus(jet_branches[i], 1);

        jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi);
        jt->SetBranchAddress("truth_const_px", &truth_const_px);
        jt->SetBranchAddress("truth_const_py", &truth_const_py);
        jt->SetBranchAddress("truth_const_pz", &truth_const_pz);
        jt->SetBranchAddress("truth_const_E", &truth_const_E);

        Int_t Cal_npts = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        Int_t Cal_pt_detid[MAX_PTS];

        st->SetBranchStatus("*", 0);
        const char *simple_branches[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y",
                                         "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++)
            st->SetBranchStatus(simple_branches[i], 1);

        st->SetBranchAddress("Cal_npts", &Cal_npts);
        st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x);
        st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z);
        st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);

        Long64_t n_entries = jt->GetEntries();

        std::vector<Deposit> jets, pts, truth;

        for (Long64_t ev = 0; ev < n_entries; ev++)
        {
            if (jt->GetEntry(ev) <= 0 || st->GetEntry(ev) <= 0)
                continue;
            n_events_read++;

            // The tag requirement needs at least one reco jet
            if (reco_eta->empty())
                continue;
            n_events_with_jet++;

            jets.clear();
            for (size_t j = 0; j < reco_eta->size(); j++)
            {
                Deposit d = {reco_eta->at(j), reco_phi->at(j), 0.0};
                jets.push_back(d);
            }

            // Points: same selection and same (x,y,z)-direction convention
            // as the reco jets' own inputs (JetMatcher.cpp, use_points)
            pts.clear();
            for (int i = 0; i < Cal_npts; i++)
            {
                if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1)
                    continue;
                if (!pass_point_energy_cut(Cal_pt_energy[i], kSqrtS_MainRun22_GeV))
                    continue;

                double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                Deposit d = {std::asinh(Cal_pt_z[i] / rho), std::atan2(Cal_pt_y[i], Cal_pt_x[i]), Cal_pt_energy[i]};
                pts.push_back(d);
            }

            // Truth-jet constituents (of the truth jets JetMatcher selected)
            truth.clear();
            for (size_t i = 0; i < truth_const_E->size(); i++)
            {
                double px = truth_const_px->at(i), py = truth_const_py->at(i), pz = truth_const_pz->at(i);
                double pt = std::sqrt(px * px + py * py);
                if (!(pt > 0))
                    continue;
                Deposit d = {std::asinh(pz / pt), std::atan2(py, px), truth_const_E->at(i)};
                truth.push_back(d);
            }

            for (size_t ip = 0; ip < points.size(); ip++)
            {
                GridPoint &p = points[ip];

                // mirror = 1 is the left-right mirrored copy of the event,
                // evaluated as the mirrored grid point in the original event
                for (int mirror = 0; mirror < 2; mirror++)
                {
                    const double p_phi = mirror ? M_PI - p.phi : p.phi;

                    bool vetoed = false, tagged = false;
                    for (size_t j = 0; j < jets.size(); j++)
                    {
                        double d_p = DeltaR(jets[j].eta, jets[j].phi, p.eta, p_phi);
                        if (d_p < ue_min_dr_to_point)
                        {
                            vetoed = true;
                            break;
                        }
                        if (std::min(d_p, DeltaR(jets[j].eta, jets[j].phi, p.eta, p_phi + M_PI)) >= p.min_tag_distance)
                            tagged = true;
                    }
                    if (vetoed || !tagged)
                        continue;

                    double e_ue = 0.0;
                    for (size_t i = 0; i < pts.size(); i++)
                        if (DeltaR(pts[i].eta, pts[i].phi, p.eta, p_phi) < cone_radius)
                            e_ue += pts[i].E;

                    double e_contam = 0.0;
                    for (size_t i = 0; i < truth.size(); i++)
                        if (DeltaR(truth[i].eta, truth[i].phi, p.eta, p_phi) < cone_radius)
                            e_contam += truth[i].E;

                    p.n_events++;
                    p.ue.Fill(e_ue, w);
                    if (e_ue <= 0.0)
                        p.sw_ue_zero += w;
                    if (e_contam > 0.0)
                        p.contam.Fill(e_contam, w);
                }
            }
        }

        fj->Close();
        delete fj;
        fs->Close();
        delete fs;
    }

    std::cout << "Read " << n_events_read << " events from " << n_pairs_used << " file pairs ("
              << n_pairs_skipped << " skipped); " << n_events_with_jet << " have at least one reco jet" << std::endl;

    // ---------------- output ----------------

    TFile fout(out_name.c_str(), "RECREATE");

    // Bins centred on the grid points
    const double hx_lo = x_min - 0.5 * grid_spacing, hx_hi = x_min + (n_grid_x - 0.5) * grid_spacing;
    const double hy_lo = y_min - 0.5 * grid_spacing, hy_hi = y_min + (n_grid_y - 0.5) * grid_spacing;

    const int n_grids = 12;
    const char *grid_names[n_grids] = {"ue_mean", "ue_sigma", "ue_method", "ue_zero_fraction",
                                       "contam_fraction", "contam_mean", "contam_sigma", "contam_method",
                                       "n_events", "n_eff", "best_distance_in_R", "contam_mean_all"};
    const char *grid_titles[n_grids] = {
        "UE mean energy;x (cm);y (cm);GeV",
        "UE energy sigma;x (cm);y (cm);GeV",
        "UE method (1 = Gaussian fit, 0 = numerical);x (cm);y (cm)",
        "Fraction of selected events with E_{UE} = 0;x (cm);y (cm)",
        "Fraction of selected events with truth-jet contamination;x (cm);y (cm)",
        "Mean contamination energy, contaminated events only;x (cm);y (cm);GeV",
        "Contamination energy sigma, contaminated events only;x (cm);y (cm);GeV",
        "Contamination method (1 = Gaussian fit, 0 = numerical);x (cm);y (cm)",
        "Selected events (raw entries, both mirrors);x (cm);y (cm)",
        "Selected events (effective, cross-section weighted);x (cm);y (cm)",
        "Best achievable jet distance / R;x (cm);y (cm)",
        "Mean contamination energy over all selected events;x (cm);y (cm);GeV"};

    TH2D *grid[n_grids];
    for (int g = 0; g < n_grids; g++)
        grid[g] = new TH2D(grid_names[g], grid_titles[g], n_grid_x, hx_lo, hx_hi, n_grid_y, hy_lo, hy_hi);

    // Which cells of each grid actually hold a value: a stored 0 can be a
    // real result (a method flag, a zero fraction), so the plots below need
    // this to tell it apart from "nothing here".
    std::vector<std::vector<char> > grid_set(n_grids, std::vector<char>((n_grid_x + 2) * (n_grid_y + 2), 0));
    auto set_cell = [&](int g, int bx, int by, double value)
    {
        grid[g]->SetBinContent(bx, by, value);
        grid_set[g][bx * (n_grid_y + 2) + by] = 1;
    };

    TTree t("points", "UE map, one row per grid point");

    double b_x, b_y, b_eta, b_phi, b_best_distance;
    int b_n_events, b_ue_method, b_contam_method, b_contam_n;
    double b_n_eff, b_ue_zero_fraction;
    double b_ue_mean, b_ue_sigma, b_ue_num_mean, b_ue_num_std, b_ue_fit_mean, b_ue_fit_sigma, b_ue_fit_chi2_ndf;
    double b_contam_fraction, b_contam_n_eff, b_contam_mean_all;
    double b_contam_mean, b_contam_sigma, b_contam_num_mean, b_contam_num_std;
    double b_contam_fit_mean, b_contam_fit_sigma, b_contam_fit_chi2_ndf;

    t.Branch("x", &b_x);
    t.Branch("y", &b_y);
    t.Branch("eta", &b_eta);
    t.Branch("phi", &b_phi);
    t.Branch("best_distance", &b_best_distance);
    t.Branch("n_events", &b_n_events);
    t.Branch("n_eff", &b_n_eff);
    t.Branch("ue_mean", &b_ue_mean);
    t.Branch("ue_sigma", &b_ue_sigma);
    t.Branch("ue_method", &b_ue_method);
    t.Branch("ue_zero_fraction", &b_ue_zero_fraction);
    t.Branch("ue_num_mean", &b_ue_num_mean);
    t.Branch("ue_num_std", &b_ue_num_std);
    t.Branch("ue_fit_mean", &b_ue_fit_mean);
    t.Branch("ue_fit_sigma", &b_ue_fit_sigma);
    t.Branch("ue_fit_chi2_ndf", &b_ue_fit_chi2_ndf);
    t.Branch("contam_fraction", &b_contam_fraction);
    t.Branch("contam_n", &b_contam_n);
    t.Branch("contam_n_eff", &b_contam_n_eff);
    t.Branch("contam_mean", &b_contam_mean);
    t.Branch("contam_sigma", &b_contam_sigma);
    t.Branch("contam_method", &b_contam_method);
    t.Branch("contam_num_mean", &b_contam_num_mean);
    t.Branch("contam_num_std", &b_contam_num_std);
    t.Branch("contam_fit_mean", &b_contam_fit_mean);
    t.Branch("contam_fit_sigma", &b_contam_fit_sigma);
    t.Branch("contam_fit_chi2_ndf", &b_contam_fit_chi2_ndf);
    t.Branch("contam_mean_all", &b_contam_mean_all);

    TDirectory *hist_dir = fout.mkdir("hists");

    int n_filled = 0, n_empty = 0;
    int n_ue_fit = 0, n_ue_num = 0, n_contam_fit = 0, n_contam_num = 0, n_contam_none = 0;

    for (size_t ip = 0; ip < points.size(); ip++)
    {
        GridPoint &p = points[ip];
        const int bx = p.ix + 1, by = p.iy + 1;

        Summary ue = Summarize(p.ue);
        Summary ct = Summarize(p.contam);

        b_x = p.x;
        b_y = p.y;
        b_eta = p.eta;
        b_phi = p.phi;
        b_best_distance = p.best_distance;
        b_n_events = p.n_events;
        b_n_eff = p.ue.NEff();

        b_ue_mean = ue.mean;
        b_ue_sigma = ue.sigma;
        b_ue_method = ue.method;
        b_ue_zero_fraction = p.ue.sw > 0 ? p.sw_ue_zero / p.ue.sw : NAN;
        b_ue_num_mean = ue.num_mean;
        b_ue_num_std = ue.num_std;
        b_ue_fit_mean = ue.fit_mean;
        b_ue_fit_sigma = ue.fit_sigma;
        b_ue_fit_chi2_ndf = ue.fit_chi2_ndf;

        // contam holds the contaminated events only; its weight relative to
        // all selected events is the contaminated fraction
        b_contam_fraction = p.ue.sw > 0 ? p.contam.sw / p.ue.sw : NAN;
        b_contam_n = p.contam.n;
        b_contam_n_eff = p.contam.NEff();
        b_contam_mean = ct.mean;
        b_contam_sigma = ct.sigma;
        b_contam_method = ct.method;
        b_contam_num_mean = ct.num_mean;
        b_contam_num_std = ct.num_std;
        b_contam_fit_mean = ct.fit_mean;
        b_contam_fit_sigma = ct.fit_sigma;
        b_contam_fit_chi2_ndf = ct.fit_chi2_ndf;
        b_contam_mean_all = p.ue.sw > 0 ? p.contam.swx / p.ue.sw : NAN;

        t.Fill();

        set_cell(10, bx, by, p.best_distance / R);

        if (p.n_events == 0)
        {
            n_empty++;
            continue;
        }
        n_filled++;

        set_cell(0, bx, by, ue.mean);
        set_cell(1, bx, by, ue.sigma);
        set_cell(2, bx, by, ue.method);
        set_cell(3, bx, by, b_ue_zero_fraction);
        set_cell(4, bx, by, b_contam_fraction);
        set_cell(8, bx, by, p.n_events);
        set_cell(9, bx, by, b_n_eff);
        set_cell(11, bx, by, b_contam_mean_all);

        (ue.method == 1 ? n_ue_fit : n_ue_num)++;

        if (ct.method < 0)
            n_contam_none++;
        else
        {
            set_cell(5, bx, by, ct.mean);
            set_cell(6, bx, by, ct.sigma);
            set_cell(7, bx, by, ct.method);
            (ct.method == 1 ? n_contam_fit : n_contam_num)++;
        }

        hist_dir->cd();
        p.ue.hist->Write();
        p.contam.hist->Write();
    }

    fout.cd();
    t.Write();
    for (int g = 0; g < n_grids; g++)
        grid[g]->Write();

    std::cout << "Grid points: " << points.size() << " in the fiducial region, " << n_filled
              << " with selected events, " << n_empty << " with none" << std::endl;
    std::cout << "UE mean/sigma:            Gaussian fit at " << n_ue_fit << " points, numerical at "
              << n_ue_num << std::endl;
    std::cout << "Contamination mean/sigma: Gaussian fit at " << n_contam_fit << " points, numerical at "
              << n_contam_num << ", no contaminated events at " << n_contam_none << std::endl;
    if (n_ue_num > 0 || n_contam_num > 0)
        std::cout << "FLAG: numerical mean/std used where the Gaussian fit was not trusted (fewer than "
                  << min_fit_entries << " entries or " << min_fit_n_eff << " effective entries, a failed fit, or chi2/ndf > "
                  << max_fit_chi2_ndf << ") -- see the *_method grids" << std::endl;

    // ---------------- plots ----------------

    gStyle->SetOptStat(0);
    gStyle->SetPalette(53);
    gStyle->SetNumberContours(100);

    std::string pdf_name = out_name;
    size_t dot = pdf_name.rfind(".root");
    if (dot != std::string::npos)
        pdf_name.erase(dot);
    pdf_name += ".pdf";

    TCanvas c("c", "", 1800, 1350);
    c.Divide(4, 3);
    for (int g = 0; g < n_grids; g++)
    {
        c.cd(g + 1);
        gPad->SetRightMargin(0.16);
        gPad->SetLeftMargin(0.12);

        // Draw a copy in which cells without a value sit far below the colour
        // scale (and so are left blank), with the scale starting just below
        // the smallest real value so that one is not blank too.
        double lo = 1e300, hi = -1e300;
        for (int bx = 1; bx <= n_grid_x; bx++)
        for (int by = 1; by <= n_grid_y; by++)
        {
            if (!grid_set[g][bx * (n_grid_y + 2) + by])
                continue;
            lo = std::min(lo, grid[g]->GetBinContent(bx, by));
            hi = std::max(hi, grid[g]->GetBinContent(bx, by));
        }
        if (lo > hi)
        {
            lo = 0.0;
            hi = 1.0;
        }
        double span = hi > lo ? hi - lo : 1.0;

        TH2D *shown = (TH2D *)grid[g]->Clone(Form("%s_shown", grid_names[g]));
        for (int bx = 1; bx <= n_grid_x; bx++)
        for (int by = 1; by <= n_grid_y; by++)
            if (!grid_set[g][bx * (n_grid_y + 2) + by])
                shown->SetBinContent(bx, by, lo - 100.0 * span);

        shown->SetMinimum(lo - 0.02 * span);
        shown->SetMaximum(hi > lo ? hi : lo + span);
        shown->Draw("colz");
    }
    c.SaveAs(pdf_name.c_str());

    std::string png_name = pdf_name.substr(0, pdf_name.size() - 4) + ".png";
    c.SaveAs(png_name.c_str());

    fout.Close();

    std::cout << "Output written to " << out_name << ", " << pdf_name << std::endl;
    return 0;
}
