// Diagnostics-only, standalone tool: for a small hardcoded list of target
// (x_center, y_center, E_low, E_high) bins (chosen externally -- see
// query used to pick them), rescans a JetTrees folder collecting just
// those bins' raw response/scale values (same match_radius/reflection
// convention as JetEnergyScaleFineGrid.cpp), reproduces gaussian_fit()'s
// exact windowing/fit logic, and saves a histogram + fit-window plot per
// bin. Lets us see directly whether the [mean-1.5sigma, mean+1.5sigma]
// window is well-placed on the true peak or clipping/missing it,
// especially at low E where JetEnergyScaleFineGrid.cpp's per-bin fits
// were found to converge far less often (85.6% vs 99.3% for ECAL+HCAL)
// and sometimes land on physically nonsensical values (e.g. negative
// response_mu).
//
// Usage: ./plot_gaussian_peaks <jettrees_dir> <output_dir>
// (target bin list is hardcoded in main() -- see TARGETS below)

#include <TFile.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include <algorithm>

static const bool reflect_to_positive_half = true;
static const double match_radius = 5.0;

static const int n_hist_bins = 50;
static const int min_total_entries = 5;
static const double n_sigma = 1.5;
static const double asym = 1.0;

// Below this truth_EMF, "reco_E / (truth_E*truth_EMF)" divides by
// near-zero and blows up -- these jets have essentially no true EM content
// to calibrate an EM-calorimeter response against anyway, so they're
// dropped from the EMF-corrected scale (not from the plain reco/truth one).
static const double min_truth_emf_for_emf_scale = 0.05;

struct Target
{
    std::string label;
    double x_center, y_center, E_low, E_high;
};

struct RawBin
{
    std::vector<double> resp_vals;
    std::vector<double> scale_vals;
    std::vector<double> scale_emf_vals; // reco_E / (truth_E * truth_EMF) -- see min_truth_emf_for_emf_scale
    std::vector<int> nconst_vals; // unfiltered -- every matched jet in this bin, for the multiplicity histogram
    long n_single_em = 0, n_not_single_em = 0; // unfiltered counts, regardless of truth_select
    long n_emf_too_small = 0; // dropped from scale_emf_vals due to min_truth_emf_for_emf_scale
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <jettrees_dir> <output_dir> [min_nconst] [truth_select]" << std::endl;
        std::cerr << "  truth_select: 0=all truth jets (default), 1=exclude single-constituent EMF=1 truth jets"
                      " (likely bare electron/photon, not a real hadronic jet), 2=ONLY those (isolate the suspected contamination)" << std::endl;
        return 1;
    }

    std::string input_dir = argv[1];
    if (!input_dir.empty() && input_dir.back() != '/') input_dir += "/";
    std::string out_dir = argv[2];
    int min_nconst = (argc > 3) ? std::atoi(argv[3]) : 0;
    int truth_select = (argc > 4) ? std::atoi(argv[4]) : 0;
    std::cout << "min_nconst cut: " << min_nconst << (min_nconst > 0 ? "" : " (no cut)") << std::endl;
    std::cout << "truth_select mode: " << truth_select << std::endl;
    { std::string cmd = "mkdir -p " + out_dir; system(cmd.c_str()); }

    // Chosen from the 8 highest-n_entries bins at the lowest energy bin
    // (E_low=10) in em_only_jet_calibration_diag.root: the one with a
    // plausible fitted scale (x=28.0,y=0.78) plus three of the ones that
    // fit to physically nonsensical negative scale_mu/response_mu.
    std::vector<Target> TARGETS = {
        {"x28_y0.8_sane",     28.036, 0.775612, 10.0, 12.1243},
        {"x28_y5.8_negative", 28.036, 5.77561,  10.0, 12.1243},
        {"x28_y-4.2_negative",28.036, -4.22439, 10.0, 12.1243},
        {"x25_y0.8_negative", 25.536, 0.775612, 10.0, 12.1243},
    };

    std::vector<RawBin> accum(TARGETS.size());

    // 2D truth_E vs reco_E diagnostic, filled for ALL truth energies (not
    // just each target's narrow E_low/E_high slice) at each target's (x,y)
    // position -- one pair of histograms (all jets / excluding single-EM
    // truth jets) per target, to see the two populations across the full
    // energy range at once. h2d_emf uses (truth_E*truth_EMF) on the x-axis
    // instead of raw truth_E -- see the EMF-scale motivation in the
    // conversation: this should recover a diagonal ridge if reco_E is
    // really tracking the jet's true EM energy, not its full energy.
    std::vector<TH2D> h2d_all, h2d_excl, h2d_emf;
    for (size_t k = 0; k < TARGETS.size(); k++)
    {
        h2d_all.emplace_back(Form("h2d_all_%zu", k), "", 100, 0, 100, 100, 0, 100);
        h2d_excl.emplace_back(Form("h2d_excl_%zu", k), "", 100, 0, 100, 100, 0, 100);
        h2d_emf.emplace_back(Form("h2d_emf_%zu", k), "", 100, 0, 100, 100, 0, 100);
    }

    // --- Discover input files ---
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
    std::cout << "Found " << root_files.size() << " ROOT files" << std::endl;

    // --- Scan, collecting only the target bins' raw values ---
    for (size_t fidx = 0; fidx < root_files.size(); fidx++)
    {
        TFile *f = TFile::Open(root_files[fidx].c_str());
        if (!f || f->IsZombie()) { if (f) delete f; continue; }
        TTree *tree = (TTree *)f->Get("jetTree");
        if (!tree) { f->Close(); delete f; continue; }

        std::vector<float> *reco_E = nullptr, *reco_x = nullptr, *reco_y = nullptr, *truth_E = nullptr, *truth_EMF = nullptr;
        std::vector<int> *match_truth_idx = nullptr, *match_reco_idx = nullptr, *reco_nconst = nullptr, *truth_nconst = nullptr;
        tree->SetBranchAddress("reco_E", &reco_E);
        tree->SetBranchAddress("reco_x", &reco_x);
        tree->SetBranchAddress("reco_y", &reco_y);
        tree->SetBranchAddress("truth_E", &truth_E);
        tree->SetBranchAddress("truth_EMF", &truth_EMF);
        tree->SetBranchAddress("match_truth_idx", &match_truth_idx);
        tree->SetBranchAddress("match_reco_idx", &match_reco_idx);
        tree->SetBranchAddress("reco_nconst", &reco_nconst);
        tree->SetBranchAddress("truth_nconst", &truth_nconst);

        Long64_t n_events = tree->GetEntries();
        for (Long64_t ev = 0; ev < n_events; ev++)
        {
            tree->GetEntry(ev);
            for (size_t m = 0; m < match_reco_idx->size(); m++)
            {
                int ri = match_reco_idx->at(m);
                int ti = match_truth_idx->at(m);
                double truth = truth_E->at(ti);
                double reco = reco_E->at(ri);
                double emf = truth_EMF->at(ti);
                double truth_em = truth * emf; // true EM-only energy in the jet (photons/e-/e+)
                int nconst = reco_nconst->at(ri);
                // Bare electron/photon suspect: truth jet built from a single
                // constituent that's 100% EM (see BumpReason-style reasoning
                // in the conversation -- a real hadronic jet at this energy
                // should rarely be a single pure-EM particle).
                bool is_single_em = (truth_nconst->at(ti) == 1) && (emf >= 0.999f);
                double x = reflect_to_positive_half ? std::fabs(reco_x->at(ri)) : reco_x->at(ri);
                double y = reco_y->at(ri);

                for (size_t k = 0; k < TARGETS.size(); k++)
                {
                    const Target &t = TARGETS[k];
                    if (truth < t.E_low || truth >= t.E_high) continue;
                    double dx = x - t.x_center, dy = y - t.y_center;
                    if (dx * dx + dy * dy > match_radius * match_radius) continue;
                    accum[k].nconst_vals.push_back(nconst);
                    if (is_single_em) accum[k].n_single_em++; else accum[k].n_not_single_em++;
                    if (nconst < min_nconst) continue;
                    if (truth_select == 1 && is_single_em) continue;
                    if (truth_select == 2 && !is_single_em) continue;
                    accum[k].resp_vals.push_back(reco);
                    accum[k].scale_vals.push_back(reco / truth);
                    if (emf < min_truth_emf_for_emf_scale) { accum[k].n_emf_too_small++; continue; }
                    accum[k].scale_emf_vals.push_back(reco / truth_em);
                }

                // 2D diagnostic: same (x,y) match, but no truth_E restriction --
                // see h2d_all/h2d_excl/h2d_emf declaration above.
                for (size_t k = 0; k < TARGETS.size(); k++)
                {
                    const Target &t = TARGETS[k];
                    double dx = x - t.x_center, dy = y - t.y_center;
                    if (dx * dx + dy * dy > match_radius * match_radius) continue;
                    h2d_all[k].Fill(truth, reco);
                    if (!is_single_em) h2d_excl[k].Fill(truth, reco);
                    if (emf >= min_truth_emf_for_emf_scale) h2d_emf[k].Fill(truth_em, reco);
                }
            }
        }
        f->Close();
        delete f;
    }

    // --- Save the 2D truth_E vs reco_E diagnostics ---
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1); // rainbow -- kBird isn't available in ROOT 5.34
    for (size_t k = 0; k < TARGETS.size(); k++)
    {
        const Target &t = TARGETS[k];
        for (int pass = 0; pass < 3; pass++)
        {
            TH2D &h = pass == 0 ? h2d_all[k] : (pass == 1 ? h2d_excl[k] : h2d_emf[k]);
            std::string tag = pass == 0 ? "all_jets" : (pass == 1 ? "excl_single_EM" : "emf_corrected");
            std::string xlabel = pass == 2 ? "truth_E*truth_EMF (GeV)" : "truth_E (GeV)";
            if (h.GetEntries() < 1) continue;

            TCanvas c2("c2", "", 900, 800);
            c2.SetLogz();
            h.SetTitle(Form("%s (x=%.1f y=%.1f), %s;%s;reco_E (GeV)",
                             t.label.c_str(), t.x_center, t.y_center, tag.c_str(), xlabel.c_str()));
            h.Draw("COLZ");

            TLine diag(0, 0, 100, 100);
            diag.SetLineColor(kRed);
            diag.SetLineStyle(2);
            diag.SetLineWidth(2);
            diag.Draw("SAME");

            std::string outpath = out_dir + "/peak2d_" + t.label + "_" + tag + ".png";
            c2.Print(outpath.c_str());
            std::cout << "2D hist (" << tag << ") for " << t.label << ": " << h.GetEntries()
                      << " entries -> " << outpath << std::endl;
        }
    }

    for (size_t k = 0; k < TARGETS.size(); k++)
    {
        const std::vector<int> &nc = accum[k].nconst_vals;
        long n1 = std::count(nc.begin(), nc.end(), 1);
        long n2 = std::count(nc.begin(), nc.end(), 2);
        double nc_mean = nc.empty() ? NAN : std::accumulate(nc.begin(), nc.end(), 0.0) / nc.size();
        long n_em_tot = accum[k].n_single_em + accum[k].n_not_single_em;
        std::cout << "Target " << TARGETS[k].label << ": " << nc.size() << " matched jets total ("
                  << accum[k].resp_vals.size() << " pass cuts: min_nconst=" << min_nconst << " truth_select=" << truth_select << ")"
                  << " | reco nconst: mean=" << nc_mean << " n==1: " << n1 << " (" << (100.0 * n1 / std::max((size_t)1, nc.size())) << "%)"
                  << " n==2: " << n2 << " (" << (100.0 * n2 / std::max((size_t)1, nc.size())) << "%)"
                  << " | single-EM truth jets: " << accum[k].n_single_em << " (" << (100.0 * accum[k].n_single_em / std::max(1L, n_em_tot)) << "%)"
                  << " | scale_emf_vals: " << accum[k].scale_emf_vals.size() << " (dropped " << accum[k].n_emf_too_small
                  << " with truth_EMF<" << min_truth_emf_for_emf_scale << ")" << std::endl;
    }

    // --- Reproduce gaussian_fit()'s exact logic, but keep the histogram/fit for plotting ---
    auto plot_one = [&](const std::string &varname, const std::vector<double> &v, const std::string &outpath, const std::string &title)
    {
        if (v.size() < (size_t)min_total_entries) { std::cout << "  (too few entries to fit " << varname << ")" << std::endl; return; }

        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double var = 0.0;
        for (double x : v) var += (x - mean) * (x - mean);
        double sigma = std::sqrt(var / std::max<size_t>(1, v.size() - 1));

        std::cout << "  " << varname << ": raw mean=" << mean << " raw sigma=" << sigma << std::endl;
        if (!std::isfinite(sigma) || sigma <= 0.0) { std::cout << "    (sigma non-finite/non-positive, skipping)" << std::endl; return; }

        TH1D h("h", "", n_hist_bins, mean - 5 * sigma, mean + 5 * sigma);
        for (double x : v) h.Fill(x);

        double xmin = mean - n_sigma * sigma;
        double xmax = mean + n_sigma * sigma * asym;

        TF1 fitf("f", "gaus", xmin, xmax);
        fitf.SetParameters(h.GetMaximum(), mean, sigma);
        h.Fit(&fitf, "RQ");

        double fit_mu = fitf.GetParameter(1), fit_sigma = fitf.GetParameter(2);
        double chi2_ndf = fitf.GetNDF() > 0 ? fitf.GetChisquare() / fitf.GetNDF() : NAN;
        std::cout << "    fit mu=" << fit_mu << " fit sigma=" << fit_sigma
                  << " chi2/ndf=" << chi2_ndf << " window=[" << xmin << "," << xmax << "]" << std::endl;

        TCanvas c("c", "", 900, 650);
        h.SetTitle(Form("%s;%s;entries", title.c_str(), varname.c_str()));
        h.SetLineColor(kBlack);
        h.Draw("HIST");

        TF1 fdraw("fdraw", "gaus", h.GetXaxis()->GetXmin(), h.GetXaxis()->GetXmax());
        fdraw.SetParameters(fitf.GetParameter(0), fit_mu, fit_sigma);
        fdraw.SetLineColor(kBlue);
        fdraw.SetLineWidth(2);
        fdraw.Draw("SAME");

        double ymax = h.GetMaximum() * 1.15;
        TLine lmin(xmin, 0, xmin, ymax), lmax(xmax, 0, xmax, ymax);
        lmin.SetLineColor(kRed); lmin.SetLineStyle(2); lmin.SetLineWidth(2);
        lmax.SetLineColor(kRed); lmax.SetLineStyle(2); lmax.SetLineWidth(2);
        lmin.Draw("SAME");
        lmax.Draw("SAME");

        TLegend leg(0.58, 0.65, 0.88, 0.88);
        leg.AddEntry(&h, "raw values", "l");
        leg.AddEntry(&fdraw, Form("windowed fit (#mu=%.3g, #sigma=%.3g)", fit_mu, fit_sigma), "l");
        leg.AddEntry(&lmin, Form("fit window (raw mean %.3g #pm %.1f#sigma_{raw}=%.3g)", mean, n_sigma, sigma), "l");
        leg.Draw();

        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.03);
        lat.DrawLatex(0.58, 0.60, Form("raw mean=%.3g  raw #sigma=%.3g", mean, sigma));
        lat.DrawLatex(0.58, 0.56, Form("chi2/ndf=%.3g  n=%d", chi2_ndf, (int)v.size()));

        c.Print(outpath.c_str());
        std::cout << "    -> " << outpath << std::endl;
    };

    for (size_t k = 0; k < TARGETS.size(); k++)
    {
        const Target &t = TARGETS[k];
        std::cout << "=== " << t.label << " (x=" << t.x_center << " y=" << t.y_center
                  << " E=[" << t.E_low << "," << t.E_high << ")) ===" << std::endl;

        std::string title = Form("%s (min_nconst=%d, truth_select=%d): x=%.1f y=%.1f E=[%.1f,%.1f)",
                                  t.label.c_str(), min_nconst, truth_select, t.x_center, t.y_center, t.E_low, t.E_high);
        std::string suffix = Form("_nconst%d_truthsel%d", min_nconst, truth_select);

        plot_one("response_mu (GeV)", accum[k].resp_vals,
                  out_dir + "/peak_" + t.label + suffix + "_response.png", title);
        plot_one("scale_mu = reco_E/truth_E", accum[k].scale_vals,
                  out_dir + "/peak_" + t.label + suffix + "_scale.png", title);
        plot_one("scale_mu = reco_E/(truth_E*truth_EMF)", accum[k].scale_emf_vals,
                  out_dir + "/peak_" + t.label + suffix + "_scale_emf.png", title);
    }

    return 0;
}
