// Does the EM-only JES model (CalibrationMap_em_only.cpp: A + B/E + C/E^2 plus
// an asymmetric skew Gaussian) still describe the response now that reco jets
// are built from FCS points and the underlying event (UE) is subtracted first?
//
// Reads sim_to_jet JetTrees that carry the UE branches of
// underlying_event/UnderlyingEvent.cpp, and for matched reco/truth jets forms
// the median scale in truth-energy bins (the same bins and the same weighted
// median as JetEnergyScaleFineGrid_em_only.cpp) for three selections:
//   raw                    reco_E / truth_E, every matched jet
//   reco UE subtracted     (reco_E - reco_ue_E) / truth_E, every matched jet
//                          with a reco UE estimate
//   UE subtracted, both    (reco_E - reco_ue_E) / (truth_E - truth_ue_E),
//                          every matched jet with both estimates
// The bins are always in the truth energy of the ratio's denominator.
// in four position regions: everywhere, and three bands of the jet's distance
// from the beam with equal numbers of jets. Each curve is fitted with the
// background alone and with the full model, as a function of the median
// reco energy of the bin (the "response", which the correction is looked up
// with). The sample is far too small for the 2.5 cm position grid of the real
// calibration; this only tests the shape.
//
// reco jets are only kept above reco_jet_pt_min_em, so a truth-energy bin near
// that threshold holds only jets that fluctuated upwards and its median scale
// is biased high. truth_pt_min (GeV) keeps only matched jets whose TRUTH pT is
// above it, which removes those bins without biasing the rest.
//
// Usage: ./points_ue_form_check <jettrees_dir/> [output_base] [unweighted|weighted] [truth_pt_min]
// Build: g++ -m64 -O2 points_ue_form_check.cpp ../shared/VectorDict.cxx -o points_ue_form_check
//            -I../shared -I../jet_scale `root-config --cflags --libs` -std=c++11

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"

// same truth-energy binning as JetEnergyScaleFineGrid_em_only.cpp
const double E0_bins = 3, Emax_bins = 100, energy_curvature = 1.5;
const int n_bins_E = 20;

// same parameter domain as CalibrationMap_em_only.cpp
const double D_MIN = 0.0, D_MAX = 0.5;
const double E0_MIN = 1.0, E0_MAX = 40.0;
const double SIGMA_MIN = 1.0, SIGMA_MAX = 35.0;

const double min_n_eff = 10.0;       // a bin needs this effective count to give a point
const int n_region = 4, n_sel = 3;

struct Jet
{
    float reco, sub, truth, truth_sub, r, w;
    bool has_ue, has_truth_ue;
};

double median_of(const std::vector<double> &v, const std::vector<double> &w)
{
    std::vector<size_t> order(v.size());
    for (size_t i = 0; i < order.size(); i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&v](size_t a, size_t b) { return v[a] < v[b]; });
    double half = 0.5 * std::accumulate(w.begin(), w.end(), 0.0);
    double cum = 0.0;
    for (size_t i = 0; i < order.size(); i++)
    {
        cum += w[order[i]];
        if (cum > half) return v[order[i]];
        if (cum == half && i + 1 < order.size()) return 0.5 * (v[order[i]] + v[order[i + 1]]);
    }
    return v[order.back()];
}

// median and its approximate standard error (as median_fit in
// JetEnergyScaleFineGrid_em_only.cpp); false if the bin is unusable
bool median_with_error(const std::vector<double> &v, const std::vector<double> &w, double &med, double &err)
{
    if (v.size() < 5 || EffectiveEntries(w) < min_n_eff) return false;
    med = median_of(v, w);
    std::vector<double> absdev(v.size());
    for (size_t i = 0; i < v.size(); i++) absdev[i] = std::fabs(v[i] - med);
    double mad_sigma = 1.4826 * median_of(absdev, w);
    if (!std::isfinite(mad_sigma) || mad_sigma <= 0.0) return false;
    err = 1.2533 * mad_sigma / std::sqrt(EffectiveEntries(w));
    return true;
}

double bg_fn(double *x, double *p)
{
    return p[0] + p[1] / x[0] + p[2] / (x[0] * x[0]);
}

double full_fn(double *x, double *p)
{
    double sigma = (x[0] < p[4]) ? p[5] : p[6];
    return bg_fn(x, p) + p[3] * std::exp(-(x[0] - p[4]) * (x[0] - p[4]) / (2 * sigma * sigma));
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    if (argc < 2)
    {
        std::cerr << "usage: ./points_ue_form_check <jettrees_dir/> [output_base] [unweighted]" << std::endl;
        return 1;
    }
    std::string input_dir = argv[1];
    std::string out_base = argc > 2 ? argv[2] : "points_ue_form_check";
    bool unweighted = argc > 3 && std::string(argv[3]) == "unweighted";
    double truth_pt_min = argc > 4 ? std::atof(argv[4]) : 0.0;

    std::vector<std::string> root_files;
    TSystemDirectory dir("input_dir", input_dir.c_str());
    TIter next(dir.GetListOfFiles());
    while (TSystemFile *file = (TSystemFile *)next())
    {
        std::string fname = file->GetName();
        if (!file->IsDirectory() && fname.find("jet_output") == 0 && fname.find(".root") != std::string::npos)
            root_files.push_back(input_dir + fname);
    }
    std::sort(root_files.begin(), root_files.end());

    // a job that ended early leaves sigma_job/N_job unset; such a file cannot
    // be weighted (ComputeXSecWeights stops on it), so it is left out
    std::vector<std::string> usable;
    for (size_t i = 0; i < root_files.size(); i++)
    {
        TFile *f = TFile::Open(root_files[i].c_str());
        TTree *tree = f && !f->IsZombie() ? (TTree *)f->Get("jetTree") : 0;
        bool ok = tree && tree->GetEntries() > 0 && tree->GetLeaf("N_job");
        if (ok) { tree->GetEntry(0); ok = tree->GetLeaf("N_job")->GetValue() > 0; }
        if (ok) usable.push_back(root_files[i]);
        if (f) { f->Close(); delete f; }
    }
    std::cout << root_files.size() - usable.size() << " of " << root_files.size()
              << " files left out (empty, or no cross section stored)" << std::endl;
    root_files = usable;

    std::vector<float> file_weight;
    if (!ComputeXSecWeights(root_files, "jetTree", file_weight)) return 1;

    std::vector<Jet> jets;
    for (size_t fidx = 0; fidx < root_files.size(); fidx++)
    {
        if (file_weight[fidx] < 0) continue;
        TFile *f = TFile::Open(root_files[fidx].c_str());
        if (!f || f->IsZombie()) continue;
        TTree *tree = (TTree *)f->Get("jetTree");
        if (!tree || !tree->GetBranch("reco_ue_E") || !tree->GetBranch("truth_ue_E")) { f->Close(); delete f; continue; }

        std::vector<float> *reco_E = nullptr, *reco_x = nullptr, *reco_y = nullptr, *truth_E = nullptr, *truth_pt = nullptr, *reco_ue_E = nullptr, *truth_ue_E = nullptr;
        std::vector<int> *match_truth_idx = nullptr, *match_reco_idx = nullptr, *reco_ue_nproxy = nullptr, *truth_ue_nproxy = nullptr;
        tree->SetBranchAddress("reco_E", &reco_E);
        tree->SetBranchAddress("reco_x", &reco_x);
        tree->SetBranchAddress("reco_y", &reco_y);
        tree->SetBranchAddress("truth_E", &truth_E);
        tree->SetBranchAddress("truth_pt", &truth_pt);
        tree->SetBranchAddress("match_truth_idx", &match_truth_idx);
        tree->SetBranchAddress("match_reco_idx", &match_reco_idx);
        tree->SetBranchAddress("reco_ue_E", &reco_ue_E);
        tree->SetBranchAddress("reco_ue_nproxy", &reco_ue_nproxy);
        tree->SetBranchAddress("truth_ue_E", &truth_ue_E);
        tree->SetBranchAddress("truth_ue_nproxy", &truth_ue_nproxy);
        tree->SetBranchStatus("*", 0);
        const char *used[] = {"reco_E", "reco_x", "reco_y", "truth_E", "truth_pt", "match_truth_idx", "match_reco_idx",
                              "reco_ue_E", "reco_ue_nproxy", "truth_ue_E", "truth_ue_nproxy"};
        for (int k = 0; k < 11; k++) tree->SetBranchStatus(used[k], 1);

        Long64_t n_events = tree->GetEntries();
        for (Long64_t ev = 0; ev < n_events; ev++)
        {
            tree->GetEntry(ev);
            for (size_t m = 0; m < match_reco_idx->size(); m++)
            {
                int ri = match_reco_idx->at(m), ti = match_truth_idx->at(m);
                float x = reco_x->at(ri), y = reco_y->at(ri);
                if (!pass_fiducial_cut(x, y, reco_fiducial_buffer_ecal, kFiducialRectEcal)) continue;
                if (truth_pt->at(ti) < truth_pt_min) continue;
                Jet j;
                j.reco = reco_E->at(ri);
                j.truth = truth_E->at(ti);
                j.r = std::sqrt(x * x + y * y);
                j.w = unweighted ? 1.0f : file_weight[fidx];
                j.has_ue = reco_ue_nproxy->at(ri) > 0;
                j.has_truth_ue = truth_ue_nproxy->at(ti) > 0;
                j.truth_sub = j.has_truth_ue ? j.truth - truth_ue_E->at(ti) : 0.0f;
                j.sub = j.has_ue ? j.reco - reco_ue_E->at(ri) : 0.0f;
                jets.push_back(j);
            }
        }
        f->Close();
        delete f;
    }

    // region boundaries: thirds of the jets that have a UE estimate
    std::vector<float> radii;
    for (size_t i = 0; i < jets.size(); i++) if (jets[i].has_ue) radii.push_back(jets[i].r);
    std::sort(radii.begin(), radii.end());
    if (radii.size() < 100) { std::cerr << "too few jets with a UE estimate" << std::endl; return 1; }
    double r1 = radii[radii.size() / 3], r2 = radii[2 * radii.size() / 3];
    std::cout << jets.size() << " matched fiducial jets, " << radii.size() << " with a UE estimate; region boundaries r = "
              << r1 << ", " << r2 << " cm" << (unweighted ? "; UNWEIGHTED" : "; cross-section weighted") << "; truth pT >= " << truth_pt_min << " GeV" << std::endl;

    std::vector<double> e_bins(n_bins_E + 1);
    for (int i = 0; i <= n_bins_E; i++)
        e_bins[i] = E0_bins + (Emax_bins - E0_bins) * std::pow((double)i / n_bins_E, energy_curvature);

    const char *region_name[n_region] = {"all positions", "inner third", "middle third", "outer third"};
    const char *sel_name[n_sel] = {"raw, all matched jets", "UE subtracted from reco only", "UE subtracted from reco and truth"};
    int sel_color[n_sel] = {kGray + 1, kRed + 1, kAzure + 2};
    int sel_marker[n_sel] = {24, 20, 25};

    TCanvas c("c", "", 1400, 1100);
    c.Divide(2, 2, 0.001, 0.001);
    std::vector<TObject *> keep;
    TFile fout((out_base + ".root").c_str(), "RECREATE");

    for (int ir = 0; ir < n_region; ir++)
    {
        c.cd(ir + 1);
        gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.11);
        gPad->SetLogx();
        TH1F *frame = gPad->DrawFrame(2.0, 0.5, 110.0, 1.6);
        frame->GetXaxis()->SetTitle("median reco jet energy of the bin (GeV)");
        frame->GetYaxis()->SetTitle("median reco / truth energy");
        frame->GetXaxis()->SetMoreLogLabels(); frame->GetXaxis()->SetNoExponent();
        frame->GetXaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetTitleSize(0.045);
        TLine *one = new TLine(2.0, 1.0, 110.0, 1.0);
        one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
        TLegend *leg = new TLegend(0.36, 0.70, 0.96, 0.91);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);

        std::cout << "\n== " << region_name[ir];
        if (ir == 1) std::cout << " (r < " << r1 << " cm)";
        if (ir == 2) std::cout << " (" << r1 << " <= r < " << r2 << " cm)";
        if (ir == 3) std::cout << " (r >= " << r2 << " cm)";
        std::cout << " ==" << std::endl;

        for (int is = 0; is < n_sel; is++)
        {
            std::vector<std::vector<double> > val(n_bins_E), sc(n_bins_E), w(n_bins_E);
            long n_used = 0;
            for (size_t i = 0; i < jets.size(); i++)
            {
                const Jet &j = jets[i];
                if (is > 0 && !j.has_ue) continue;
                if (is == 2 && !j.has_truth_ue) continue;
                const double truth = is == 2 ? j.truth_sub : j.truth;
                if (ir == 1 && !(j.r < r1)) continue;
                if (ir == 2 && !(j.r >= r1 && j.r < r2)) continue;
                if (ir == 3 && !(j.r >= r2)) continue;
                if (truth < e_bins.front() || truth >= e_bins.back()) continue;
                int eb = (int)(std::upper_bound(e_bins.begin(), e_bins.end(), truth) - e_bins.begin()) - 1;
                double reco = is > 0 ? j.sub : j.reco;
                val[eb].push_back(reco); sc[eb].push_back(reco / truth); w[eb].push_back(j.w);
                n_used++;
            }

            TGraphErrors *g = new TGraphErrors();
            for (int eb = 0; eb < n_bins_E; eb++)
            {
                double resp, resp_err, scale, scale_err;
                if (!median_with_error(val[eb], w[eb], resp, resp_err)) continue;
                if (!median_with_error(sc[eb], w[eb], scale, scale_err)) continue;
                if (resp <= 0) continue;
                int n = g->GetN();
                g->SetPoint(n, resp, scale);
                g->SetPointError(n, 0.0, scale_err);
            }
            g->SetMarkerStyle(sel_marker[is]); g->SetMarkerColor(sel_color[is]); g->SetLineColor(sel_color[is]);
            g->SetMarkerSize(1.0);
            g->SetName(Form("g_region%d_sel%d", ir, is));

            std::cout << "  " << sel_name[is] << ": " << n_used << " jets, " << g->GetN() << " points";
            if (g->GetN() < 9) { std::cout << " -- too few to fit" << std::endl; g->Draw("P same"); keep.push_back(g); continue; }
            double x_lo = g->GetX()[0] * 0.95, x_hi = g->GetX()[g->GetN() - 1] * 1.05;

            TF1 *bg = new TF1(Form("bg_%d_%d", ir, is), bg_fn, x_lo, x_hi, 3);
            bg->SetParameters(1.0, 0.0, 0.0);
            g->Fit(bg, "QN0R");
            double chi2_bg = bg->GetChisquare();
            int ndf_bg = g->GetN() - 3;

            // full model: scan seeds of the bump position and widths, keep the best
            TF1 *best = 0;
            double e0_seed[5] = {3, 6, 10, 18, 30}, s_seed[3] = {3, 8, 18};
            for (int a = 0; a < 5; a++) for (int b = 0; b < 3; b++) for (int d = 0; d < 3; d++)
            {
                TF1 *full = new TF1("full_try", full_fn, x_lo, x_hi, 7);
                full->SetParameters(bg->GetParameter(0), bg->GetParameter(1), bg->GetParameter(2), 0.1, e0_seed[a], s_seed[b], s_seed[d]);
                full->SetParLimits(3, D_MIN, D_MAX);
                full->SetParLimits(4, E0_MIN, E0_MAX);
                full->SetParLimits(5, SIGMA_MIN, SIGMA_MAX);
                full->SetParLimits(6, SIGMA_MIN, SIGMA_MAX);
                g->Fit(full, "QN0R");
                if (!best || full->GetChisquare() < best->GetChisquare()) { delete best; best = full; }
                else delete full;
            }
            best->SetName(Form("full_%d_%d", ir, is));
            int ndf_full = g->GetN() - 7;
            std::cout << "; background only chi2/ndf = " << chi2_bg << "/" << ndf_bg << " = " << chi2_bg / ndf_bg
                      << "; full model " << best->GetChisquare() << "/" << ndf_full << " = " << best->GetChisquare() / ndf_full << std::endl;
            std::cout << "      A=" << best->GetParameter(0) << " B=" << best->GetParameter(1) << " C=" << best->GetParameter(2)
                      << " D=" << best->GetParameter(3) << " E0=" << best->GetParameter(4)
                      << " sigmaL=" << best->GetParameter(5) << " sigmaR=" << best->GetParameter(6) << std::endl;
            double worst = 0, worst_E = 0, sum_abs = 0;
            for (int i = 0; i < g->GetN(); i++)
            {
                double d = g->GetY()[i] - best->Eval(g->GetX()[i]);
                sum_abs += std::fabs(d);
                if (std::fabs(d) > std::fabs(worst)) { worst = d; worst_E = g->GetX()[i]; }
            }
            std::cout << "      mean |point - fit| = " << sum_abs / g->GetN() << ", largest " << worst << " at " << worst_E << " GeV" << std::endl;

            best->SetLineColor(sel_color[is]); best->SetLineWidth(2); best->SetNpx(400);
            best->Draw("same");
            g->Draw("P same");
            leg->AddEntry(g, Form("%s  (#chi^{2}/ndf %.1f)", sel_name[is], best->GetChisquare() / ndf_full), "pl");
            g->Write(); best->Write();
            keep.push_back(g); keep.push_back(best);
        }
        leg->Draw();
        TLatex *t = new TLatex();
        t->SetNDC(); t->SetTextSize(0.05); t->SetTextFont(62);
        if (ir == 0) t->DrawLatex(0.11, 0.935, region_name[ir]);
        if (ir == 1) t->DrawLatex(0.11, 0.935, Form("%s of jets, r < %.0f cm", region_name[ir], r1));
        if (ir == 2) t->DrawLatex(0.11, 0.935, Form("%s, %.0f #leq r < %.0f cm", region_name[ir], r1, r2));
        if (ir == 3) t->DrawLatex(0.11, 0.935, Form("%s, r #geq %.0f cm", region_name[ir], r2));
    }

    c.SaveAs((out_base + ".png").c_str());
    c.SaveAs((out_base + ".pdf").c_str());
    fout.Close();
    return 0;
}
