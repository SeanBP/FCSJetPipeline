// How much does the mean UE estimate depend on the ptHat bin, at fixed jet pT
// and position? In simulation UnderlyingEvent.cpp draws a jet's proxy events
// only from its own ptHat bin (events are unweighted within a bin), while in
// data the proxies come from the natural mix of events. If the UE at fixed
// jet pT and position varies with ptHat, the simulation subtraction differs
// from the data one.
//
// For matched, fiducial jets with an estimate, in bins of jet pT and in three
// radial bands (thirds of the jets), takes the mean reco_ue_E (and
// truth_ue_E, binned in truth pT) per ptHat bin, UE_k, and the
// cross-section-weighted mean over bins, UE_w = sum w_k n_k UE_k / sum w_k n_k,
// which is what a pool of all bins with cross-section-weighted proxies would
// give on average. The difference UE_k - UE_w is the shift of a jet's
// subtracted energy caused by pooling per bin; its cross-section-weighted RMS
// over jets is compared with the mean jet energy.
//
// Usage: ./ue_pthat_dependence <file list> [output_base]
// Build: g++ -m64 -O2 ue_pthat_dependence.cpp ../shared/VectorDict.cxx -o ue_pthat_dependence
//            -I../shared -I../jet_scale `root-config --cflags --libs` -std=c++11

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"

const int n_pt = 7;
const double pt_edges[n_pt + 1] = {0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 10.0};
const int n_region = 3, n_level = 2;

struct Acc
{
    double n, sum_ue, sum_E;
    Acc() : n(0), sum_ue(0), sum_E(0) {}
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    if (argc < 2) { std::cerr << "usage: ./ue_pthat_dependence <file list> [output_base]" << std::endl; return 1; }
    std::string out_base = argc > 2 ? argv[2] : "ue_pthat_dependence";

    std::vector<std::string> files;
    {
        std::ifstream in(argv[1]);
        std::string s;
        while (in >> s) files.push_back(s);
    }
    std::vector<float> file_weight;
    if (!ComputeXSecWeights(files, "jetTree", file_weight)) return 1;

    // first pass would be needed for exact thirds; use the boundaries of the
    // whole sample measured by the form check instead
    const double r1 = 66.35, r2 = 87.53;

    // acc[level][region][pt][pthat] (unweighted counts and sums), and the bin weight
    std::map<float, double> bin_weight;
    std::vector<std::vector<std::vector<std::map<float, Acc> > > > acc(
        n_level, std::vector<std::vector<std::map<float, Acc> > >(n_region, std::vector<std::map<float, Acc> >(n_pt)));

    for (size_t fidx = 0; fidx < files.size(); fidx++)
    {
        if (file_weight[fidx] < 0) continue;
        TFile *f = TFile::Open(files[fidx].c_str());
        TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("jetTree") : 0;
        if (!t || !t->GetBranch("truth_ue_E")) { if (f) { f->Close(); delete f; } continue; }
        std::vector<float> *reco_pt = 0, *reco_E = 0, *reco_x = 0, *reco_y = 0, *reco_ue_E = 0;
        std::vector<float> *truth_pt = 0, *truth_E = 0, *truth_ue_E = 0;
        std::vector<int> *reco_ue_n = 0, *truth_ue_n = 0, *mr = 0, *mt = 0;
        Float_t pthat = 0;
        t->SetBranchStatus("*", 0);
        const char *used[] = {"reco_pt", "reco_E", "reco_x", "reco_y", "reco_ue_E", "reco_ue_nproxy", "truth_pt", "truth_E",
                              "truth_ue_E", "truth_ue_nproxy", "match_reco_idx", "match_truth_idx", "mc_pthatmin_gev"};
        for (int k = 0; k < 13; k++) t->SetBranchStatus(used[k], 1);
        t->SetBranchAddress("reco_pt", &reco_pt);
        t->SetBranchAddress("reco_E", &reco_E);
        t->SetBranchAddress("reco_x", &reco_x);
        t->SetBranchAddress("reco_y", &reco_y);
        t->SetBranchAddress("reco_ue_E", &reco_ue_E);
        t->SetBranchAddress("reco_ue_nproxy", &reco_ue_n);
        t->SetBranchAddress("truth_pt", &truth_pt);
        t->SetBranchAddress("truth_E", &truth_E);
        t->SetBranchAddress("truth_ue_E", &truth_ue_E);
        t->SetBranchAddress("truth_ue_nproxy", &truth_ue_n);
        t->SetBranchAddress("match_reco_idx", &mr);
        t->SetBranchAddress("match_truth_idx", &mt);
        t->GetEntry(0);
        pthat = (Float_t)t->GetLeaf("mc_pthatmin_gev")->GetValue();
        bin_weight[pthat] = file_weight[fidx];

        for (Long64_t ev = 0; ev < t->GetEntries(); ev++)
        {
            t->GetEntry(ev);
            for (size_t m = 0; m < mr->size(); m++)
            {
                const int ri = mr->at(m), ti = mt->at(m);
                const float x = reco_x->at(ri), y = reco_y->at(ri);
                if (!pass_fiducial_cut(x, y, reco_fiducial_buffer_ecal, kFiducialRectEcal)) continue;
                const double r = std::sqrt(x * x + y * y);
                const int ir = r < r1 ? 0 : (r < r2 ? 1 : 2);
                for (int lv = 0; lv < n_level; lv++)
                {
                    const bool has = lv == 0 ? reco_ue_n->at(ri) > 0 : truth_ue_n->at(ti) > 0;
                    if (!has) continue;
                    const double pt = lv == 0 ? reco_pt->at(ri) : truth_pt->at(ti);
                    const int ip = (int)(std::upper_bound(pt_edges, pt_edges + n_pt + 1, pt) - pt_edges) - 1;
                    if (ip < 0 || ip >= n_pt) continue;
                    Acc &a = acc[lv][ir][ip][pthat];
                    a.n += 1;
                    a.sum_ue += lv == 0 ? reco_ue_E->at(ri) : truth_ue_E->at(ti);
                    a.sum_E += lv == 0 ? reco_E->at(ri) : truth_E->at(ti);
                }
            }
        }
        f->Close();
        delete f;
    }

    const char *level_name[n_level] = {"reco", "truth"};
    const char *region_name[n_region] = {"inner", "middle", "outer"};
    TCanvas c("c", "", 1500, 900);
    c.Divide(3, 2, 0.001, 0.001);
    std::vector<TObject *> keep;
    int colors[n_pt] = {kGray + 2, kViolet + 1, kAzure + 2, kTeal + 2, kGreen + 2, kOrange + 1, kRed + 1};

    for (int lv = 0; lv < n_level; lv++)
    {
        std::cout << "\n===== " << level_name[lv] << " level (UE estimate " << level_name[lv] << "_ue_E, binned in " << level_name[lv]
                  << " pT) =====" << std::endl;
        for (int ir = 0; ir < n_region; ir++)
        {
            c.cd(lv * 3 + ir + 1);
            gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.12);
            gPad->SetLogy();
            TH1F *frame = gPad->DrawFrame(-3, 0.05, 105, 30);
            frame->GetXaxis()->SetTitle("ptHat bin lower edge (GeV)");
            frame->GetYaxis()->SetTitle(Form("mean %s UE energy in cone (GeV)", level_name[lv]));
            frame->GetXaxis()->SetTitleSize(0.05); frame->GetYaxis()->SetTitleSize(0.05);
            TLegend *leg = new TLegend(0.55, 0.62, 0.97, 0.91);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
            TLatex *tl = new TLatex();
            tl->SetNDC(); tl->SetTextSize(0.055); tl->SetTextFont(62);
            tl->DrawLatex(0.13, 0.935, Form("%s level, %s third", level_name[lv], region_name[ir]));
            keep.push_back(leg); keep.push_back(tl);

            std::cout << "-- " << region_name[ir] << " third --" << std::endl;
            for (int ip = 0; ip < n_pt; ip++)
            {
                const std::map<float, Acc> &m = acc[lv][ir][ip];
                double sw = 0, sw_ue = 0, sw_E = 0, sw_jets = 0;
                for (std::map<float, Acc>::const_iterator it = m.begin(); it != m.end(); ++it)
                {
                    const double w = bin_weight[it->first];
                    sw += w * it->second.n; sw_ue += w * it->second.sum_ue; sw_E += w * it->second.sum_E;
                    sw_jets += it->second.n;
                }
                if (sw <= 0) continue;
                const double ue_w = sw_ue / sw, E_w = sw_E / sw;
                double var = 0, lo = 1e9, hi = -1e9, dominant = 0;
                float dominant_bin = 0;
                TGraph *g = new TGraph();
                for (std::map<float, Acc>::const_iterator it = m.begin(); it != m.end(); ++it)
                {
                    if (it->second.n < 20) continue;
                    const double w = bin_weight[it->first];
                    const double ue_k = it->second.sum_ue / it->second.n;
                    var += w * it->second.n * (ue_k - ue_w) * (ue_k - ue_w);
                    lo = std::min(lo, ue_k); hi = std::max(hi, ue_k);
                    if (w * it->second.n > dominant) { dominant = w * it->second.n; dominant_bin = it->first; }
                    g->SetPoint(g->GetN(), it->first, ue_k);
                }
                const double rms = std::sqrt(var / sw);
                std::printf("  pT %5.2f-%5.2f GeV: %8.0f jets, weighted mean UE %6.3f GeV, jet E %6.2f GeV; per-bin UE %6.3f to %6.3f;"
                            " weighted RMS of (UE_k - UE_w) %6.3f GeV = %5.2f %% of jet E; %4.1f %% of the weight in ptHat %g-\n",
                            pt_edges[ip], pt_edges[ip + 1], sw_jets, ue_w, E_w, lo, hi, rms, 100.0 * rms / E_w,
                            100.0 * dominant / sw, dominant_bin);
                g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->SetMarkerColor(colors[ip]); g->SetLineColor(colors[ip]);
                g->Draw("PL same");
                TLine *l = new TLine(-3, ue_w, 105, ue_w);
                l->SetLineColor(colors[ip]); l->SetLineStyle(2); l->Draw();
                leg->AddEntry(g, Form("pT %.2g-%.2g GeV", pt_edges[ip], pt_edges[ip + 1]), "pl");
                keep.push_back(g); keep.push_back(l);
            }
            leg->Draw();
        }
    }
    c.SaveAs((out_base + ".png").c_str());
    c.SaveAs((out_base + ".pdf").c_str());
    return 0;
}
