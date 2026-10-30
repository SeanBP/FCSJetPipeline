// How much of the FCS pixel ring a cut around the recoil side removes, per
// jet: the ring (in_event_geo.h BuildPixelRing; |eta - eta_jet| < R, >= 2R
// from the jet, no other-jet exclusion) with cells at |dphi| above a cap left
// out, caps kCaps (180 = no cut). Geometry only: the jets' positions (tree
// "jets" of in_event_ue outputs, jet_eta, jet_phi). Per cap: the live ring
// area kept (fraction of the uncut ring's), the share of the jet's live area
// in eta bins the ring reaches (cov; 1 = every eta bin matched), the share of
// jets with cov below 1, 0.9, 0.5 and with no ring at all (the template
// fallback), and the largest cell weight. Also a plot of the distributions.
//
// Build and run from this folder:
//   g++ -m64 -O2 ring_cut_check.cpp -o ring_cut_check -I../../shared -I../../underlying_event `root-config --cflags --libs` -std=c++11
//   ./ring_cut_check "<in_event outputs>" <out.png>

#include <TROOT.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "in_event_geo.h"

static const int kNC = 4;
static const double kCaps[kNC] = {180, 175, 150, 130}; // deg

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::fprintf(stderr, "usage: ./ring_cut_check \"<in_event outputs>\" <out.png>\n"); return 1; }
    TChain ch("jets");
    std::vector<std::string> f = Split(argv[1]);
    for (size_t i = 0; i < f.size(); i++) ch.Add(f[i].c_str());
    Float_t eta, phi;
    ch.SetBranchStatus("*", 0); ch.SetBranchStatus("jet_eta", 1); ch.SetBranchStatus("jet_phi", 1);
    ch.SetBranchAddress("jet_eta", &eta); ch.SetBranchAddress("jet_phi", &phi);
    Geometry geo; geo.Build();
    const std::vector<std::pair<float, float> > none;

    TH1D *h_keep[kNC], *h_cov[kNC];
    double n = 0, s_keep[kNC] = {0}, s_cov[kNC] = {0}, n_c1[kNC] = {0}, n_c09[kNC] = {0}, n_c05[kNC] = {0}, n_none[kNC] = {0}, s_wmax[kNC] = {0}, n_k05[kNC] = {0};
    std::vector<double> keeps[kNC];
    double worst[kNC], worst_eta[kNC], worst_phi[kNC];
    for (int c = 0; c < kNC; c++) worst[c] = 2;
    for (int c = 0; c < kNC; c++)
    {
        h_keep[c] = new TH1D(Form("hk%d", c), ";live ring area kept / uncut ring;jets (normalised)", 50, 0, 1.0001);
        h_cov[c] = new TH1D(Form("hc%d", c), ";share of the jet's live area in eta bins the ring reaches;jets (normalised)", 50, 0, 1.0001);
    }
    PixelRing pr;
    for (Long64_t i = 0; i < ch.GetEntries(); i++)
    {
        ch.GetEntry(i);
        n++;
        double area0 = 0;
        for (int c = 0; c < kNC; c++)
        {
            BuildPixelRing(geo, eta, phi, none, kAnnJetExcl, pr, kCaps[c] * M_PI / 180 + 1e-9);
            if (c == 0) area0 = pr.area_ring;
            const double keep = area0 > 0 ? pr.area_ring / area0 : 0, cov = pr.area_jet > 0 ? pr.den / pr.area_jet : 0;
            if (keep < worst[c]) { worst[c] = keep; worst_eta[c] = eta; worst_phi[c] = phi; }
            s_keep[c] += keep; s_cov[c] += cov; s_wmax[c] += pr.w_max; keeps[c].push_back(keep);
            if (cov < 0.999) n_c1[c]++;
            if (cov < 0.9) n_c09[c]++;
            if (cov < 0.5) n_c05[c]++;
            if (!(pr.den > 0)) n_none[c]++;
            if (keep < 0.5) n_k05[c]++;
            h_keep[c]->Fill(keep); h_cov[c]->Fill(cov);
        }
    }
    std::printf("ring_cut_check: %.0f jets\n  cut (|dphi| <=)   area kept: mean  min    5%%-ile  <50%%     cov: mean   <1      <0.9    <0.5    no ring   mean w_max\n", n);
    for (int c = 0; c < kNC; c++)
    {
        std::sort(keeps[c].begin(), keeps[c].end());
        std::printf("  %3.0f deg            %.3f  %.3f  %.3f  %.4f      %.4f  %.4f  %.4f  %.4f  %.4f    %.3f\n", kCaps[c], s_keep[c] / n, keeps[c].front(), keeps[c][(size_t)(0.05 * keeps[c].size())],
                    n_k05[c] / n, s_cov[c] / n, n_c1[c] / n, n_c09[c] / n, n_c05[c] / n, n_none[c] / n, s_wmax[c] / n);
    }

    for (int c = 1; c < kNC; c++)
        std::printf("  worst jet at |dphi| <= %.0f deg: keeps %.3f, eta %.3f, phi %.1f deg\n", kCaps[c], worst[c], worst_eta[c], worst_phi[c] * 180 / M_PI);

    gStyle->SetOptStat(0);
    const int col[kNC] = {TColor::GetColor("#222222"), TColor::GetColor("#4292c6"), TColor::GetColor("#d9480f"), TColor::GetColor("#2f9e44")};
    TCanvas cv("cv", "", 1500, 620);
    cv.Divide(2, 1, 0.002, 0.002);
    for (int p = 0; p < 2; p++)
    {
        cv.cd(p + 1);
        gPad->SetMargin(0.12, 0.03, 0.12, 0.08); gPad->SetLogy();
        TLegend *lg = new TLegend(0.15, 0.70, 0.6, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
        for (int c = (p == 0 ? 1 : 0); c < kNC; c++)
        {
            TH1D *h = p ? h_cov[c] : h_keep[c];
            h->Scale(1.0 / h->Integral());
            h->SetLineColor(col[c]); h->SetLineWidth(2); h->SetMinimum(1e-5); h->SetMaximum(3);
            h->GetYaxis()->SetTitleOffset(1.3);
            h->Draw(c == (p == 0 ? 1 : 0) ? "HIST" : "HIST SAME");
            lg->AddEntry(h, kCaps[c] == 180 ? "no cut" : Form("ring |#Delta#phi| #leq %.0f#circ", kCaps[c]), "l");
        }
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.036);
        t.DrawLatex(0.12, 0.94, Form("FCS pixel ring (2R, no other-jet exclusion), %.0f jets: %s", n, p ? "eta coverage" : "area kept by the cut"));
    }
    cv.SaveAs(argv[2]);
    return 0;
}
