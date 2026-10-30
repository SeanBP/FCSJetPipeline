// Plot of the recoil_truth outputs: per jet-energy bin, (top) the EM energy
// density in the jet's eta band against |dphi| with and without multi-parton
// interactions, and their difference; (bottom) where the recoiling hard
// parton goes in eta relative to the jet. ptHat bins are combined with
// weights sigma / N (all jobs of a bin pooled).
//
// Build and run (plain ROOT):
//   g++ -O2 recoil_plot.cpp -o recoil_plot `root-config --cflags --libs` -std=c++11
//   ./recoil_plot <dir with rt_mpi{0,1}_pt*_*.root> <out.png>

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

static const int kNB = 4;
static const double kEdge[kNB + 1] = {10, 20, 35, 60, 120};
static const double R = 0.7;

static double jk_a[2][4][11], jk_c[2][4][11], jk_nm[4][11], jk_nr[4][11];

struct Sum
{
    TH1D *band[kNB], *rdeta[kNB], *rwhere[kNB], *njet;
    TProfile *cone[kNB];
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::fprintf(stderr, "usage: ./recoil_plot <dir> <out.png>\n"); return 1; }
    const std::string dir = argv[1];
    // files per (mpi, ptHat bin)
    std::map<std::pair<int, int>, std::vector<std::string> > files;
    void *d = gSystem->OpenDirectory(dir.c_str());
    const char *e;
    while ((e = gSystem->GetDirEntry(d)))
    {
        int mpi, pt, j;
        if (std::sscanf(e, "rt_mpi%d_pt%d_%d.root", &mpi, &pt, &j) == 3) files[std::make_pair(mpi, pt)].push_back(dir + "/" + e);
    }
    gSystem->FreeDirectory(d);

    Sum S[2];
    for (int m = 0; m < 2; m++)
    {
        S[m].njet = 0;
        for (int b = 0; b < kNB; b++) { S[m].band[b] = 0; S[m].rdeta[b] = 0; S[m].rwhere[b] = 0; S[m].cone[b] = 0; }
    }
    for (std::map<std::pair<int, int>, std::vector<std::string> >::iterator it = files.begin(); it != files.end(); ++it)
    {
        const int m = it->first.first;
        // pooled weight of this ptHat bin: mean sigma / total events
        double sig = 0, nev = 0;
        for (size_t f = 0; f < it->second.size(); f++)
        {
            TFile in(it->second[f].c_str());
            TTree *t = (TTree *)in.Get("info");
            if (!t) continue;
            double s, n; t->SetBranchAddress("sigma", &s); t->SetBranchAddress("events", &n); t->GetEntry(0);
            sig += s * n; nev += n;
        }
        if (!(nev > 0)) continue;
        const double w = (sig / nev) / nev;
        std::printf("mpi %d ptHat %d: %zu files, %.0f events, sigma %.4g mb\n", m, it->first.second, it->second.size(), nev, sig / nev);
        for (size_t f = 0; f < it->second.size(); f++)
        {
            TFile in(it->second[f].c_str());
            if (in.IsZombie() || !in.Get("info")) continue;
            gROOT->cd();
            TH1D *nj = (TH1D *)in.Get("h_njet");
            if (!S[m].njet) { S[m].njet = (TH1D *)nj->Clone(Form("nj_%d", m)); S[m].njet->Reset(); }
            S[m].njet->Add(nj, w);
            for (int b = 0; b < kNB; b++)
            {
                TH1D *hb = (TH1D *)in.Get(Form("h_band_%d", b)), *hr = (TH1D *)in.Get(Form("h_rec_deta_%d", b)), *hw = (TH1D *)in.Get(Form("h_rec_where_%d", b));
                TProfile *pc = (TProfile *)in.Get(Form("p_cone_%d", b));
                if (!S[m].band[b])
                {
                    S[m].band[b] = (TH1D *)hb->Clone(Form("band_%d_%d", m, b)); S[m].band[b]->Reset();
                    S[m].rdeta[b] = (TH1D *)hr->Clone(Form("rdeta_%d_%d", m, b)); S[m].rdeta[b]->Reset();
                    S[m].rwhere[b] = (TH1D *)hw->Clone(Form("rwhere_%d_%d", m, b)); S[m].rwhere[b]->Reset();
                    S[m].cone[b] = (TProfile *)pc->Clone(Form("cone_%d_%d", m, b)); S[m].cone[b]->Reset();
                }
                S[m].band[b]->Add(hb, w); S[m].rdeta[b]->Add(hr, w); S[m].rwhere[b]->Add(hw, w);
                S[m].cone[b]->Add(pc, w);
            }
        }
    }
    if (!S[1].njet || !S[0].njet) { std::fprintf(stderr, "ERROR: need both mpi 0 and 1 outputs\n"); return 1; }

    // jackknife over the job index (job j left out of every ptHat bin):
    // UE180/UE90 for MPI on and off, and the recoil fraction inside the 180 cone
    const int kJ = 10;
    double jk_val[3][kNB][kJ + 1];
    {
        double a[2][kNB][kJ + 1] = {}, c[2][kNB][kJ + 1] = {}, nm[kNB][kJ + 1] = {}, nr[kNB][kJ + 1] = {};
        for (std::map<std::pair<int, int>, std::vector<std::string> >::iterator it = files.begin(); it != files.end(); ++it)
        {
            const int m = it->first.first;
            double sig = 0, nev = 0;
            for (size_t f = 0; f < it->second.size(); f++)
            {
                int mm, pt, j;
                const std::string base = it->second[f].substr(it->second[f].find_last_of('/') + 1);
                std::sscanf(base.c_str(), "rt_mpi%d_pt%d_%d.root", &mm, &pt, &j);
                TFile in(it->second[f].c_str());
                TTree *t = (TTree *)in.Get("info");
                if (!t) continue;
                double sg, n; t->SetBranchAddress("sigma", &sg); t->SetBranchAddress("events", &n); t->GetEntry(0);
                sig += sg * n; nev += n;
                for (int b = 0; b < kNB; b++)
                {
                    TProfile *p = (TProfile *)in.Get(Form("p_cone_%d", b));
                    TH1D *hw = (TH1D *)in.Get(Form("h_rec_where_%d", b));
                    const double x180 = p->GetBinContent(1) * p->GetBinEntries(1);
                    const double x90 = 0.5 * (p->GetBinContent(4) * p->GetBinEntries(4) + p->GetBinContent(5) * p->GetBinEntries(5));
                    for (int k = 0; k <= kJ; k++)
                    {
                        if (k == j) continue;
                        a[m][b][k] += x180; c[m][b][k] += x90;
                        if (m == 1) { nm[b][k] += hw->GetBinContent(1); nr[b][k] += hw->GetBinContent(3); }
                    }
                }
            }
            // weights: sigma / events per bin; a leave-one-out sum lacks one job's events
            const double w_all = (sig / nev) / nev, w_loo = (sig / nev) / (nev * (kJ - 1) / kJ);
            // fold this bin into the weighted totals
            for (int b = 0; b < kNB; b++)
                for (int k = 0; k <= kJ; k++)
                {
                    const double w = k == kJ ? w_all : w_loo;
                    jk_a[m][b][k] += w * a[m][b][k]; jk_c[m][b][k] += w * c[m][b][k];
                    if (m == 1) { jk_nm[b][k] += w * nm[b][k]; jk_nr[b][k] += w * nr[b][k]; }
                    a[m][b][k] = c[m][b][k] = 0;
                    if (m == 1) nm[b][k] = nr[b][k] = 0;
                }
        }
        for (int b = 0; b < kNB; b++)
            for (int k = 0; k <= kJ; k++)
            {
                jk_val[0][b][k] = jk_a[1][b][k] / jk_c[1][b][k];
                jk_val[1][b][k] = jk_a[0][b][k] / jk_c[0][b][k];
                jk_val[2][b][k] = jk_nr[b][k] / jk_nm[b][k];
            }
    }
    double jk_err[3][kNB];
    for (int q = 0; q < 3; q++)
        for (int b = 0; b < kNB; b++)
        {
            double mu = 0, v = 0;
            for (int k = 0; k < kJ; k++) mu += jk_val[q][b][k] / kJ;
            for (int k = 0; k < kJ; k++) v += (jk_val[q][b][k] - mu) * (jk_val[q][b][k] - mu);
            jk_err[q][b] = std::sqrt(v * (kJ - 1) / kJ);
        }

    gStyle->SetOptStat(0);
    const int kFull = TColor::GetColor("#222222"), kHard = TColor::GetColor("#d9480f"), kMpi = TColor::GetColor("#1c6bb0");
    TCanvas cv("cv", "", 1800, 900);
    cv.Divide(kNB, 2, 0.002, 0.004);
    TLatex tx; tx.SetNDC(); tx.SetTextFont(42);
    std::printf("\n jet E      UE180/UE90 MPI on (jackknife)   MPI off            recoil parton in band, in 180 cone, in 90 cones\n");
    std::vector<TObject *> keep;
    for (int b = 0; b < kNB; b++)
    {
        // density per jet per unit eta-phi area: band height 2R, |dphi| folds both sides
        TH1D *h[3];
        for (int m = 0; m < 2; m++)
        {
            h[m] = (TH1D *)S[m].band[b]->Clone(Form("dens_%d_%d", m, b));
            h[m]->Rebin(2);
            const double nj = S[m].njet->GetBinContent(b + 1);
            h[m]->Scale(1.0 / (nj * 2 * R * 2 * h[m]->GetBinWidth(1)));
        }
        // MPI on is index 1, off index 0
        h[2] = (TH1D *)h[1]->Clone(Form("diff_%d", b)); h[2]->Add(h[0], -1);
        cv.cd(b + 1);
        gPad->SetMargin(0.15, 0.03, 0.13, 0.10); gPad->SetLogy();
        h[1]->SetLineColor(kFull); h[1]->SetLineWidth(2);
        h[0]->SetLineColor(kHard); h[0]->SetLineWidth(2);
        h[2]->SetLineColor(kMpi); h[2]->SetLineWidth(2); h[2]->SetLineStyle(2);
        h[1]->SetTitle(";|#Delta#phi| to the jet (rad);EM dE/d#etad#phi per jet (GeV), |#Delta#eta| < R");
        h[1]->GetYaxis()->SetTitleSize(0.05); h[1]->GetXaxis()->SetTitleSize(0.05);
        h[1]->GetYaxis()->SetLabelSize(0.045); h[1]->GetXaxis()->SetLabelSize(0.045);
        h[1]->SetMinimum(0.1); h[1]->SetMaximum(h[1]->GetMaximum() * 3);
        h[1]->Draw("HIST");
        TBox *bx = new TBox(0, 0.1, R, h[1]->GetMaximum()); bx->SetFillColor(TColor::GetColor("#eeeeee")); bx->Draw();
        keep.push_back(bx);
        for (int q = 0; q < 2; q++)
        {
            const double x = q == 0 ? M_PI / 2 : M_PI;
            TBox *c = new TBox(x - R, 0.1, std::min(x + R, M_PI), h[1]->GetMaximum());
            c->SetFillColor(TColor::GetColor(q == 0 ? "#e7f0fa" : "#fbeee6")); c->Draw(); keep.push_back(c);
        }
        h[1]->Draw("HIST SAME"); h[0]->Draw("HIST SAME"); h[2]->Draw("HIST SAME");
        gPad->RedrawAxis();
        tx.SetTextSize(0.055); tx.SetTextColor(kFull);
        tx.DrawLatex(0.17, 0.92, Form("truth jets %.0f-%.0f GeV, 2.5 < #eta < 4", kEdge[b], kEdge[b + 1]));
        tx.SetTextSize(0.04);
        tx.DrawLatex(0.20, 0.835, "jet");
        tx.DrawLatex(0.43, 0.835, "90#circ cone");
        tx.DrawLatex(0.78, 0.835, "180#circ cone");
        if (b == 0)
        {
            TLegend *lg = new TLegend(0.30, 0.57, 0.97, 0.80);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
            lg->AddEntry(h[1], "Pythia6 FMS tune (MPI on)", "l");
            lg->AddEntry(h[0], "MPI off (hard scatter, showers)", "l");
            lg->AddEntry(h[2], "difference: MPI", "l");
            lg->Draw(); keep.push_back(lg);
        }
        // cone ratios
        const double c180[2] = {S[0].cone[b]->GetBinContent(1), S[1].cone[b]->GetBinContent(1)};
        const double c90[2] = {0.5 * (S[0].cone[b]->GetBinContent(4) + S[0].cone[b]->GetBinContent(5)), 0.5 * (S[1].cone[b]->GetBinContent(4) + S[1].cone[b]->GetBinContent(5))};
        const double hard_ex = (c180[0] - c90[0]) / c90[1];
        tx.SetTextSize(0.045);
        tx.DrawLatex(0.20, 0.21, Form("UE_{180}/UE_{90} = %.2f #pm %.2f (MPI on)", jk_val[0][b][kJ], jk_err[0][b]));
        tx.SetTextColor(kHard);
        tx.DrawLatex(0.20, 0.155, Form("MPI off: %.2f #pm %.2f", jk_val[1][b][kJ], jk_err[1][b]));
        tx.SetTextColor(kFull);

        cv.cd(kNB + b + 1);
        gPad->SetMargin(0.15, 0.03, 0.13, 0.06);
        TH1D *rd = S[1].rdeta[b];
        rd->Rebin(4);
        const double nm = S[1].rwhere[b]->GetBinContent(1);
        rd->Scale(1.0 / (nm * rd->GetBinWidth(1)));
        rd->SetTitle(";#eta_{recoil parton} - #eta_{jet};recoil partons per jet per unit #Delta#eta");
        rd->GetYaxis()->SetTitleSize(0.05); rd->GetXaxis()->SetTitleSize(0.05);
        rd->GetYaxis()->SetLabelSize(0.045); rd->GetXaxis()->SetLabelSize(0.045);
        rd->SetLineColor(kHard); rd->SetFillColor(TColor::GetColor("#f7c59f")); rd->SetLineWidth(2);
        rd->SetMaximum(rd->GetMaximum() * 1.35);
        rd->Draw("HIST");
        TBox *bb = new TBox(-R, 0, R, rd->GetMaximum()); bb->SetFillColorAlpha(TColor::GetColor("#1c6bb0"), 0.0);
        bb->SetFillStyle(3354); bb->SetFillColor(TColor::GetColor("#1c6bb0")); bb->Draw(); keep.push_back(bb);
        rd->Draw("HIST SAME"); gPad->RedrawAxis();
        const double f_band = S[1].rwhere[b]->GetBinContent(2) / nm, f180 = S[1].rwhere[b]->GetBinContent(3) / nm, f90 = S[1].rwhere[b]->GetBinContent(4) / nm;
        tx.SetTextSize(0.045);
        tx.DrawLatex(0.18, 0.86, "recoil parton (jet within R of a hard parton)");
        tx.DrawLatex(0.18, 0.80, Form("in the jet's #eta band (hatched): %.1f%%", 100 * f_band));
        tx.DrawLatex(0.18, 0.74, Form("inside the 180#circ cone: %.1f #pm %.1f%%", 100 * f180, 100 * jk_err[2][b]));
        tx.DrawLatex(0.18, 0.68, Form("inside a #pm90#circ cone: %.1f%%", 100 * f90));
        std::printf("%5.0f-%-5.0f  %.3f +- %.3f              %.3f +- %.3f      %.3f  %.3f  %.3f\n", kEdge[b], kEdge[b + 1],
                    jk_val[0][b][kJ], jk_err[0][b], jk_val[1][b][kJ], jk_err[1][b], f_band, f180, f90);
        (void)hard_ex;
    }
    cv.SaveAs(argv[2]);
    return 0;
}
