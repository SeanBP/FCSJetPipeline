// Run 22 st_fwd vertex summary (EM/JP-triggered events, SimpleTree vertex
// branches from vertex_test_v1): per-run fractions of events with each vertex,
// per-run mean and RMS of z, and the z distributions/correlations.
// x axis of the per-run plots: day of 2022 from the run number (23DDDRRR),
// plus the run's sequence within the day / 1000.
// Usage: root -b -q 'vertex_summary.C+("out/vtx_*.root,out2/vtx_*.root", "plots")' (comma-separated patterns)
#ifndef __CINT__
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <map>
#include <vector>
#include <cmath>
#include <cstdio>

namespace {
const int kNV = 5;
// vertex definitions: name, "valid" condition, z expression
const char *vName[kNV]  = {"BBC TAC (any z)", "BBC TAC |z| < 150 cm", "EPD (any z)", "EPD |z| < 150 cm", "old BBC zVertex (|z| < 80 cm)"};
const int   vCol[kNV]   = {kAzure + 2, kAzure - 4, kOrange + 7, kOrange - 3, kGray + 1};
const int   vMk[kNV]    = {20, 24, 21, 25, 22};
struct RunAcc { double n = 0, pass[kNV] = {0}, comb = 0; double s1[2] = {0}, s2[2] = {0}, nz[2] = {0}; };
double DayX(int run) { return (run / 1000) % 1000 + (run % 1000) / 1000.0; }
void Frame(const char *name, const char *title, double xlo, double xhi, double ylo, double yhi)
{
    TH1D *f = new TH1D(name, title, 100, xlo, xhi);
    f->SetMinimum(ylo); f->SetMaximum(yhi); f->SetStats(0);
    f->GetXaxis()->SetTitleSize(0.04); f->GetYaxis()->SetTitleSize(0.045); f->GetYaxis()->SetTitleOffset(0.75); f->GetXaxis()->SetTitleOffset(1.0);
    f->Draw("axis");
}
void Note(const char *txt) { TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.03); l.DrawLatex(0.01, 0.965, txt); }
}

#endif
void vertex_summary(const char *pat = "out/vtx_*.root,out2/vtx_*.root", const char *outdir = "plots")
#ifndef __CINT__
{
    gStyle->SetOptStat(0); gStyle->SetPalette(55); gStyle->SetPadGridY(1);
    gSystem->mkdir(outdir, true);
    TChain c("data");
    { TString p(pat); TObjArray *t = p.Tokenize(","); for (int i = 0; i < t->GetEntries(); i++) c.Add(((TObjString*)t->At(i))->String()); delete t; }
    int run; Float_t bbc, epd, det; Int_t src;
    c.SetBranchStatus("*", 0);
    for (const char *b : {"event_run", "bbc_vz_tac", "epd_vz", "bbc_vz_det", "event_vertex_source"}) c.SetBranchStatus(b, 1);
    c.SetBranchAddress("event_run", &run); c.SetBranchAddress("bbc_vz_tac", &bbc); c.SetBranchAddress("epd_vz", &epd);
    c.SetBranchAddress("bbc_vz_det", &det); c.SetBranchAddress("event_vertex_source", &src);
    TH1D hB("hB", "", 120, -600, 600), hE("hE", "", 120, -600, 600), hO("hO", "", 120, -600, 600);
    TH1D hBe("hBe", "", 120, -600, 600), hBl("hBl", "", 120, -600, 600);  // BBC before / from run 23050000
    TH2D hBE("hBE", ";EPD z (cm);BBC TAC z (cm)", 60, -300, 300, 60, -300, 300);
    TH1D hD("hD", ";BBC TAC z #minus EPD z (cm);events", 80, -400, 400);
    std::map<int, RunAcc> acc;
    const Long64_t N = c.GetEntries();
    for (Long64_t i = 0; i < N; i++)
    {
        c.GetEntry(i);
        RunAcc &a = acc[run]; a.n++;
        const bool vb = bbc > -9000, ve = epd > -9000;
        const bool p[kNV] = {vb, vb && fabs(bbc) < 150, ve, ve && fabs(epd) < 150, src == 3};
        for (int k = 0; k < kNV; k++) a.pass[k] += p[k];
        if ((ve && fabs(epd) < 150) || (!ve && vb && fabs(bbc) < 150)) a.comb++;
        if (vb) { hB.Fill(bbc); (run < 23050000 ? hBe : hBl).Fill(bbc); a.s1[0] += bbc; a.s2[0] += bbc * bbc; a.nz[0]++; }
        if (ve) { hE.Fill(epd); a.s1[1] += epd; a.s2[1] += epd * epd; a.nz[1]++; }
        if (fabs(det) < 9000) hO.Fill(det);
        if (vb && ve) { hBE.Fill(epd, bbc); hD.Fill(bbc - epd); }
    }
    printf("%lld events, %zu runs\n", N, acc.size());
    double xlo = 1e9, xhi = -1e9;
    for (auto &r : acc) { xlo = std::min(xlo, DayX(r.first)); xhi = std::max(xhi, DayX(r.first)); }
    xlo = floor(xlo) - 2; xhi = ceil(xhi) + 2;
    const char *xt = "Day of 2022 (run number 23DDDsss: DDD + sss/1000)";
    TString sample = Form("Run 22 st_fwd P24ia, EM/JP-triggered events: %zu runs, 1 MuDst each, %lld events", acc.size(), N);

    // ---- page 1: fraction of events with each vertex vs run
    TCanvas c1("c1", "", 1600, 900);
    c1.SetLeftMargin(0.07); c1.SetRightMargin(0.02); c1.SetTopMargin(0.17); c1.SetBottomMargin(0.1);
    Frame("f1", Form(";%s;fraction of events with a vertex", xt), xlo, xhi, 0, 1.05);
    TLegend lg1(0.07, 0.835, 0.98, 0.94); lg1.SetNColumns(3); lg1.SetBorderSize(0); lg1.SetFillStyle(0); lg1.SetTextSize(0.026);
    std::vector<TGraphErrors *> g1;
    for (int k = 0; k <= kNV; k++)
    {
        TGraphErrors *g = new TGraphErrors();
        for (auto &r : acc)
        {
            const double f = (k < kNV ? r.second.pass[k] : r.second.comb) / r.second.n;
            const int ip = g->GetN(); g->SetPoint(ip, DayX(r.first), f); g->SetPointError(ip, 0, sqrt(f * (1 - f) / r.second.n));
        }
        g->SetMarkerStyle(k < kNV ? vMk[k] : 33); g->SetMarkerSize(k < kNV ? 1.0 : 1.5);
        g->SetMarkerColor(k < kNV ? vCol[k] : kBlack); g->SetLineColor(k < kNV ? vCol[k] : kBlack);
        g->Draw("P SAME"); g1.push_back(g);
        lg1.AddEntry(g, k < kNV ? vName[k] : "combined: EPD if present, else BBC TAC (|z| < 150 cm)", "p");
    }
    lg1.Draw();
    TLine l50(50.0, 0, 50.0, 1.05); l50.SetLineStyle(2); l50.Draw();
    Note(sample + "; dashed: day 50 (EPD vertex from here)");
    c1.SaveAs(Form("%s/vertex_fraction_vs_run.png", outdir)); c1.SaveAs(Form("%s/vertex_fraction_vs_run.pdf", outdir));

    // ---- page 2: per-run mean and RMS of z, BBC TAC and EPD
    TCanvas c2("c2", "", 1600, 1000); c2.Divide(1, 2, 0.001, 0.001);
    const char *ztit[2] = {"mean z (cm)", "RMS of z (cm)"};
    const double zlo[2] = {-80, 0}, zhi[2] = {80, 160};
    TLegend lg2(0.08, 0.80, 0.6, 0.93); lg2.SetNColumns(2); lg2.SetBorderSize(0); lg2.SetFillStyle(0); lg2.SetTextSize(0.05);
    for (int pnl = 0; pnl < 2; pnl++)
    {
        c2.cd(pnl + 1); gPad->SetLeftMargin(0.07); gPad->SetRightMargin(0.02); gPad->SetTopMargin(pnl ? 0.04 : 0.1); gPad->SetBottomMargin(0.14);
        Frame(Form("f2%d", pnl), Form(";%s;%s", xt, ztit[pnl]), xlo, xhi, zlo[pnl], zhi[pnl]);
        for (int d = 0; d < 2; d++)
        {
            TGraphErrors *g = new TGraphErrors();
            for (auto &r : acc)
            {
                const RunAcc &a = r.second; if (a.nz[d] < 20) continue;
                const double m = a.s1[d] / a.nz[d], rms = sqrt(std::max(0.0, a.s2[d] / a.nz[d] - m * m));
                const int ip = g->GetN();
                g->SetPoint(ip, DayX(r.first), pnl ? rms : m);
                g->SetPointError(ip, 0, pnl ? rms / sqrt(2 * a.nz[d]) : rms / sqrt(a.nz[d]));
            }
            const int k = d ? 2 : 0;
            g->SetMarkerStyle(vMk[k]); g->SetMarkerColor(vCol[k]); g->SetLineColor(vCol[k]); g->Draw("P SAME");
            if (pnl == 0) lg2.AddEntry(g, d ? "EPD" : "BBC TAC", "p");
        }
        if (pnl == 0) { lg2.Draw(); TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.045); l.DrawLatex(0.01, 0.94, sample + " (all z; runs with #geq 20 events)"); }
    }
    c2.SaveAs(Form("%s/vertex_z_vs_run.png", outdir)); c2.SaveAs(Form("%s/vertex_z_vs_run.pdf", outdir));

    // ---- page 3: z distributions and BBC-EPD correlation
    TCanvas c3("c3", "", 1800, 650); c3.Divide(3, 1, 0.001, 0.001);
    c3.cd(1); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.14);
    TH1D *hs[3] = {&hB, &hE, &hO};
    const char *hl[3] = {"BBC TAC", "EPD, days #geq 50 only", "old BBC zVertex"};
    const int hc[3] = {kAzure + 2, kOrange + 7, kGray + 1};
    double ymax = 0;
    for (int k = 0; k < 3; k++) { if (hs[k]->Integral() > 0) hs[k]->Scale(1.0 / hs[k]->Integral()); ymax = std::max(ymax, hs[k]->GetMaximum()); }
    TLegend lg3(0.16, 0.71, 0.97, 0.85); lg3.SetBorderSize(0); lg3.SetFillStyle(1001); lg3.SetFillColor(kWhite); lg3.SetTextSize(0.034);
    for (int k = 0; k < 3; k++)
    {
        hs[k]->SetLineColor(hc[k]); hs[k]->SetLineWidth(2); hs[k]->SetMaximum(1.4 * ymax);
        hs[k]->SetTitle(";z (cm);fraction of events per 10 cm"); hs[k]->GetYaxis()->SetTitleOffset(1.5); hs[k]->GetXaxis()->SetRangeUser(-500, 500);
        hs[k]->Draw(k ? "hist same" : "hist");
        lg3.AddEntry(hs[k], Form("%s (RMS %.0f cm)", hl[k], hs[k]->GetRMS()), "l");
    }
    lg3.Draw();
    { TLatex l; l.SetNDC(); l.SetTextFont(62); l.SetTextSize(0.042); l.DrawLatex(0.13, 0.9, "z of each vertex (unit area)"); }
    c3.cd(2); gPad->SetLogz(); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.13); gPad->SetTopMargin(0.14);
    hBE.GetYaxis()->SetTitleOffset(1.5); hBE.Draw("colz");
    { TLatex l; l.SetNDC(); l.SetTextFont(62); l.SetTextSize(0.042); l.DrawLatex(0.13, 0.9, Form("BBC TAC vs EPD z (correlation %.2f)", hBE.GetCorrelationFactor())); }
    c3.cd(3); gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.14);
    hD.GetYaxis()->SetTitleOffset(2.0); hD.SetMaximum(1.35 * hD.GetMaximum()); hD.SetLineColor(kBlack); hD.SetLineWidth(2); hD.Fit("gaus", "Q", "", -150, 150); hD.Draw("hist same");
    { TLatex l; l.SetNDC(); l.SetTextFont(62); l.SetTextSize(0.042); l.DrawLatex(0.13, 0.9, "BBC TAC z #minus EPD z");
      TLegend *bx = new TLegend(0.5, 0.74, 0.96, 0.85); bx->SetBorderSize(0); bx->SetFillStyle(1001); bx->SetFillColor(kWhite); bx->SetTextSize(0.034);
      bx->AddEntry((TObject*)0, Form("RMS %.0f cm", hD.GetRMS()), "");
      bx->AddEntry(hD.GetFunction("gaus"), Form("core Gaussian #sigma %.0f cm", hD.GetFunction("gaus")->GetParameter(2)), "l"); bx->Draw(); }
    c3.cd(0); Note(sample);
    c3.SaveAs(Form("%s/vertex_z_distributions.png", outdir)); c3.SaveAs(Form("%s/vertex_z_distributions.pdf", outdir));
}
#else
;
#endif
