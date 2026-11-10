#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TMath.h"
#include <cstdio>
// Vertex test summary: availability per detector, the new chain, correlations.
void analyze(const char *pat = "out/vtx_*.root", const char *png = "vertex_test.png")
{
    gStyle->SetOptStat(0); gStyle->SetPalette(55);
    TChain c("data"); c.Add(pat);
    const double N = c.GetEntries();
    printf("events (EM/JP filtered): %.0f in %d files\n", N, c.GetNtrees());
    struct S { const char *name, *cut; } s[] = {
        {"old event_vertex (|z|<80)", "event_vertex_source>0"},
        {"old BBC zVertex any z", "abs(bbc_vz_det)<9000"},
        {"BBC TAC valid", "bbc_vz_tac>-9000"},
        {"BBC TAC |z|<150", "abs(bbc_vz_tac)<150"},
        {"BBC TAC |z|<200", "abs(bbc_vz_tac)<200"},
        {"EPD valid", "epd_vz>-9000"},
        {"EPD |z|<150", "abs(epd_vz)<150"},
        {"VPD valid", "vpd_vz>-9000"},
        {"EPD or BBC-TAC valid", "epd_vz>-9000||bbc_vz_tac>-9000"},
        {"chain EPD>BBC |z|<150", "(epd_vz>-9000&&abs(epd_vz)<150)||(epd_vz<-9000&&abs(bbc_vz_tac)<150)"}};
    for (auto &x : s) printf("  %-28s %6.1f%%\n", x.name, 100. * c.GetEntries(x.cut) / N);
    TCanvas cv("cv", "", 1800, 1000); cv.Divide(4, 2);
    cv.cd(1); c.Draw("bbc_tac_diff>>h0(400,-1000,1000)", "bbc_tac_diff>-4000"); ((TH1*)gROOT->FindObject("h0"))->SetTitle("BBC TAC diff (-4096 offset);TAC counts");
    cv.cd(2); c.Draw("bbc_vz_tac>>h1(100,-300,300)", "bbc_vz_tac>-9000"); c.Draw("epd_vz>>h2(100,-300,300)", "epd_vz>-9000", "same"); ((TH1*)gROOT->FindObject("h2"))->SetLineColor(kRed);
    ((TH1*)gROOT->FindObject("h1"))->SetTitle("z: BBC TAC (blue), EPD (red);z (cm)");
    cv.cd(3); c.Draw("bbc_vz_det>>h3(100,-300,300)", "abs(bbc_vz_det)<9000"); ((TH1*)gROOT->FindObject("h3"))->SetTitle("old BBC zVertex;z (cm)");
    cv.cd(4); c.Draw("vpd_vz>>h4(100,-300,300)", "vpd_vz>-9000"); ((TH1*)gROOT->FindObject("h4"))->SetTitle("VPD z;z (cm)");
    cv.cd(5)->SetLogz(); c.Draw("bbc_vz_tac:epd_vz>>c1(60,-300,300,60,-300,300)", "bbc_vz_tac>-9000&&epd_vz>-9000", "colz"); ((TH1*)gROOT->FindObject("c1"))->SetTitle("BBC TAC vs EPD;EPD z;BBC TAC z");
    cv.cd(6)->SetLogz(); c.Draw("bbc_vz_tac:bbc_vz_det>>c2(60,-300,300,60,-300,300)", "bbc_vz_tac>-9000&&abs(bbc_vz_det)<9000", "colz"); ((TH1*)gROOT->FindObject("c2"))->SetTitle("BBC TAC vs old zVertex;old zVertex;BBC TAC z");
    cv.cd(7)->SetLogz(); c.Draw("bbc_vz_tac:vpd_vz>>c3(60,-300,300,60,-300,300)", "bbc_vz_tac>-9000&&vpd_vz>-9000", "colz"); ((TH1*)gROOT->FindObject("c3"))->SetTitle("BBC TAC vs VPD;VPD z;BBC TAC z");
    cv.cd(8)->SetLogz(); c.Draw("epd_vz:vpd_vz>>c4(60,-300,300,60,-300,300)", "epd_vz>-9000&&vpd_vz>-9000", "colz"); ((TH1*)gROOT->FindObject("c4"))->SetTitle("EPD vs VPD;VPD z;EPD z");
    cv.SaveAs(png);
    const char *pr[][3] = {{"bbc_vz_tac","epd_vz","bbc_vz_tac>-9000&&epd_vz>-9000&&abs(bbc_vz_tac-epd_vz)<100"},
                           {"bbc_vz_tac","vpd_vz","bbc_vz_tac>-9000&&vpd_vz>-9000&&abs(bbc_vz_tac-vpd_vz)<100"},
                           {"epd_vz","vpd_vz","epd_vz>-9000&&vpd_vz>-9000&&abs(epd_vz-vpd_vz)<100"},
                           {"bbc_vz_tac","bbc_vz_det","bbc_vz_tac>-9000&&abs(bbc_vz_det)<9000"}};
    for (auto &p : pr) {
        c.Draw(Form("%s-%s>>d(200,-100,100)", p[0], p[1]), p[2], "goff");
        TH1 *d = (TH1*)gROOT->FindObject("d");
        double n = d->GetEntries();
        c.Draw(Form("%s:%s>>cc(120,-300,300,120,-300,300)", p[0], p[1]), p[2], "goff");
        double r = ((TH2*)gROOT->FindObject("cc"))->GetCorrelationFactor();
        printf("  %s - %s: n %.0f mean %.1f rms %.1f cm, corr %.2f\n", p[0], p[1], n, d->GetMean(), d->GetRMS(), r);
    }
    for (const char *v : {"bbc_vz_tac", "epd_vz"}) {
        c.Draw(Form("%s>>m(300,-300,300)", v), Form("%s>-9000", v), "goff");
        TH1 *m = (TH1*)gROOT->FindObject("m"); printf("  %s mean %.1f rms %.1f\n", v, m->GetMean(), m->GetRMS());
    }
}
