// Plots a sys_slices.py output (2026-11-17): for one (level, class, beam) row,
// top: A_N of the tight sample, the loose sample and the slice (loose minus
// tight) vs <x_F>, on the common A_N axes (AN_YLO/AN_YHI/AN_XR as in
// sivers/replot_report.sh); bottom: delta = A_loose - A_tight with sig_diff
// error bars, and the Barlow systematic as a shaded band where |delta| >
// sig_diff. The row's chi^2/ndf is written on the plot. dlim: the lower
// panel's range +-dlim; 0 = from the bins with sig_diff < 0.02 (a few
// low-statistics bins would otherwise flatten the rest; theirs run off scale).
// Usage: root -b -q -l 'plot_sys_slices.C("trig_115_vs_130.txt","out.png","all","all","blue","trigger 1.15 vs 1.3")'
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdio>

void plot_sys_slices(const char *in, const char *out_png, const char *level = "all", const char *cls = "all",
                     const char *beam = "blue", const char *label = "", double ylo = -0.02, double yhi = 0.08, double xr = 0.6,
                     double dlim = 0)
{
    TGraphErrors *gt = new TGraphErrors(), *gl = new TGraphErrors(), *gs = new TGraphErrors();
    TGraphErrors *gd = new TGraphErrors(), *gsys = new TGraphErrors();
    double chi2 = 0, dmax = 0;
    bool disj = false;
    int ndf = 0;
    std::ifstream f(in);
    std::string line;
    while (std::getline(f, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string lv, cl, bm, flag;
        double lo, hi, xf, al, el, at, et, w, as, es, d, sd, pull, sys;
        long nl, nt;
        if (!(ss >> lv >> cl >> bm >> lo >> hi >> xf >> al >> el >> at >> et >> nl >> nt >> w >> as >> es >> d >> sd >> pull >> sys >> flag)) continue;
        if (lv != level || cl != cls || bm != beam) continue;
        disj |= flag == "disjoint";
        const double x = std::fabs(xf), dx = 0.004;
        int n = gt->GetN();
        gt->SetPoint(n, x, at); gt->SetPointError(n, 0, et);
        gl->SetPoint(n, x - dx, al); gl->SetPointError(n, 0, el);
        if (es > 0) { int m = gs->GetN(); gs->SetPoint(m, x + dx, as); gs->SetPointError(m, 0, es); }
        gd->SetPoint(n, x, d); gd->SetPointError(n, 0, sd);
        if (sys > 0) { int m = gsys->GetN(); gsys->SetPoint(m, x, 0); gsys->SetPointError(m, 0.5 * (hi - lo), sys); }
        if (sd > 0) { chi2 += pull * pull; ndf++; }
        if (sd < 0.02) dmax = std::max(dmax, std::fabs(d) + sd);
    }
    if (gt->GetN() == 0) { printf("plot_sys_slices: no %s %s %s bins in %s\n", level, cls, beam, in); return; }

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    TPad *p1 = new TPad("p1", "", 0, 0.38, 1, 1), *p2 = new TPad("p2", "", 0, 0, 1, 0.38);
    p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.13); p1->SetTopMargin(0.08);
    p2->SetTopMargin(0.02); p2->SetBottomMargin(0.25); p2->SetLeftMargin(0.13);
    p1->Draw(); p2->Draw();

    p1->cd(); p1->SetGridy();
    TH1F *fr = new TH1F("fr", Form("%s, %s, %s beam%s%s;;A_{N}", level, cls, beam, label[0] ? ": " : "", label), 10, 0, xr);
    fr->SetMinimum(ylo); fr->SetMaximum(yhi); fr->SetStats(0);
    fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.1);
    fr->Draw("AXIS");
    {
        TLatex tt; tt.SetNDC(); tt.SetTextFont(62); tt.SetTextSize(0.05);
        tt.DrawLatex(0.13, 0.94, Form("%s, %s, %s beam%s%s", level, cls, beam, label[0] ? ": " : "", label));
    }
    gt->SetMarkerStyle(20); gt->SetMarkerColor(kBlue + 1); gt->SetLineColor(kBlue + 1);
    gl->SetMarkerStyle(24); gl->SetMarkerColor(kGray + 2); gl->SetLineColor(kGray + 2);
    gs->SetMarkerStyle(21); gs->SetMarkerColor(kOrange + 1); gs->SetLineColor(kOrange + 1);
    gl->Draw("P same"); gs->Draw("P same"); gt->Draw("P same");
    TLegend *lg = new TLegend(0.16, 0.66, 0.55, 0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
    if (disj)
    {
        lg->AddEntry(gl, "first sample", "p");
        lg->AddEntry(gt, "second sample", "p");
    }
    else
    {
        lg->AddEntry(gt, "tight (nominal comparison)", "p");
        lg->AddEntry(gl, "loose", "p");
        lg->AddEntry(gs, "slice = loose minus tight", "p");
    }
    lg->Draw();

    p2->cd(); p2->SetGridy();
    const double r = dlim > 0 ? dlim : dmax > 0 ? 1.3 * dmax : 0.01;
    TH1F *fd = new TH1F("fd", ";x_{F};A_{loose} #minus A_{tight}", 10, 0, xr);
    fd->SetMinimum(-r); fd->SetMaximum(r); fd->SetStats(0);
    fd->GetXaxis()->SetTitleSize(0.09); fd->GetXaxis()->SetLabelSize(0.075);
    fd->GetYaxis()->SetTitleSize(0.075); fd->GetYaxis()->SetLabelSize(0.065); fd->GetYaxis()->SetTitleOffset(0.8); fd->GetYaxis()->SetNdivisions(505);
    if (disj) fd->GetYaxis()->SetTitle("A_{1} #minus A_{2}");
    fd->Draw("AXIS");
    gsys->SetFillColor(kOrange - 9); gsys->SetFillStyle(1001); gsys->SetLineColor(kOrange - 9);
    if (gsys->GetN() && !disj) gsys->Draw("2 same");
    TLine *l0 = new TLine(0, 0, xr, 0); l0->SetLineStyle(2); l0->Draw();
    gd->SetMarkerStyle(20); gd->Draw("P same");
    TLatex t; t.SetNDC(); t.SetTextSize(0.07); t.SetTextFont(42);
    t.DrawLatex(0.16, 0.86, disj ? Form("#chi^{2}/ndf = %.1f/%d (independent samples)", chi2, ndf)
                                 : Form("#chi^{2}/ndf = %.1f/%d;  band: Barlow sys where |#Delta| > #sigma_{diff}", chi2, ndf));
    c->SaveAs(out_png);
}
