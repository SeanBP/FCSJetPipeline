// Final A_N plots from a combine_sys.py table (Sean 2026-11-17: statistical
// only; systematic only, the total and each component; combined), one row
// (level, class). Blue beam at +x_F, yellow at -x_F (as stored by AnAsymmetry), as plot_an.C; every plot
// on the common A_N axes (AN_YLO / AN_YHI / AN_XR of smoke_report.sh:
// -0.02, 0.08, 0.6 by default).
// Writes
//   <prefix>_stat.png            A_N with statistical bars only
//   <prefix>_sys.png             A_N with the total systematic boxes only
//   <prefix>_sys_components.png  each systematic component and the total vs x_F
//   <prefix>_combined.png        statistical bars + systematic boxes (box
//                                width: the x_F uncertainty, xF_err)
//   <prefix>_total.png           bars = statistical (+) systematic in quadrature
// The polarization scale uncertainty is quoted separately (a note on the plots).
// Usage: root -b -q 'plot_an_final.C+("final_jets.txt","final_all_all","all","all","EM + JP union, all photon classes")'
#include "TBox.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct FinalRow { int beam; double xf, an, stat, sys, tot, xferr; std::vector<double> comp; };

void plot_an_final(const char *table, const char *prefix, const char *level = "all", const char *cls = "all", const char *title = "",
                   double ylo = -0.02, double yhi = 0.08, double xr = 0.6, const char *pol_note = "polarization scale uncertainty not included")
{
    gStyle->SetOptStat(0);
    std::ifstream in(table);
    if (!in) { printf("plot_an_final: cannot open %s\n", table); return; }
    std::vector<std::string> cname;
    std::vector<FinalRow> rows;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.compare(0, 8, "# level ") == 0)
        {
            std::istringstream h(line.substr(2));
            std::string w;
            int i = 0;
            while (h >> w) { if (i++ >= 11) cname.push_back(w.substr(4)); }  // sys_<name>
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        std::istringstream s(line);
        std::string l, c, b;
        double lo, hi;
        FinalRow r;
        if (!(s >> l >> c >> b >> lo >> hi >> r.xf >> r.an >> r.stat >> r.sys >> r.tot >> r.xferr)) continue;
        if (l != level || c != cls) continue;
        r.beam = b == "yellow" ? 1 : 0;
        double v;
        while (s >> v) r.comp.push_back(v);
        rows.push_back(r);
    }
    if (rows.empty()) { printf("plot_an_final: no %s %s rows in %s\n", level, cls, table); return; }
    const int ncomp = cname.size();
    const int col[2] = {kBlue + 1, kOrange + 1}, mk[2] = {20, 21}, bcol[2] = {kBlue - 9, kOrange - 9};
    const int ccol[8] = {kRed + 1, kGreen + 2, kMagenta + 1, kCyan + 2, kViolet + 1, kGray + 2, kSpring - 6, kPink + 6};
    const TString head = Form("Run 22 FCS EM-jets, p^{#uparrow}+p #sqrt{s} = 508 GeV%s%s", TString(title) != "" ? ", " : "", title);

    auto frame = [&](const char *name, double y0, double y1, const char *ytit)
    {
        TH1D *fr = new TH1D(name, Form(";x_{F};%s", ytit), 100, -xr, xr);
        fr->SetMinimum(y0); fr->SetMaximum(y1);
        fr->GetYaxis()->SetTitleOffset(1.4);
        fr->Draw();
        TLine *z = new TLine(-xr, 0, xr, 0); z->SetLineStyle(2); z->Draw();
    };
    auto header = [&](const char *what)
    {
        TLatex t; t.SetNDC(); t.SetTextSize(0.033);
        t.DrawLatex(0.15, 0.86, head);
        t.DrawLatex(0.15, 0.81, what);
        t.SetTextSize(0.026); t.SetTextColor(kGray + 2);
        t.DrawLatex(0.15, 0.17, pol_note);
        t.SetTextColor(col[0]); t.DrawLatex(0.62, 0.81, "blue beam (x_{F} > 0)");
        t.SetTextColor(col[1]); t.DrawLatex(0.62, 0.77, "yellow beam (x_{F} < 0)");
    };
    // points with bars from err(r) (0 = none) and optional boxes from box(r)
    auto points = [&](int mode)
    {
        // mode 0 stat, 1 sys boxes only, 2 stat + boxes, 3 total
        for (int b = 0; b < 2; b++)
        {
            TGraphAsymmErrors *g = new TGraphAsymmErrors();
            for (size_t i = 0; i < rows.size(); i++)
            {
                const FinalRow &r = rows[i];
                if (r.beam != b) continue;
                const double x = r.xf;  // AnAsymmetry tables: yellow already at -x_F
                if (mode >= 1 && mode <= 2)
                {
                    const double w = std::max(std::fabs(r.xferr), 0.006);  // the x_F uncertainty; at least visible
                    // filled box (light beam colour) with an outline, clipped to the frame
                    const double y0 = std::max(r.an - r.sys, ylo), y1 = std::min(r.an + r.sys, yhi);
                    if (y1 > y0 && r.sys > 0)
                    {
                        TBox *bf = new TBox(x - w, y0, x + w, y1);
                        bf->SetFillColor(bcol[b]); bf->SetFillStyle(1001); bf->Draw();
                        TBox *bo = new TBox(x - w, y0, x + w, y1);
                        bo->SetFillStyle(0); bo->SetLineColor(col[b]); bo->Draw();
                    }
                }
                const int n = g->GetN();
                g->SetPoint(n, x, r.an);
                const double e = mode == 0 || mode == 2 ? r.stat : mode == 3 ? r.tot : 0;
                g->SetPointError(n, 0, 0, e, e);
            }
            g->SetMarkerStyle(mk[b]); g->SetMarkerColor(col[b]); g->SetLineColor(col[b]); g->SetMarkerSize(1.0);
            g->Draw("P same");
        }
    };
    const char *what[4] = {"statistical uncertainties only", "systematic uncertainties only (boxes: total)",
                           "bars: statistical; boxes: systematic", "bars: statistical #oplus systematic"};
    const char *suf[4] = {"_stat", "_sys", "_combined", "_total"};
    for (int m = 0; m < 4; m++)
    {
        TCanvas c(Form("c%d", m), "", 900, 700);
        c.SetLeftMargin(0.13);
        frame(Form("fr%d", m), ylo, yhi, "A_{N}");
        points(m);
        header(what[m]);
        c.SaveAs(Form("%s%s.png", prefix, suf[m]));
    }
    // each component and the total vs x_F (both beams, yellow at -x_F)
    {
        TCanvas c("ccomp", "", 900, 700);
        c.SetLeftMargin(0.13);
        double smax = 0;
        for (size_t i = 0; i < rows.size(); i++) smax = std::max(smax, rows[i].sys);
        frame("frc", 0, smax > 0 ? 1.6 * smax : 0.01, "systematic uncertainty on A_{N}");
        TLegend *lg = new TLegend(0.15, 0.55, 0.55, 0.76);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.028);
        for (int k = -1; k < ncomp; k++)
            for (int b = 0; b < 2; b++)  // one line per beam
            {
                TGraph *g = new TGraph();
                for (size_t i = 0; i < rows.size(); i++)
                {
                    const FinalRow &r = rows[i];
                    if (r.beam != b) continue;
                    const double y = k < 0 ? r.sys : (k < (int)r.comp.size() ? r.comp[k] : 0);
                    g->SetPoint(g->GetN(), r.xf, y);
                }
                if (g->GetN() == 0) continue;
                g->Sort();
                g->SetMarkerStyle(k < 0 ? 20 : 24 + (k % 4)); g->SetMarkerSize(k < 0 ? 1.1 : 0.9);
                const int cc = k < 0 ? kBlack : ccol[k % 8];
                g->SetMarkerColor(cc); g->SetLineColor(cc); g->SetLineWidth(k < 0 ? 2 : 1);
                g->Draw("PL same");
                if (b == 0) lg->AddEntry(g, k < 0 ? "total (quadrature)" : cname[k].c_str(), "pl");
            }
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.033);
        t.DrawLatex(0.15, 0.86, head);
        t.DrawLatex(0.15, 0.81, "systematic components (yellow beam at x_{F} < 0)");
        c.SaveAs(Form("%s_sys_components.png", prefix));
    }
    printf("plot_an_final: %s %s, %zu points, %d components -> %s_{stat,sys,sys_components,combined,total}.png\n", level, cls, rows.size(), ncomp, prefix);
}
