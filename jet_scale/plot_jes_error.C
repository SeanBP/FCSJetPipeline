// JES map error vs energy from jes_map_compare.py error (Sean 2026-11-18):
//   <prefix>_sys.png    systematic per alternative map (noise removed; band:
//                       the spatial spread of the smoothed, noise-subtracted
//                       difference) and their combination (black)
//   <prefix>_total.png  the standard map's statistical error (band: spatial
//                       spread), the combined systematic and the total
//                       stat (+) sys that AnHistograms applies (_vjesup / _vjesdn)
// All in percent of F(E), vs jet energy (linear axis, or log with logx = true).
// Usage: root -b -q 'plot_jes_error.C+("err.txt","jes_error")'
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void plot_jes_error(const char *table, const char *prefix = "jes_error", bool logx = false)
{
    gStyle->SetOptStat(0);
    std::ifstream in(table);
    if (!in) { printf("plot_jes_error: cannot open %s\n", table); return; }
    std::vector<std::string> head, notes;
    std::vector<std::vector<double> > rows;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.compare(0, 8, "# E_GeV ") == 0)
        {
            std::istringstream h(line.substr(2));
            std::string w;
            while (h >> w) head.push_back(w);
            continue;
        }
        if (line.empty() || line[0] == '#') { notes.push_back(line); continue; }
        std::istringstream s(line);
        std::vector<double> r;
        double v;
        while (s >> v) r.push_back(v);
        if (r.size() == head.size()) rows.push_back(r);
    }
    if (rows.empty() || head.size() < 5) { printf("plot_jes_error: no rows in %s\n", table); return; }
    auto col = [&](const std::string &n) { for (size_t i = 0; i < head.size(); i++) if (head[i] == n) return (int)i; return -1; };
    std::vector<std::string> alts;
    for (size_t i = 5; i < head.size(); i++)
        if (head[i].size() > 4 && head[i].compare(head[i].size() - 4, 4, "_sys") == 0) alts.push_back(head[i].substr(0, head[i].size() - 4));
    const int n = rows.size();
    const double emin = logx ? rows.front()[0] * 0.8 : 0, emax = rows.back()[0] * (logx ? 1.25 : 1.04);
    // the standard map's sample, from the first header line ("... standard map <path> (sample <name>, ...")
    TString sample = "standard map";
    for (const auto &nt : notes)
    {
        const size_t p = nt.find("(sample ");
        if (p != std::string::npos) { sample = nt.substr(p + 8, nt.find(',', p) - p - 8).c_str(); break; }
    }
    const int acol[8] = {kRed + 1, kGreen + 2, kMagenta + 1, kCyan + 2, kOrange + 7, kViolet + 1, kSpring - 6, kPink + 6};
    auto frame = [&](const char *name, double ymax, const char *ytit)
    {
        TH1D *fr = new TH1D(name, Form(";jet energy (GeV);%s", ytit), 100, emin, emax);
        fr->SetMinimum(0); fr->SetMaximum(ymax);
        if (logx) { fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent(); }
        fr->GetYaxis()->SetTitleOffset(1.2);
        fr->Draw();
    };
    // systematic per alternative
    {
        TCanvas c("cs", "", 1000, 650);
        if (logx) c.SetLogx();
        c.SetLeftMargin(0.09); c.SetRightMargin(0.03); c.SetTopMargin(0.14);
        double ymax = 0;
        for (const auto &r : rows)
        {
            ymax = std::max(ymax, 100 * r[col("sys")]);
            for (const auto &a : alts) ymax = std::max(ymax, 100 * (r[col(a + "_loc_mean")] + r[col(a + "_loc_sd")]));
        }
        frame("frs", 1.25 * std::max(ymax, 0.5), "systematic error of F(E) (%)");
        TLegend *lg = new TLegend(0.55, 0.62, 0.96, 0.85);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
        for (size_t k = 0; k < alts.size(); k++)
        {
            const int cs = col(alts[k] + "_sys"), cm = col(alts[k] + "_loc_mean"), cd = col(alts[k] + "_loc_sd");
            TGraphAsymmErrors *band = new TGraphAsymmErrors(n);
            TGraph *g = new TGraph(n);
            for (int i = 0; i < n; i++)
            {
                const double m = 100 * rows[i][cm], sd = 100 * rows[i][cd];
                band->SetPoint(i, rows[i][0], m);
                band->SetPointError(i, 0, 0, std::min(sd, m), sd);
                g->SetPoint(i, rows[i][0], 100 * rows[i][cs]);
            }
            band->SetFillColorAlpha(acol[k % 8], 0.15); band->SetLineColor(acol[k % 8]); band->Draw("3 same");
            g->SetLineColor(acol[k % 8]); g->SetMarkerColor(acol[k % 8]); g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->SetLineWidth(2);
            g->Draw("PL same");
            lg->AddEntry(g, alts[k].c_str(), "pl");
        }
        if (alts.size() > 1)  // one alternative: it is the systematic
        {
            TGraph *gc = new TGraph(n);
            for (int i = 0; i < n; i++) gc->SetPoint(i, rows[i][0], 100 * rows[i][col("sys")]);
            gc->SetLineWidth(3); gc->SetMarkerStyle(20); gc->Draw("PL same");
            lg->AddEntry(gc, "combined systematic", "pl");
        }
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.033);
        t.DrawLatex(0.09, 0.93, "Systematic error from tune / generator choice vs energy (vs. " + sample + ")");
        t.SetTextSize(0.026); t.SetTextColor(kGray + 2);
        t.DrawLatex(0.09, 0.89, "statistical noise of both maps removed; lines: RMS over the grid; bands: spatial std. dev.");
        c.SaveAs(Form("%s_sys.png", prefix));
    }
    // statistical, systematic, total
    {
        TCanvas c("ct", "", 1000, 650);
        if (logx) c.SetLogx();
        c.SetLeftMargin(0.09); c.SetRightMargin(0.03); c.SetTopMargin(0.14);
        double ymax = 0;
        for (const auto &r : rows) ymax = std::max(ymax, 100 * std::max(r[col("total")], r[col("stat")] + r[col("stat_sp")]));
        frame("frt", 1.25 * std::max(ymax, 0.5), "error of F(E) (%)");
        TGraphAsymmErrors *sb = new TGraphAsymmErrors(n);
        TGraph *gs = new TGraph(n), *gy = new TGraph(n), *gt = new TGraph(n);
        for (int i = 0; i < n; i++)
        {
            const double st = 100 * rows[i][col("stat")], sp = 100 * rows[i][col("stat_sp")];
            sb->SetPoint(i, rows[i][0], st); sb->SetPointError(i, 0, 0, std::min(sp, st), sp);
            gs->SetPoint(i, rows[i][0], st);
            gy->SetPoint(i, rows[i][0], 100 * rows[i][col("sys")]);
            gt->SetPoint(i, rows[i][0], 100 * rows[i][col("total")]);
        }
        sb->SetFillColorAlpha(kBlue + 1, 0.15); sb->Draw("3 same");
        gs->SetLineColor(kBlue + 1); gs->SetMarkerColor(kBlue + 1); gs->SetMarkerStyle(20); gs->SetMarkerSize(0.8); gs->SetLineWidth(2); gs->Draw("PL same");
        gy->SetLineColor(kBlack); gy->SetMarkerStyle(21); gy->SetMarkerSize(0.8); gy->SetLineWidth(2); gy->Draw("PL same");
        gt->SetLineColor(kRed + 1); gt->SetMarkerColor(kRed + 1); gt->SetMarkerStyle(22); gt->SetLineWidth(3); gt->SetLineStyle(2); gt->Draw("PL same");
        TLegend *lg = new TLegend(0.40, 0.66, 0.96, 0.85);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
        lg->AddEntry(gs, "statistical, standard map (band: spatial spread)", "pl");
        lg->AddEntry(gy, "systematic (combined)", "pl");
        lg->AddEntry(gt, "total, applied as #pm#delta(E)", "pl");
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.033);
        t.DrawLatex(0.09, 0.93, "JES map error vs energy: statistical #oplus systematic (" + sample + ")");
        c.SaveAs(Form("%s_total.png", prefix));
    }
    printf("plot_jes_error: %d energies, %zu alternatives -> %s_{sys,total}.png\n", n, alts.size(), prefix);
}
