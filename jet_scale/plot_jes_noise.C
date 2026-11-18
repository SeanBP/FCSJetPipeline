// Statistical error of the JES map from the subset maps (jes_map_compare.py
// noise output <prefix>.txt; Sean 2026-11-18, after his earlier PYTHIA 8 study:
// four 2M maps -> ~2% per map, ~1% expected for 8M):
//   <out>_stat_vs_E.png     expected statistical error of the full-sample map
//                           vs jet energy: spatial mean over the grid nodes,
//                           band = spatial std. dev.; thin dashed lines: the
//                           same from each division k alone (they agree if the
//                           error scales as 1/sqrt(N))
//   <out>_stat_vs_N.png     average statistical error vs number of simulated
//                           events: the subset maps measured at N/k (points,
//                           bars = spatial std. dev.), the 1/sqrt(N) curve
//                           through them with the spread as a band, and the
//                           full map at N (extrapolated)
// Errors are relative, in percent of F(E). Energies averaged for the second
// plot: [emin, emax] GeV (spatial mean per energy, then the mean over energy).
// Usage: root -b -q 'plot_jes_noise.C+("noise.txt", 9.86, "jes_noise", 10, 150)'
//   9.86: events of the full sample in millions (the subsets have N/k)
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

struct NoiseRow { int k, dof; double e, sub_mean, sub_sd; };

void plot_jes_noise(const char *table, double nfull_m, const char *out = "jes_noise", double emin = 0, double emax = 1e9, bool logx = false)
{
    gStyle->SetOptStat(0);
    std::ifstream in(table);
    if (!in) { printf("plot_jes_noise: cannot open %s\n", table); return; }
    std::vector<std::string> head;
    std::vector<NoiseRow> rows;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.compare(0, 8, "# k E_GeV") == 0 || line.compare(0, 9, "# k E_GeV") == 0)
        {
            std::istringstream h(line.substr(2));
            std::string w;
            while (h >> w) head.push_back(w);
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        std::istringstream s(line);
        std::vector<double> v;
        double x;
        while (s >> x) v.push_back(x);
        if (v.size() != head.size() || head.size() < 10) continue;
        auto c = [&](const char *n) -> double { for (size_t i = 0; i < head.size(); i++) if (head[i] == n) return v[i]; return std::nan(""); };
        NoiseRow r;
        r.k = (int)c("k"); r.e = c("E_GeV"); r.dof = (int)c("dof"); r.sub_mean = c("sub_mean"); r.sub_sd = c("sub_sd");
        if (r.sub_mean == r.sub_mean) rows.push_back(r);
    }
    if (rows.empty()) { printf("plot_jes_noise: no rows with sub_mean in %s (rerun jes_map_compare.py noise)\n", table); return; }
    std::vector<int> ks;
    std::vector<double> es;
    for (const auto &r : rows)
    {
        if (std::find(ks.begin(), ks.end(), r.k) == ks.end()) ks.push_back(r.k);
        if (std::find(es.begin(), es.end(), r.e) == es.end()) es.push_back(r.e);
    }
    std::sort(ks.begin(), ks.end());
    std::sort(es.begin(), es.end());
    // full-map error per energy: each division's subset error / sqrt(k), weighted by its degrees of freedom
    std::map<double, double> fm, fs;
    std::map<int, std::map<double, double> > perk;
    for (double e : es)
    {
        double w = 0, m = 0, sd = 0;
        for (const auto &r : rows)
            if (r.e == e)
            {
                w += r.dof; m += r.dof * r.sub_mean / std::sqrt((double)r.k); sd += r.dof * r.sub_sd / std::sqrt((double)r.k);
                perk[r.k][e] = r.sub_mean / std::sqrt((double)r.k);
            }
        fm[e] = m / w; fs[e] = sd / w;
    }
    const TString nfull = Form("%.3g M events", nfull_m);
    // 1. statistical error vs energy
    {
        TCanvas c("c1", "", 1000, 620);
        c.SetLeftMargin(0.09); c.SetRightMargin(0.03); c.SetTopMargin(0.13);
        if (logx) c.SetLogx();
        double ymax = 0;
        for (double e : es) ymax = std::max(ymax, 100 * (fm[e] + fs[e]));
        TH1D *fr = new TH1D("fr1", ";jet energy (GeV);statistical error in F(E) (%)", 100, logx ? es.front() * 0.8 : 0, es.back() * (logx ? 1.25 : 1.04));
        fr->SetMinimum(0); fr->SetMaximum(1.15 * ymax);
        if (logx) { fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent(); }
        fr->Draw();
        const int n = es.size();
        TGraphAsymmErrors *band = new TGraphAsymmErrors(n);
        TGraph *g = new TGraph(n);
        for (int i = 0; i < n; i++)
        {
            const double m = 100 * fm[es[i]], sd = 100 * fs[es[i]];
            band->SetPoint(i, es[i], m); band->SetPointError(i, 0, 0, std::min(sd, m), sd);
            g->SetPoint(i, es[i], m);
        }
        band->SetFillColorAlpha(kAzure + 1, 0.25); band->Draw("3 same");
        TLegend *lg = new TLegend(0.47, 0.64, 0.96, 0.85);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.029);
        const int kc[6] = {kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1};
        for (size_t q = 0; q < ks.size(); q++)
        {
            TGraph *gk = new TGraph();
            for (double e : es) if (perk[ks[q]].count(e)) gk->SetPoint(gk->GetN(), e, 100 * perk[ks[q]][e]);
            gk->SetLineColor(kc[q % 6]); gk->SetLineStyle(2); gk->Draw("L same");
            lg->AddEntry(gk, Form("from 1/%d subsets (%.3g M each) / #sqrt{%d}", ks[q], nfull_m / ks[q], ks[q]), "l");
        }
        g->SetLineColor(kAzure + 2); g->SetMarkerColor(kAzure + 2); g->SetMarkerStyle(20); g->SetMarkerSize(0.7); g->SetLineWidth(3);
        g->Draw("PL same");
        lg->AddEntry(g, "standard map (" + nfull + "): spatial mean", "pl");
        lg->AddEntry(band, "spatial std. dev.", "f");
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.034);
        t.DrawLatex(0.09, 0.93, "Expected statistical error of the standard JES map vs energy (" + nfull + ")");
        c.SaveAs(Form("%s_stat_vs_E.png", out));
    }
    // 2. average statistical error vs number of simulated events
    {
        std::vector<double> pn, pm, ps;
        double wsum = 0, full_m = 0, full_s = 0;
        for (int k : ks)
        {
            double m = 0, sd = 0;
            int ne = 0, dof = 0;
            for (const auto &r : rows)
                if (r.k == k && r.e >= emin && r.e <= emax) { m += r.sub_mean; sd += r.sub_sd; ne++; dof = r.dof; }
            if (!ne) continue;
            m /= ne; sd /= ne;
            pn.push_back(nfull_m / k); pm.push_back(100 * m); ps.push_back(100 * sd);
            wsum += dof; full_m += dof * 100 * m / std::sqrt((double)k); full_s += dof * 100 * sd / std::sqrt((double)k);
        }
        if (pn.empty()) { printf("plot_jes_noise: no energies in [%g, %g]\n", emin, emax); return; }
        full_m /= wsum; full_s /= wsum;
        TCanvas c("c2", "", 850, 650);
        c.SetLeftMargin(0.11); c.SetRightMargin(0.03); c.SetTopMargin(0.12);
        double ymax = 0;
        for (size_t i = 0; i < pn.size(); i++) ymax = std::max(ymax, pm[i] + ps[i]);
        TH1D *fr = new TH1D("fr2", ";number of simulated events (millions);average statistical error on the JES map (%)", 100, 0, nfull_m * 1.04);
        fr->SetMinimum(0); fr->SetMaximum(std::min(1.5 * ymax, 4 * ymax));
        fr->GetYaxis()->SetTitleOffset(1.3);
        fr->Draw();
        const int nc = 400;
        TGraphAsymmErrors *band = new TGraphAsymmErrors(nc);
        TGraph *curve = new TGraph(nc);
        for (int i = 0; i < nc; i++)
        {
            const double nn = nfull_m * (0.02 + 0.98 * i / (nc - 1.0)), f = std::sqrt(nfull_m / nn);
            curve->SetPoint(i, nn, full_m * f);
            band->SetPoint(i, nn, full_m * f); band->SetPointError(i, 0, 0, std::min(full_s * f, full_m * f), full_s * f);
        }
        band->SetFillColorAlpha(kAzure + 1, 0.25); band->Draw("3 same");
        curve->SetLineColor(kAzure + 2); curve->SetLineWidth(3); curve->Draw("L same");
        TGraphErrors *meas = new TGraphErrors(pn.size());
        for (size_t i = 0; i < pn.size(); i++) { meas->SetPoint(i, pn[i], pm[i]); meas->SetPointError(i, 0, ps[i]); }
        meas->SetMarkerStyle(21); meas->SetMarkerColor(kRed + 1); meas->SetLineColor(kRed + 1); meas->SetMarkerSize(1.5); meas->SetLineWidth(2);
        meas->Draw("P same");
        TGraphErrors *full = new TGraphErrors(1);
        full->SetPoint(0, nfull_m, full_m); full->SetPointError(0, 0, full_s);
        full->SetMarkerStyle(20); full->SetMarkerSize(1.6); full->SetLineWidth(2);
        full->Draw("P same");
        TLegend *lg = new TLegend(0.40, 0.64, 0.96, 0.86);
        lg->SetBorderSize(1); lg->SetFillColor(kWhite); lg->SetTextSize(0.032);
        lg->AddEntry(curve, "expected error #propto 1/#sqrt{N}", "l");
        lg->AddEntry(band, "expected spatial spread (#pm1#sigma)", "f");
        TString ksz;
        for (size_t i = 0; i < pn.size(); i++) ksz += Form("%s%.3g", i ? ", " : "", pn[i]);
        lg->AddEntry(meas, "subset maps (" + ksz + " M, measured)", "pe");
        lg->AddEntry(full, Form("standard map (%.3g M, extrapolated): %.2f%%", nfull_m, full_m), "pe");
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.033);
        t.DrawLatex(0.11, 0.935, "Expected statistical precision of the JES map vs sample size");
        t.SetTextSize(0.026); t.SetTextColor(kGray + 2);
        t.DrawLatex(0.11, 0.895, Form("bars / band: spatial std. dev.; energies %g-%g GeV averaged", std::max(emin, es.front()), std::min(emax, es.back())));
        c.SaveAs(Form("%s_stat_vs_N.png", out));
        printf("plot_jes_noise: standard map (%.3g M) average statistical error %.2f%% (spatial sd %.2f%%), energies %g-%g GeV\n",
               nfull_m, full_m, full_s, std::max(emin, es.front()), std::min(emax, es.back()));
    }
    printf("plot_jes_noise: %zu divisions, %zu energies -> %s_stat_vs_E.png, %s_stat_vs_N.png\n", ks.size(), es.size(), out, out);
}
