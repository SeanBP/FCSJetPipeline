// Combines ue_truth_decomposition outputs (one per ptHat bin and build) with
// cross-section weights (a sampled jet of a bin stands for
// sigma / trials x n_estimable / n_sampled) and compares the variants:
//   sep_veto, sep (no veto; the pipeline method), nosep (no separation, no veto)
// and, for reference, the jet's own same-event truth cones (the truth UE
// jet_scale subtracts).
//
// Prints, per jet-pT and jet-eta bin, the weighted mean truth energy in the
// UE cone and its parts (proxy jet's truth jet, other truth jets >= 2 GeV,
// 1-2 GeV, below 1 GeV), the method's reco cone energy and the proxy
// events' mean |dphi|; and the same against the band half-width each jet's
// position gives (the window is sized to the 5% area, so wide windows are
// where the band reaches far from 90 deg; the own-event reference is placed
// at the pipeline variant's width). Pages <base>_pt, <base>_eta,
// <base>_width: one panel per part, plus the total and the reco cone, each
// variant a curve (own event dashed). Page <base>_vs_truth: each variant's
// truth cone energy against the truth UE of the same jets (own-event cones),
// as the ratio of the weighted sums and the mean difference, vs pT, eta and
// the variant's own band half-width.
//
// Usage: ./truth_decomp_compare <file list> <output base>
// Build: g++ -m64 -O2 truth_decomp_compare.cpp -o truth_decomp_compare `root-config --cflags --libs` -std=c++11

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static const int kNV = 4; // 0 sep_veto, 1 sep, 2 nosep, 3 own event
static const char *kName[kNV] = {"separation + veto", "separation, no veto (pipeline)", "no separation, no veto", "jet's own event, #pm90#circ cones (truth UE)"};
static const char *kShort[kNV] = {"sep_veto", "sep", "nosep", "own"};
static const char *kColor[kNV] = {"#2a78d6", "#8a5cd1", "#1baf7a", "#3a3a36"};
static const int kVMarker[kNV] = {20, 22, 33, 24};
static const int kNQ = 7; // quantities
static const char *kQName[kNQ] = {"truth cone E, all (GeV)", "proxy jet's own truth jet (GeV)", "other truth jets, p_{T} #geq 2 GeV (GeV)",
                                  "other truth jets, 1-2 GeV (GeV)", "truth jets below 1 GeV (GeV)", "reco cone E (GeV)", "proxy jet |#Delta#phi| (deg)"};

struct Acc
{
    std::vector<double> sw, sw2;
    std::vector<double> s[kNQ], s2[kNQ];
    void Init(int nb)
    {
        sw.assign(nb, 0); sw2.assign(nb, 0);
        for (int q = 0; q < kNQ; q++) { s[q].assign(nb, 0); s2[q].assign(nb, 0); }
    }
    void Add(int b, double w, const double *x)
    {
        sw[b] += w; sw2[b] += w * w;
        for (int q = 0; q < kNQ; q++) { s[q][b] += w * x[q]; s2[q][b] += w * x[q] * x[q]; }
    }
    double Mean(int q, int b) const { return sw[b] > 0 ? s[q][b] / sw[b] : NAN; }
    double Err(int q, int b) const
    {
        if (!(sw[b] > 0)) return 0;
        const double m = Mean(q, b), var = std::max(0.0, s2[q][b] / sw[b] - m * m), neff = sw[b] * sw[b] / sw2[b];
        return neff > 1 ? std::sqrt(var / neff) : 0;
    }
};

// A variant's truth cone energy a against the same jets' truth UE b (own event)
struct RAcc
{
    std::vector<double> sw, sw2, sa, sb, swa2, swb2, swab, sd2;
    void Init(int nb) { sw.assign(nb, 0); sw2 = sa = sb = swa2 = swb2 = swab = sd2 = sw; }
    void Add(int k, double w, double a, double b)
    {
        sw[k] += w; sw2[k] += w * w; sa[k] += w * a; sb[k] += w * b;
        swa2[k] += w * w * a * a; swb2[k] += w * w * b * b; swab[k] += w * w * a * b; sd2[k] += w * (a - b) * (a - b);
    }
    double Ratio(int k) const { return sb[k] > 0 ? sa[k] / sb[k] : NAN; }
    // linearized: the spread of a - r b over the jets
    double RatioErr(int k) const
    {
        const double r = Ratio(k);
        return sb[k] > 0 ? std::sqrt(std::max(0.0, swa2[k] - 2 * r * swab[k] + r * r * swb2[k])) / sb[k] : 0;
    }
    double Diff(int k) const { return sw[k] > 0 ? (sa[k] - sb[k]) / sw[k] : NAN; }
    double DiffErr(int k) const
    {
        if (!(sw[k] > 0)) return 0;
        const double m = Diff(k), var = std::max(0.0, sd2[k] / sw[k] - m * m), neff = sw[k] * sw[k] / sw2[k];
        return neff > 1 ? std::sqrt(var / neff) : 0;
    }
};

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "usage: ./truth_decomp_compare <file list> <output base>" << std::endl; return 1; }
    std::vector<std::string> files;
    {
        std::ifstream in(argv[1]);
        std::string f;
        while (in >> f) files.push_back(f);
    }
    const std::string out = argv[2];
    const int kNVar = 3;
    const char *vtitle[kNVar] = {"jet p_{T} (GeV)", "jet #eta", "band half-width at the jet (deg)"};
    const char *vtag[kNVar] = {"pt", "eta", "width"};
    std::vector<double> edges[kNVar];
    edges[0] = {2, 2.5, 3, 3.5, 4, 5, 6, 8, 12};
    edges[1] = {2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.6};
    edges[2] = {0, 10, 15, 20, 25, 30, 35, 40, 46};
    Acc acc[kNVar][kNV];
    RAcc racc[kNVar][3];
    for (int k = 0; k < kNVar; k++)
    {
        for (int v = 0; v < kNV; v++) acc[k][v].Init((int)edges[k].size() - 1);
        for (int v = 0; v < 3; v++) racc[k][v].Init((int)edges[k].size() - 1);
    }

    long n_rows = 0;
    for (size_t i = 0; i < files.size(); i++)
    {
        TFile f(files[i].c_str());
        TTree *ti = (TTree *)f.Get("info");
        TTree *tj = (TTree *)f.Get("jets");
        if (!ti || !tj) { std::cerr << "skip " << files[i] << std::endl; continue; }
        Double_t sigma, trials;
        Long64_t n_est, n_samp;
        ti->SetBranchAddress("sigma_mb", &sigma);
        ti->SetBranchAddress("trials", &trials);
        ti->SetBranchAddress("n_estimable", &n_est);
        ti->SetBranchAddress("n_sampled", &n_samp);
        ti->GetEntry(0);
        if (!(sigma > 0) || !(trials > 0) || n_samp <= 0) { std::cerr << "no weight in " << files[i] << std::endl; continue; }
        const double w = sigma / trials * (double)n_est / n_samp;
        Int_t variant;
        Float_t pt, eta, rc, ta, tp, t2, t1, ts, dphi, oa, oj, o2, o1, os, hw;
        tj->SetBranchAddress("variant", &variant);
        tj->SetBranchAddress("jet_pt", &pt);
        tj->SetBranchAddress("jet_eta", &eta);
        tj->SetBranchAddress("reco_cone", &rc);
        tj->SetBranchAddress("t_all", &ta);
        tj->SetBranchAddress("t_proxy", &tp);
        tj->SetBranchAddress("t_other2", &t2);
        tj->SetBranchAddress("t_other1", &t1);
        tj->SetBranchAddress("t_soft", &ts);
        tj->SetBranchAddress("proxy_dphi", &dphi);
        tj->SetBranchAddress("half_width", &hw);
        tj->SetBranchAddress("own_all", &oa);
        tj->SetBranchAddress("own_jet", &oj);
        tj->SetBranchAddress("own_other2", &o2);
        tj->SetBranchAddress("own_other1", &o1);
        tj->SetBranchAddress("own_soft", &os);
        for (Long64_t e = 0; e < tj->GetEntries(); e++)
        {
            tj->GetEntry(e);
            n_rows++;
            const double x[kNQ] = {ta, tp, t2, t1, ts, rc, dphi};
            // The own-event reference once per jet: from the variant-0 or -2 rows (one per build)
            const bool own = variant == 0 || variant == 2;
            const double xo[kNQ] = {oa, oj, o2, o1, os, NAN, NAN};
            const double xv[kNVar] = {pt, eta, hw};
            for (int k = 0; k < kNVar; k++)
            {
                const std::vector<double> &ed = edges[k];
                if (xv[k] < ed.front() || xv[k] >= ed.back()) continue;
                const int b = (int)(std::upper_bound(ed.begin(), ed.end(), xv[k]) - ed.begin()) - 1;
                acc[k][variant].Add(b, w, x);
                racc[k][variant].Add(b, w, ta, oa);
                if (variant == 0 && own) acc[k][3].Add(b, w, xo);
            }
        }
    }
    printf("truth_decomp_compare: %zu files, %ld rows\n", files.size(), n_rows);

    for (int k = 0; k < kNVar; k++)
    {
        printf("\n  %s: weighted mean per jet, GeV: truth cone all = proxy jet + other >=2 + other 1-2 + soft | reco cone | proxy |dphi|\n", vtitle[k]);
        for (int b = 0; b + 1 < (int)edges[k].size(); b++)
        {
            printf("  %5.2f-%-5.2f\n", edges[k][b], edges[k][b + 1]);
            for (int v = 0; v < kNV; v++)
            {
                const Acc &a = acc[k][v];
                if (!(a.sw[b] > 0)) continue;
                printf("    %-9s %6.3f = %6.3f + %6.3f + %6.3f + %6.3f", kShort[v], a.Mean(0, b), a.Mean(1, b), a.Mean(2, b), a.Mean(3, b), a.Mean(4, b));
                if (v < 3) printf("  | %6.3f | %5.1f", a.Mean(5, b), a.Mean(6, b));
                printf("   (n_eff %.0f)\n", a.sw[b] * a.sw[b] / a.sw2[b]);
            }
        }
    }

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(5);
    const int ink = TColor::GetColor("#3a3a36");
    std::vector<TObject *> keep;
    for (int k = 0; k < kNVar; k++)
    {
        TCanvas c(Form("c%d", k), "", 2100, 1150);
        c.Divide(4, 2, 0.004, 0.004);
        for (int q = 0; q < kNQ; q++)
        {
            c.cd(q + 1);
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.05); gPad->SetRightMargin(0.03);
            double ymax = 0;
            TGraphErrors *g[kNV];
            for (int v = 0; v < kNV; v++)
            {
                g[v] = new TGraphErrors();
                keep.push_back(g[v]);
                if (v == 3 && q >= 5) continue;
                for (int b = 0; b + 1 < (int)edges[k].size(); b++)
                {
                    const double m = acc[k][v].Mean(q, b);
                    if (!std::isfinite(m)) continue;
                    const double x = 0.5 * (edges[k][b] + edges[k][b + 1]) + (v - 1.5) * 0.012 * (edges[k].back() - edges[k].front());
                    const int n = g[v]->GetN();
                    g[v]->SetPoint(n, x, m);
                    g[v]->SetPointError(n, 0, acc[k][v].Err(q, b));
                    ymax = std::max(ymax, m + acc[k][v].Err(q, b));
                }
            }
            TH1D *fr = new TH1D(Form("fr%d_%d", k, q), "", 1, edges[k].front(), edges[k].back());
            keep.push_back(fr);
            fr->SetMinimum(0); fr->SetMaximum(q == 6 ? 180 : 1.35 * ymax + 1e-3);
            fr->GetXaxis()->SetTitle(vtitle[k]); fr->GetYaxis()->SetTitle(kQName[q]);
            fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleOffset(1.4);
            fr->Draw("AXIS");
            for (int v = 0; v < kNV; v++)
            {
                if (g[v]->GetN() == 0) continue;
                const int col = TColor::GetColor(kColor[v]);
                g[v]->SetMarkerStyle(kVMarker[v]); g[v]->SetMarkerSize(kVMarker[v] == 33 ? 1.6 : 1.2);
                g[v]->SetMarkerColor(col); g[v]->SetLineColor(col); g[v]->SetLineWidth(2);
                if (v == 3) g[v]->SetLineStyle(2);
                g[v]->Draw("PL same");
            }
        }
        c.cd(8);
        TLegend *leg = new TLegend(0.03, 0.25, 0.97, 0.8);
        keep.push_back(leg);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.055); leg->SetTextColor(ink);
        leg->SetHeader("UE cone (R) of the proxy events; simulation, #sigma-weighted");
        TGraphErrors *dummy[kNV];
        for (int v = 0; v < kNV; v++)
        {
            dummy[v] = new TGraphErrors();
            keep.push_back(dummy[v]);
            dummy[v]->SetMarkerStyle(kVMarker[v]); dummy[v]->SetMarkerColor(TColor::GetColor(kColor[v]));
            dummy[v]->SetLineColor(TColor::GetColor(kColor[v])); dummy[v]->SetLineWidth(2);
            if (v == 3) dummy[v]->SetLineStyle(2);
            leg->AddEntry(dummy[v], kName[v], "lp");
        }
        leg->Draw();
        c.SaveAs((out + "_" + vtag[k] + ".png").c_str());
        c.SaveAs((out + "_" + vtag[k] + ".pdf").c_str());
    }

    printf("\n  each variant's truth cone energy vs the same jets' truth UE (own-event cones): ratio of sums | mean difference (GeV)\n");
    for (int k = 0; k < kNVar; k++)
    {
        printf("  %s\n", vtitle[k]);
        for (int b = 0; b + 1 < (int)edges[k].size(); b++)
        {
            printf("  %5.2f-%-5.2f", edges[k][b], edges[k][b + 1]);
            for (int v = 0; v < 3; v++)
            {
                const RAcc &r = racc[k][v];
                if (r.sw[b] > 0) printf("   %-8s %.3f +- %.3f | %+.3f +- %.3f", kShort[v], r.Ratio(b), r.RatioErr(b), r.Diff(b), r.DiffErr(b));
                else printf("   %-8s %31s", kShort[v], "");
            }
            printf("\n");
        }
    }

    TCanvas cr("cr", "", 2100, 1150);
    cr.Divide(4, 2, 0.004, 0.004);
    for (int row = 0; row < 2; row++)
        for (int k = 0; k < kNVar; k++)
        {
            cr.cd(row * 4 + k + 1);
            gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.05); gPad->SetRightMargin(0.03);
            TH1D *fr = new TH1D(Form("frr%d_%d", row, k), "", 1, edges[k].front(), edges[k].back());
            keep.push_back(fr);
            fr->SetMinimum(row == 0 ? 0.5 : -0.8); fr->SetMaximum(row == 0 ? 1.5 : 0.8);
            fr->GetXaxis()->SetTitle(k == 2 ? "band half-width of the method at the jet (deg)" : vtitle[k]);
            fr->GetYaxis()->SetTitle(row == 0 ? "method / truth UE (truth cone E)" : "method #minus truth UE (GeV)");
            fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleOffset(1.4);
            fr->Draw("AXIS");
            TLine *l = new TLine(edges[k].front(), row == 0 ? 1 : 0, edges[k].back(), row == 0 ? 1 : 0);
            keep.push_back(l);
            l->SetLineColor(ink); l->SetLineStyle(2); l->Draw();
            for (int v = 0; v < 3; v++)
            {
                TGraphErrors *g = new TGraphErrors();
                keep.push_back(g);
                for (int b = 0; b + 1 < (int)edges[k].size(); b++)
                {
                    const RAcc &r = racc[k][v];
                    // too few jets for a stable mean (the narrow edge bins of the separation variants)
                    if (!(r.sw[b] > 0) || r.sw[b] * r.sw[b] / r.sw2[b] < 20) continue;
                    const double x = 0.5 * (edges[k][b] + edges[k][b + 1]) + (v - 1) * 0.012 * (edges[k].back() - edges[k].front());
                    const int n = g->GetN();
                    g->SetPoint(n, x, row == 0 ? r.Ratio(b) : r.Diff(b));
                    g->SetPointError(n, 0, row == 0 ? r.RatioErr(b) : r.DiffErr(b));
                }
                if (g->GetN() == 0) continue;
                const int col = TColor::GetColor(kColor[v]);
                g->SetMarkerStyle(kVMarker[v]); g->SetMarkerSize(kVMarker[v] == 33 ? 1.6 : 1.2);
                g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2);
                g->Draw("PL same");
            }
        }
    cr.cd(4);
    TLegend *lr = new TLegend(0.02, 0.2, 0.98, 0.85);
    keep.push_back(lr);
    lr->SetBorderSize(0); lr->SetFillStyle(0); lr->SetTextSize(0.05); lr->SetTextColor(ink);
    lr->SetHeader("truth E in the method's UE cones");
    for (int v = 0; v < 3; v++)
    {
        TGraphErrors *d = new TGraphErrors();
        keep.push_back(d);
        d->SetMarkerStyle(kVMarker[v]); d->SetMarkerColor(TColor::GetColor(kColor[v]));
        d->SetLineColor(TColor::GetColor(kColor[v])); d->SetLineWidth(2);
        lr->AddEntry(d, kName[v], "lp");
    }
    lr->AddEntry((TObject *)0, "truth UE: own event, #pm90#circ cones", "");
    lr->AddEntry((TObject *)0, "same jets; simulation, #sigma-weighted", "");
    lr->Draw();
    cr.SaveAs((out + "_vs_truth.png").c_str());
    cr.SaveAs((out + "_vs_truth.pdf").c_str());
    return 0;
}
