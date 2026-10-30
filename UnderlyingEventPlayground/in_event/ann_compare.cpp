// The eta-weighted annulus estimate (in_event_ue v3 outputs) against the cone
// estimates, per jet-energy bin:
//   top     simulation closure, estimate / truth target, from truth photons on
//           live cells: annulus whole ring, transverse, away side; 180 deg cone
//           and the +-90 deg pair (area scaling, sides >= 30% live)
//   middle  per-jet fluctuation in data: RMS / mean of the estimate for the
//           180 cone, the +-90 pair, the transverse annulus, the whole ring
//   bottom  away / transverse annulus ratio of the means: data reco,
//           simulation reco, simulation truth photons on live cells
// Simulation means use the per-jet weights. Also prints the tables.
//
// Build and run from this folder:
//   g++ -m64 -O2 ann_compare.cpp -o ann_compare `root-config --cflags --libs` -std=c++11
//   ./ann_compare "<data files>" "<sim files>" <out.png>

#include <TROOT.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TColor.h>

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

static const int kNB = 4, kNC = 5, kNA = 6;  // annulus: ring, transverse, away; with, then without the other-jet exclusion
static const double kE[kNB + 1] = {10, 20, 35, 60, 120};
static const double kMinLive = 0.3;
static const double kMaxFill = 1 / 0.7;  // annulus valid if its covered strips carry >= 70% of the jet template

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}
struct Acc
{
    double w, wx, wxx, w2;
    Acc() : w(0), wx(0), wxx(0), w2(0) {}
    void Add(double x, double ww) { w += ww; wx += ww * x; wxx += ww * x * x; w2 += ww * ww; }
    double Mean() const { return w > 0 ? wx / w : NAN; }
    double Rms() const { const double m = Mean(); return w > 0 ? std::sqrt(std::max(0.0, wxx / w - m * m)) : NAN; }
    double Err() const { return w > 0 ? Rms() * std::sqrt(w2) / w : NAN; }
};
// ratio of means with the error of a ratio of paired sums (per-jet covariance)
struct Ratio
{
    double w, w2, sx, sy, sxx, syy, sxy;
    Ratio() : w(0), w2(0), sx(0), sy(0), sxx(0), syy(0), sxy(0) {}
    void Add(double x, double y, double ww) { w += ww; w2 += ww * ww; sx += ww * x; sy += ww * y; sxx += ww * x * x; syy += ww * y * y; sxy += ww * x * y; }
    double Val() const { return sy != 0 ? sx / sy : NAN; }
    double Err() const
    {
        if (!(w > 0) || sy == 0) return NAN;
        const double mx = sx / w, my = sy / w, r = mx / my;
        const double vx = sxx / w - mx * mx, vy = syy / w - my * my, cxy = sxy / w - mx * my;
        const double neff = w * w / w2;
        return std::fabs(r) * std::sqrt(std::max(0.0, vx / (mx * mx) + vy / (my * my) - 2 * cxy / (mx * my)) / neff);
    }
};

enum { kRingX, kTransX, kAwayX, kRing, kTrans, kAway, kC180, kC90, kNEst };
static const char *kName[kNEst] = {"ring, other jets excl.", "transverse, other jets excl.", "away, other jets excl.",
                                   "annulus, whole ring", "annulus, transverse", "annulus, away side", "180#circ cone", "#pm90#circ cone pair"};

struct Res
{
    Ratio clos[kNB][kNEst];   // truth estimate / truth target (simulation)
    Acc est[kNB][kNEst];      // reco estimate
    Ratio away_trans[2][kNB], away_trans_truth[2][kNB];  // [0] other jets excluded, [1] not
};

static void Fill(const std::vector<std::string> &files, bool sim, Res &r)
{
    TChain ch("jets");
    for (size_t i = 0; i < files.size(); i++) ch.Add(files[i].c_str());
    Float_t E, w, a_jet, e[kNC], a[kNC], tf[kNC], tla[kNC], an_e[kNA], an_fill[kNA], an_tl[kNA], an_tf[kNA];
    ch.SetBranchAddress("jet_E", &E); ch.SetBranchAddress("weight", &w); ch.SetBranchAddress("a_jet", &a_jet);
    ch.SetBranchAddress("e", e); ch.SetBranchAddress("a", a); ch.SetBranchAddress("tf", tf); ch.SetBranchAddress("tla", tla);
    ch.SetBranchAddress("an_e", an_e); ch.SetBranchAddress("an_fill", an_fill); ch.SetBranchAddress("an_tl", an_tl); ch.SetBranchAddress("an_tf", an_tf);
    const Long64_t n = ch.GetEntries();
    for (Long64_t i = 0; i < n; i++)
    {
        ch.GetEntry(i);
        int b = -1;
        for (int q = 0; q < kNB; q++) if (E >= kE[q] && E < kE[q + 1]) b = q;
        if (b < 0) continue;
        const double ww = sim ? w : 1.0;
        for (int v = 0; v < kNA; v++)
            if (an_fill[v] > 0 && an_fill[v] <= kMaxFill)
            {
                r.est[b][v].Add(an_e[v], ww);
                if (sim) r.clos[b][v].Add(an_tl[v], an_tf[v], ww);
            }
        for (int x = 0; x < 2; x++)
        {
            const int it = x ? kTrans : kTransX, ia = x ? kAway : kAwayX;
            if (an_fill[it] > 0 && an_fill[it] <= kMaxFill && an_fill[ia] > 0 && an_fill[ia] <= kMaxFill)
            {
                r.away_trans[x][b].Add(an_e[ia], an_e[it], ww);
                if (sim) r.away_trans_truth[x][b].Add(an_tl[ia], an_tl[it], ww);
            }
        }
        if (a[0] >= kMinLive)
        {
            r.est[b][kC180].Add(e[0] * a_jet / a[0], ww);
            if (sim) r.clos[b][kC180].Add(tla[0], tf[0], ww);
        }
        double s = 0, st = 0, sf = 0; int ns = 0;
        for (int k = 3; k <= 4; k++)
            if (a[k] >= kMinLive) { s += e[k] * a_jet / a[k]; st += tla[k]; sf += tf[k]; ns++; }
        if (ns > 0)
        {
            r.est[b][kC90].Add(s / ns, ww);
            if (sim) r.clos[b][kC90].Add(st / ns, sf / ns, ww);
        }
    }
    std::printf("%s: %lld jets\n", sim ? "simulation" : "data", n);
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 4) { std::fprintf(stderr, "usage: ./ann_compare \"<data files>\" \"<sim files>\" <out.png>\n"); return 1; }
    Res D, S;
    Fill(Split(argv[1]), false, D);
    Fill(Split(argv[2]), true, S);

    std::printf("\nUE estimate per jet (GeV), data | sim; per-jet RMS/mean in data\n");
    for (int b = 0; b < kNB; b++)
    {
        std::printf("%3.0f-%-3.0f", kE[b], kE[b + 1]);
        for (int k = 0; k < kNEst; k++) std::printf("  %5.2f|%5.2f (%.2f)", D.est[b][k].Mean(), S.est[b][k].Mean(), D.est[b][k].Rms() / D.est[b][k].Mean());
        std::printf("\n");
    }
    std::printf("columns: ring, transverse, away (other jets excluded); the same not excluded; 180 cone, +-90 pair\n\nsimulation closure (truth estimate / truth target)\n");
    for (int b = 0; b < kNB; b++)
    {
        std::printf("%3.0f-%-3.0f", kE[b], kE[b + 1]);
        for (int k = 0; k < kNEst; k++) std::printf("  %.3f+-%.3f", S.clos[b][k].Val(), S.clos[b][k].Err());
        std::printf("\n");
    }
    for (int x = 0; x < 2; x++)
    {
        std::printf("\naway / transverse annulus, other jets %s: data reco | sim reco | sim truth on live cells\n", x ? "not excluded" : "excluded");
        for (int b = 0; b < kNB; b++)
            std::printf("%3.0f-%-3.0f  %.3f+-%.3f | %.3f+-%.3f | %.3f+-%.3f\n", kE[b], kE[b + 1], D.away_trans[x][b].Val(), D.away_trans[x][b].Err(),
                        S.away_trans[x][b].Val(), S.away_trans[x][b].Err(), S.away_trans_truth[x][b].Val(), S.away_trans_truth[x][b].Err());
    }

    gStyle->SetOptStat(0);
    const int col[kNEst] = {TColor::GetColor("#1c4f8a"), TColor::GetColor("#2f9e44"), TColor::GetColor("#d9480f"),
                            TColor::GetColor("#1c4f8a"), TColor::GetColor("#2f9e44"), TColor::GetColor("#d9480f"), TColor::GetColor("#868e96"), TColor::GetColor("#9c36b5")};
    const int mk[kNEst] = {24, 25, 26, 20, 21, 22, 27, 33};
    const std::vector<int> show[2] = {{kRing, kTrans, kAway, kC180, kC90, kRingX}, {kRing, kTrans, kC180, kC90}};
    TCanvas cv("cv", "", 1800, 600);
    cv.Divide(3, 1);
    std::vector<TObject *> keep;
    for (int p = 0; p < 3; p++)
    {
        cv.cd(p + 1);
        gPad->SetMargin(0.14, 0.03, 0.12, 0.08); gPad->SetLogx();
        TH1D *fr = new TH1D(Form("fr%d", p), "", 1, 10, 120); keep.push_back(fr);
        const char *yt[3] = {"closure: truth estimate / truth target", "data: per-jet RMS / mean of the estimate", "away / transverse annulus"};
        const double lo[3] = {0.8, 0.0, 0.9}, hi[3] = {1.25, 3.0, 1.6};
        fr->SetTitle(Form(";jet E (GeV);%s", yt[p])); fr->SetMinimum(lo[p]); fr->SetMaximum(hi[p]);
        fr->GetXaxis()->SetMoreLogLabels(); fr->GetYaxis()->SetTitleOffset(1.3);
        fr->Draw();
        TLegend *lg = new TLegend(0.15, 0.70, 0.97, 0.91); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.031);
        keep.push_back(lg);
        if (p < 2)
        {
            if (p == 0) { TLine *l = new TLine(10, 1, 120, 1); l->SetLineStyle(2); l->Draw(); keep.push_back(l); }
            for (size_t ik = 0; ik < show[p].size(); ik++)
            {
                const int k = show[p][ik];
                TGraphErrors *g = new TGraphErrors();
                for (int b = 0; b < kNB; b++)
                {
                    const double x = std::sqrt(kE[b] * kE[b + 1]) * (1 + 0.025 * ((int)ik - 2));
                    const double y = p == 0 ? S.clos[b][k].Val() : D.est[b][k].Rms() / D.est[b][k].Mean();
                    const double ey = p == 0 ? S.clos[b][k].Err() : 0;
                    g->SetPoint(b, x, y); g->SetPointError(b, 0, ey);
                }
                g->SetMarkerColor(col[k]); g->SetLineColor(col[k]); g->SetMarkerStyle(mk[k]); g->SetMarkerSize(1.3);
                if (k < 3) g->SetLineStyle(2);
                g->Draw("PL SAME"); keep.push_back(g);
                lg->AddEntry(g, kName[k], "pl");
            }
            if (p == 0) lg->SetNColumns(2);
            if (p == 1) lg->SetNColumns(2);
        }
        else
        {
            TLine *l = new TLine(10, 1, 120, 1); l->SetLineStyle(2); l->Draw(); keep.push_back(l);
            const char *nm[3] = {"data, reco", "sim, reco", "sim, truth on live cells"};
            const int cc[3] = {TColor::GetColor("#222222"), TColor::GetColor("#d9480f"), TColor::GetColor("#1c6bb0")};
            lg->SetNColumns(2);
            for (int s = 0; s < 3; s++)
                for (int x = 1; x >= 0; x--)
                {
                    TGraphErrors *g = new TGraphErrors();
                    for (int b = 0; b < kNB; b++)
                    {
                        const Ratio &r = s == 0 ? D.away_trans[x][b] : s == 1 ? S.away_trans[x][b] : S.away_trans_truth[x][b];
                        g->SetPoint(b, std::sqrt(kE[b] * kE[b + 1]) * (1 + 0.025 * (2 * s + x - 2.5)), r.Val()); g->SetPointError(b, 0, r.Err());
                    }
                    g->SetMarkerColor(cc[s]); g->SetLineColor(cc[s]); g->SetMarkerStyle(x ? 20 + s : 24 + s); g->SetMarkerSize(1.3);
                    if (!x) g->SetLineStyle(2);
                    g->Draw("PL SAME"); keep.push_back(g);
                    lg->AddEntry(g, Form("%s%s", nm[s], x ? "" : ", other jets excl."), "pl");
                }
        }
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.04);
        const char *ti[3] = {"Truth closure of the coverage / #eta weighting", "Fluctuation per jet (data)", "Away-side excess (recoil + correlations)"};
        t.DrawLatex(0.14, 0.94, ti[p]);
    }
    cv.SaveAs(argv[3]);
    return 0;
}
