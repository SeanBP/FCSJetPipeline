// Data against simulation for in_event_ue outputs.
//   page 1  UE estimate per jet against jet energy: in-event 180 deg
//           (template and area scaling), +-135 and +-90 deg (template), and,
//           if given, the proxy method (proxy_xf_study config 0, the pipeline)
//   page 2  recoil test: <UE(180)> / <UE(90)> and <UE(135)> / <UE(90)> on the
//           same jets against jet energy; reco (template scaling) and, in
//           simulation, truth photons in full cones (no acceptance)
//   page 3  the UE density map against eta, data and simulation
//   page 4  simulation closure of the scaling: estimate / truth target, where
//           the target is the truth energy in the jet cone's live footprint
//           moved to the rotated cone and the estimate scales the truth on
//           the rotated cone's live cells (area, template, eta strips)
// Rotated cones count when their live area is >= 30% of the cone; +-90 and
// +-135 average the valid sides. Simulation means use the per-jet weights.
// Also prints the tables.
//
// Build and run from this folder:
//   g++ -m64 -O2 in_event_compare.cpp -o in_event_compare `root-config --cflags --libs` -std=c++11
//   ./in_event_compare "<data files>" "<sim files>" <out.pdf> [sim label] ["<proxy_xf data files>" "<proxy_xf sim files>"] [-u]
//   -u  simulation unweighted (one ptHat bin)

#include <TROOT.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

static const int kNB = 4, kNC = 5;
static const double kE[kNB + 1] = {10, 20, 35, 60, 120};
static const double kMinLive = 0.3;

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}

// weighted mean with its error
struct Acc
{
    double w, wx, wxx, w2;
    Acc() : w(0), wx(0), wxx(0), w2(0) {}
    void Add(double x, double ww) { w += ww; wx += ww * x; wxx += ww * x * x; w2 += ww * ww; }
    double Mean() const { return w > 0 ? wx / w : NAN; }
    double Err() const
    {
        if (!(w > 0)) return NAN;
        const double m = wx / w, var = wxx / w - m * m, neff = w * w / w2;
        return neff > 1 ? std::sqrt(std::max(0.0, var) / (neff - 1)) : NAN;
    }
};

// ratio of weighted sums sum w a / sum w b with its linearised error
struct Ratio
{
    double sa, sb, saa, sbb, sab;
    Ratio() : sa(0), sb(0), saa(0), sbb(0), sab(0) {}
    void Add(double a, double b, double w) { sa += w * a; sb += w * b; saa += w * w * a * a; sbb += w * w * b * b; sab += w * w * a * b; }
    double R() const { return sb > 0 ? sa / sb : NAN; }
    double Err() const
    {
        if (!(sb > 0)) return NAN;
        const double r = sa / sb;
        return std::sqrt(std::max(0.0, saa - 2 * r * sab + r * r * sbb)) / sb;
    }
};

enum { kT180, kA180, kS180, kT135, kT90, kS90, kProxy, kNEst };
static const char *kEstName[kNEst] = {"in-event 180#circ, template", "in-event 180#circ, area", "in-event 180#circ, #eta strips", "in-event #pm135#circ, template",
                                      "in-event #pm90#circ, template", "in-event #pm90#circ, #eta strips", "proxy method (pipeline)"};
static const char *kEstTxt[kNEst] = {"180 tmpl", "180 area", "180 strip", "135 tmpl", "90 tmpl", "90 strip", "proxy"};
// closure: cone groups and scaling methods
static const int kNG = 5, kNM = 3;
static const char *kGName[kNG] = {"180#circ", "#pm135#circ", "#pm90#circ", "+90#circ", "-90#circ"};
static const char *kGTxt[kNG] = {"180", "+-135", "+-90", "+90", "-90"};
static const char *kMName[kNM] = {"area", "template", "#eta strips"};
enum { kR180Reco, kR135Reco, kR180Truth, kR135Truth, kNRat };

struct Sample
{
    Acc est[kNB][kNEst];
    Ratio rat[kNB][kNRat];
    Ratio clo[kNB][kNG][kNM];
    TH1D *rho_eta;
    bool truth;
};

static int Bin(double E) { for (int b = 0; b < kNB; b++) if (E >= kE[b] && E < kE[b + 1]) return b; return -1; }

static void Fill(const std::vector<std::string> &files, bool unweighted, Sample &s, const char *tag)
{
    TChain t("jets");
    for (size_t f = 0; f < files.size(); f++) t.Add(files[f].c_str());
    Float_t E, w, t_jet, e[kNC], a[kNC], tt[kNC], te[kNC], a_jet, sv[kNC], tf[kNC], tla[kNC], tlt[kNC], tls[kNC];
    const bool has_s = t.GetBranch("s") != 0;
    if (has_s)
    {
        t.SetBranchAddress("s", sv); t.SetBranchAddress("tf", tf);
        t.SetBranchAddress("tla", tla); t.SetBranchAddress("tlt", tlt); t.SetBranchAddress("tls", tls);
    }
    t.SetBranchAddress("jet_E", &E); t.SetBranchAddress("weight", &w); t.SetBranchAddress("t_jet", &t_jet); t.SetBranchAddress("a_jet", &a_jet);
    t.SetBranchAddress("e", e); t.SetBranchAddress("a", a); t.SetBranchAddress("t", tt); t.SetBranchAddress("te", te);
    s.truth = false;
    for (Long64_t i = 0; i < t.GetEntries(); i++)
    {
        t.GetEntry(i);
        const int b = Bin(E);
        if (b < 0) continue;
        const double ww = unweighted ? 1.0 : w;
        if (!(ww > 0)) continue;
        double u[kNC];
        bool ok[kNC];
        for (int k = 0; k < kNC; k++) { ok[k] = a[k] >= kMinLive && tt[k] > 0; u[k] = ok[k] ? e[k] * t_jet / tt[k] : 0; }
        if (ok[0])
        {
            s.est[b][kT180].Add(u[0], ww);
            s.est[b][kA180].Add(e[0] * a_jet / a[0], ww);
            if (has_s) s.est[b][kS180].Add(sv[0], ww);
        }
        if (has_s)
        {
            const int n90s = ok[3] + ok[4];
            if (n90s) s.est[b][kS90].Add(((ok[3] ? sv[3] : 0) + (ok[4] ? sv[4] : 0)) / n90s, ww);
            // closure per cone, grouped
            const int grp[kNC][2] = {{0, -1}, {1, -1}, {1, -1}, {2, 3}, {2, 4}};
            for (int k = 0; k < kNC; k++)
            {
                if (!ok[k] || !(tf[k] >= 0)) continue;
                const double est[kNM] = {tla[k], tlt[k], tls[k]};
                for (int gq = 0; gq < 2; gq++)
                {
                    const int g = grp[k][gq];
                    if (g < 0) continue;
                    for (int m = 0; m < kNM; m++) s.clo[b][g][m].Add(est[m], tf[k], ww);
                }
            }
        }
        const int n135 = ok[1] + ok[2], n90 = ok[3] + ok[4];
        const double u135 = n135 ? (u[1] + u[2]) / n135 : 0, u90 = n90 ? (u[3] + u[4]) / n90 : 0;
        if (n135) s.est[b][kT135].Add(u135, ww);
        if (n90) s.est[b][kT90].Add(u90, ww);
        if (ok[0] && n90) s.rat[b][kR180Reco].Add(u[0], u90, ww);
        if (n135 && n90) s.rat[b][kR135Reco].Add(u135, u90, ww);
        if (te[0] + te[3] + te[4] > 0) s.truth = true;
        s.rat[b][kR180Truth].Add(te[0], 0.5 * (te[3] + te[4]), ww);
        s.rat[b][kR135Truth].Add(0.5 * (te[1] + te[2]), 0.5 * (te[3] + te[4]), ww);
    }
    s.rho_eta = 0;
    for (size_t f = 0; f < files.size(); f++)
    {
        TFile tf(files[f].c_str());
        TH1D *h = (TH1D *)tf.Get("h_rho_eta");
        if (!h) continue;
        if (!s.rho_eta) { gROOT->cd(); s.rho_eta = (TH1D *)h->Clone(Form("rho_%s", tag)); s.rho_eta->SetDirectory(0); s.rho_eta->Reset(); }
        // unweighted mean over jobs of the maps (each job one run or one ptHat bin)
        s.rho_eta->Add(h, 1.0 / files.size());
    }
}

static void FillProxy(const std::vector<std::string> &files, bool unweighted, Sample &s)
{
    if (files.empty()) return;
    TChain t("jets");
    for (size_t f = 0; f < files.size(); f++) t.Add(files[f].c_str());
    Float_t E, w, ue[11];
    t.SetBranchAddress("jet_E", &E); t.SetBranchAddress("weight", &w); t.SetBranchAddress("ue", ue);
    for (Long64_t i = 0; i < t.GetEntries(); i++)
    {
        t.GetEntry(i);
        const int b = Bin(E);
        const double ww = unweighted ? 1.0 : w;
        if (b < 0 || !(ww > 0) || !(ue[0] >= 0)) continue;
        s.est[b][kProxy].Add(ue[0], ww);
    }
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;
    if (argc < 4) { std::printf("usage: ./in_event_compare \"<data files>\" \"<sim files>\" <out.pdf> [sim label] [\"<proxy data>\" \"<proxy sim>\"] [-u]\n"); return 1; }
    bool unweighted = false;
    std::vector<std::string> args;
    for (int i = 1; i < argc; i++) { if (std::string(argv[i]) == "-u") unweighted = true; else args.push_back(argv[i]); }
    const TString pdf = args[2].c_str();
    const std::string slab = args.size() > 3 ? args[3] : "simulation (EM-triggered)";
    TString png = pdf; png.ReplaceAll(".pdf", "");
    static Sample S[2];
    Fill(Split(args[0].c_str()), false, S[0], "d");
    Fill(Split(args[1].c_str()), unweighted, S[1], "s");
    if (args.size() > 5) { FillProxy(Split(args[4].c_str()), false, S[0]); FillProxy(Split(args[5].c_str()), unweighted, S[1]); }
    const char *lab[2] = {"data (run 23007048)", slab.c_str()};

    std::printf("UE estimate per jet (GeV): data | sim\n  %-10s", "jet E");
    for (int k = 0; k < kNEst; k++) std::printf("  %-27s", kEstTxt[k]);
    std::printf("\n");
    for (int b = 0; b < kNB; b++)
    {
        std::printf("  %3.0f-%-6.0f", kE[b], kE[b + 1]);
        for (int k = 0; k < kNEst; k++) std::printf("  %5.2f+-%.2f | %5.2f+-%.2f  ", S[0].est[b][k].Mean(), S[0].est[b][k].Err(), S[1].est[b][k].Mean(), S[1].est[b][k].Err());
        std::printf("\n");
    }
    if (S[1].truth)
    {
        std::printf("\nSimulation closure of the scaling: estimate / truth target; area | template | eta strips\n");
        for (int b = 0; b < kNB; b++)
        {
            std::printf("  %3.0f-%-6.0f", kE[b], kE[b + 1]);
            for (int g = 0; g < kNG; g++)
                std::printf("  %s: %.3f+-%.3f %.3f+-%.3f %.3f+-%.3f", kGTxt[g], S[1].clo[b][g][0].R(), S[1].clo[b][g][0].Err(), S[1].clo[b][g][1].R(), S[1].clo[b][g][1].Err(), S[1].clo[b][g][2].R(), S[1].clo[b][g][2].Err());
            std::printf("\n");
        }
    }
    std::printf("\nRecoil test, same jets: <UE(180)>/<UE(90)> and <UE(135)>/<UE(90)>; data reco | sim reco | sim truth (full cones)\n");
    for (int b = 0; b < kNB; b++)
        std::printf("  %3.0f-%-6.0f 180/90: %.3f+-%.3f | %.3f+-%.3f | %.3f+-%.3f    135/90: %.3f+-%.3f | %.3f+-%.3f | %.3f+-%.3f\n", kE[b], kE[b + 1],
                    S[0].rat[b][kR180Reco].R(), S[0].rat[b][kR180Reco].Err(), S[1].rat[b][kR180Reco].R(), S[1].rat[b][kR180Reco].Err(),
                    S[1].rat[b][kR180Truth].R(), S[1].rat[b][kR180Truth].Err(),
                    S[0].rat[b][kR135Reco].R(), S[0].rat[b][kR135Reco].Err(), S[1].rat[b][kR135Reco].R(), S[1].rat[b][kR135Reco].Err(),
                    S[1].rat[b][kR135Truth].R(), S[1].rat[b][kR135Truth].Err());

    gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.045, "XY"); gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetPadLeftMargin(0.16); gStyle->SetPadBottomMargin(0.14);
    const int ecol[kNEst] = {kBlack, kGray + 2, kOrange + 7, kAzure + 2, kGreen + 2, kGreen + 3, kRed + 1}, emk[kNEst] = {20, 24, 34, 21, 22, 23, 33};

    // page 1: UE estimates
    {
        TCanvas c("c1", "", 1400, 650);
        TPad top("top", "", 0, 0.88, 1, 1); top.Draw();
        TPad body("body", "", 0, 0, 1, 0.88); body.Draw(); body.Divide(2, 1);
        TGraphErrors *g[2][kNEst];
        double hi = 0;
        for (int m = 0; m < 2; m++)
            for (int k = 0; k < kNEst; k++)
            {
                g[m][k] = new TGraphErrors();
                for (int b = 0; b < kNB; b++)
                {
                    const double v = S[m].est[b][k].Mean();
                    if (!(v == v)) continue;
                    const int n = g[m][k]->GetN();
                    g[m][k]->SetPoint(n, std::sqrt(kE[b] * kE[b + 1]) * (1 + 0.012 * (k - 2)), v);
                    g[m][k]->SetPointError(n, 0, S[m].est[b][k].Err());
                    hi = std::max(hi, v);
                }
                g[m][k]->SetMarkerColor(ecol[k]); g[m][k]->SetLineColor(ecol[k]); g[m][k]->SetMarkerStyle(emk[k]); g[m][k]->SetMarkerSize(1.3); g[m][k]->SetLineWidth(2);
            }
        for (int m = 0; m < 2; m++)
        {
            body.cd(m + 1); gPad->SetLogx();
            TH1D *fr = new TH1D(Form("f1%d", m), ";reco jet E (GeV);UE in the jet cone (GeV)", 1, 10, 120);
            fr->SetDirectory(0); fr->SetMinimum(0); fr->SetMaximum(1.25 * hi); fr->GetXaxis()->SetMoreLogLabels(); fr->Draw("AXIS");
            for (int k = 0; k < kNEst; k++) if (g[m][k]->GetN()) g[m][k]->Draw("P SAME");
            TLatex t; t.SetNDC(); t.SetTextSize(0.05); t.SetTextAlign(21); t.DrawLatex(0.56, 0.93, lab[m]);
        }
        top.cd();
        TLegend *lg = new TLegend(0.02, 0.05, 0.98, 0.95); lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.3);
        for (int k = 0; k < kNEst; k++) if (g[0][k]->GetN() || g[1][k]->GetN()) lg->AddEntry(g[0][k]->GetN() ? g[0][k] : g[1][k], kEstName[k], "p");
        lg->Draw();
        c.Print(png + "_1.png"); c.Print(pdf + "(");
    }
    // page 2: recoil test
    {
        TCanvas c("c2", "", 1400, 650);
        TPad top("top", "", 0, 0.88, 1, 1); top.Draw();
        TPad body("body", "", 0, 0, 1, 0.88); body.Draw(); body.Divide(2, 1);
        const char *rt[2] = {"#LTUE(180#circ)#GT / #LTUE(90#circ)#GT", "#LTUE(135#circ)#GT / #LTUE(90#circ)#GT"};
        const int rc[3] = {kBlack, kAzure + 2, kRed + 1}, rm[3] = {20, 24, 25};
        const char *rl[3] = {"data, reco", "", "sim, truth photons (full cones)"};
        TGraphErrors *keep[3] = {0, 0, 0};
        for (int p = 0; p < 2; p++)
        {
            body.cd(p + 1); gPad->SetLogx();
            TGraphErrors *g[3];
            double lo = 1e9, hi = -1e9;
            for (int q = 0; q < 3; q++)
            {
                const Sample &s = q == 0 ? S[0] : S[1];
                const int ri = q == 2 ? (p ? kR135Truth : kR180Truth) : (p ? kR135Reco : kR180Reco);
                g[q] = new TGraphErrors();
                if (q == 2 && !S[1].truth) continue;
                for (int b = 0; b < kNB; b++)
                {
                    const double v = s.rat[b][ri].R(), er = s.rat[b][ri].Err();
                    if (!(v == v) || !(er < 0.2)) continue;
                    const int n = g[q]->GetN();
                    g[q]->SetPoint(n, std::sqrt(kE[b] * kE[b + 1]) * (1 + 0.02 * (q - 1)), v); g[q]->SetPointError(n, 0, er);
                    lo = std::min(lo, v - er); hi = std::max(hi, v + er);
                }
                g[q]->SetMarkerColor(rc[q]); g[q]->SetLineColor(rc[q]); g[q]->SetMarkerStyle(rm[q]); g[q]->SetMarkerSize(1.3); g[q]->SetLineWidth(2);
                if (!keep[q]) keep[q] = g[q];
            }
            lo = std::min(lo, 0.95); hi = std::max(hi, 1.05);
            const double pd = 0.15 * (hi - lo);
            TH1D *fr = new TH1D(Form("f2%d", p), Form(";reco jet E (GeV);%s", rt[p]), 1, 10, 120);
            fr->SetDirectory(0); fr->SetMinimum(lo - pd); fr->SetMaximum(hi + pd); fr->GetXaxis()->SetMoreLogLabels(); fr->Draw("AXIS");
            TLine *one = new TLine(10, 1, 120, 1); one->SetLineStyle(2); one->Draw();
            for (int q = 0; q < 3; q++) if (g[q]->GetN()) g[q]->Draw("P SAME");
        }
        top.cd();
        TLegend *lg = new TLegend(0.05, 0.05, 0.95, 0.95); lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.3);
        const std::string sreco = slab + ", reco";
        for (int q = 0; q < 3; q++) if (keep[q] && keep[q]->GetN()) lg->AddEntry(keep[q], q == 1 ? sreco.c_str() : rl[q], "p");
        lg->Draw();
        c.Print(png + "_2.png"); c.Print(pdf);
    }
    // page 3: UE density against eta
    {
        TCanvas c("c3", "", 800, 600);
        gPad->SetLogy();
        double hi = 0;
        for (int m = 0; m < 2; m++) if (S[m].rho_eta) hi = std::max(hi, S[m].rho_eta->GetMaximum());
        TH1D *fr = new TH1D("f3", ";#eta;mean UE density away from jets (GeV per unit #eta-#phi area)", 1, 2.2, 4.4);
        fr->SetDirectory(0); fr->SetMinimum(0.05); fr->SetMaximum(2 * hi); fr->Draw("AXIS");
        TLegend *lg = new TLegend(0.20, 0.75, 0.70, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
        const int col[2] = {kBlack, kAzure + 2};
        for (int m = 0; m < 2; m++)
        {
            if (!S[m].rho_eta) continue;
            S[m].rho_eta->SetLineColor(col[m]); S[m].rho_eta->SetMarkerColor(col[m]); S[m].rho_eta->SetMarkerStyle(m ? 24 : 20);
            S[m].rho_eta->Draw("P SAME");
            lg->AddEntry(S[m].rho_eta, lab[m], "p");
        }
        lg->Draw();
        c.Print(png + "_3.png"); c.Print(S[1].truth ? pdf : pdf + ")");
    }
    // page 4: closure
    if (S[1].truth)
    {
        TCanvas c("c4", "", 1500, 600);
        TPad top("top", "", 0, 0.88, 1, 1); top.Draw();
        TPad body("body", "", 0, 0, 1, 0.88); body.Draw(); body.Divide(3, 1);
        const int mc[kNM] = {kGray + 2, kBlack, kOrange + 7}, mm[kNM] = {24, 20, 21};
        const int gsh[3] = {0, 3, 4};
        TGraphErrors *keep[kNM] = {0, 0, 0};
        for (int pq = 0; pq < 3; pq++)
        {
            const int g = gsh[pq];
            body.cd(pq + 1); gPad->SetLogx();
            TH1D *fr = new TH1D(Form("f4%d", pq), ";reco jet E (GeV);estimate / truth in the jet-cone footprint", 1, 10, 120);
            fr->SetDirectory(0); fr->SetMinimum(0.7); fr->SetMaximum(1.3); fr->GetXaxis()->SetMoreLogLabels(); fr->Draw("AXIS");
            TLine *one = new TLine(10, 1, 120, 1); one->SetLineStyle(2); one->Draw();
            for (int m = 0; m < kNM; m++)
            {
                TGraphErrors *gr = new TGraphErrors();
                for (int b = 0; b < kNB; b++)
                {
                    const double v = S[1].clo[b][g][m].R(), er = S[1].clo[b][g][m].Err();
                    if (!(v == v) || !(er < 0.2)) continue;
                    const int n = gr->GetN();
                    gr->SetPoint(n, std::sqrt(kE[b] * kE[b + 1]) * (1 + 0.02 * (m - 1)), v); gr->SetPointError(n, 0, er);
                }
                gr->SetMarkerColor(mc[m]); gr->SetLineColor(mc[m]); gr->SetMarkerStyle(mm[m]); gr->SetMarkerSize(1.3); gr->SetLineWidth(2);
                gr->Draw("P SAME");
                if (!keep[m]) keep[m] = gr;
            }
            TLatex t; t.SetNDC(); t.SetTextSize(0.055); t.SetTextAlign(21); t.DrawLatex(0.56, 0.93, Form("cone at %s (%s)", kGName[g], slab.c_str()));
        }
        top.cd();
        TLegend *lg = new TLegend(0.1, 0.05, 0.9, 0.95); lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.4);
        for (int m = 0; m < kNM; m++) lg->AddEntry(keep[m], Form("%s scaling", kMName[m]), "p");
        lg->Draw();
        c.Print(png + "_4.png"); c.Print(pdf + ")");
    }
    return 0;
}
