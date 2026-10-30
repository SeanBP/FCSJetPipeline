// Checks of the pixel ring (in_event_ue v7+ outputs), every jet used (no cuts):
//   per jet-energy bin: simulation closure of the pixel ring and of the strip
//   ring (both without the other-jet exclusion; target an_tf[3]); template
//   fallbacks; jets with eta bins the ring does not reach; mean largest weight; mean UE
//   data | sim with other jets excluded
//   per live ring area bin (20-120 GeV): share of jets, closure, data per-jet
//   RMS / mean
// Also writes a figure: closure against jet energy and against ring area.
//
// Build and run from this folder:
//   g++ -m64 -O2 px_check.cpp -o px_check `root-config --cflags --libs` -std=c++11
//   ./px_check "<data files>" "<sim files>" <out.png>

#include <TROOT.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

static const int kNA = 6, kNE = 4, kNBA = 6;
static const double kE[kNE + 1] = {10, 20, 35, 60, 120};
static const double kAEdge[kNBA + 1] = {0, 0.75, 1.0, 1.25, 1.5, 2.0, 99};
static const double kCone = M_PI * 0.49, kMaxFill = 1 / 0.7;

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}
struct Ratio
{
    double w = 0, w2 = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    void Add(double x, double y, double ww) { w += ww; w2 += ww * ww; sx += ww * x; sy += ww * y; sxx += ww * x * x; syy += ww * y * y; sxy += ww * x * y; }
    double Val() const { return sy != 0 ? sx / sy : NAN; }
    double Err() const
    {
        if (!(w > 0) || sy == 0) return NAN;
        const double mx = sx / w, my = sy / w, r = mx / my;
        const double vx = sxx / w - mx * mx, vy = syy / w - my * my, cxy = sxy / w - mx * my;
        return std::fabs(r) * std::sqrt(std::max(0.0, vx / (mx * mx) + vy / (my * my) - 2 * cxy / (mx * my)) / (w * w / w2));
    }
};
struct Acc { double w = 0, wx = 0, wxx = 0; void Add(double x, double ww) { w += ww; wx += ww * x; wxx += ww * x * x; } double M() const { return w > 0 ? wx / w : NAN; } double Rms() const { const double m = M(); return w > 0 ? std::sqrt(std::max(0.0, wxx / w - m * m)) : NAN; } };

struct Res
{
    Ratio cl_px[kNE], cl_strip[kNE], cl_px_a[kNBA], cl_strip_a[kNBA];
    double n[kNE] = {}, nfb[kNE] = {}, nwide[kNE] = {}, na[kNBA] = {}, ntot = 0;
    Acc wmax[kNE], ue[kNE], est_a[kNBA];
};

static void Fill(const std::vector<std::string> &files, bool sim, Res &r)
{
    TChain ch("jets");
    for (size_t i = 0; i < files.size(); i++) ch.Add(files[i].c_str());
    Float_t E, w, an_tl[kNA], an_tf[kNA], an_fill[kNA], px_e[2], px_tl[2], px_area[2], px_cov[2], px_wmax[2];
    Int_t px_fb[2];
    ch.SetBranchAddress("jet_E", &E); ch.SetBranchAddress("weight", &w);
    ch.SetBranchAddress("an_tl", an_tl); ch.SetBranchAddress("an_tf", an_tf); ch.SetBranchAddress("an_fill", an_fill);
    ch.SetBranchAddress("px_e", px_e); ch.SetBranchAddress("px_tl", px_tl); ch.SetBranchAddress("px_area", px_area);
    ch.SetBranchAddress("px_cov", px_cov); ch.SetBranchAddress("px_wmax", px_wmax); ch.SetBranchAddress("px_fb", px_fb);
    for (Long64_t i = 0; i < ch.GetEntries(); i++)
    {
        ch.GetEntry(i);
        int b = -1;
        for (int q = 0; q < kNE; q++) if (E >= kE[q] && E < kE[q + 1]) b = q;
        if (b < 0) continue;
        const double ww = sim ? w : 1.0;
        r.n[b] += ww; if (px_fb[0]) r.nfb[b] += ww; if (px_cov[0] < 0.999) r.nwide[b] += ww;
        r.wmax[b].Add(px_wmax[0], ww); r.ue[b].Add(px_e[0], ww);
        if (sim)
        {
            r.cl_px[b].Add(px_tl[1], an_tf[3], ww);
            const double st = an_fill[3] > 0 ? an_tl[3] : 0;  // strip ring, no cut: a jet with nothing covered counts as 0
            r.cl_strip[b].Add(st, an_tf[3], ww);
        }
        if (b == 0) continue;  // area table: 20-120 GeV
        int a = -1, a1 = -1;
        for (int q = 0; q < kNBA; q++) { if (px_area[0] / kCone >= kAEdge[q] && px_area[0] / kCone < kAEdge[q + 1]) a = q; if (px_area[1] / kCone >= kAEdge[q] && px_area[1] / kCone < kAEdge[q + 1]) a1 = q; }
        if (a >= 0) { r.na[a] += ww; r.ntot += ww; if (!sim) r.est_a[a].Add(px_e[0], 1.0); }
        if (sim && a1 >= 0) { r.cl_px_a[a1].Add(px_tl[1], an_tf[3], ww); r.cl_strip_a[a1].Add(an_fill[3] > 0 ? an_tl[3] : 0, an_tf[3], ww); }
    }
    std::printf("%s: %lld jets\n", sim ? "simulation" : "data", ch.GetEntries());
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 4) { std::fprintf(stderr, "usage: ./px_check \"<data files>\" \"<sim files>\" <out.png>\n"); return 1; }
    Res D, S;
    Fill(Split(argv[1]), false, D);
    Fill(Split(argv[2]), true, S);
    std::printf("\nall jets, no cuts\n jet E     closure pixel    closure strip    fallback data|sim  ring short of jet eta data|sim  <W_max> data  UE data | sim (other jets excl.)\n");
    for (int b = 0; b < kNE; b++)
        std::printf("%3.0f-%-3.0f   %.3f+-%.3f    %.3f+-%.3f     %.4f|%.4f        %.3f|%.3f             %.2f         %.2f | %.2f\n", kE[b], kE[b + 1],
                    S.cl_px[b].Val(), S.cl_px[b].Err(), S.cl_strip[b].Val(), S.cl_strip[b].Err(), D.nfb[b] / D.n[b], S.nfb[b] / S.n[b],
                    D.nwide[b] / D.n[b], S.nwide[b] / S.n[b], D.wmax[b].M(), D.ue[b].M(), S.ue[b].M());
    std::printf("\n20-120 GeV, by live ring area (pi R^2)\n area        share data  sim   closure pixel    closure strip    data RMS/mean\n");
    for (int a = 0; a < kNBA; a++)
        std::printf(" %4.2f-%-5.2f  %5.3f %5.3f   %.3f+-%.3f    %.3f+-%.3f    %.2f\n", kAEdge[a], kAEdge[a + 1] > 50 ? 9.99 : kAEdge[a + 1],
                    D.na[a] / D.ntot, S.na[a] / S.ntot, S.cl_px_a[a].Val(), S.cl_px_a[a].Err(), S.cl_strip_a[a].Val(), S.cl_strip_a[a].Err(),
                    D.est_a[a].Rms() / D.est_a[a].M());

    gStyle->SetOptStat(0);
    TCanvas cv("cv", "", 1400, 560);
    cv.Divide(2, 1);
    const int cpx = TColor::GetColor("#1c4f8a"), cst = TColor::GetColor("#868e96");
    std::vector<TObject *> keep;
    for (int p = 0; p < 2; p++)
    {
        cv.cd(p + 1);
        gPad->SetMargin(0.13, 0.03, 0.13, 0.09);
        if (p == 0) gPad->SetLogx();
        TH1D *fr = new TH1D(Form("fr%d", p), "", 1, p == 0 ? 10 : 0.5, p == 0 ? 120 : 3.0);
        fr->SetTitle(p == 0 ? ";jet E (GeV);truth estimate / truth target" : ";live ring area (#piR^{2});truth estimate / truth target");
        fr->SetMinimum(0.8); fr->SetMaximum(1.3); fr->GetYaxis()->SetTitleOffset(1.2);
        if (p == 0) fr->GetXaxis()->SetMoreLogLabels();
        fr->Draw(); keep.push_back(fr);
        TLine *l = new TLine(fr->GetXaxis()->GetXmin(), 1, fr->GetXaxis()->GetXmax(), 1); l->SetLineStyle(2); l->Draw(); keep.push_back(l);
        TLegend *lg = new TLegend(0.16, 0.74, 0.95, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04); keep.push_back(lg);
        for (int k = 0; k < 2; k++)
        {
            TGraphErrors *g = new TGraphErrors();
            const int n = p == 0 ? kNE : kNBA;
            for (int q = 0; q < n; q++)
            {
                const Ratio &rr = p == 0 ? (k == 0 ? S.cl_px[q] : S.cl_strip[q]) : (k == 0 ? S.cl_px_a[q] : S.cl_strip_a[q]);
                if (!(rr.w > 0)) continue;
                double x = p == 0 ? std::sqrt(kE[q] * kE[q + 1]) : 0.5 * (kAEdge[q] + std::min(kAEdge[q + 1], 3.0));
                x *= p == 0 ? (1 + 0.02 * (k ? 1 : -1)) : 1;
                x += p == 1 ? 0.02 * (k ? 1 : -1) : 0;
                const int i = g->GetN();
                g->SetPoint(i, x, rr.Val()); g->SetPointError(i, 0, rr.Err());
            }
            g->SetMarkerStyle(k == 0 ? 20 : 24); g->SetMarkerColor(k == 0 ? cpx : cst); g->SetLineColor(k == 0 ? cpx : cst); g->SetMarkerSize(1.3);
            g->Draw("PL SAME"); keep.push_back(g);
            lg->AddEntry(g, k == 0 ? "pixel ring (per-cell #eta weights, every jet)" : "strip ring (every jet, no validity cut)", "pl");
        }
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
        t.DrawLatex(0.13, 0.93, p == 0 ? "Simulation closure of the whole ring (#geq 2R), all jets" : "against the live ring area, 20-120 GeV jets");
    }
    cv.SaveAs(argv[3]);
    return 0;
}
