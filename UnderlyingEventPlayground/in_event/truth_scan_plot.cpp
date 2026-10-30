// The truth-level azimuthal UE scans of truth_scan.cpp, laid out as the FMS
// test (FMSJetPipeline/compare/wedge_density.cpp): left, the ring wedges from
// 0 to 180 deg (no jet exclusion), vertical lines at the jet edge (R) and the
// 2R separation, horizontal lines at the mean of the wedges outside 2R and
// outside R; right, the two cones moved from 0 to 180 deg, with the same
// lines.
// Means over jets weighted by sigma / N of the ptHat bin. "-t": events with
// an EM bit (the simulation's own emulation) only; "-1": both scans on one
// |dphi| axis (wedge centres and cone placements), left the whole range on a
// log scale, right a linear zoom on the UE from kZoomFrom deg.
//
// Build and run from this folder:
//   g++ -m64 -O2 truth_scan_plot.cpp -o truth_scan_plot `root-config --cflags --libs` -std=c++11
//   ./truth_scan_plot <truth_scan.root> <out.png> [-t] [-1]

#include <TROOT.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

static const double kR = 0.7;
static const int kNW = 18;
static const double kEdges[kNW + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180};
static const int kNOff = 19;
static const double kOff[kNOff] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180};
static const double kOutsideFrom = 80, kOutsideJetFrom = 40; // deg: wedges from 2R and from R
static const double kCutFrom = 150; // deg: the candidate cut, 30 deg around 180 left out

struct WMean
{
    double w = 0, w2 = 0, s = 0, ss = 0;
    void Add(double v, double ww) { w += ww; w2 += ww * ww; s += ww * v; ss += ww * v * v; }
    double Val() const { return w > 0 ? s / w : NAN; }
    double Err() const { const double m = Val(); return w > 0 ? std::sqrt(std::max(0.0, ss / w - m * m) * w2) / w : NAN; }
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::fprintf(stderr, "usage: ./truth_scan_plot <truth_scan.root> <out.png> [-t]\n"); return 1; }
    bool trig = false, one = false;
    for (int a = 3; a < argc; a++) { if (std::string(argv[a]) == "-t") trig = true; if (std::string(argv[a]) == "-1") one = true; }
    TFile f(argv[1]);
    TTree *t = (TTree *)f.Get("jets");
    TH1D *hc = (TH1D *)f.Get("h_wcov");
    if (!t || !hc) { std::fprintf(stderr, "no jets tree / h_wcov\n"); return 1; }
    Float_t w_rho[kNW], c_rho[kNOff], w;
    Int_t em;
    t->SetBranchAddress("w_rho", w_rho); t->SetBranchAddress("c_rho", c_rho); t->SetBranchAddress("weight", &w); t->SetBranchAddress("em", &em);
    WMean mw[kNW], mc[kNOff];
    long nj = 0;
    for (Long64_t i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);
        if (trig && em < 0) continue;
        nj++;
        for (int k = 0; k < kNW; k++) if (w_rho[k] >= 0) mw[k].Add(w_rho[k], w);
        for (int k = 0; k < kNOff; k++) mc[k].Add(c_rho[k], w);
    }
    const double sep = 2 * kR * 180 / M_PI;
    double sout = 0, sall = 0; int nout = 0, nall = 0;
    for (int k = 0; k < kNW; k++)
    {
        if (!(mw[k].w > 0)) continue;
        sall += mw[k].Val(); nall++;
        if (kEdges[k] >= kOutsideFrom) { sout += mw[k].Val(); nout++; }
    }
    sall = 0; nall = 0;
    for (int k = 0; k < kNW; k++) if (mw[k].w > 0 && kEdges[k] >= kOutsideJetFrom) { sall += mw[k].Val(); nall++; }
    const double aout = sout / nout, aall = sall / nall;
    // the candidate method: wedges between the 2R separation and the cut
    double scut = 0; int ncut = 0;
    for (int k = 0; k < kNW; k++) if (mw[k].w > 0 && kEdges[k] >= kOutsideFrom && kEdges[k + 1] <= kCutFrom) { scut += mw[k].Val(); ncut++; }
    const double acut = scut / ncut;
    std::printf("mean of the wedges from 2R to %.0f deg: %.4f\n", kCutFrom, acut);
    std::printf("truth_scan_plot: %ld truth jets%s\n", nj, trig ? " (EM-triggered events)" : "");
    std::printf("mean of the wedges outside the 2R separation (from %.0f deg): %.4f; outside the jet radius (from %.0f deg): %.4f GeV per unit eta-phi area\n", kOutsideFrom, aout, kOutsideJetFrom, aall);
    std::printf("  wedge |dphi| (deg)   rho              coverage      cones at +-dphi   rho\n");
    for (int k = 0; k < kNW; k++)
    {
        std::printf("  %3.0f-%-3.0f              %.4f+-%.4f   %.2f", kEdges[k], kEdges[k + 1], mw[k].Val(), mw[k].Err(), hc->GetBinContent(k + 1));
        const int c = k;
        if (c < kNOff) std::printf("          %3.0f            %.4f+-%.4f", kOff[c], mc[c].Val(), mc[c].Err());
        std::printf("\n");
    }

    gStyle->SetOptStat(0);
    const int cRing = TColor::GetColor("#d9480f"), cCone = TColor::GetColor("#1c4f8a"), cInk = TColor::GetColor("#222222"), cAll = TColor::GetColor("#868e96");
    double ymax = 0, ymin = 1e9;
    for (int k = 0; k < kNW; k++) if (mw[k].w > 0) { ymax = std::max(ymax, mw[k].Val() + mw[k].Err()); ymin = std::min(ymin, mw[k].Val() - mw[k].Err()); }
    for (int k = 0; k < kNOff; k++) { ymax = std::max(ymax, mc[k].Val() + mc[k].Err()); ymin = std::min(ymin, mc[k].Val() - mc[k].Err()); }
    const double y0 = 0.6 * ymin, y1 = 4 * ymax; // log y: the jet's own density near 0 deg is far above the UE's
    // linear zoom on the UE: the points from kZoomFrom deg
    const double kZoomFrom = 60;
    double zmin = 1e9, zmax = 0;
    for (int k = 0; k < kNW; k++) if (mw[k].w > 0 && kEdges[k] >= kZoomFrom) { zmin = std::min(zmin, mw[k].Val() - mw[k].Err()); zmax = std::max(zmax, mw[k].Val() + mw[k].Err()); }
    for (int k = 0; k < kNOff; k++) if (kOff[k] >= kZoomFrom) { zmin = std::min(zmin, mc[k].Val() - mc[k].Err()); zmax = std::max(zmax, mc[k].Val() + mc[k].Err()); }
    // the two scans as graphs
    TGraphErrors gw, gc;
    for (int k = 0; k < kNW; k++)
    {
        if (!(mw[k].w > 0)) continue;
        const int i = gw.GetN();
        gw.SetPoint(i, 0.5 * (kEdges[k] + kEdges[k + 1]), mw[k].Val()); gw.SetPointError(i, 0, mw[k].Err());
    }
    for (int k = 0; k < kNOff; k++) { gc.SetPoint(k, kOff[k], mc[k].Val()); gc.SetPointError(k, 0, mc[k].Err()); }
    gw.SetMarkerStyle(21); gw.SetMarkerColor(cRing); gw.SetLineColor(cRing);
    gc.SetMarkerStyle(20); gc.SetMarkerColor(cCone); gc.SetLineColor(cCone);
    for (TGraphErrors *g : {&gw, &gc}) { g->SetMarkerSize(1.3); g->SetLineWidth(2); }
    // the FMS standard placement, cones at +-90 deg: a red star
    const int cStd = TColor::GetColor("#e03131");
    TGraphErrors gstd;
    for (int k = 0; k < kNOff; k++) if (kOff[k] == 90) { gstd.SetPoint(0, 90, mc[k].Val()); gstd.SetPointError(0, 0, mc[k].Err()); }
    gstd.SetMarkerStyle(29); gstd.SetMarkerSize(3.6); gstd.SetMarkerColor(cStd); gstd.SetLineColor(cStd); gstd.SetLineWidth(2);
    TGraphErrors gstd_o(gstd); // a black outline under it, so it stands out over the curves
    gstd_o.SetMarkerStyle(29); gstd_o.SetMarkerSize(4.6); gstd_o.SetMarkerColor(kBlack); gstd_o.SetLineColor(kBlack);

    TCanvas cv("cv", "", 1500, one ? 720 : 620);
    // "-1": a strip across the top for the title and legend, the two panels below
    TPad *pads[2];
    if (one)
    {
        pads[0] = new TPad("pa", "", 0, 0, 0.5, 0.80); pads[1] = new TPad("pb", "", 0.5, 0, 1, 0.80);
        pads[0]->Draw(); pads[1]->Draw();
    }
    else cv.Divide(2, 1, 0.002, 0.002);
    std::vector<TObject *> keep;
    for (int p = 0; p < 2; p++)
    {
        if (one) pads[p]->cd(); else cv.cd(p + 1);
        gPad->SetMargin(0.12, 0.03, 0.12, 0.04);
        const bool zoom = one && p == 1, cones_only = !one && p == 1;
        if (!zoom) gPad->SetLogy();
        const double x0 = zoom ? kZoomFrom - 5 : 0, x1 = 185;
        const double ylo = zoom ? zmin - 0.25 * (zmax - zmin) : y0, yhi = zoom ? zmax + 0.25 * (zmax - zmin) : (one ? 1.5 * ymax : y1);
        TH1D *fr = new TH1D(Form("fr%d", p), "", 1, x0, x1); keep.push_back(fr);
        fr->SetTitle(one ? ";|#Delta#phi| from the jet: wedge centre or cone placement (deg);truth energy density per unit #eta-#phi area (GeV)"
                         : cones_only ? ";cone placement |#Delta#phi| from the jet (deg);truth energy density per unit #eta-#phi area (GeV)"
                                      : ";|#Delta#phi| from the jet, wedge centre (deg);truth energy density per unit #eta-#phi area (GeV)");
        fr->SetMinimum(ylo); fr->SetMaximum(yhi); fr->GetYaxis()->SetTitleOffset(1.3); fr->Draw();
        TLine *lo = new TLine(x0, aout, x1, aout), *la = new TLine(x0, aall, x1, aall);
        keep.push_back(lo); keep.push_back(la);
        lo->SetLineColor(cInk); lo->SetLineStyle(2); lo->SetLineWidth(2); lo->Draw();
        la->SetLineColor(cAll); la->SetLineStyle(3); la->SetLineWidth(3); // not drawn
        const int cCut = TColor::GetColor("#2f9e44");
        TLine *lc = new TLine(x0, acut, x1, acut); keep.push_back(lc);
        lc->SetLineColor(cCut); lc->SetLineWidth(2); lc->Draw();
        {
            TLatex ts; ts.SetTextFont(42); ts.SetTextSize(0.032); ts.SetTextColor(cInk); ts.SetTextAngle(90); ts.SetTextAlign(33); // at the top of the line, reading down to it, clear of the points
            for (int v = 0; v < 3; v++)
            {
                const double xv = v == 2 ? kCutFrom : (v + 1) * kR * 180 / M_PI;
                if (xv < x0) continue;
                TLine *ls = new TLine(xv, ylo, xv, yhi); keep.push_back(ls);
                ls->SetLineColor(v == 2 ? cCut : cInk); ls->SetLineStyle(v ? 7 : 1); ls->SetLineWidth(v ? 2 : 1); ls->Draw();
                ts.SetTextColor(v == 2 ? cCut : cInk);
                ts.DrawLatex(xv + 1.5, zoom ? yhi - 0.03 * (yhi - ylo) : yhi / 1.15, v == 2 ? "180#circ #minus 30#circ" : v ? "2R separation" : "jet edge (R)");
            }
        }
        if (!cones_only)
        {
            gw.Draw("PL SAME");
        }
        if (one || cones_only) { gc.Draw("PL SAME"); gstd_o.Draw("P SAME"); gstd.Draw("P SAME"); }
        if (one && p == 1) continue;
        if (one) cv.cd();
        TLegend *lg = one ? new TLegend(0.06, 0.81, 0.98, 0.93) : new TLegend(0.15, 0.74, 0.97, 0.91); keep.push_back(lg);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(one ? 0.028 : 0.034);
        if (one) lg->SetNColumns(2);
        if (!cones_only) lg->AddEntry(&gw, "pixel ring: average energy density per wedge choice", "pl");
        if (one || cones_only) { lg->AddEntry(&gc, "two cones at #pm#Delta#phi: average energy density per cone placement", "pl"); lg->AddEntry(&gstd, "two cones at #pm90#circ: the FMS standard", "p"); }
        lg->AddEntry(lo, "ring wedges outside the 2R separation, average", "l");
        lg->AddEntry(lc, "ring wedges from 2R to 150#circ, average", "l");
        lg->Draw();
        TLatex tt; tt.SetNDC(); tt.SetTextFont(42); tt.SetTextSize(0.036);
        if (one) tt.SetTextSize(0.03);
        tt.DrawLatex(one ? 0.06 : 0.12, one ? 0.955 : 0.94, cones_only ? "Two-cone method, cones moved from 0#circ to 180#circ"
                                            : Form("FCS sim truth (EM photons, full #phi), %ld jets%s%s", nj, trig ? ", EM-triggered" : "", one ? ": left, full range (log); right, the UE region (linear)" : ""));
    }
    cv.SaveAs(argv[2]);
    return 0;
}
