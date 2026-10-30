// Diagram of the pixel ring UE estimate (in_event_geo.h BuildPixelRing) for
// one example jet on the FCS ECal: live cells, the jet cone and its live
// footprint, the ring cells (>= 2R from the jet, >= R from other jets) shaded
// by their cell weight w_i (a picture; the estimate weights each point by
// the ratio of its own eta bin, which the weighted-ring curve shows);
// beside it the area-per-eta of the jet cone, of the raw ring and of the
// weighted ring (0.05 bins, from sub-pixels), and the per-bin ratio
// h_jet / h_ring that the cell weights are built from.
//
// Build and run from this folder:
//   g++ -O2 annulus_diagram.cpp -o annulus_diagram -I../../shared -I../../underlying_event `root-config --cflags --libs` -std=c++11
//   ./annulus_diagram <out.png> [jet eta = 3.0] [jet phi deg = -62] [other eta = 3.3 | none] [other phi deg = 150] [|dphi| cap deg = 180]
// ("none": no other jet in the event; a cap below 180 leaves out the ring
// cells at |dphi| above it, drawn as two lines, with the share of the uncut
// ring's live area kept)

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "in_event_geo.h"

// detector-plane curve of the points at eta-phi distance d from (eta, phi)
static TGraph *Circle(double eta, double phi, double d)
{
    TGraph *g = new TGraph();
    for (int i = 0; i <= 360; i++)
    {
        const double t = i * M_PI / 180.0, e = eta + d * std::cos(t), p = phi + d * std::sin(t);
        const double r = z_proj / std::sinh(e);
        g->SetPoint(i, r * std::cos(p), r * std::sin(p));
    }
    return g;
}
// the same, only where |eta - band_eta| < R (one graph per contiguous piece)
static std::vector<TGraph *> CircleInBand(double eta, double phi, double d, double band_eta)
{
    std::vector<TGraph *> out;
    TGraph *g = 0;
    for (int i = 0; i <= 720; i++)
    {
        const double t = i * M_PI / 360.0, e = eta + d * std::cos(t), p = phi + d * std::sin(t);
        if (std::fabs(e - band_eta) >= R) { g = 0; continue; }
        if (!g) { g = new TGraph(); out.push_back(g); }
        const double r = z_proj / std::sinh(e);
        g->SetPoint(g->GetN(), r * std::cos(p), r * std::sin(p));
    }
    return out;
}
static TGraph *Ring(double eta)
{
    TGraph *g = new TGraph();
    const double r = z_proj / std::sinh(eta);
    for (int i = 0; i <= 360; i++) g->SetPoint(i, r * std::cos(i * M_PI / 180.0), r * std::sin(i * M_PI / 180.0));
    return g;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2) { std::fprintf(stderr, "usage: ./annulus_diagram <out.png> [jet eta] [jet phi deg] [other eta] [other phi deg]\n"); return 1; }
    const double jeta = argc > 2 ? std::atof(argv[2]) : 3.0, jphi = (argc > 3 ? std::atof(argv[3]) : -62.0) * M_PI / 180.0;
    const bool has_other = !(argc > 4 && std::string(argv[4]) == "none");
    const double oeta = argc > 4 && has_other ? std::atof(argv[4]) : 3.3, ophi = (argc > 5 ? std::atof(argv[5]) : 150.0) * M_PI / 180.0;

    Geometry geo; geo.Build();
    std::vector<int> cl;
    geo.Near(jeta, jphi, R, cl);
    std::vector<char> in_jet(geo.cells.size(), 0);
    for (size_t q = 0; q < cl.size(); q++) in_jet[cl[q]] = 1;
    std::vector<std::pair<float, float> > others;
    if (has_other) others.push_back(std::make_pair((float)oeta, (float)ophi));
    PixelRing pr;
    const double cap = (argc > 6 ? std::atof(argv[6]) : 180.0) * M_PI / 180.0;
    BuildPixelRing(geo, jeta, jphi, others, kAnnJetExcl, pr);
    const double area_uncut = pr.area_ring;
    BuildPixelRing(geo, jeta, jphi, others, kAnnJetExcl, pr, cap + 1e-9);
    const double kept = area_uncut > 0 ? pr.area_ring / area_uncut : 0;
    std::printf("|dphi| cap %.0f deg: ring live area %.3f of the uncut %.3f (%.3f)\n", cap * 180 / M_PI, pr.area_ring, area_uncut, kept);

    // eta distributions per strip
    // area per eta bin from the sub-pixels: jet cone, raw ring, weighted ring
    TH1D hj("hj", "", kPixNB, pr.eta0, pr.eta0 + kPixNB * kPixBin), ha("ha", "", kPixNB, pr.eta0, pr.eta0 + kPixNB * kPixBin);
    TH1D hw("hw", "", kPixNB, pr.eta0, pr.eta0 + kPixNB * kPixBin), hs("hs", "", kPixNB, pr.eta0, pr.eta0 + kPixNB * kPixBin);
    double wmax = 0;
    std::vector<char> in_ring(geo.cells.size(), 0);
    const int ns = kSub * kSub;
    for (size_t c = 0; c < geo.cells.size(); c++)
    {
        const Cell &cc = geo.cells[c];
        in_ring[c] = std::fabs(cc.eta - jeta) <= R + 0.2 && InPixelRing(cc.eta, cc.phi, jeta, jphi, kAnnJetExcl, others, cap + 1e-9);
        for (int k = 0; k < ns; k++)
        {
            const Pix &p = geo.pix[c * ns + k];
            if (ProxyDeltaR(p.eta, p.phi, jeta, jphi) < R) hj.Fill(p.eta, p.area);
            if (in_ring[c]) { ha.Fill(p.eta, p.area); hw.Fill(p.eta, p.area * pr.Weight((int)c, p.eta)); }
        }
        wmax = std::max(wmax, (double)pr.w[c]);
    }
    // the ratio the estimate uses (0 where the ring covers under kPixMinRatio of the jet's area)
    for (int b = 1; b <= kPixNB; b++) hs.SetBinContent(b, pr.ratio[b - 1]);
    std::printf("jet live area %.3f (fraction %.2f), ring live area %.3f, weighted ring %.3f, matched share %.3f\n",
                pr.area_jet, pr.area_jet / (M_PI * R * R), pr.area_ring, pr.den, pr.den / pr.area_jet);

    // sequential palette, one hue light to dark
    const int NC = 50;
    double stops[3] = {0.0, 0.5, 1.0}, red[3] = {0.87, 0.42, 0.03}, green[3] = {0.92, 0.61, 0.19}, blue[3] = {0.97, 0.84, 0.42};
    const int c0 = TColor::CreateGradientColorTable(3, stops, red, green, blue, NC);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(NC, 0);
    {
        static int pal[NC];
        for (int i = 0; i < NC; i++) pal[i] = c0 + i;
        gStyle->SetPalette(NC, pal);
    }
    const int kLive = TColor::GetColor("#e3e3e3"), kJet = TColor::GetColor("#d9480f"), kInk = TColor::GetColor("#222222");

    TCanvas cv("cv", "", 1600, 860);
    TPad pm("pm", "", 0.0, 0.0, 0.58, 1.0), pt("pt", "", 0.58, 0.5, 1.0, 1.0), pb("pb", "", 0.58, 0.0, 1.0, 0.5);
    pm.SetMargin(0.11, 0.17, 0.09, 0.07); pt.SetMargin(0.15, 0.05, 0.14, 0.12); pb.SetMargin(0.15, 0.05, 0.16, 0.10);
    pm.Draw(); pt.Draw(); pb.Draw();

    pm.cd();
    const double L = 160;
    TH2D fr("fr", ";x (cm);y (cm);cell weight w_{i}", 10, -L, L, 10, -110, 165);
    fr.SetMinimum(0); fr.SetMaximum(wmax);
    fr.SetBinContent(1, 1, -1);  // keep the empty frame, show the palette axis
    fr.GetZaxis()->SetTitleOffset(1.6); fr.GetYaxis()->SetTitleOffset(1.2);
    fr.Draw("COLZ");
    const double dx = (kFiducialRectEcal.x_outer - kFiducialRectEcal.x_inner) / geo.nx, dy = (kFiducialRectEcal.y_max - kFiducialRectEcal.y_min) / geo.ny;
    std::vector<TBox *> boxes;
    for (size_t c = 0; c < geo.cells.size(); c++)
    {
        const Cell &cc = geo.cells[c];
        TBox *b = new TBox(cc.x - dx / 2, cc.y - dy / 2, cc.x + dx / 2, cc.y + dy / 2);
        int col = kLive;
        if (in_ring[c]) col = c0 + std::min(NC - 1, (int)(pr.w[c] / wmax * (NC - 1) + 0.5));
        b->SetFillColor(col); b->SetLineWidth(0);
        b->Draw();
        boxes.push_back(b);
    }
    // jet live footprint, outlined cell by cell
    for (size_t c = 0; c < geo.cells.size(); c++)
        if (in_jet[c])
        {
            const Cell &cc = geo.cells[c];
            TBox *b = new TBox(cc.x - dx / 2, cc.y - dy / 2, cc.x + dx / 2, cc.y + dy / 2);
            b->SetFillColor(TColor::GetColor("#f7c59f")); b->SetLineWidth(0); b->Draw();
        }
    TGraph *gj = Circle(jeta, jphi, R); gj->SetLineColor(kJet); gj->SetLineWidth(3); gj->Draw("L");
    std::vector<TGraph *> gxs = CircleInBand(jeta, jphi, kAnnJetExcl, jeta);
    for (size_t q = 0; q < gxs.size(); q++) { gxs[q]->SetLineColor(kJet); gxs[q]->SetLineWidth(2); gxs[q]->SetLineStyle(2); gxs[q]->Draw("L"); }
    TGraph *gx = gxs.empty() ? new TGraph() : gxs[0];
    TGraph *go = Circle(oeta, ophi, R); go->SetLineColor(kInk); go->SetLineWidth(2); go->SetLineStyle(2); if (has_other) go->Draw("L");
    TGraph *g1 = Ring(jeta - R), *g2 = Ring(jeta + R);
    g1->SetLineColor(kInk); g2->SetLineColor(kInk); g1->SetLineWidth(2); g2->SetLineWidth(2);
    g1->Draw("L"); g2->Draw("L");
    // the cut around the recoil side: lines at |dphi| = cap across the band
    TLine *lcut = new TLine();
    lcut->SetLineColor(TColor::GetColor("#c2255c")); lcut->SetLineWidth(3); lcut->SetLineStyle(2);
    if (cap < M_PI - 1e-6)
        for (int s = -1; s <= 1; s += 2)
        {
            const double ph = jphi + s * cap, r0 = z_proj / std::sinh(jeta + R), r1 = z_proj / std::sinh(jeta - R);
            lcut->DrawLine(r0 * std::cos(ph), r0 * std::sin(ph), r1 * std::cos(ph), r1 * std::sin(ph));
        }
    TLatex tx; tx.SetTextFont(42); tx.SetTextSize(0.026); tx.SetTextColor(kInk);
    const double rj = z_proj / std::sinh(jeta);
    tx.SetTextColor(kJet); tx.SetTextAlign(22);
    tx.DrawLatex(rj * std::cos(jphi), rj * std::sin(jphi), "jet");
    tx.SetTextColor(kInk);
    const double ro = z_proj / std::sinh(oeta);
    if (has_other) tx.DrawLatex(ro * std::cos(ophi), ro * std::sin(ophi), "other jet");
    tx.SetTextAlign(12); tx.SetTextSize(0.024);
    tx.DrawLatexNDC(0.12, 0.965, Form("Example jet #eta = %.2f, #phi = %.0f#circ: ring cells |#Delta#eta| < R, #geq 2R from the jet%s", jeta, jphi * 180 / M_PI, has_other ? ", #geq R from other jets" : ""));
    TLegend lg(0.12, 0.765, 0.82, 0.925);
    lg.SetNColumns(2);
    lg.SetBorderSize(0); lg.SetFillStyle(1001); lg.SetFillColor(0); lg.SetTextSize(0.022);
    TBox bl(0, 0, 1, 1); bl.SetFillColor(kLive);
    TBox bf(0, 0, 1, 1); bf.SetFillColor(TColor::GetColor("#f7c59f"));
    lg.AddEntry(gj, "jet cone (R = 0.7)", "l");
    lg.AddEntry(&bf, "jet cone, live cells", "f");
    lg.AddEntry(g1, "band #eta_{jet} #pm R", "l");
    lg.AddEntry(gx, "2R jet exclusion (in the band)", "l");
    if (has_other) lg.AddEntry(go, "other jet, R exclusion", "l");
    if (cap < M_PI - 1e-6) lg.AddEntry(lcut, Form("cut |#Delta#phi| > %.0f#circ (%.0f%% of ring area kept)", cap * 180 / M_PI, 100 * kept), "l");
    lg.AddEntry(&bl, "live ECal, not sampled", "f");
    lg.Draw();
    pm.RedrawAxis();

    pt.cd();
    TH1D hjn(hj), han(ha), hwn(hw);
    han.Scale(hj.Integral() / ha.Integral());
    hjn.SetLineColor(kJet); hjn.SetLineWidth(3); hjn.SetFillColorAlpha(kJet, 0.0);
    han.SetLineColor(TColor::GetColor("#868e96")); han.SetLineWidth(2); han.SetLineStyle(2);
    hwn.SetLineColor(TColor::GetColor("#1c4f8a")); hwn.SetLineWidth(2); hwn.SetMarkerColor(TColor::GetColor("#1c4f8a")); hwn.SetMarkerStyle(20); hwn.SetMarkerSize(1.1);
    hjn.SetTitle(";#eta;live #eta-#phi area per 0.05 in #eta");
    hjn.SetMaximum(1.7 * std::max(hjn.GetMaximum(), han.GetMaximum())); hjn.SetMinimum(0);
    hjn.GetXaxis()->SetTitleSize(0.055); hjn.GetYaxis()->SetTitleSize(0.055); hjn.GetXaxis()->SetLabelSize(0.05); hjn.GetYaxis()->SetLabelSize(0.05);
    hjn.Draw("HIST"); han.Draw("HIST SAME"); hwn.Draw("P SAME");
    TLegend l2(0.18, 0.66, 0.95, 0.87); l2.SetBorderSize(0); l2.SetFillStyle(0); l2.SetTextSize(0.048);
    l2.AddEntry(&hjn, "jet cone, live", "l");
    l2.AddEntry(&han, "ring, live (unweighted, scaled to the jet)", "l");
    l2.AddEntry(&hwn, "ring, live, each point weighted by its #eta bin", "p");
    l2.Draw();
    tx.SetTextSize(0.05); tx.DrawLatexNDC(0.15, 0.93, "#eta distribution: the weights match the ring to the jet cone");

    pb.cd();
    hs.SetTitle(";#eta;h_{jet}(#eta) / h_{ring}(#eta)");
    hs.SetLineColor(TColor::GetColor("#1c4f8a")); hs.SetLineWidth(2);
    hs.SetFillColor(TColor::GetColor("#9ab8d8"));
    hs.GetXaxis()->SetTitleSize(0.055); hs.GetYaxis()->SetTitleSize(0.055); hs.GetXaxis()->SetLabelSize(0.05); hs.GetYaxis()->SetLabelSize(0.05);
    hs.SetMinimum(0); hs.SetMaximum(1.25 * hs.GetMaximum());
    hs.Draw("HIST");
    tx.DrawLatexNDC(0.15, 0.94, Form("point weight r(#eta) per #eta bin, 0 where unreached (jet %.0f%% live)", 100 * pr.area_jet / (M_PI * R * R)));
    TLatex t2; t2.SetTextFont(42); t2.SetTextSize(0.046); t2.SetTextColor(kInk);
    t2.DrawLatexNDC(0.18, 0.80, "#rho = #Sigma_{points} E_{p} r(#eta_{p}) / #Sigma_{bins} r h_{ring},   UE = #rho A_{jet}");

    cv.SaveAs(argv[1]);
    return 0;
}
