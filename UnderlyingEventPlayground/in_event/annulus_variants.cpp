// The six annulus variants (in_event_geo.h) for one example jet on the FCS
// ECal: columns whole ring / transverse / away side, rows with and without the
// R exclusion around other reco jets. Each panel shades the sampled cells by
// their point weight (common scale), marks the jet cone's live cells whose eta
// strip the annulus does not cover (their UE comes from the template
// scale-up), and states the covered share of the jet template, the scale-up
// an_fill and whether the variant passes the validity cut (share >= 70%).
//
// Build and run from this folder:
//   g++ -O2 annulus_variants.cpp -o annulus_variants -I../../shared -I../../underlying_event `root-config --cflags --libs` -std=c++11
//   ./annulus_variants <out.png> [jet eta = 3.0] [jet phi deg = 0] [other eta = 3.2] [other phi deg = 100]

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "in_event_geo.h"

static const double kMaxFill = 1 / 0.7;

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
    if (argc < 2) { std::fprintf(stderr, "usage: ./annulus_variants <out.png> [jet eta] [jet phi deg] [other eta] [other phi deg]\n"); return 1; }
    const double jeta = argc > 2 ? std::atof(argv[2]) : 3.0, jphi = (argc > 3 ? std::atof(argv[3]) : 0.0) * M_PI / 180.0;
    const double oeta = argc > 4 ? std::atof(argv[4]) : 3.2, ophi = (argc > 5 ? std::atof(argv[5]) : 100.0) * M_PI / 180.0;

    Geometry geo; geo.Build();
    std::vector<double> rho(geo.cells.size(), 1.0);  // flat template: the covered share is an area share
    std::vector<int> cl;
    geo.Near(jeta, jphi, R, cl);
    std::vector<char> in_jet(geo.cells.size(), 0);
    double Ajs[kNS] = {0}, Tjs[kNS] = {0}, T = 0;
    for (size_t q = 0; q < cl.size(); q++)
    {
        const Cell &c = geo.cells[cl[q]];
        in_jet[cl[q]] = 1;
        const int s = StripOf(c.eta - jeta);
        if (s >= 0) { Ajs[s] += c.area; Tjs[s] += c.area; T += c.area; }
    }
    std::vector<std::pair<float, float> > others(1, std::make_pair((float)oeta, (float)ophi)), none;
    Annulus an[2][kNA];
    double wmax = 0;
    for (int x = 0; x < 2; x++)
        for (int v = 0; v < kNA; v++)
        {
            BuildAnnulus(geo, rho, jeta, jphi, x == 0 ? others : none, Ajs, T, Tjs, v, an[x][v]);
            if (an[x][v].fill > 0 && an[x][v].fill <= kMaxFill)  // scale set by the valid variants; others saturate
                for (size_t c = 0; c < geo.cells.size(); c++) wmax = std::max(wmax, (double)an[x][v].w[c]);
        }

    const int NC = 50;
    double stops[3] = {0.0, 0.5, 1.0}, red[3] = {0.87, 0.42, 0.03}, green[3] = {0.92, 0.61, 0.19}, blue[3] = {0.97, 0.84, 0.42};
    const int c0 = TColor::CreateGradientColorTable(3, stops, red, green, blue, NC);
    gStyle->SetOptStat(0);
    {
        static int pal[NC];
        for (int i = 0; i < NC; i++) pal[i] = c0 + i;
        gStyle->SetPalette(NC, pal);
    }
    const int kLive = TColor::GetColor("#e3e3e3"), kJet = TColor::GetColor("#d9480f"), kInk = TColor::GetColor("#222222");
    const int kJetCov = TColor::GetColor("#f7c59f"), kJetUncov = TColor::GetColor("#c92a2a");
    const int kPass = TColor::GetColor("#2b8a3e"), kFail = TColor::GetColor("#c92a2a");

    if (!(wmax > 0)) wmax = 1;
    TCanvas cv("cv", "", 1950, 1400);
    TPad top("top", "", 0, 0.93, 1, 1), grid("grid", "", 0, 0.08, 0.93, 0.93), side("side", "", 0.93, 0.08, 1, 0.93), bot("bot", "", 0, 0, 1, 0.08);
    top.Draw(); grid.Draw(); side.Draw(); bot.Draw();
    grid.cd(); grid.Divide(3, 2, 0.002, 0.002);
    const double dx = (kFiducialRectEcal.x_outer - kFiducialRectEcal.x_inner) / geo.nx, dy = (kFiducialRectEcal.y_max - kFiducialRectEcal.y_min) / geo.ny;
    const char *vname[kNA] = {"whole ring", "transverse (60-120#circ)", "away side (> 120#circ)"};
    const char *xname[2] = {"other jets excluded (R)", "other jets not excluded"};
    const double r_in = z_proj / std::sinh(jeta + R), r_out = z_proj / std::sinh(jeta - R);
    std::vector<TObject *> keep;
    for (int x = 0; x < 2; x++)
        for (int v = 0; v < kNA; v++)
        {
            grid.cd(x * kNA + v + 1);
            gPad->SetMargin(0.09, 0.02, 0.07, 0.17);
            const Annulus &A = an[x][v];
            TH2D *fr = new TH2D(Form("fr%d%d", x, v), ";x (cm);y (cm)", 10, -150, 150, 10, -105, 100);
            fr->GetXaxis()->SetLabelSize(0.04); fr->GetYaxis()->SetLabelSize(0.04);
            fr->GetXaxis()->SetTitleSize(0.04); fr->GetYaxis()->SetTitleSize(0.04); fr->GetYaxis()->SetTitleOffset(1.0); fr->GetXaxis()->SetTitleOffset(0.8);
            fr->Draw(); keep.push_back(fr);
            double cov_share = 0;
            for (int s = 0; s < kNS; s++) if (A.cov[s]) cov_share += Tjs[s];
            cov_share = T > 0 ? cov_share / T : 0;
            for (size_t c = 0; c < geo.cells.size(); c++)
            {
                const Cell &cc = geo.cells[c];
                int col = kLive;
                if (A.w[c] > 0) col = c0 + std::min(NC - 1, (int)(A.w[c] / wmax * (NC - 1) + 0.5));
                if (in_jet[c])
                {
                    const int s = StripOf(cc.eta - jeta);
                    col = (s >= 0 && A.cov[s]) ? kJetCov : kJetUncov;
                }
                TBox *b = new TBox(cc.x - dx / 2, cc.y - dy / 2, cc.x + dx / 2, cc.y + dy / 2);
                b->SetFillColor(col); b->SetLineWidth(0); b->Draw(); keep.push_back(b);
            }
            TGraph *gj = Circle(jeta, jphi, R); gj->SetLineColor(kJet); gj->SetLineWidth(3); gj->Draw("L"); keep.push_back(gj);
            std::vector<TGraph *> gx = CircleInBand(jeta, jphi, kAnnJetExcl, jeta);
            for (size_t q = 0; q < gx.size(); q++) { gx[q]->SetLineColor(kJet); gx[q]->SetLineWidth(2); gx[q]->SetLineStyle(2); gx[q]->Draw("L"); }
            TGraph *go = Circle(oeta, ophi, R); go->SetLineColor(kInk); go->SetLineWidth(2); go->SetLineStyle(x == 0 ? 2 : 3); go->Draw("L"); keep.push_back(go);
            TGraph *g1 = Ring(jeta - R), *g2 = Ring(jeta + R);
            g1->SetLineColor(kInk); g2->SetLineColor(kInk); g1->Draw("L"); g2->Draw("L");
            for (int sgn = -1; sgn <= 1; sgn += 2)
                for (int k = 1; k <= 2; k++)
                {
                    const double p = jphi + sgn * k * M_PI / 3;
                    TLine *l = new TLine(r_in * std::cos(p), r_in * std::sin(p), r_out * std::cos(p), r_out * std::sin(p));
                    l->SetLineColor(kInk); l->SetLineStyle(7); l->Draw(); keep.push_back(l);
                }
            TLatex t; t.SetTextFont(42); t.SetTextAlign(22); t.SetTextSize(0.045);
            const double rj = z_proj / std::sinh(jeta), ro = z_proj / std::sinh(oeta);
            t.SetTextColor(kJet); t.DrawLatex(rj * std::cos(jphi), rj * std::sin(jphi) + 8, "jet");
            t.SetTextColor(kInk); t.DrawLatex(ro * std::cos(ophi), ro * std::sin(ophi), "other jet");
            t.SetNDC(); t.SetTextAlign(12); t.SetTextSize(0.05);
            t.DrawLatex(0.09, 0.95, Form("#bf{%s}, %s", vname[v], xname[x]));
            const bool pass = A.fill > 0 && A.fill <= kMaxFill;
            t.SetTextSize(0.042);
            t.DrawLatex(0.09, 0.895, Form("covered strips: %.0f%% of the jet cone, an_fill = %.2f", 100 * cov_share, A.fill));
            t.SetTextColor(pass ? kPass : kFail);
            t.DrawLatex(0.09, 0.85, pass ? "#bf{valid} (#geq 70%)" : "#bf{rejected} by the validity cut (< 70%)");
            gPad->RedrawAxis();
            std::printf("%-26s %-26s covered %.2f fill %.3f live area %.3f %s\n", vname[v], xname[x], cov_share, A.fill, A.area_live, pass ? "valid" : "rejected");
        }

    side.cd();
    side.SetMargin(0.05, 0.62, 0.10, 0.10);
    TH2D *pz = new TH2D("pz", "", 1, 0, 1, NC, 0, wmax);
    for (int i = 0; i < NC; i++) pz->SetBinContent(1, i + 1, (i + 0.5) * wmax / NC);
    pz->SetMinimum(0); pz->SetMaximum(wmax);
    pz->GetXaxis()->SetLabelSize(0); pz->GetXaxis()->SetTickLength(0);
    pz->GetYaxis()->SetLabelSize(0);
    pz->Draw("COL"); keep.push_back(pz);
    TGaxis *ax = new TGaxis(1, 0, 1, wmax, 0, wmax, 510, "+L");
    ax->SetLabelSize(0.20); ax->SetLabelFont(42); ax->SetTitle("point weight w_{s} (saturates)"); ax->SetTitleSize(0.20); ax->SetTitleFont(42); ax->SetTitleOffset(1.45);
    ax->Draw();

    top.cd();
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.30); t.SetTextColor(kInk);
    t.DrawLatex(0.01, 0.60, Form("Annulus variants for an example jet at #eta = %.1f, #phi = %.0f#circ, other jet at #eta = %.1f, #phi = %.0f#circ", jeta, jphi * 180 / M_PI, oeta, ophi * 180 / M_PI));
    t.SetTextSize(0.24);
    t.DrawLatex(0.01, 0.18, "Inside the band the 2R jet exclusion needs |#Delta#phi| #geq #sqrt{4R^{2} - #Delta#eta^{2}} #geq 69#circ, so the transverse range is effectively ~80-120#circ; red jet cells: strips the annulus does not cover, extrapolated with an_fill");

    bot.cd();
    TLegend lg(0.01, 0.05, 0.99, 0.95);
    lg.SetNColumns(4); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.24);
    TBox b1(0, 0, 1, 1); b1.SetFillColor(kJetCov);
    TBox b2(0, 0, 1, 1); b2.SetFillColor(kJetUncov);
    TBox b3(0, 0, 1, 1); b3.SetFillColor(kLive);
    TGraph lj; lj.SetLineColor(kJet); lj.SetLineWidth(3);
    TGraph lx; lx.SetLineColor(kJet); lx.SetLineWidth(2); lx.SetLineStyle(2);
    TGraph lo; lo.SetLineColor(kInk); lo.SetLineWidth(2); lo.SetLineStyle(2);
    TGraph lb; lb.SetLineColor(kInk);
    TGraph lr; lr.SetLineColor(kInk); lr.SetLineStyle(7);
    lg.AddEntry(&lj, "jet cone (R = 0.7)", "l");
    lg.AddEntry(&b1, "jet cone, live, strip covered", "f");
    lg.AddEntry(&b2, "jet cone, live, strip not covered", "f");
    lg.AddEntry(&b3, "live ECal, not sampled", "f");
    lg.AddEntry(&lb, "band #eta_{jet} #pm R", "l");
    lg.AddEntry(&lx, "2R jet exclusion", "l");
    lg.AddEntry(&lo, "other jet, R (dashed = excluded)", "l");
    lg.AddEntry(&lr, "#Delta#phi = 60#circ, 120#circ", "l");
    lg.Draw();

    cv.SaveAs(argv[1]);
    return 0;
}
