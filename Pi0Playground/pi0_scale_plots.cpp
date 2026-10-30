// Plots motivating the ECal photon energy scale k (data photons sit below the
// simulation). Uses pi0_pairs' histograms (ECal gain calibration's selection,
// default clustering, common 45 MeV tower floor), the same inputs and peak fit
// (Gaussian + quadratic) as pi0_overlap_study:
//   1. pi0 mass spectra, data against simulation, in pair-energy bins
//   2. pi0 + eta mass spectra for lower photon energy 3-20 GeV (eta photons
//      are ~4x further apart: little cluster overlap)
//   3. data/sim peak ratio against energy: pi0 vs pair energy (all events, EM
//      bit), pi0 vs single-photon energy, eta vs lower photon energy, with the
//      adopted k = 0.93 +- 0.04 band
//   4. pi0 ratio against pair energy per ECal ring
// The ratio points are also printed as a table.
//
// Usage: ./pi0_scale_plots <data_all.root> <data_em.root> <sim.root> <out.pdf>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TBox.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TError.h>
#include <cmath>
#include <cstdio>

struct Peak { double m, em, s; bool ok; };

static Peak FitG(TH1D *h, double lo, double hi, double mlo, double mhi, double slo, double shi, double seedlo, double seedhi, TF1 *keep = 0)
{
    Peak p = {0, 0, 0, false};
    if (!h || h->GetEffectiveEntries() < 150) return p;
    double best = -1, at = 0.5 * (seedlo + seedhi);
    for (int i = h->FindBin(seedlo); i <= h->FindBin(seedhi); i++) if (h->GetBinContent(i) > best) { best = h->GetBinContent(i); at = h->GetBinCenter(i); }
    TF1 f("f", "[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])^2)", lo, hi);
    f.SetParameters(0.3 * best, 0, 0, 0.5 * best, at, 0.5 * (slo + shi));
    f.SetParLimits(3, 0, 3 * best); f.SetParLimits(4, mlo, mhi); f.SetParLimits(5, slo, shi);
    if (h->Fit(&f, "QRN0") != 0) return p;
    p.m = f.GetParameter(4); p.em = f.GetParError(4); p.s = f.GetParameter(5);
    p.ok = p.em > 0 && p.em < 0.02 && f.GetParameter(3) > 3 * f.GetParError(3) && p.m > mlo + 0.002 && p.m < mhi - 0.002 && p.s > slo * 1.05 && p.s < shi * 0.98;
    if (keep) { keep->SetRange(lo, hi); for (int i = 0; i < 6; i++) keep->SetParameter(i, f.GetParameter(i)); }
    return p;
}
static Peak FitPi0(TH1D *h, TF1 *keep = 0) { return FitG(h, 0.05, 0.30, 0.09, 0.19, 0.006, 0.045, 0.08, 0.20, keep); }
static Peak FitEta(TH1D *h, TF1 *keep = 0) { return FitG(h, 0.35, 0.78, 0.45, 0.65, 0.015, 0.09, 0.45, 0.65, keep); }

static TH1D *P3(TFile &f, const char *name, double y0, double y1, double z0, double z1)
{
    TH3D *h = (TH3D *)f.Get(name);
    if (!h) return 0;
    return h->ProjectionX(Form("%s_%p_%g_%g_%g_%g", name, (void *)&f, y0, y1, z0, z1), h->GetYaxis()->FindBin(y0 + 1e-6), h->GetYaxis()->FindBin(y1 - 1e-6),
                          h->GetZaxis()->FindBin(z0 + 1e-6), h->GetZaxis()->FindBin(z1 - 1e-6), "e");
}
static TH1D *P2(TFile &f, const char *name, double y0, double y1)
{
    TH2D *h = (TH2D *)f.Get(name);
    if (!h) return 0;
    return h->ProjectionX(Form("%s_%p_%g_%g", name, (void *)&f, y0, y1), h->GetYaxis()->FindBin(y0 + 1e-6), h->GetYaxis()->FindBin(y1 - 1e-6), "e");
}

// data/sim ratio point at x (bin [x0, x1]); returns false if either fit failed
static bool AddRatio(TGraphErrors *g, double x0, double x1, const Peak &d, const Peak &s, const char *label)
{
    if (!(d.ok && s.ok)) { std::printf("  %-28s %5.1f-%-5.1f           -\n", label, x0, x1); return false; }
    const double r = d.m / s.m, er = r * std::hypot(d.em / d.m, s.em / s.m);
    const int n = g->GetN();
    g->SetPoint(n, 0.5 * (x0 + x1), r);
    g->SetPointError(n, 0.5 * (x1 - x0), er);
    std::printf("  %-28s %5.1f-%-5.1f  %6.1f | %6.1f MeV  %6.3f +- %5.3f\n", label, x0, x1, 1000 * d.m, 1000 * s.m, r, er);
    return true;
}

// Area-normalise data and sim to the same integral over [lo, hi] and draw
static void Overlay(TH1D *d, TH1D *s, double lo, double hi, const char *title, const Peak &pd, const Peak &ps, bool eta = false, const Peak *ed = 0, const Peak *es = 0)
{
    const double nd = d->Integral(d->FindBin(lo), d->FindBin(hi)), ns = s->Integral(s->FindBin(lo), s->FindBin(hi));
    if (nd > 0) d->Scale(1.0 / nd);
    if (ns > 0) s->Scale(1.0 / ns);
    d->SetTitle(title);
    d->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    d->GetYaxis()->SetTitle("normalised counts");
    d->GetXaxis()->SetRangeUser(lo, hi);
    s->GetXaxis()->SetRangeUser(lo, hi); // GetMaximum below then ignores bins outside the window
    d->SetMarkerStyle(20); d->SetMarkerSize(0.6); d->SetLineColor(kBlack); d->SetMarkerColor(kBlack);
    s->SetLineColor(kRed + 1); s->SetLineWidth(2);
    d->SetMaximum(1.25 * std::max(d->GetMaximum(), s->GetMaximum()));
    d->SetMinimum(0);
    d->Draw("E1");
    s->Draw("HIST SAME");
    const double ymax = d->GetMaximum();
    TLine *l;
    if (pd.ok) { l = new TLine(pd.m, 0, pd.m, ymax); l->SetLineStyle(2); l->SetLineColor(kBlack); l->Draw(); }
    if (ps.ok) { l = new TLine(ps.m, 0, ps.m, ymax); l->SetLineStyle(2); l->SetLineColor(kRed + 1); l->Draw(); }
    TLatex t; t.SetNDC(); t.SetTextSize(0.042);
    if (pd.ok && ps.ok) t.DrawLatex(0.42, 0.82, Form("#pi^{0}: %.1f / %.1f MeV = %.3f", 1000 * pd.m, 1000 * ps.m, pd.m / ps.m));
    if (eta && ed && es && ed->ok && es->ok)
    {
        l = new TLine(ed->m, 0, ed->m, ymax); l->SetLineStyle(2); l->SetLineColor(kBlack); l->Draw();
        l = new TLine(es->m, 0, es->m, ymax); l->SetLineStyle(2); l->SetLineColor(kRed + 1); l->Draw();
        t.DrawLatex(0.42, 0.75, Form("#eta: %.1f / %.1f MeV = %.3f", 1000 * ed->m, 1000 * es->m, ed->m / es->m));
    }
}

static void Style(TGraphErrors *g, int col, int mk)
{
    g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(1.1); g->SetLineWidth(2);
}

// frame + k band + unity line for a ratio pad
static void RatioFrame(double xlo, double xhi, const char *xt, const char *title)
{
    TH1D *fr = new TH1D(Form("fr_%p", (void *)xt), title, 1, xlo, xhi);
    fr->SetDirectory(0);
    fr->SetMinimum(0.82); fr->SetMaximum(1.08);
    fr->GetXaxis()->SetTitle(xt);
    fr->GetYaxis()->SetTitle("peak position data / sim");
    fr->SetStats(0);
    fr->Draw("AXIS");
    TBox *b = new TBox(xlo, 0.89, xhi, 0.97);
    b->SetFillColor(kAzure - 9); // solid: ROOT 5 PNG output drops alpha
    b->SetFillStyle(1001);
    b->Draw();
    TLine *l = new TLine(xlo, 0.93, xhi, 0.93); l->SetLineColor(kAzure + 2); l->SetLineWidth(2); l->Draw();
    l = new TLine(xlo, 1.0, xhi, 1.0); l->SetLineStyle(2); l->Draw();
    fr->Draw("AXIS SAME");
    TLatex t; t.SetNDC(); t.SetTextSize(0.05); t.SetTextAlign(21);
    t.DrawLatex(0.5, 0.93, title);
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError;
    if (argc < 5) { std::printf("usage: ./pi0_scale_plots <data_all.root> <data_em.root> <sim.root> <out.pdf>\n"); return 1; }
    TFile fa(argv[1]), fe(argv[2]), fs(argv[3]);
    const TString pdf = argv[4];
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(1.25, "Y");
    gStyle->SetPaperSize(TStyle::kUSLetter);
    TString png = pdf; png.ReplaceAll(".pdf", "");

    // 1. pi0 spectra in pair-energy bins, all events, all rings
    const int nP = 6;
    const double pE[nP + 1] = {4, 6, 8, 10, 14, 20, 40};
    TCanvas c1("c1", "", 1200, 800);
    c1.Divide(3, 2);
    for (int i = 0; i < nP; i++)
    {
        c1.cd(i + 1);
        TH1D *d = P3(fa, "h3_liangf_m_E_ring_all", pE[i], pE[i + 1], 0, 3), *s = P3(fs, "h3_liangf_m_E_ring_all", pE[i], pE[i + 1], 0, 3);
        const Peak pd = FitPi0(d), ps = FitPi0(s);
        Overlay(d, s, 0.0, 0.34, Form("pair energy %.0f-%.0f GeV", pE[i], pE[i + 1]), pd, ps);
        if (i == 0)
        {
            TLegend *lg = new TLegend(0.42, 0.58, 0.88, 0.74);
            lg->SetBorderSize(0);
            lg->AddEntry(d, "data (run 23007048, all events)", "pe");
            lg->AddEntry(s, "simulation (FMS tune)", "l");
            lg->Draw();
        }
    }
    c1.Print(png + "_1.png");
    c1.Print(pdf + "(");

    // 2. pi0 + eta, lower photon energy 3-20 GeV
    TCanvas c2("c2", "", 1200, 500);
    c2.Divide(2, 1);
    const char *cls[2] = {"all", "em"};
    TFile *fd[2] = {&fa, &fe};
    const char *ctitle[2] = {"all events", "EM-triggered"};
    for (int c = 0; c < 2; c++)
    {
        c2.cd(c + 1);
        TH1D *d = P2(*fd[c], Form("h2_leta_m_Eg_%s", cls[c]), 3, 20), *s = P2(fs, Form("h2_leta_m_Eg_%s", cls[c]), 3, 20);
        const Peak pd = FitPi0(d), ps = FitPi0(s), ed = FitEta(d), es = FitEta(s);
        d->Rebin(2); s->Rebin(2);
        Overlay(d, s, 0.0, 0.78, Form("lower photon energy 3-20 GeV, %s", ctitle[c]), pd, ps, true, &ed, &es);
        gPad->SetLogy();
        d->SetMinimum(1e-4);
    }
    c2.Print(png + "_2.png");
    c2.Print(pdf);

    // 3. ratio summary
    std::printf("data/sim peak ratio (data | sim peak, ratio)\n");
    TGraphErrors *gPairAll = new TGraphErrors(), *gPairEm = new TGraphErrors(), *gGam = new TGraphErrors(), *gEta = new TGraphErrors();
    const int nR = 9;
    const double rE[nR + 1] = {4, 6, 8, 10, 12, 14, 16, 20, 30, 40};
    for (int i = 0; i < nR; i++)
    {
        AddRatio(gPairAll, rE[i], rE[i + 1], FitPi0(P3(fa, "h3_liangf_m_E_ring_all", rE[i], rE[i + 1], 0, 3)), FitPi0(P3(fs, "h3_liangf_m_E_ring_all", rE[i], rE[i + 1], 0, 3)), "pi0 vs pair E, all");
    }
    for (int i = 0; i < nR; i++)
        AddRatio(gPairEm, rE[i], rE[i + 1], FitPi0(P3(fe, "h3_liangf_m_E_ring_em", rE[i], rE[i + 1], 0, 3)), FitPi0(P3(fs, "h3_liangf_m_E_ring_em", rE[i], rE[i + 1], 0, 3)), "pi0 vs pair E, EM bit");
    const int nG = 8;
    const double gE[nG + 1] = {1, 2, 3, 4, 5, 6, 8, 12, 20};
    for (int i = 0; i < nG; i++)
        AddRatio(gGam, gE[i], gE[i + 1], FitPi0(P3(fa, "h3_lph_m_Eg_Z_all", gE[i], gE[i + 1], 0, 0.7)), FitPi0(P3(fs, "h3_lph_m_Eg_Z_all", gE[i], gE[i + 1], 0, 0.7)), "pi0 vs photon E, all");
    const int nL = 3;
    const double lE[nL + 1] = {3, 5, 8, 20};
    for (int i = 0; i < nL; i++)
        AddRatio(gEta, lE[i], lE[i + 1], FitEta(P2(fa, "h2_leta_m_Eg_all", lE[i], lE[i + 1])), FitEta(P2(fs, "h2_leta_m_Eg_all", lE[i], lE[i + 1])), "eta vs lower photon E, all");
    Style(gPairAll, kBlack, 20); Style(gPairEm, kGray + 1, 24); Style(gGam, kRed + 1, 21); Style(gEta, kGreen + 2, 33);
    gEta->SetMarkerSize(1.8);

    TCanvas c3("c3", "", 1200, 500);
    c3.Divide(2, 1);
    c3.cd(1);
    RatioFrame(0, 40, "pair energy (GeV)", "#pi^{0} peak against pair energy");
    gPairAll->Draw("P SAME"); gPairEm->Draw("P SAME");
    TLegend *l3 = new TLegend(0.40, 0.15, 0.88, 0.32);
    l3->SetBorderSize(0); l3->SetFillStyle(0);
    l3->AddEntry(gPairAll, "#pi^{0}, all events", "pe");
    l3->AddEntry(gPairEm, "#pi^{0}, EM-triggered", "pe");
    l3->AddEntry((TObject *)0, "band: k = 0.93 #pm 0.04", "");
    l3->Draw();
    c3.cd(2);
    RatioFrame(0, 20, "photon energy (GeV)", "against single-photon energy");
    gGam->Draw("P SAME"); gEta->Draw("P SAME");
    TLegend *l4 = new TLegend(0.40, 0.15, 0.88, 0.32);
    l4->SetBorderSize(0); l4->SetFillStyle(0);
    l4->AddEntry(gGam, "#pi^{0} vs each photon's energy", "pe");
    l4->AddEntry(gEta, "#eta vs lower photon energy", "pe");
    l4->Draw();
    c3.Print(png + "_3.png");
    c3.Print(pdf);

    // 4. per ring
    TCanvas c4("c4", "", 800, 600);
    RatioFrame(0, 40, "pair energy (GeV)", "#pi^{0} peak per ECal ring, all events");
    const char *rname[3] = {"inner", "middle", "outer"};
    const int rcol[3] = {kBlue + 1, kOrange + 7, kMagenta + 1}, rmk[3] = {20, 21, 22};
    TLegend *l5 = new TLegend(0.55, 0.15, 0.88, 0.32);
    l5->SetBorderSize(0); l5->SetFillStyle(0);
    const int nQ = 6;
    const double qE[nQ + 1] = {4, 8, 12, 16, 20, 30, 40};
    for (int r = 0; r < 3; r++)
    {
        TGraphErrors *g = new TGraphErrors();
        for (int i = 0; i < nQ; i++)
            AddRatio(g, qE[i], qE[i + 1], FitPi0(P3(fa, "h3_liangf_m_E_ring_all", qE[i], qE[i + 1], r, r + 1)), FitPi0(P3(fs, "h3_liangf_m_E_ring_all", qE[i], qE[i + 1], r, r + 1)), Form("pi0 vs pair E, %s ring", rname[r]));
        for (int i = 0; i < g->GetN(); i++) g->SetPoint(i, g->GetX()[i] + 0.6 * (r - 1), g->GetY()[i]);
        Style(g, rcol[r], rmk[r]);
        g->Draw("P SAME");
        l5->AddEntry(g, rname[r], "pe");
    }
    l5->Draw();
    c4.Print(png + "_4.png");
    c4.Print(pdf + ")");
    return 0;
}
