// Fits the pi0 peak of pi0_pairs outputs, data against simulation, and splits
// a data/simulation difference into energy scale and two-photon separation.
//
// Peaks: Gaussian + quadratic background on [kFitLo, kFitHi]; the peak
// position and width are quoted with their fit errors.
// Tables:
//   1. per pair energy bin (Z < 0.7, all d): data, data with the event vertex,
//      simulation, tagged simulation pi0s; data / simulation - 1
//   2. per pair energy bin and point distance d: data / simulation - 1 of the
//      peak (does the difference follow E or d?)
//   3. per pair energy bin and asymmetry Z: the same; 3b. by cluster class
//      (both points from one cluster, or from two); 3c. the eta peak
//   4. tagged simulation pi0s: peak of the reco mass, of the mass with the
//      true photon energies (positions only) and with the true positions
//      (energies only), each / m_pi0; point E / photon E (median) per photon
//      energy; reco - true two-photon distance per true distance
// Writes <base>.pdf: the mass spectra of data and simulation per energy bin,
// and per energy and d bin.
//
// Usage: ./pi0_compare <data.root> <sim.root> <base>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TError.h>
#include <cmath>
#include <cstdio>
#include <string>

static const int kNE = 7;
static const double kEEdges[kNE + 1] = {5, 10, 15, 20, 30, 45, 70, 100};
static const int kND = 7;
static const double kDEdges[kND + 1] = {0, 6, 9, 12, 16, 22, 30, 60};
static const int kNZ = 5;
static const double kZEdges[kNZ + 1] = {0, 0.2, 0.4, 0.55, 0.7, 1.0};
static const double kMPi0 = 0.13498;
static const double kFitLo = 0.04, kFitHi = 0.30;
static const double kMEta = 0.54786, kEtaLo = 0.36, kEtaHi = 0.78;

struct Peak { double m, em, s, n; bool ok; };

static Peak Fit(TH1D *h)
{
    Peak p = {0, 0, 0, 0, false};
    if (!h) return p;
    double neff = h->GetEffectiveEntries();
    p.n = neff;
    if (neff < 50) return p;
    TF1 f("f", "gaus(0)+pol2(3)", kFitLo, kFitHi);
    const int b = h->FindBin(kMPi0);
    f.SetParameters(h->GetBinContent(b), 0.13, 0.02, 0, 0, 0);
    f.SetParLimits(0, 0, 10 * h->GetMaximum());
    f.SetParLimits(1, 0.07, 0.20);
    f.SetParLimits(2, 0.006, 0.06);
    if (h->Fit(&f, "QRN0") != 0) return p;
    p.m = f.GetParameter(1); p.em = f.GetParError(1); p.s = f.GetParameter(2);
    p.ok = p.em > 0 && p.em < 0.02;
    return p;
}

static Peak FitEta(TH1D *h)
{
    Peak p = {0, 0, 0, 0, false};
    if (!h || h->GetEffectiveEntries() < 200) return p;
    TF1 f("fe", "gaus(0)+pol2(3)", kEtaLo, kEtaHi);
    f.SetParameters(0.2 * h->GetBinContent(h->FindBin(kMEta)), 0.52, 0.04, h->GetBinContent(h->FindBin(kMEta)), 0, 0);
    f.SetParLimits(0, 0, 10 * h->GetMaximum());
    f.SetParLimits(1, 0.42, 0.66);
    f.SetParLimits(2, 0.015, 0.09);
    if (h->Fit(&f, "QRN0") != 0) return p;
    p.m = f.GetParameter(1); p.em = f.GetParError(1); p.s = f.GetParameter(2);
    p.ok = p.em > 0 && p.em < 0.03 && f.GetParameter(0) > 2 * f.GetParError(0);
    return p;
}

// A. Ogawa's massfit.C model (FCSJetAnalysis/reference/massfit.C): Gaussian
// on A/tau^2 x^p exp(-x/tau), background fitted alone first, on [0, 0.4]
static double AkioBg(double *x, double *p) { return x[0] <= 0 ? 0 : p[0] / p[1] / p[1] * std::pow(x[0], p[2]) * std::exp(-x[0] / p[1]); }
static double AkioAll(double *x, double *p) { return AkioBg(x, p) + p[3] * std::exp(-0.5 * std::pow((x[0] - p[4]) / p[5], 2)); }
static Peak FitAkio(TH1D *h)
{
    Peak p = {0, 0, 0, 0, false};
    if (!h || h->GetEffectiveEntries() < 150) return p;
    const double mx = h->GetMaximum();
    double peak = h->GetBinCenter(h->GetMaximumBin());
    if (peak < 0.05 || peak > 0.20) peak = 0.135;
    TF1 b("akbg", AkioBg, 0, 0.4, 3);
    b.SetParameters(mx / 2, 0.2, 1.0);
    h->Fit(&b, "QRN0");
    TF1 f("akall", AkioAll, 0, 0.4, 6);
    f.SetParameters(b.GetParameter(0), b.GetParameter(1), b.GetParameter(2), mx / 2, peak, 0.04);
    f.SetParLimits(3, 0, mx); f.SetParLimits(5, 0.005, 0.1); f.SetParLimits(4, 0.06, 0.25);
    if (h->Fit(&f, "QRN0") != 0) return p;
    p.m = f.GetParameter(4); p.em = f.GetParError(4); p.s = f.GetParameter(5);
    p.ok = p.em > 0 && p.em < 0.02;
    return p;
}

static TH1D *Proj(TFile &f, const char *name, int ie, int j0, int j1, const char *tag)
{
    TH3D *h = (TH3D *)f.Get(name);
    if (!h) return 0;
    TH1D *p = h->ProjectionX(Form("%s_%s_%d_%d_%d", name, tag, ie, j0, j1), ie + 1, ie + 1, j0 + 1, j1 + 1, "e");
    return p;
}

static void Cell(const Peak &d, const Peak &s)
{
    if (d.ok && s.ok)
    {
        const double r = d.m / s.m - 1, er = (d.m / s.m) * std::hypot(d.em / d.m, s.em / s.m);
        std::printf("  %+6.1f +- %4.1f", 100 * r, 100 * er);
    }
    else std::printf("  %14s", "-");
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);
    if (argc < 4) { std::printf("usage: ./pi0_compare <data.root> <sim.root> <base>\n"); return 1; }
    TFile fd(argv[1]), fs(argv[2]);
    const std::string base = argv[3];

    std::printf("1. Peak per pair energy (Z < 0.7, all d): position (MeV) +- error [width], and data/sim - 1 (%%)\n");
    std::printf("  %-9s %22s %22s %22s %22s %16s %16s\n", "E (GeV)", "data", "data, vertex", "sim", "sim tagged pi0", "data/sim (%)", "sim/m_pi0 (%)");
    TCanvas cv("cv", "", 1600, 900);
    cv.Print((base + ".pdf[").c_str());
    cv.Divide(4, 2);
    for (int e = 0; e < kNE; e++)
    {
        TH1D *hd = Proj(fd, "h3_m_E_d", e, 0, kND - 1, "d"), *hv = Proj(fd, "h3_m_E_d_vtx", e, 0, kND - 1, "v");
        TH1D *hs = Proj(fs, "h3_m_E_d", e, 0, kND - 1, "s"), *ht = Proj(fs, "h3_m_E_d_tag", e, 0, kND - 1, "t");
        Peak pd = Fit(hd), pv = Fit(hv), ps = Fit(hs), pt = Fit(ht);
        std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
        const Peak *all[4] = {&pd, &pv, &ps, &pt};
        for (int k = 0; k < 4; k++)
            if (all[k]->ok) std::printf("   %6.1f +- %4.1f [%4.1f]", 1000 * all[k]->m, 1000 * all[k]->em, 1000 * all[k]->s);
            else std::printf("   %22s", "-");
        Cell(pd, ps);
        if (ps.ok) std::printf("  %+6.1f", 100 * (ps.m / kMPi0 - 1)); else std::printf("  %6s", "-");
        std::printf("\n");

        cv.cd(e + 1);
        if (hd && hs && hd->Integral() > 0 && hs->Integral() > 0)
        {
            hd->Scale(1.0 / hd->Integral()); hs->Scale(1.0 / hs->Integral());
            hs->SetLineColor(kRed + 1); hs->SetFillColorAlpha(kRed - 9, 0.4);
            hd->SetMarkerStyle(20); hd->SetMarkerSize(0.6);
            hs->SetMaximum(1.25 * std::max(hs->GetMaximum(), hd->GetMaximum())); hs->SetMinimum(0);
            hs->SetTitle(Form("pair E %.0f-%.0f GeV, Z < 0.7;m_{#gamma#gamma} (GeV);fraction", kEEdges[e], kEEdges[e + 1]));
            hs->Draw("hist e2"); hd->Draw("e same");
            if (e == 0) { TLegend *l = new TLegend(0.55, 0.75, 0.88, 0.88); l->AddEntry(hd, "data", "p"); l->AddEntry(hs, "trig. sim", "f"); l->SetBorderSize(0); l->Draw(); }
        }
    }
    cv.Print((base + ".pdf").c_str());

    std::printf("\n2. data/sim - 1 of the peak (%%), per pair energy (rows) and point distance d at the ECal (columns, cm)\n");
    std::printf("  %-9s", "E (GeV)");
    for (int d = 0; d < kND; d++) std::printf("  %14s", Form("d %.0f-%.0f", kDEdges[d], kDEdges[d + 1]));
    std::printf("\n");
    cv.Clear(); cv.Divide(kND, kNE, 0.001, 0.001);
    for (int e = 0; e < kNE; e++)
    {
        std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
        for (int d = 0; d < kND; d++)
        {
            TH1D *hd = Proj(fd, "h3_m_E_d", e, d, d, "dd"), *hs = Proj(fs, "h3_m_E_d", e, d, d, "sd");
            Cell(Fit(hd), Fit(hs));
            cv.cd(e * kND + d + 1);
            if (hd && hs && hd->Integral() > 0 && hs->Integral() > 0)
            {
                hd->Scale(1.0 / hd->Integral()); hs->Scale(1.0 / hs->Integral());
                hs->SetLineColor(kRed + 1); hs->SetFillColorAlpha(kRed - 9, 0.4); hd->SetMarkerStyle(20); hd->SetMarkerSize(0.3);
                hs->SetMaximum(1.25 * std::max(hs->GetMaximum(), hd->GetMaximum())); hs->SetMinimum(0);
                hs->SetTitle(""); hs->Draw("hist"); hd->Draw("e same");
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.09); tx.DrawLatex(0.35, 0.82, Form("E %.0f-%.0f, d %.0f-%.0f", kEEdges[e], kEEdges[e + 1], kDEdges[d], kDEdges[d + 1]));
            }
        }
        std::printf("\n");
    }
    cv.Print((base + ".pdf").c_str());
    cv.Print((base + ".pdf]").c_str());

    std::printf("\n3. data/sim - 1 of the peak (%%), per pair energy (rows) and asymmetry Z (columns)\n");
    std::printf("  %-9s", "E (GeV)");
    for (int z = 0; z < kNZ; z++) std::printf("  %14s", Form("Z %.2f-%.2f", kZEdges[z], kZEdges[z + 1]));
    std::printf("\n");
    for (int e = 0; e < kNE; e++)
    {
        std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
        for (int z = 0; z < kNZ; z++) Cell(Fit(Proj(fd, "h3_m_E_Z", e, z, z, "dz")), Fit(Proj(fs, "h3_m_E_Z", e, z, z, "sz")));
        std::printf("\n");
    }

    std::printf("\n3b. Peak by cluster class (Z < 0.7): position (MeV), and data/sim - 1 (%%); one cluster = split by the shower-shape fit\n");
    std::printf("  %-9s %16s %16s %16s   %16s %16s %16s\n", "E (GeV)", "1 cl: data", "sim", "data/sim", "2 cl: data", "sim", "data/sim");
    for (int e = 0; e < kNE; e++)
    {
        std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
        for (int c = 0; c < 2; c++)
        {
            Peak pd = Fit(Proj(fd, "h3_m_E_c", e, c, c, "dc")), ps = Fit(Proj(fs, "h3_m_E_c", e, c, c, "sc"));
            if (pd.ok) std::printf("   %6.1f +- %4.1f", 1000 * pd.m, 1000 * pd.em); else std::printf("   %14s", "-");
            if (ps.ok) std::printf("   %6.1f +- %4.1f", 1000 * ps.m, 1000 * ps.em); else std::printf("   %14s", "-");
            Cell(pd, ps);
            std::printf(" ");
        }
        std::printf("\n");
    }

    std::printf("\n3d. With the dropped shower tails added back (hits %.2f-%.2f GeV within %.0f cm): peak (MeV) and data/sim - 1 (%%); mean tail fraction of point E\n", 0.05, 0.4, 8.0);
    for (int e = 0; e < kNE; e++)
    {
        Peak pd = Fit(Proj(fd, "h3_m_E_tail", e, 0, kND - 1, "dt")), ps = Fit(Proj(fs, "h3_m_E_tail", e, 0, kND - 1, "st"));
        std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
        if (pd.ok) std::printf("   data %6.1f +- %4.1f", 1000 * pd.m, 1000 * pd.em); else std::printf("   data %14s", "-");
        if (ps.ok) std::printf("   sim %6.1f +- %4.1f", 1000 * ps.m, 1000 * ps.em); else std::printf("   sim %14s", "-");
        Cell(pd, ps);
        std::printf("\n");
    }
    {
        TH2D *td = (TH2D *)fd.Get("h2_tail_E"), *ts = (TH2D *)fs.Get("h2_tail_E");
        if (td && ts)
        {
            std::printf("   mean tail / point E per point energy (GeV), data | sim:");
            const double pb[6] = {1, 3, 5, 10, 20, 40};
            for (int k = 0; k < 5; k++)
            {
                const int b0 = td->GetXaxis()->FindBin(pb[k] + 1e-6), b1 = td->GetXaxis()->FindBin(pb[k + 1] - 1e-6);
                TH1D *a1 = td->ProjectionY(Form("tda%d", k), b0, b1), *a2 = ts->ProjectionY(Form("tsa%d", k), b0, b1);
                std::printf("  %.0f-%.0f: %.3f | %.3f", pb[k], pb[k + 1], a1->GetMean(), a2->GetMean());
            }
            std::printf("\n");
        }
    }

    std::printf("\n3f. Clusters in place of points (every pair, Z < 0.7, same cuts): peak (MeV), data/sim - 1 (%%), and the point result for comparison\n");
    for (int e = 0; e < kNE; e++)
    {
        Peak pd = Fit(Proj(fd, "h3_mclu_E_d", e, 0, kND - 1, "dcl")), ps = Fit(Proj(fs, "h3_mclu_E_d", e, 0, kND - 1, "scl"));
        Peak qd = Fit(Proj(fd, "h3_m_E_d", e, 0, kND - 1, "dpt")), qs = Fit(Proj(fs, "h3_m_E_d", e, 0, kND - 1, "spt"));
        std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
        if (pd.ok) std::printf("   clusters: data %6.1f", 1000 * pd.m); else std::printf("   clusters: data %6s", "-");
        if (ps.ok) std::printf("  sim %6.1f", 1000 * ps.m); else std::printf("  sim %6s", "-");
        Cell(pd, ps);
        std::printf("   | points:");
        Cell(qd, qs);
        std::printf("\n");
    }

    std::printf("\n3e. Gain-calibration selection (StFcsPi0FinderForEcal: best cluster pair, clusters > 1 GeV, Z < 0.7, <= 40 hits > 1 GeV), fitted with A. Ogawa's massfit model\n");
    {
        const double eb[6] = {0, 5, 10, 15, 20, 40};
        for (int k = 0; k < 6; k++)
        {
            const double lo = k < 5 ? eb[k] : 0, hi = k < 5 ? eb[k + 1] : 40;
            Peak pk[2];
            for (int f = 0; f < 2; f++)
            {
                TH2D *h = (TH2D *)(f == 0 ? fd : fs).Get("h2_liang_m_E");
                pk[f] = h ? FitAkio(h->ProjectionX(Form("li%d_%d", f, k), h->GetYaxis()->FindBin(lo + 1e-6), h->GetYaxis()->FindBin(hi - 1e-6), "e")) : Peak();
            }
            std::printf("  pair E %4.0f-%-4.0f", lo, hi);
            for (int f = 0; f < 2; f++)
                if (pk[f].ok) std::printf("   %s %6.1f +- %4.1f [%4.1f] (%+5.1f%% of m_pi0)", f ? "sim" : "data", 1000 * pk[f].m, 1000 * pk[f].em, 1000 * pk[f].s, 100 * (pk[f].m / kMPi0 - 1));
                else std::printf("   %s %37s", f ? "sim" : "data", "-");
            Cell(pk[0], pk[1]);
            std::printf("\n");
        }
    }

    std::printf("\n3c. eta peak (Z < 0.7, all d, fit %.2f-%.2f GeV): position (MeV) +- error [width], and data/sim - 1 (%%)\n", kEtaLo, kEtaHi);
    for (int e = 0; e < kNE; e++)
    {
        Peak pd = FitEta(Proj(fd, "h3_m_E_d", e, 0, kND - 1, "de")), ps = FitEta(Proj(fs, "h3_m_E_d", e, 0, kND - 1, "se"));
        std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
        if (pd.ok) std::printf("   data %6.1f +- %4.1f [%4.1f]", 1000 * pd.m, 1000 * pd.em, 1000 * pd.s); else std::printf("   data %22s", "-");
        if (ps.ok) std::printf("   sim %6.1f +- %4.1f [%4.1f]", 1000 * ps.m, 1000 * ps.em, 1000 * ps.s); else std::printf("   sim %22s", "-");
        Cell(pd, ps);
        std::printf("\n");
    }

    std::printf("\n4. Tagged simulation pi0s, peak / m_pi0 - 1 (%%): reco mass; true energies at reco positions (positions only); reco energies at true positions (energies only)\n");
    std::printf("  %-9s %16s %16s %16s\n", "E (GeV)", "reco", "positions only", "energies only");
    for (int e = 0; e < kNE; e++)
    {
        std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
        const char *nm[3] = {"h3_m_E_d_tag", "h3_mt_E_d_tag", "h3_mp_E_d_tag"};
        for (int k = 0; k < 3; k++)
        {
            Peak p = Fit(Proj(fs, nm[k], e, 0, kND - 1, "tg"));
            if (p.ok) std::printf("   %+6.1f +- %4.1f", 100 * (p.m / kMPi0 - 1), 100 * p.em / kMPi0); else std::printf("   %14s", "-");
        }
        std::printf("\n");
    }
    TH2D *hr = (TH2D *)fs.Get("h2_resp_E");
    if (hr)
    {
        std::printf("\n   simulation point E / photon E, median per photon energy (GeV):");
        const double eb[8] = {2, 5, 10, 15, 20, 30, 45, 70};
        for (int k = 0; k < 7; k++)
        {
            TH1D *p = hr->ProjectionY(Form("r%d", k), hr->GetXaxis()->FindBin(eb[k] + 1e-6), hr->GetXaxis()->FindBin(eb[k + 1] - 1e-6));
            double q = 0.5, med = 0;
            if (p->Integral() > 0) { p->GetQuantiles(1, &med, &q); std::printf("  %.0f-%.0f: %.3f", eb[k], eb[k + 1], med); }
        }
        std::printf("\n");
    }
    TH3D *hsep = (TH3D *)fs.Get("h3_dsep_E_d");
    if (hsep)
    {
        std::printf("\n   simulation reco - true two-photon distance (cm, mean), per pair energy (rows) and true distance (columns)\n");
        for (int e = 0; e < kNE; e++)
        {
            std::printf("  %4.0f-%-4.0f", kEEdges[e], kEEdges[e + 1]);
            for (int d = 0; d < kND; d++)
            {
                TH1D *p = hsep->ProjectionX(Form("sp_%d_%d", e, d), e + 1, e + 1, d + 1, d + 1);
                if (p->GetEffectiveEntries() > 20) std::printf("  %+6.2f", p->GetMean()); else std::printf("  %6s", "-");
            }
            std::printf("\n");
        }
    }
    return 0;
}
