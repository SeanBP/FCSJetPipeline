// The ECal gain calibration's pi0 selection (pi0_pairs' h3_liang_m_E_ring_*:
// StFcsPi0FinderForEcal's highest-energy cluster pair, clusters > 1 GeV,
// Z < 0.7, <= 40 ECal hits above 1 GeV), data against simulation, per pair
// energy bin, ring and trigger class, fitted with A. Ogawa's massfit.C model.
//
// Prints, for each trigger class (any EM, fcsEM2) and ring (all, inner,
// middle, outer): per pair-energy bin the data and simulation peaks with fit
// errors, their ratio, and the entries; then the data peaks of each further
// data file (other runs) in the same bins, with the ratio to the first data
// file. Peaks: A. Ogawa's model with the mean held to 90-190 MeV and the width
// to 6-45 MeV (fits at a limit or without a significant peak are dropped);
// the data / sim ratio is also given for a Gaussian + quadratic fit on
// 50-300 MeV as a fit-model check. Writes <base>.pdf with the fitted spectra
// (all rings, any EM trigger).
//
// Usage: ./pi0_calib_compare <sim.root> <base> <label>:<data.root> [<label>:<data.root> ...]
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TError.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

static const double kMPi0 = 0.13498;
// pi0_pairs puts masses above 0.4 GeV in the last bin (0.395-0.4): fits and
// plots stop below it
static const double kFitMax = 0.39;
static const int kNB = 8;
static const double kBEdges[kNB + 1] = {2, 4, 6, 8, 10, 14, 20, 30, 40};
static const char *kRingName[4] = {"all rings", "inner (r < 66 cm)", "middle (66-88 cm)", "outer (r > 88 cm)"};
static const char *kClass[2] = {"em", "em2"};  // h3_liang_m_E_ring_<class>
static const char *kClassName[2] = {"any EM trigger", "fcsEM2"};

struct Peak { double m, em, s, n; bool ok; };

// A. Ogawa's massfit.C model: Gaussian on A/tau^2 x^p exp(-x/tau) (x0 = 0),
// background fitted alone first, on [0, kFitMax]
static double AkioBg(double *x, double *p) { return x[0] <= 0 ? 0 : p[0] / p[1] / p[1] * std::pow(x[0], p[2]) * std::exp(-x[0] / p[1]); }
static double AkioAll(double *x, double *p) { return AkioBg(x, p) + p[3] * std::exp(-0.5 * std::pow((x[0] - p[4]) / p[5], 2)); }

static double PeakSeed(TH1D *h)
{
    // highest bin between 80 and 200 MeV
    double best = -1, at = 0.135;
    for (int i = h->FindBin(0.08); i <= h->FindBin(0.2); i++)
        if (h->GetBinContent(i) > best) { best = h->GetBinContent(i); at = h->GetBinCenter(i); }
    return at;
}

// Mean limited to 90-190 MeV and width to 6-45 MeV; a fit that hits a limit
// or has an insignificant peak is not used
static Peak Judge(TF1 *f, double lo_m, double hi_m)
{
    Peak p = {f->GetParameter(4), f->GetParError(4), f->GetParameter(5), 0, false};
    p.ok = p.em > 0 && p.em < 0.015 && f->GetParameter(3) > 3 * f->GetParError(3) &&
           p.m > lo_m + 0.002 && p.m < hi_m - 0.002 && p.s > 0.0065 && p.s < 0.0445;
    return p;
}

static Peak FitAkio(TH1D *h, TF1 **keep = 0)
{
    Peak p = {0, 0, 0, 0, false};
    if (!h) return p;
    const double n = h->GetEffectiveEntries();
    p.n = n;
    if (n < 150) return p;
    const double peak = PeakSeed(h), hp = h->GetBinContent(h->FindBin(peak));
    TF1 b("akbg", AkioBg, 0, kFitMax, 3);
    b.SetParameters(h->GetMaximum() / 2, 0.2, 1.0);
    h->Fit(&b, "QRN0");
    TF1 *f = new TF1(Form("ak_%s", h->GetName()), AkioAll, 0, kFitMax, 6);
    f->SetParameters(b.GetParameter(0), b.GetParameter(1), b.GetParameter(2), 0.5 * hp, peak, 0.015);
    f->SetParLimits(3, 0, 3 * hp); f->SetParLimits(4, 0.09, 0.19); f->SetParLimits(5, 0.006, 0.045);
    if (h->Fit(f, "QRN0") != 0) { delete f; return p; }
    p = Judge(f, 0.09, 0.19);
    p.n = n;
    if (keep) *keep = f; else delete f;
    return p;
}

// Cross-check model: Gaussian + quadratic on [0.05, 0.30]
static double GPol(double *x, double *p) { return p[0] + p[1] * x[0] + p[2] * x[0] * x[0] + p[3] * std::exp(-0.5 * std::pow((x[0] - p[4]) / p[5], 2)); }
static Peak FitGPol(TH1D *h)
{
    Peak p = {0, 0, 0, 0, false};
    if (!h || h->GetEffectiveEntries() < 150) return p;
    const double peak = PeakSeed(h), hp = h->GetBinContent(h->FindBin(peak));
    TF1 f("gp", GPol, 0.05, 0.30, 6);
    f.SetParameters(0.3 * hp, 0, 0, 0.5 * hp, peak, 0.015);
    f.SetParLimits(3, 0, 3 * hp); f.SetParLimits(4, 0.09, 0.19); f.SetParLimits(5, 0.006, 0.045);
    if (h->Fit(&f, "QRN0") != 0) return p;
    p = Judge(&f, 0.09, 0.19);
    p.n = h->GetEffectiveEntries();
    return p;
}

static TH1D *Slice(TFile &f, int cls, int ib, int ring)
{
    TH3D *h = (TH3D *)f.Get(Form("h3_liang_m_E_ring_%s", kClass[cls]));
    if (!h) return 0;
    const int y0 = h->GetYaxis()->FindBin(kBEdges[ib] + 1e-6), y1 = h->GetYaxis()->FindBin(kBEdges[ib + 1] - 1e-6);
    const int z0 = ring == 0 ? 1 : ring, z1 = ring == 0 ? 3 : ring;
    return h->ProjectionX(Form("s_%s_%d_%d_%d_%p", kClass[cls], ib, ring, cls, (void *)&f), y0, y1, z0, z1, "e");
}

static void PrintPeak(const Peak &p)
{
    if (p.ok) std::printf("  %6.1f +- %4.1f [%4.1f]", 1000 * p.m, 1000 * p.em, 1000 * p.s);
    else std::printf("  %22s", "-");
}

static void PrintRatio(const Peak &a, const Peak &b)
{
    if (a.ok && b.ok) std::printf("  %6.3f +- %5.3f", a.m / b.m, (a.m / b.m) * std::hypot(a.em / a.m, b.em / b.m));
    else std::printf("  %15s", "-");
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);
    if (argc < 4) { std::printf("usage: ./pi0_calib_compare <sim.root> <base> <label>:<data.root> ...\n"); return 1; }
    TFile fs(argv[1]);
    const std::string base = argv[2];
    std::vector<std::string> label;
    std::vector<TFile *> fd;
    for (int i = 3; i < argc; i++)
    {
        const std::string a = argv[i];
        const size_t c = a.find(':');
        label.push_back(a.substr(0, c));
        fd.push_back(new TFile(a.substr(c + 1).c_str()));
    }

    for (int cls = 0; cls < 2; cls++)
        for (int ring = 0; ring < 4; ring++)
        {
            std::printf("\n== %s, %s: peak (MeV) +- fit error [sigma]; data / sim; data entries | sim effective entries\n", kClassName[cls], kRingName[ring]);
            std::printf("  %-9s  %22s  %22s  %15s  %15s  %16s\n", "pair E", ("data " + label[0]).c_str(), "sim", "data / sim", "(gaus+pol2)", "entries");
            for (int ib = 0; ib < kNB; ib++)
            {
                TH1D *hd = Slice(*fd[0], cls, ib, ring), *hs = Slice(fs, cls, ib, ring);
                Peak pd = FitAkio(hd), ps = FitAkio(hs);
                std::printf("  %4.0f-%-4.0f", kBEdges[ib], kBEdges[ib + 1]);
                PrintPeak(pd); PrintPeak(ps); PrintRatio(pd, ps);
                PrintRatio(FitGPol(hd), FitGPol(hs));
                std::printf("  %7.0f | %7.0f\n", pd.n, ps.n);
            }
        }

    if (fd.size() > 1)
        for (int cls = 0; cls < 2; cls++)
        {
            std::printf("\n== Data runs, %s, all rings: peak (MeV) +- fit error, and ratio to %s\n", kClassName[cls], label[0].c_str());
            std::printf("  %-9s", "pair E");
            for (size_t k = 0; k < fd.size(); k++) std::printf("  %22s", label[k].c_str());
            std::printf("\n");
            for (int ib = 0; ib < kNB; ib++)
            {
                std::printf("  %4.0f-%-4.0f", kBEdges[ib], kBEdges[ib + 1]);
                std::vector<Peak> pk;
                for (size_t k = 0; k < fd.size(); k++) pk.push_back(FitAkio(Slice(*fd[k], cls, ib, 0)));
                for (size_t k = 0; k < fd.size(); k++)
                {
                    if (!pk[k].ok) { std::printf("  %22s", "-"); continue; }
                    if (k == 0 || !pk[0].ok) std::printf("  %8.1f +- %4.1f      ", 1000 * pk[k].m, 1000 * pk[k].em);
                    else std::printf("  %6.1f +- %3.1f (%5.3f)", 1000 * pk[k].m, 1000 * pk[k].em, pk[k].m / pk[0].m);
                }
                std::printf("\n");
            }
        }

    // Spectra with fits: data (first file) and simulation, any EM trigger, all rings
    TCanvas cv("cv", "", 1600, 800);
    cv.Divide(4, 2);
    for (int ib = 0; ib < kNB; ib++)
    {
        cv.cd(ib + 1);
        TH1D *hd = Slice(*fd[0], 0, ib, 0), *hs = Slice(fs, 0, ib, 0);
        if (!hd || !hs || hd->Integral() <= 0 || hs->Integral() <= 0) continue;
        TF1 *f1 = 0, *f2 = 0;
        FitAkio(hd, &f1); FitAkio(hs, &f2);
        const double nd = hd->Integral(), ns = hs->Integral();
        hd->Scale(1 / nd); hs->Scale(1 / ns);
        hs->SetLineColor(kRed + 1); hd->SetMarkerStyle(20); hd->SetMarkerSize(0.6);
        hd->GetXaxis()->SetRangeUser(0, kFitMax); hs->GetXaxis()->SetRangeUser(0, kFitMax);
        hs->SetMaximum(1.3 * std::max(hd->GetMaximum(), hs->GetMaximum())); hs->SetMinimum(0);
        hs->SetTitle(Form("pair E %.0f-%.0f GeV;m_{#gamma#gamma} (GeV);fraction", kBEdges[ib], kBEdges[ib + 1]));
        hs->Draw("hist"); hd->Draw("e same");
        if (f1) { f1->SetParameter(0, f1->GetParameter(0) / nd); f1->SetParameter(3, f1->GetParameter(3) / nd); f1->SetLineColor(kBlack); f1->Draw("same"); }
        if (f2) { f2->SetParameter(0, f2->GetParameter(0) / ns); f2->SetParameter(3, f2->GetParameter(3) / ns); f2->SetLineColor(kRed + 1); f2->SetLineStyle(2); f2->Draw("same"); }
        TLatex t; t.SetNDC(); t.SetTextSize(0.05);
        t.DrawLatex(0.5, 0.82, Form("data (black), sim (red)"));
    }
    cv.SaveAs((base + ".pdf").c_str());
    return 0;
}
