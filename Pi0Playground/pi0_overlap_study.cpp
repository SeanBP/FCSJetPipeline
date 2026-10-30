// Is the data/simulation pi0 peak ratio a real energy non-linearity (it would
// follow single-photon energy) or a two-cluster overlap artifact (it would
// follow cluster separation)? Uses pi0_pairs' overlap histograms (ECal gain
// calibration's selection, 45 MeV tower floor), data against simulation:
//   A. ratio against cluster separation in towers, 4-20 GeV pairs, and per
//      pair-energy bin
//   B. ratio against single-photon energy, all Z and per Z range
//   C. eta and pi0 peaks against the lower photon energy (eta photons are ~4x
//      further apart: little overlap)
//   D. cluster width (energy-weighted RMS radius) against cluster energy
//   E. MIP peak of isolated single ECal towers
// Peaks: Gaussian + quadratic (pi0 on [0.05, 0.30], eta on [0.35, 0.78]).
// Each table gives the "all" class (every event, both sides) and "em" (EM bit).
//
// Usage: ./pi0_overlap_study <data_all.root> <data_em.root> <sim.root>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TROOT.h>
#include <TError.h>
#include <cmath>
#include <cstdio>
#include <string>

struct Peak { double m, em, s; bool ok; };

static Peak FitG(TH1D *h, double lo, double hi, double mlo, double mhi, double slo, double shi, double seedlo, double seedhi)
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
    return p;
}
static Peak FitPi0(TH1D *h) { return FitG(h, 0.05, 0.30, 0.09, 0.19, 0.006, 0.045, 0.08, 0.20); }
static Peak FitEta(TH1D *h) { return FitG(h, 0.35, 0.78, 0.45, 0.65, 0.015, 0.09, 0.45, 0.65); }

static void Ratio(const Peak &a, const Peak &b)
{
    if (a.ok && b.ok) std::printf("  %6.3f +- %5.3f", a.m / b.m, (a.m / b.m) * std::hypot(a.em / a.m, b.em / b.m));
    else std::printf("  %15s", "-");
}

// 1D mass projection of a TH3 over y in [y0, y1) and z in [z0, z1)
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

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError;
    if (argc < 4) { std::printf("usage: ./pi0_overlap_study <data_all.root> <data_em.root> <sim.root>\n"); return 1; }
    TFile fa(argv[1]), fe(argv[2]), fs(argv[3]);
    TFile *fd[2] = {&fa, &fe};
    const char *cls[2] = {"all", "em"};

    std::printf("A. data/sim pi0 peak against cluster separation (towers), pairs 4-20 GeV; columns: all events | EM bit\n");
    const int nS = 8;
    const double sE[nS + 1] = {0, 2, 2.5, 3, 3.5, 4, 5, 6.5, 10};
    for (int i = 0; i < nS; i++)
    {
        std::printf("  sep %4.1f-%-4.1f", sE[i], sE[i + 1]);
        for (int c = 0; c < 2; c++)
            Ratio(FitPi0(P3(*fd[c], Form("h3_lsep_m_E_sep_%s", cls[c]), 4, 20, sE[i], sE[i + 1])), FitPi0(P3(fs, Form("h3_lsep_m_E_sep_%s", cls[c]), 4, 20, sE[i], sE[i + 1])));
        std::printf("\n");
    }
    std::printf("\n   per pair-energy bin (rows) and separation (columns), all events:\n  %-9s", "pair E");
    const int nE = 4;
    const double eE[nE + 1] = {4, 8, 10, 14, 20};
    for (int i = 0; i < nS; i++) std::printf("  %15s", Form("%.1f-%.1f", sE[i], sE[i + 1]));
    std::printf("\n");
    for (int e = 0; e < nE; e++)
    {
        std::printf("  %4.0f-%-4.0f", eE[e], eE[e + 1]);
        for (int i = 0; i < nS; i++)
            Ratio(FitPi0(P3(fa, "h3_lsep_m_E_sep_all", eE[e], eE[e + 1], sE[i], sE[i + 1])), FitPi0(P3(fs, "h3_lsep_m_E_sep_all", eE[e], eE[e + 1], sE[i], sE[i + 1])));
        std::printf("\n");
    }

    std::printf("\nB. data/sim pi0 peak against single-photon energy (each photon of the pair), all Z; columns: all events | EM bit\n");
    const int nG = 8;
    const double gE[nG + 1] = {1, 2, 3, 4, 5, 6, 8, 12, 20};
    for (int i = 0; i < nG; i++)
    {
        std::printf("  E_g %4.0f-%-4.0f", gE[i], gE[i + 1]);
        for (int c = 0; c < 2; c++)
            Ratio(FitPi0(P3(*fd[c], Form("h3_lph_m_Eg_Z_%s", cls[c]), gE[i], gE[i + 1], 0, 0.7)), FitPi0(P3(fs, Form("h3_lph_m_Eg_Z_%s", cls[c]), gE[i], gE[i + 1], 0, 0.7)));
        std::printf("\n");
    }
    std::printf("\n   per photon energy (rows) and Z (columns), all events:\n  %-11s", "E_g");
    const int nZ = 3;
    const double zE[nZ + 1] = {0, 0.3, 0.5, 0.7};
    for (int z = 0; z < nZ; z++) std::printf("  %15s", Form("Z %.1f-%.1f", zE[z], zE[z + 1]));
    std::printf("\n");
    for (int i = 0; i < nG; i++)
    {
        std::printf("  %4.0f-%-5.0f", gE[i], gE[i + 1]);
        for (int z = 0; z < nZ; z++)
            Ratio(FitPi0(P3(fa, "h3_lph_m_Eg_Z_all", gE[i], gE[i + 1], zE[z], zE[z + 1])), FitPi0(P3(fs, "h3_lph_m_Eg_Z_all", gE[i], gE[i + 1], zE[z], zE[z + 1])));
        std::printf("\n");
    }

    std::printf("\nC. eta and pi0 against the lower photon energy of the pair: peak (MeV) data | sim, and data/sim; all events, then EM bit\n");
    const int nL = 5;
    const double lE[nL + 1] = {1, 2, 3, 5, 8, 20};
    for (int c = 0; c < 2; c++)
    {
        std::printf("  [%s]\n", cls[c]);
        for (int i = 0; i < nL; i++)
        {
            TH1D *dh = P2(*fd[c], Form("h2_leta_m_Eg_%s", cls[c]), lE[i], lE[i + 1]), *sh = P2(fs, Form("h2_leta_m_Eg_%s", cls[c]), lE[i], lE[i + 1]);
            Peak dp = FitPi0(dh), sp = FitPi0(sh), de = FitEta(dh), se = FitEta(sh);
            std::printf("  E_low %4.0f-%-4.0f pi0:", lE[i], lE[i + 1]);
            if (dp.ok && sp.ok) std::printf(" %6.1f | %6.1f", 1000 * dp.m, 1000 * sp.m); else std::printf(" %15s", "-");
            Ratio(dp, sp);
            std::printf("   eta:");
            if (de.ok && se.ok) std::printf(" %6.1f | %6.1f", 1000 * de.m, 1000 * se.m); else std::printf(" %15s", "-");
            Ratio(de, se);
            std::printf("\n");
        }
    }

    std::printf("\nD. cluster width: mean energy-weighted RMS radius (cm) of the pair clusters, data (all | EM) and simulation, against cluster energy\n");
    TH2D *wd[2] = {(TH2D *)fa.Get("h2_clwidth_E"), (TH2D *)fe.Get("h2_clwidth_E")}, *ws = (TH2D *)fs.Get("h2_clwidth_E");
    for (int i = 0; i < nG; i++)
    {
        std::printf("  E %4.0f-%-4.0f", gE[i], gE[i + 1]);
        double m[3];
        TH2D *hh[3] = {wd[0], wd[1], ws};
        for (int k = 0; k < 3; k++)
        {
            TH1D *p = hh[k]->ProjectionY(Form("w%d_%d", k, i), hh[k]->GetXaxis()->FindBin(gE[i] + 1e-6), hh[k]->GetXaxis()->FindBin(gE[i + 1] - 1e-6));
            m[k] = p->GetMean();
            std::printf("  %5.2f", m[k]);
        }
        std::printf("   data/sim %.3f | %.3f\n", m[0] / m[2], m[1] / m[2]);
    }

    std::printf("\nE. MIP peak of isolated single ECal towers (Gaussian on 0.1-0.6 GeV): data all | data EM | sim\n");
    {
        TH1D *h[3] = {(TH1D *)fa.Get("h_isohit_E"), (TH1D *)fe.Get("h_isohit_E"), (TH1D *)fs.Get("h_isohit_E")};
        Peak p[3];
        for (int k = 0; k < 3; k++) p[k] = FitG(h[k], 0.10, 0.60, 0.12, 0.50, 0.02, 0.20, 0.12, 0.45);
        for (int k = 0; k < 3; k++)
            if (p[k].ok) std::printf("  %6.1f +- %4.1f MeV [%4.1f]", 1000 * p[k].m, 1000 * p[k].em, 1000 * p[k].s); else std::printf("  %22s", "-");
        std::printf("\n  data/sim:"); Ratio(p[0], p[2]); Ratio(p[1], p[2]); std::printf("\n");
    }
    return 0;
}
