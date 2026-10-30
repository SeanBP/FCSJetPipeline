// Trigger dependence of the pi0 peak in the ECal gain calibration's selection
// (pi0_pairs' h3_liang_m_E_ring_<class>, all rings), data against simulation.
//
//   1. Peak per pair-energy bin, Gaussian + quadratic on [0.05, 0.30] GeV:
//      data (all events, EM-triggered, no EM bit) and simulation (all,
//      emulated EM trigger)
//   2. Trigger turn-on: fraction of events with each EM bit against the
//      event's highest ECal cluster E_T, data (all events) and simulation;
//      E_T where each reaches 50% of its plateau
//   3. EM-triggered events, pair inside / outside the offline firing-patch
//      stand-in, data and simulation
//   4. Pooled triggered data against triggered simulation over 4-20 GeV: one
//      fit of the summed 4-20 GeV spectra, and the error-weighted mean of the
//      per-bin ratios, for four fit models (spread = model uncertainty)
//   5. Widths: sigma_data, sigma_sim per bin; extra relative smearing of the
//      pair mass sqrt(sd^2 - ss^2) / m_data; per cluster (x sqrt 2); per tower
//      (/ sqrt <sum f_i^2>, the data's tower energy sharing)
//
// An optional 4th argument "floor" uses the common-45-MeV-tower-floor
// histograms (h3_liangf_*).
//
// Usage: ./pi0_trigger_study <data_all.root> <data_em.root> <sim.root> [floor]
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TError.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

static const int kNB = 7;
static const double kBEdges[kNB + 1] = {4, 6, 8, 10, 14, 20, 30, 40};
static const double kMPi0 = 0.13498;

struct Peak { double m, em, s, es; bool ok; };

static double PeakSeed(TH1D *h)
{
    double best = -1, at = 0.135;
    for (int i = h->FindBin(0.08); i <= h->FindBin(0.2); i++)
        if (h->GetBinContent(i) > best) { best = h->GetBinContent(i); at = h->GetBinCenter(i); }
    return at;
}

// model 0: Gaussian + quadratic on [0.05, 0.30]; 1: Gaussian + cubic on
// [0.05, 0.30]; 2: Gaussian + quadratic on [0.07, 0.22]; 3: A. Ogawa's model
// (Gaussian on A/tau^2 x^p exp(-x/tau)) on [0, 0.39]
static double Akio(double *x, double *p) { return (x[0] <= 0 ? 0 : p[0] / p[1] / p[1] * std::pow(x[0], p[2]) * std::exp(-x[0] / p[1])) + p[3] * std::exp(-0.5 * std::pow((x[0] - p[4]) / p[5], 2)); }
static Peak Fit(TH1D *h, int model = 0)
{
    Peak p = {0, 0, 0, 0, false};
    if (!h || h->GetEffectiveEntries() < 150) return p;
    const double seed = PeakSeed(h), hp = h->GetBinContent(h->FindBin(seed));
    TF1 *f = 0;
    if (model == 3)
    {
        TF1 b("akb", "[0]/[1]/[1]*pow(x,[2])*exp(-x/[1])", 0.001, 0.39);
        b.SetParameters(h->GetMaximum() / 2, 0.2, 1.0);
        h->Fit(&b, "QRN0");
        f = new TF1("fa", Akio, 0, 0.39, 6);
        f->SetParameters(b.GetParameter(0), b.GetParameter(1), b.GetParameter(2), 0.5 * hp, seed, 0.015);
    }
    else
    {
        const double lo = model == 2 ? 0.07 : 0.05, hi = model == 2 ? 0.22 : 0.30;
        f = new TF1("fg", model == 1 ? "[0]+[1]*x+[2]*x*x+[6]*x*x*x+[3]*exp(-0.5*((x-[4])/[5])^2)" : "[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])^2)", lo, hi);
        f->SetParameters(0.3 * hp, 0, 0, 0.5 * hp, seed, 0.015);
        if (model == 1) f->SetParameter(6, 0);
    }
    f->SetParLimits(3, 0, 3 * hp); f->SetParLimits(4, 0.09, 0.19); f->SetParLimits(5, 0.006, 0.045);
    if (h->Fit(f, "QRN0") == 0)
    {
        p.m = f->GetParameter(4); p.em = f->GetParError(4); p.s = f->GetParameter(5); p.es = f->GetParError(5);
        p.ok = p.em > 0 && p.em < 0.015 && f->GetParameter(3) > 3 * f->GetParError(3) &&
               p.m > 0.092 && p.m < 0.188 && p.s > 0.0065 && p.s < 0.0445;
    }
    delete f;
    return p;
}

static const char *gPrefix = "h3_liang";
static TH1D *Slice(TFile &f, const char *cls, double e0, double e1, const char *tag)
{
    TH3D *h = (TH3D *)f.Get(Form("%s_m_E_ring_%s", gPrefix, cls));
    if (!h) return 0;
    return h->ProjectionX(Form("%s_%s_%.0f_%.0f_%p", tag, cls, e0, e1, (void *)&f), h->GetYaxis()->FindBin(e0 + 1e-6), h->GetYaxis()->FindBin(e1 - 1e-6), 1, 3, "e");
}

static void Cell(const Peak &p) { if (p.ok) std::printf("  %6.1f +- %4.1f", 1000 * p.m, 1000 * p.em); else std::printf("  %14s", "-"); }
static void Ratio(const Peak &a, const Peak &b)
{
    if (a.ok && b.ok) std::printf("  %6.3f +- %5.3f", a.m / b.m, (a.m / b.m) * std::hypot(a.em / a.m, b.em / b.m));
    else std::printf("  %15s", "-");
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError;
    if (argc < 4) { std::printf("usage: ./pi0_trigger_study <data_all.root> <data_em.root> <sim.root>\n"); return 1; }
    TFile fa(argv[1]), fe(argv[2]), fs(argv[3]);
    if (argc > 4 && std::string(argv[4]) == "floor") gPrefix = "h3_liangf";
    std::printf("Histograms: %s_m_E_ring_<class>\n\n", gPrefix);

    std::printf("1. Peak (MeV) +- fit error per pair energy, Gaussian + quadratic; ratios to the simulation of the same class\n");
    std::printf("  %-9s %16s %16s %16s %16s %16s %15s %15s %15s\n", "pair E", "data all", "data EM", "data no EM", "sim all", "sim EM",
                "data/sim all", "data/sim EM", "data EM/all");
    for (int b = 0; b < kNB; b++)
    {
        Peak da = Fit(Slice(fa, "all", kBEdges[b], kBEdges[b + 1], "da")), de = Fit(Slice(fa, "em", kBEdges[b], kBEdges[b + 1], "de"));
        Peak dn = Fit(Slice(fa, "noem", kBEdges[b], kBEdges[b + 1], "dn"));
        Peak sa = Fit(Slice(fs, "all", kBEdges[b], kBEdges[b + 1], "sa")), se = Fit(Slice(fs, "em", kBEdges[b], kBEdges[b + 1], "se"));
        std::printf("  %4.0f-%-4.0f", kBEdges[b], kBEdges[b + 1]);
        Cell(da); Cell(de); Cell(dn); Cell(sa); Cell(se); Ratio(da, sa); Ratio(de, se); Ratio(de, da);
        std::printf("\n");
    }

    std::printf("\n1b. Offline-trigger-matched: EM bit and highest ECal cluster E_T above 2 / 3 / 4 GeV, both sides; data/sim per class\n");
    std::printf("  %-9s %15s %15s %15s %15s\n", "pair E", "EM bit only", "E_T > 2", "E_T > 3", "E_T > 4");
    {
        const char *oc[4] = {"em", "emet2", "emet3", "emet4"};
        for (int b = 0; b < kNB; b++)
        {
            std::printf("  %4.0f-%-4.0f", kBEdges[b], kBEdges[b + 1]);
            for (int k = 0; k < 4; k++) Ratio(Fit(Slice(fe, oc[k], kBEdges[b], kBEdges[b + 1], "od")), Fit(Slice(fs, oc[k], kBEdges[b], kBEdges[b + 1], "os")));
            std::printf("\n");
        }
        std::printf("  4-20 GeV in one fit:");
        for (int k = 0; k < 4; k++) Ratio(Fit(Slice(fe, oc[k], 4, 20, "o1")), Fit(Slice(fs, oc[k], 4, 20, "o2")));
        std::printf("\n");
    }

    std::printf("\n2. Trigger turn-on: fraction of events with the EM bit against the highest ECal cluster E_T (GeV)\n");
    std::printf("  %-8s", "E_T");
    const char *lev[3] = {"em0", "em1", "em2"};
    for (int k = 0; k < 3; k++) std::printf("   data %-4s  sim %-4s", lev[k], lev[k]);
    std::printf("\n");
    TH1D *ta[2] = {(TH1D *)fa.Get("h_maxclet_all"), (TH1D *)fs.Get("h_maxclet_all")};
    std::vector<std::vector<double> > frac(6);
    for (int i = 1; i <= ta[0]->GetNbinsX(); i += 2)
    {
        std::printf("  %4.1f-%-3.1f", ta[0]->GetBinLowEdge(i), ta[0]->GetBinLowEdge(i + 2));
        for (int k = 0; k < 3; k++)
            for (int f = 0; f < 2; f++)
            {
                TH1D *h = (TH1D *)(f == 0 ? fa : fs).Get(Form("h_maxclet_%s", lev[k]));
                const double n = ta[f]->GetBinContent(i) + ta[f]->GetBinContent(i + 1), m = h->GetBinContent(i) + h->GetBinContent(i + 1);
                std::printf("  %8.3f", n > 0 ? m / n : -1);
            }
        std::printf("\n");
    }
    std::printf("  E_T at half of the 6-8 GeV plateau:");
    for (int k = 0; k < 3; k++)
        for (int f = 0; f < 2; f++)
        {
            TH1D *h = (TH1D *)(f == 0 ? fa : fs).Get(Form("h_maxclet_%s", lev[k]));
            TH1D *r = (TH1D *)h->Clone(Form("r%d%d", k, f)); r->Divide(ta[f]);
            double plat = 0, np = 0;
            for (int i = r->FindBin(6.0); i <= r->FindBin(7.9); i++) { plat += r->GetBinContent(i); np++; }
            plat /= np;
            double x50 = -1;
            for (int i = 1; i < r->GetNbinsX(); i++)
                if (r->GetBinContent(i) < 0.5 * plat && r->GetBinContent(i + 1) >= 0.5 * plat) { x50 = r->GetBinLowEdge(i + 1); }
            std::printf("  %s %s %.2f (plateau %.2f)", f ? "sim" : "data", lev[k], x50, plat);
        }
    std::printf("\n");

    std::printf("\n3. EM-triggered events, pair inside / outside the offline firing-patch stand-in: peak (MeV)\n");
    std::printf("  %-9s %16s %16s %16s %16s %15s %15s\n", "pair E", "data inside", "data outside", "sim inside", "sim outside", "data/sim in", "data/sim out");
    for (int b = 0; b < kNB; b++)
    {
        Peak di = Fit(Slice(fe, "inpatch", kBEdges[b], kBEdges[b + 1], "di")), dout = Fit(Slice(fe, "outpatch", kBEdges[b], kBEdges[b + 1], "do"));
        Peak si = Fit(Slice(fs, "inpatch", kBEdges[b], kBEdges[b + 1], "si")), so = Fit(Slice(fs, "outpatch", kBEdges[b], kBEdges[b + 1], "so"));
        TH1D *ci = Slice(fe, "inpatch", kBEdges[b], kBEdges[b + 1], "ci"), *co = Slice(fe, "outpatch", kBEdges[b], kBEdges[b + 1], "co");
        std::printf("  %4.0f-%-4.0f", kBEdges[b], kBEdges[b + 1]);
        Cell(di); Cell(dout); Cell(si); Cell(so); Ratio(di, si); Ratio(dout, so);
        std::printf("   (data inside share %.2f)\n", ci && co && ci->Integral() + co->Integral() > 0 ? ci->Integral() / (ci->Integral() + co->Integral()) : -1);
    }

    std::printf("\n4. Pooled triggered data (%s) / triggered simulation over 4-20 GeV\n", argv[2]);
    const char *mname[4] = {"gaus+pol2 [0.05,0.30]", "gaus+pol3 [0.05,0.30]", "gaus+pol2 [0.07,0.22]", "Ogawa model [0,0.39]"};
    std::vector<double> one, wmean;
    for (int m = 0; m < 4; m++)
    {
        Peak d = Fit(Slice(fe, "em", 4, 20, "pd"), m), s = Fit(Slice(fs, "em", 4, 20, "ps"), m);
        double sw = 0, swr = 0;
        for (int b = 0; b < 5; b++)
        {
            Peak db = Fit(Slice(fe, "em", kBEdges[b], kBEdges[b + 1], "pb"), m), sb = Fit(Slice(fs, "em", kBEdges[b], kBEdges[b + 1], "qb"), m);
            if (!db.ok || !sb.ok) continue;
            const double r = db.m / sb.m, er = r * std::hypot(db.em / db.m, sb.em / sb.m);
            sw += 1 / (er * er); swr += r / (er * er);
        }
        std::printf("  %-24s one fit 4-20 GeV:", mname[m]);
        if (d.ok && s.ok) { std::printf(" data %.1f, sim %.1f MeV, ratio %.3f +- %.3f", 1000 * d.m, 1000 * s.m, d.m / s.m, (d.m / s.m) * std::hypot(d.em / d.m, s.em / s.m)); one.push_back(d.m / s.m); }
        else std::printf(" %41s", "-");
        if (sw > 0) { std::printf(";  weighted mean of bins %.3f +- %.3f", swr / sw, 1 / std::sqrt(sw)); wmean.push_back(swr / sw); }
        std::printf("\n");
    }
    auto spread = [](const std::vector<double> &v, const char *n) {
        if (v.empty()) return;
        double lo = v[0], hi = v[0], s = 0;
        for (double x : v) { lo = std::min(lo, x); hi = std::max(hi, x); s += x; }
        std::printf("  %s: mean %.3f, model spread %.3f-%.3f (half range %.3f)\n", n, s / v.size(), lo, hi, (hi - lo) / 2);
    };
    spread(one, "one fit");
    spread(wmean, "bin mean");

    std::printf("\n5. Widths (EM-triggered, Gaussian + quadratic): sigma (MeV); extra relative pair-mass smearing sqrt(sd^2 - ss^2)/m_d;\n"
                "   per cluster (x sqrt 2); data <sum f_i^2>; per tower (/ sqrt <sum f^2>)\n");
    TProfile *pf = (TProfile *)fe.Get("p_liang_sumf2_E");
    for (int b = 0; b < kNB; b++)
    {
        Peak d = Fit(Slice(fe, "em", kBEdges[b], kBEdges[b + 1], "wd")), s = Fit(Slice(fs, "em", kBEdges[b], kBEdges[b + 1], "ws"));
        std::printf("  %4.0f-%-4.0f", kBEdges[b], kBEdges[b + 1]);
        if (!d.ok || !s.ok) { std::printf("  -\n"); continue; }
        const double ex2 = d.s * d.s - s.s * s.s;
        double f2 = 0, nf = 0;
        for (int i = pf->GetXaxis()->FindBin(kBEdges[b] + 1e-6); i <= pf->GetXaxis()->FindBin(kBEdges[b + 1] - 1e-6); i++)
        { f2 += pf->GetBinContent(i) * pf->GetBinEntries(i); nf += pf->GetBinEntries(i); }
        f2 = nf > 0 ? f2 / nf : -1;
        std::printf("  data %4.1f +- %3.1f  sim %4.1f +- %3.1f", 1000 * d.s, 1000 * d.es, 1000 * s.s, 1000 * s.es);
        if (ex2 > 0)
        {
            const double pair = std::sqrt(ex2) / d.m, clu = std::sqrt(2.0) * pair;
            std::printf("  pair %.3f  cluster %.3f  <sum f^2> %.2f  tower %.3f\n", pair, clu, f2, f2 > 0 ? clu / std::sqrt(f2) : -1);
        }
        else std::printf("  (data not wider)\n");
    }
    return 0;
}
