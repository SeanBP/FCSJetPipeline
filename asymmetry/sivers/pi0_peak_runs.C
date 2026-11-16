// Run-by-run pi0 peak of the FMS-style pi0 candidates (Sean 2026-11-16): a
// check of the single pi0 energy scale factor (AnConfig.h kPi0EnergyScale) and
// a bad-run finder, in the spirit of D. Kapukchyan's pi0 QA list.
// Per yield file (one run; AnHistograms' "runs" entry):
//   - sums h_xFpi_* over EM levels + JP-not-EM (_jponly), sides, spin states,
//     every z_em class and x_F;
//   - fits M_gg (already on the pi0 scale) with a Gaussian + quadratic in
//     0.03-0.3 GeV.
// Flags (written in the table, the "flag" column; "-" = good), computed on the
// RAW peak and width (the fit / the factor the yields used), so they do not
// depend on the scale table:
//   LowStat   fewer than kMinCand candidates in 0-0.3 GeV, or the fit failed
//             (mean outside 0.08-0.2, sigma outside 0.008-0.08)
//   PeakShift |mean - median| > 3 sqrt(err^2 + spread^2), spread = the robust
//             run-to-run scatter (1.4826 x MAD of the good-statistics runs)
//   Wide      sigma > median sigma + 3 x its robust scatter
// Bad-run list (Sean 2026-11-16): <out>_bad_runs.txt, in D. Kapukchyan's
// format, holds the PeakShift / Wide runs. AnAsymmetry -Q masks them on top of
// his list; LowStat runs appear as comments, not masked.
// Writes <out>.txt (run, day, candidates, mean, err, sigma, err, mean/0.135,
// the scale the yields used, raw peak, flag) and <out>.png (peak and width vs
// day; the median, the +-3% band and the gain changes 23027047 / 23066056
// marked; flagged runs in red).
// Scale tables (Sean 2026-11-16: choose global, per period or per run after
// the pass, then rehist with submit_rehist.sh -s <table>). The raw peak is the
// fitted peak / the factor the yields used (TNamed "pi0_scale" of each file;
// AnConfig.h kPi0EnergyScale if missing). Each factor is 0.135 / raw peak:
//   <out>_scale_global.txt      default = 0.135 / median raw peak
//   <out>_scale_per_period.txt  one factor per period between the gain
//                               changes (median raw peak of its good runs)
//   <out>_scale_per_run.txt     each good, unflagged run its own factor; the
//                               others the global one
// Input: a list of per-run yield files, or an AnMerge output (.root; its
// pi0_runs spectra, the merged report flow).
// Usage (compiled, ACLiC): root -b -q 'pi0_peak_runs.C+("yields.list", "pi0_peak_runs")'
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TNamed.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

struct RunPeak { int run; double day, n, mean, emean, sig, esig, scale; bool ok; string flag; };

// the factor a yield file used for this run ("<source>; <run>=<factor> ..."), 0 if not recorded
static double ScaleOf(TFile &f, int run)
{
    TNamed *t = (TNamed *)f.Get("pi0_scale");
    if (!t) return 0;
    const string s = t->GetTitle();
    const string key = Form(" %d=", run);
    const size_t p = s.find(key);
    return p == string::npos ? 0 : atof(s.c_str() + p + key.size());
}

// fits one run's M_gg spectrum (Gaussian + quadratic, 0.03-0.3 GeV)
static void FitRun(TH1D *m, RunPeak &r, double min_cand, double mpi)
{
    if (r.n >= min_cand)
    {
        TF1 fg("fg", "gaus(0)+pol2(3)", 0.03, 0.3);
        fg.SetParameters(m->GetMaximum(), mpi, 0.025, 0, 0, 0);
        fg.SetParLimits(1, 0.06, 0.25);
        fg.SetParLimits(2, 0.005, 0.1);
        m->Fit(&fg, "RQN0");
        r.mean = fg.GetParameter(1); r.emean = fg.GetParError(1);
        r.sig = fabs(fg.GetParameter(2)); r.esig = fg.GetParError(2);
        r.ok = r.mean > 0.08 && r.mean < 0.2 && r.sig > 0.008 && r.sig < 0.08 && r.emean > 0;
    }
    if (!r.ok) r.flag = "LowStat";
}

// period between the gain changes 23027047 and 23066056
static int PeriodOf(int run) { return run < 23027047 ? 0 : run < 23066056 ? 1 : 2; }

static double Median(vector<double> v)
{
    if (v.empty()) return 0;
    sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// day of Run 22 (2021 days negative), from the run number
static double DayOf(int run) { return (run / 1000) % 1000 + (run / 1000000 == 22 ? -365 : 0) + (run % 1000) / 1000.0; }

void pi0_peak_runs(const char *yield_list, const char *out)
{
    gStyle->SetOptStat(0);
    const double kMinCand = 300, kMpi = 0.135;
    const char *lv[6] = {"EM0", "EM1", "EM2", "JP0", "JP1", "JP2"};
    const char *sd[2] = {"xneg", "xpos"};
    const char *sp[4] = {"BuYu", "BuYd", "BdYu", "BdYd"};
    vector<RunPeak> rp;
    int nmiss = 0, nnoscale = 0;
    // an AnMerge output (merged report flow): its per-run spectra (pi0_runs/m_<run>,
    // title = the scale factor used); otherwise a list of per-run yield files
    if (TString(yield_list).EndsWith(".root"))
    {
        TFile mf(yield_list);
        TDirectory *d = mf.IsZombie() ? 0 : (TDirectory *)mf.Get("pi0_runs");
        if (!d) { printf("pi0_peak_runs: no pi0_runs in %s\n", yield_list); return; }
        TIter nk(d->GetListOfKeys());
        while (TObject *k = nk())
        {
            TH1D *m = (TH1D *)d->Get(k->GetName());
            const int run = atoi(k->GetName() + 2);
            double sc = atof(m->GetTitle());
            if (!(sc > 0)) { sc = 0.135 / 0.1248; nnoscale++; }
            RunPeak r = {run, DayOf(run), m->Integral(), 0, 0, 0, 0, sc, false, "-"};
            FitRun(m, r, kMinCand, kMpi);
            rp.push_back(r);
        }
    }
    ifstream in(TString(yield_list).EndsWith(".root") ? "" : yield_list);
    string fn;
    while (in >> fn)
    {
        TFile f(fn.c_str());
        if (f.IsZombie()) continue;
        TNamed *rn = (TNamed *)f.Get("runs");
        int run = 0, extra = 0;
        if (rn) { istringstream is(rn->GetTitle()); if (!(is >> run) || (is >> extra)) run = 0; }
        if (run == 0) continue;
        TH1D *m = 0;
        for (int l = 0; l < 6; l++)
            for (int s = 0; s < 2; s++)
                for (int p = 0; p < 4; p++)
                {
                    TH3D *h = (TH3D *)f.Get(Form("h_xFpi_%s_%s_%s%s", lv[l], sd[s], sp[p], l >= 3 ? "_jponly" : ""));
                    if (!h) continue;
                    TH1D *py = h->ProjectionY("py_tmp");
                    if (!m) { m = (TH1D *)py->Clone("m_run"); m->SetDirectory(0); m->Reset(); }
                    m->Add(py);
                    delete py;
                }
        if (!m) { nmiss++; continue; }
        double sc = ScaleOf(f, run);
        if (!(sc > 0)) { sc = 0.135 / 0.1248; nnoscale++; }  // AnConfig.h kPi0EnergyScale of 2026-11-16
        RunPeak r = {run, DayOf(run), m->Integral(), 0, 0, 0, 0, sc, false, "-"};
        FitRun(m, r, kMinCand, kMpi);
        rp.push_back(r);
        delete m;
    }
    if (rp.empty()) { printf("pi0_peak_runs: no yield files with h_xFpi_* (%d without)\n", nmiss); return; }
    sort(rp.begin(), rp.end(), [](const RunPeak &a, const RunPeak &b) { return a.run < b.run; });
    // robust centre and scatter of the good runs: on the yields' (scaled) peaks
    // for the plot, and on the RAW peaks (fit / the factor used) for the flags,
    // so a per-run rescaling cannot hide a shifted run
    vector<double> vm, vs, vr, vrs;
    for (size_t i = 0; i < rp.size(); i++)
        if (rp[i].ok) { vm.push_back(rp[i].mean); vs.push_back(rp[i].sig); vr.push_back(rp[i].mean / rp[i].scale); vrs.push_back(rp[i].sig / rp[i].scale); }
    const double med = Median(vm), meds = Median(vs), medr = Median(vr), medrs = Median(vrs);
    vector<double> dm, ds, dr, drs;
    for (size_t i = 0; i < vm.size(); i++)
    {
        dm.push_back(fabs(vm[i] - med)); ds.push_back(fabs(vs[i] - meds));
        dr.push_back(fabs(vr[i] - medr)); drs.push_back(fabs(vrs[i] - medrs));
    }
    const double spread = 1.4826 * Median(dm), spreads = 1.4826 * Median(ds);
    const double spreadr = 1.4826 * Median(dr), spreadrs = 1.4826 * Median(drs);
    int nshift = 0, nwide = 0, nlow = 0;
    for (size_t i = 0; i < rp.size(); i++)
    {
        RunPeak &r = rp[i];
        if (!r.ok) { nlow++; continue; }
        const double rm = r.mean / r.scale, re = r.emean / r.scale, rs = r.sig / r.scale;
        string fl;
        if (vr.size() >= 5 && fabs(rm - medr) > 3 * sqrt(re * re + spreadr * spreadr)) { fl = "PeakShift"; nshift++; }
        if (vrs.size() >= 5 && rs > medrs + 3 * spreadrs) { fl += fl.empty() ? "Wide" : ",Wide"; nwide++; }
        if (!fl.empty()) r.flag = fl;
    }
    FILE *t = fopen(Form("%s.txt", out), "w");
    fprintf(t, "# run-by-run pi0 peak of the FMS-style candidates (M_gg on the pi0 scale, AnConfig.h kPi0EnergyScale); Gaussian + quadratic, 0.03-0.3 GeV\n");
    fprintf(t, "# median mean %.4f GeV (%.4f x 0.135), robust run-to-run scatter %.4f GeV (%.1f%%); median sigma %.4f, scatter %.4f\n", med, med / kMpi, spread, med > 0 ? 100 * spread / med : 0, meds, spreads);
    fprintf(t, "# flags: LowStat < %.0f candidates or failed fit; PeakShift |mean - median| > 3 sqrt(err^2 + scatter^2); Wide sigma > median + 3 scatter\n", kMinCand);
    fprintf(t, "# run day candidates mean err sigma err mean/0.135 scale_used raw_peak flag\n");
    for (size_t i = 0; i < rp.size(); i++)
        fprintf(t, "%d %.3f %.0f %.5f %.5f %.5f %.5f %.4f %.5f %.5f %s\n", rp[i].run, rp[i].day, rp[i].n, rp[i].mean, rp[i].emean, rp[i].sig, rp[i].esig, rp[i].ok ? rp[i].mean / kMpi : 0,
                rp[i].scale, rp[i].ok ? rp[i].mean / rp[i].scale : 0, rp[i].flag.c_str());
    fclose(t);
    // bad-run list in D. Kapukchyan's format, for AnAsymmetry -Q: the flagged
    // runs (masked); low-statistics runs only as comments (unchecked, not masked)
    {
        FILE *b = fopen(Form("%s_bad_runs.txt", out), "w");
        fprintf(b, "# pi0-peak bad-run list from pi0_peak_runs.C (raw peak median %.5f GeV, run-to-run scatter %.5f; width median %.5f, scatter %.5f)\n", medr, spreadr, medrs, spreadrs);
        fprintf(b, "# masked by AnAsymmetry -Q: every line with a flag. LowStat runs (< %.0f candidates) are listed as comments only.\n", kMinCand);
        fprintf(b, "Runnumber | Flag | AdditionalNotes\n");
        for (size_t i = 0; i < rp.size(); i++)
        {
            const RunPeak &r = rp[i];
            if (r.ok && r.flag != "-")
                fprintf(b, "%d | %s | raw peak %.4f +- %.4f GeV (median %.4f), raw sigma %.4f (median %.4f), %.0f candidates\n", r.run, r.flag.c_str(),
                        r.mean / r.scale, r.emean / r.scale, medr, r.sig / r.scale, medrs, r.n);
            else if (!r.ok) fprintf(b, "# %d LowStat: %.0f candidates (not checked)\n", r.run, r.n);
        }
        fclose(b);
    }
    if (nnoscale) printf("pi0_peak_runs: WARNING: %d file(s) without a pi0_scale record; assumed 0.135/0.1248\n", nnoscale);
    // scale tables from the raw peaks of the good, unflagged runs
    {
        vector<double> raw, rawp[3];
        for (size_t i = 0; i < rp.size(); i++)
            if (rp[i].ok && rp[i].flag == "-") { const double rr = rp[i].mean / rp[i].scale; raw.push_back(rr); rawp[PeriodOf(rp[i].run)].push_back(rr); }
        if (!raw.empty())
        {
            const double glob = kMpi / Median(raw);
            const char *hdr = "# pi0 energy scale table (asymmetry/pi0_scale.txt format) from pi0_peak_runs.C";
            FILE *g1 = fopen(Form("%s_scale_global.txt", out), "w");
            fprintf(g1, "%s: global factor = 0.135 / median raw peak of %zu good runs\ndefault %.6f\n", hdr, raw.size(), glob);
            fclose(g1);
            FILE *g2 = fopen(Form("%s_scale_per_period.txt", out), "w");
            fprintf(g2, "%s: one factor per period between the gain changes 23027047 / 23066056\ndefault %.6f\n", hdr, glob);
            double pf[3];
            for (int k = 0; k < 3; k++)
            {
                pf[k] = rawp[k].empty() ? glob : kMpi / Median(rawp[k]);
                fprintf(g2, "# period %d: %zu good runs, factor %.6f%s\n", k, rawp[k].size(), pf[k], rawp[k].empty() ? " (no good runs: global)" : "");
            }
            for (size_t i = 0; i < rp.size(); i++) fprintf(g2, "%d %.6f\n", rp[i].run, pf[PeriodOf(rp[i].run)]);
            fclose(g2);
            FILE *g3 = fopen(Form("%s_scale_per_run.txt", out), "w");
            fprintf(g3, "%s: each good, unflagged run 0.135 / its raw peak; flagged and low-statistics runs the global factor\ndefault %.6f\n", hdr, glob);
            for (size_t i = 0; i < rp.size(); i++)
                if (rp[i].ok && rp[i].flag == "-") fprintf(g3, "%d %.6f\n", rp[i].run, kMpi * rp[i].scale / rp[i].mean);
            fclose(g3);
            printf("pi0_peak_runs: scale tables: global %.5f; per period %.5f / %.5f / %.5f (runs %zu / %zu / %zu); %s_scale_{global,per_period,per_run}.txt\n",
                   glob, pf[0], pf[1], pf[2], rawp[0].size(), rawp[1].size(), rawp[2].size(), out);
        }
    }
    printf("pi0_peak_runs: %zu runs, %zu with a fit; median peak %.4f GeV (%.4f x 0.135), run-to-run scatter %.1f%%; flagged: %d PeakShift, %d Wide, %d LowStat\n",
           rp.size(), vm.size(), med, med / kMpi, med > 0 ? 100 * spread / med : 0, nshift, nwide, nlow);
    // plot: peak / 0.135 and sigma vs day; good runs blue, flagged red
    TGraphErrors *g[2][2];
    for (int q = 0; q < 2; q++) for (int k = 0; k < 2; k++) g[q][k] = new TGraphErrors();
    double dlo = 1e9, dhi = -1e9;
    for (size_t i = 0; i < rp.size(); i++)
    {
        if (!rp[i].ok) continue;
        const int k = rp[i].flag == "-" ? 0 : 1;
        const double v[2] = {rp[i].mean / kMpi, rp[i].sig}, e[2] = {rp[i].emean / kMpi, rp[i].esig};
        for (int q = 0; q < 2; q++) { const int ip = g[q][k]->GetN(); g[q][k]->SetPoint(ip, rp[i].day, v[q]); g[q][k]->SetPointError(ip, 0, e[q]); }
        dlo = min(dlo, rp[i].day); dhi = max(dhi, rp[i].day);
    }
    if (dlo > dhi) { dlo = -20; dhi = 110; }
    const double xlo = dlo - 5, xhi = dhi + 5;
    TCanvas *c = new TCanvas("c", "", 1500, 1000);
    TPad *mp = new TPad("mp", "", 0, 0, 1, 0.93); mp->Draw(); mp->Divide(1, 2, 0.003, 0.003);
    const char *yt[2] = {"fitted peak / 0.135", "fitted #sigma (GeV)"};
    const double cen[2] = {med / kMpi, meds};
    for (int q = 0; q < 2; q++)
    {
        mp->cd(q + 1); gPad->SetLeftMargin(0.08); gPad->SetRightMargin(0.02); gPad->SetBottomMargin(0.12); gPad->SetGridy();
        double ylo = q == 0 ? 0.85 : 0.0, yhi = q == 0 ? 1.2 : 0.07;
        TH1D *fr = new TH1D(Form("fr%d", q), Form(";day of Run 22 (2021 days negative);%s", yt[q]), 10, xlo, xhi);
        fr->SetMinimum(ylo); fr->SetMaximum(yhi);
        fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleOffset(0.65); fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetLabelSize(0.045); fr->GetXaxis()->SetLabelSize(0.045);
        fr->Draw();
        if (q == 0 && med > 0)
        {
            TBox *b = new TBox(xlo, cen[0] * 0.97, xhi, cen[0] * 1.03); b->SetFillColor(kGray); b->SetFillStyle(1001); b->Draw();
            fr->Draw("axis same");
        }
        TLine *lm = new TLine(xlo, cen[q], xhi, cen[q]); lm->SetLineColor(kAzure + 2); lm->SetLineStyle(2); lm->SetLineWidth(2); lm->Draw();
        if (q == 0) { TLine *l1 = new TLine(xlo, 1, xhi, 1); l1->SetLineColor(kBlack); l1->SetLineWidth(1); l1->Draw(); }
        const int gc[2] = {23027047, 23066056};
        for (int k = 0; k < 2; k++)
        {
            const double d = DayOf(gc[k]);
            TLine *lg = new TLine(d, ylo, d, yhi); lg->SetLineColor(kGreen + 2); lg->SetLineStyle(7); lg->SetLineWidth(2); lg->Draw();
        }
        g[q][0]->SetMarkerStyle(20); g[q][0]->SetMarkerSize(0.8); g[q][0]->SetMarkerColor(kAzure + 2); g[q][0]->SetLineColor(kAzure + 2);
        g[q][1]->SetMarkerStyle(21); g[q][1]->SetMarkerSize(0.9); g[q][1]->SetMarkerColor(kRed + 1); g[q][1]->SetLineColor(kRed + 1);
        if (g[q][0]->GetN()) g[q][0]->Draw("P");
        if (g[q][1]->GetN()) g[q][1]->Draw("P");
        if (q == 0)
        {
            TLegend *lg = new TLegend(0.09, 0.80, 0.98, 0.89); lg->SetNColumns(4); lg->SetTextFont(42); lg->SetTextSize(0.042); lg->SetBorderSize(1); lg->SetFillColor(kWhite);
            lg->AddEntry(g[0][0], "run", "p");
            lg->AddEntry(g[0][1], "flagged run", "p");
            lg->AddEntry(lm, Form("median %.3f (#pm3%% band)", cen[0]), "l");
            TLine *dum = new TLine(); dum->SetLineColor(kGreen + 2); dum->SetLineStyle(7); dum->SetLineWidth(2);
            lg->AddEntry(dum, "gain changes", "l");
            lg->Draw();
        }
    }
    c->cd();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(42); tt.SetTextSize(0.021);
    tt.DrawLatex(0.01, 0.970, Form("#bf{#pi^{0} peak run by run}: FMS-style candidates (all z_{em}, EM + JP union), M_{#gamma#gamma} scaled by kPi0EnergyScale; Gaussian + quadratic fit, runs with #geq %.0f candidates", kMinCand));
    tt.DrawLatex(0.01, 0.943, Form("%zu runs fitted of %zu; run-to-run scatter %.1f%%; flagged: %d peak shift, %d wide (red), %d low statistics (not drawn); table: %s.txt", vm.size(), rp.size(), med > 0 ? 100 * spread / med : 0, nshift, nwide, nlow, out));
    c->SaveAs(Form("%s.png", out));
}
