// Trigger levels over Run 22 relative to the HARDWARE settings (Sean 2026-11-17):
// Akio's nominal thresholds kTrgThresholdGeV (EM0-2 2 / 3.5 / 5 GeV, JP0-2
// 4 / 6 / 8 GeV) and, for EM, the HCal veto setting r = 4 H_max / E_4x4 = 1
// (the hardware vetoes at r > 1). One factor f scales the trigger value and the
// veto buffer together: an EM jet passes at f if E_4x4 >= f x E_hw and
// r <= 1 / f, i.e. its value is min(E_4x4 / E_hw, 1 / r); a JP jet passes if its
// patch sum >= f x E_hw. The event's value is its largest jet value. Per group
// of runs and trigger, the f where
//   efficiency = N(offline and hardware) / N(hardware)  = 50%  (falls with f)
//   purity     = N(offline and hardware) / N(offline)   = 50%  (rises with f)
//   x50        = the bit is set in 50% of the events at that f (direct crossing
//                of the bit fraction in bins of 0.05, no fit; spectrum-independent)
// The median x50 of each gain period (23027047, 23066056) is drawn as a line:
// the period's nominal factor. Crossings are interpolated; errors: the range
// where the 68% Clopper-Pearson band contains 0.5. Reference events only (EM: JP,
// DiJP, Had, HHT, HTOT, HcalTot; JP: EM0/1, DY, J/psi, ELE2); hardware bit:
// fcs_tcu_bits before prescale, JP0/JP1 any region (hardware events count with a
// jet >= 0.5 x E_hw). Hot JP run 23070019 is left out of JP.
// Groups: consecutive runs within a gain period, [runs per group] each; with a
// trg_safety.txt (trg_eff_scan.cpp) its per-trigger groups instead (wider for
// EM1/EM2, which have few events near threshold).
// Output: <outdir>/trg_time_scan.txt (per trigger and group; per period medians)
// and trg_time_scan.png.
// Usage: ./trg_time_scan <jet file list> <outdir> [runs per group = 5] [trg_safety.txt]
// Build: g++ -m64 -O2 -I../shared -I../asymmetry trg_time_scan.cpp ../shared/VectorDict.cxx -o trg_time_scan
//            `root-config --cflags --libs` -std=c++11
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "VectorDict.h"
#include "AnConfig.h"

static const int kNS = 6;  // EM0-2, JP0-2
static const char *kScanName[kNS] = {"EM0", "EM1", "EM2", "JP0", "JP1", "JP2"};
static const int kPeriodEdge[2] = {23027047, 23066056};
static int PeriodOf(int run) { return run < kPeriodEdge[0] ? 0 : run < kPeriodEdge[1] ? 1 : 2; }
// cut grid: f = 0.05 ... 2.0
static const int kNF = 196;
static double GridOf(int i) { return 0.05 + 0.01 * i; }
// bit fraction per value bin (x50): bins of 0.05 in f up to 2.5
static const int kNB = 50;
static const double kDB = 0.05;

struct Tally { long hw = 0, off = 0, both = 0; };
struct RunData
{
    std::vector<Tally> t;            // [scan][grid]
    std::vector<long> tot, pass;     // [scan][bin]
    RunData() : t(kNS * kNF), tot(kNS * kNB, 0), pass(kNS * kNB, 0) {}
};

// 0.5 crossing of k/n along x in the given direction (falling: from >= 0.5 to < 0.5);
// err: the x range near it where the 68% Clopper-Pearson band contains 0.5
static bool Cross(const std::vector<double> &xs, const std::vector<long> &ks, const std::vector<long> &ns, bool falling,
                  double &x50, double &lo, double &hi, long nmin = 5)
{
    std::vector<double> x, y, ylo, yhi;
    for (size_t i = 0; i < xs.size(); i++)
    {
        if (ns[i] < nmin) continue;
        x.push_back(xs[i]);
        y.push_back(double(ks[i]) / ns[i]);
        ylo.push_back(TEfficiency::ClopperPearson(ns[i], ks[i], 0.683, false));
        yhi.push_back(TEfficiency::ClopperPearson(ns[i], ks[i], 0.683, true));
    }
    if (x.size() < 3) return false;
    x50 = -1;
    for (size_t i = 1; i < x.size(); i++)
    {
        const double a = y[i - 1] - 0.5, b = y[i] - 0.5;
        if ((falling && a >= 0 && b < 0) || (!falling && a <= 0 && b > 0)) { x50 = x[i - 1] + (x[i] - x[i - 1]) * (a / (a - b)); break; }
    }
    if (x50 < 0) return false;
    lo = hi = x50;
    for (size_t i = 0; i < x.size(); i++)
        if (ylo[i] <= 0.5 && yhi[i] >= 0.5 && std::fabs(x[i] - x50) < 0.3) { lo = std::min(lo, x[i]); hi = std::max(hi, x[i]); }
    return true;
}

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <jet file list> <outdir> [runs per group = 5] [trg_safety.txt]" << std::endl; return 1; }
    const std::string outdir = argv[2];
    const int per_group = argc > 3 ? std::max(1, std::atoi(argv[3])) : 5;
    gSystem->mkdir(outdir.c_str(), true);
    std::map<int, RunData> runs;

    std::ifstream list(argv[1]);
    std::string fn;
    int nfiles = 0;
    while (list >> fn)
    {
        TFile f(fn.c_str());
        TTree *t = f.IsZombie() ? 0 : (TTree *)f.Get("jetTree");
        if (!t || !t->GetBranch("reco_trg_em_et") || !t->GetBranch("fcs_tcu_bits")) { std::cerr << "skip " << fn << std::endl; continue; }
        int file_run = 0, seq = 0;
        const std::string base = fn.substr(fn.find_last_of('/') + 1);
        if (std::sscanf(base.c_str(), "jet_output_%d_%d.root", &file_run, &seq) != 2) file_run = 0;
        Int_t ev_run = 0;
        if (t->GetBranch("event_run")) t->SetBranchAddress("event_run", &ev_run);
        else if (file_run == 0) { std::cerr << "skip (no run) " << fn << std::endl; continue; }
        nfiles++;
        std::vector<float> *em = 0, *hmax = 0, *jp = 0;
        UInt_t tcu = 0;
        // reference triggers, as trg_cut_perf.cpp / trg_eff_scan.cpp
        const char *ref_em_names[] = {"fcsJP2", "fcsJPA1", "fcsJPBC1", "fcsJPDE1", "fcsJPA0", "fcsJPBC0", "fcsJPDE0", "fcsDiJP", "fcsDiJPAsy",
                                      "fcsHad0", "fcsHad1", "fcsHad2", "fcsHHT_N", "fcsHHT_S", "fcsHTOT_N", "fcsHTOT_S", "fcs_HcalTot"};
        const char *ref_jp_names[] = {"fcsEM0", "fcsEM1", "fcsEM0_tpc", "fcsEM1_tpc", "fcsDY", "fcsDYAsy", "fcsJPsi", "fcsELE2"};
        const int n_ref_em = sizeof(ref_em_names) / sizeof(ref_em_names[0]), n_ref_jp = sizeof(ref_jp_names) / sizeof(ref_jp_names[0]);
        std::vector<Int_t> ref_em(n_ref_em, 0), ref_jp(n_ref_jp, 0);
        for (int k = 0; k < n_ref_em; k++) if (t->GetBranch(ref_em_names[k])) t->SetBranchAddress(ref_em_names[k], &ref_em[k]);
        for (int k = 0; k < n_ref_jp; k++) if (t->GetBranch(ref_jp_names[k])) t->SetBranchAddress(ref_jp_names[k], &ref_jp[k]);
        t->SetBranchAddress("reco_trg_em_et", &em);
        t->SetBranchAddress("reco_trg_em_hmax", &hmax);
        t->SetBranchAddress("reco_trg_jp_et", &jp);
        t->SetBranchAddress("fcs_tcu_bits", &tcu);
        for (Long64_t i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);
            const int run = ev_run > 0 ? ev_run : file_run;
            if (run <= 0) continue;
            bool em_ref = false, jp_ref = false;
            for (int k = 0; k < n_ref_em; k++) em_ref |= ref_em[k] != 0;
            for (int k = 0; k < n_ref_jp; k++) jp_ref |= ref_jp[k] != 0;
            if (!em_ref && !jp_ref) continue;
            RunData &rd = runs[run];
            double ljp = -1, lE = -1;
            for (size_t j = 0; j < em->size(); j++) { ljp = std::max(ljp, (double)(*jp)[j]); lE = std::max(lE, (double)(*em)[j]); }
            auto fill = [&](int s, double value, bool hw) {
                for (int fi = 0; fi < kNF; fi++)
                {
                    const bool off = value >= GridOf(fi);
                    Tally &c = rd.t[s * kNF + fi];
                    c.hw += hw; c.off += off; c.both += hw && off;
                }
                const int b = std::min(kNB - 1, std::max(0, (int)(value / kDB)));
                rd.tot[s * kNB + b]++;
                if (hw) rd.pass[s * kNB + b]++;
            };
            if (em_ref && lE > 0)
                for (int l = 0; l < 3; l++)
                {
                    // the event's value: the largest jet min(E_4x4 / E_hw, 1 / r) (threshold and veto, one factor)
                    double v = -1;
                    for (size_t j = 0; j < em->size(); j++)
                    {
                        const double e = (*em)[j], h = (*hmax)[j];
                        if (!(e > 0) || h < 0) continue;
                        const double fe = e / kTrgThresholdGeV[l], fv = h > 0 ? e / (4 * h) : 1e9;
                        v = std::max(v, std::min(fe, fv));
                    }
                    const bool hw = (tcu >> (3 + l)) & 1;
                    // hardware events count with a jet >= 0.5 x E_hw (any veto), as JP
                    if (v > 0) fill(l, v, hw && lE >= 0.5 * kTrgThresholdGeV[l]);
                }
            if (jp_ref && run != 23070019 && ljp > 0)
                for (int l = 0; l < 3; l++)
                {
                    const int tr = 3 + l;
                    const bool hw = tr == 5 ? ((tcu >> 6) & 1) : (((tcu >> (tr == 3 ? 10 : 7)) & 7) != 0);
                    // hardware events count only with a jet near threshold (trg_cut_perf: >= 0.5 x thr)
                    fill(tr, ljp / kTrgThresholdGeV[tr], hw && ljp >= 0.5 * kTrgThresholdGeV[tr]);
                }
        }
    }
    printf("trg_time_scan: %d files, %zu runs\n", nfiles, runs.size());

    // default groups of consecutive runs within a gain period
    std::vector<std::vector<int> > groups;
    {
        std::vector<int> cur;
        int curp = -1;
        for (std::map<int, RunData>::const_iterator it = runs.begin(); it != runs.end(); ++it)
        {
            const int p = PeriodOf(it->first);
            if (p != curp && !cur.empty())
            {
                if ((int)cur.size() < per_group && !groups.empty() && PeriodOf(groups.back().front()) == curp)
                    groups.back().insert(groups.back().end(), cur.begin(), cur.end());
                else groups.push_back(cur);
                cur.clear();
            }
            curp = p;
            cur.push_back(it->first);
            if ((int)cur.size() == per_group) { groups.push_back(cur); cur.clear(); }
        }
        if (!cur.empty())
        {
            if ((int)cur.size() < per_group && !groups.empty() && PeriodOf(groups.back().front()) == curp)
                groups.back().insert(groups.back().end(), cur.begin(), cur.end());
            else groups.push_back(cur);
        }
    }
    // per-trigger groups from trg_safety.txt (a run outside every range joins the nearest group of its period)
    std::vector<std::pair<long, long> > xr[6];
    if (argc > 4)
    {
        std::ifstream sf(argv[4]);
        std::string ln;
        while (std::getline(sf, ln))
        {
            if (ln.empty() || ln[0] == '#') continue;
            std::istringstream ss(ln);
            std::string t; long a, b;
            if (!(ss >> t >> a >> b)) continue;
            for (int tr = 0; tr < 6; tr++) if (t == kScanName[tr]) xr[tr].push_back(std::make_pair(a, b));
        }
    }
    std::vector<std::vector<int> > sgroups[kNS];
    for (int tr = 0; tr < 6; tr++)
    {
        if (xr[tr].empty()) { sgroups[tr] = groups; continue; }
        std::vector<std::vector<int> > g(xr[tr].size());
        for (std::map<int, RunData>::const_iterator it = runs.begin(); it != runs.end(); ++it)
        {
            const long run = it->first;
            int best = -1; long bd = -1;
            for (size_t k = 0; k < xr[tr].size(); k++)
            {
                if (PeriodOf((int)xr[tr][k].first) != PeriodOf((int)run)) continue;
                const long d = run < xr[tr][k].first ? xr[tr][k].first - run : run > xr[tr][k].second ? run - xr[tr][k].second : 0;
                if (best < 0 || d < bd) { best = (int)k; bd = d; }
            }
            if (best >= 0) g[best].push_back((int)run);
        }
        for (size_t k = 0; k < g.size(); k++) if (!g[k].empty()) sgroups[tr].push_back(g[k]);
    }

    std::ofstream txt(outdir + "/trg_time_scan.txt");
    txt << "# per trigger and group: factor f (relative to the hardware settings; EM: E >= f E_hw and r <= 1/f) at\n"
        << "# 50% efficiency, 50% purity and x50 (bit fraction 50%), with 68% ranges\n"
        << "# groups: " << (argc > 4 ? "per trigger from the trg_safety file" : "consecutive runs") << "\n"
        << "# scan first_run last_run n_runs n_hw n_events eff50 lo hi pur50 lo hi x50 lo hi\n";
    TGraphAsymmErrors *ge[kNS], *gp[kNS], *gx[kNS];
    for (int s = 0; s < kNS; s++) { ge[s] = new TGraphAsymmErrors(); gp[s] = new TGraphAsymmErrors(); gx[s] = new TGraphAsymmErrors(); }
    size_t ngroups = 0;
    for (int s = 0; s < kNS; s++)
    {
        const bool veto = false;
        const double u = 1.0;  // results in units of the factor f
        for (size_t g = 0; g < sgroups[s].size(); g++)
        {
            const std::vector<int> &grp = sgroups[s][g];
            ngroups++;
            std::vector<double> xg(kNF), xb(kNB);
            std::vector<long> kb(kNF), nhw(kNF), noff(kNF), bt(kNB, 0), bp(kNB, 0);
            for (int fi = 0; fi < kNF; fi++) xg[fi] = GridOf(fi);
            for (int b = 0; b < kNB; b++) xb[b] = (b + 0.5) * kDB;
            for (size_t r = 0; r < grp.size(); r++)
            {
                const RunData &rd = runs[grp[r]];
                for (int fi = 0; fi < kNF; fi++) { const Tally &c = rd.t[s * kNF + fi]; kb[fi] += c.both; nhw[fi] += c.hw; noff[fi] += c.off; }
                for (int b = 0; b < kNB; b++) { bt[b] += rd.tot[s * kNB + b]; bp[b] += rd.pass[s * kNB + b]; }
            }
            long nev = 0;
            for (int b = 0; b < kNB; b++) nev += bt[b];
            double fe = -1, el = 0, eh = 0, fp = -1, pl = 0, ph = 0, fx = -1, xl = 0, xh = 0;
            // thresholds: efficiency falls, purity rises, bit fraction rises; veto: the reverse
            const bool oke = Cross(xg, kb, nhw, !veto, fe, el, eh), okp = Cross(xg, kb, noff, veto, fp, pl, ph);
            const bool okx = Cross(xb, bp, bt, veto, fx, xl, xh, 10);
            const int r0 = grp.front(), r1 = grp.back();
            const double xc = 0.5 * (r0 + r1);
            char line[400];
            std::snprintf(line, sizeof line, "%s %d %d %zu %ld %ld %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", kScanName[s], r0, r1, grp.size(), nhw[0], nev,
                          oke ? u * fe : -1, oke ? u * el : -1, oke ? u * eh : -1, okp ? u * fp : -1, okp ? u * pl : -1, okp ? u * ph : -1,
                          okx ? u * fx : -1, okx ? u * xl : -1, okx ? u * xh : -1);
            txt << line;
            if (oke) { const int n = ge[s]->GetN(); ge[s]->SetPoint(n, xc, u * fe); ge[s]->SetPointError(n, xc - r0, r1 - xc, u * (fe - el), u * (eh - fe)); }
            if (okp) { const int n = gp[s]->GetN(); gp[s]->SetPoint(n, xc, u * fp); gp[s]->SetPointError(n, xc - r0, r1 - xc, u * (fp - pl), u * (ph - fp)); }
            if (okx) { const int n = gx[s]->GetN(); gx[s]->SetPoint(n, xc, u * fx); gx[s]->SetPointError(n, xc - r0, r1 - xc, u * (fx - xl), u * (xh - fx)); }
        }
    }
    // per period: median x50 over the groups (the period's nominal factor)
    double pmed[kNS][3];
    for (int s = 0; s < kNS; s++)
        for (int p = 0; p < 3; p++)
        {
            std::vector<double> v;
            for (int i = 0; i < gx[s]->GetN(); i++)
            {
                double x, y; gx[s]->GetPoint(i, x, y);
                if (PeriodOf((int)(x + 0.5)) == p) v.push_back(y);
            }
            std::sort(v.begin(), v.end());
            pmed[s][p] = v.empty() ? -1 : (v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]));
            char line[200];
            std::snprintf(line, sizeof line, "# period_median_x50 %s period %d %.3f (%zu groups) = %.3f GeV\n", kScanName[s], p, pmed[s][p], v.size(), pmed[s][p] * kTrgThresholdGeV[s]);
            txt << line;
        }
    txt.close();

    gStyle->SetOptStat(0);
    TCanvas c("c", "c", 1800, 1000);
    c.Divide(3, 2, 0.004, 0.004);
    const double rlo = runs.empty() ? 23000000 : runs.begin()->first - 2000, rhi = runs.empty() ? 23110000 : runs.rbegin()->first + 2000;
    for (int s = 0; s < kNS; s++)
    {
        const bool em = s < 3;
        c.cd(s + 1);
        gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.12);
        const double ymax = 1.6;
        TH1F *fr = new TH1F(Form("fr%d", s), Form(";run number;factor f of the hardware setting%s", em ? "s" : ""), 10, rlo, rhi);
        fr->SetMinimum(0.0); fr->SetMaximum(ymax);
        fr->GetXaxis()->SetNoExponent(); fr->GetXaxis()->SetNdivisions(505);
        fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetLabelSize(0.045);
        fr->Draw("AXIS");
        TLatex tl; tl.SetNDC(); tl.SetTextFont(62); tl.SetTextSize(0.052);
        if (em) tl.DrawLatex(0.13, 0.93, Form("%s:  E_{4x4} #geq f #times %.1f GeV  and  r #leq 1/f", kScanName[s], kTrgThresholdGeV[s]));
        else tl.DrawLatex(0.13, 0.93, Form("%s:  patch #geq f #times %.1f GeV", kScanName[s], kTrgThresholdGeV[s]));
        TLine *thl[3];
        for (int p = 0; p < 3; p++)
        {
            const double a = p == 0 ? rlo : kPeriodEdge[p - 1], b = p == 2 ? rhi : kPeriodEdge[p];
            thl[p] = new TLine(a, pmed[s][p], b, pmed[s][p]); thl[p]->SetLineColor(kGreen + 2); thl[p]->SetLineWidth(3);
            if (pmed[s][p] > 0) thl[p]->Draw();
        }
        for (int k = 0; k < 2; k++)
        {
            TLine *l = new TLine(kPeriodEdge[k], 0.0, kPeriodEdge[k], ymax);
            l->SetLineColor(kGray + 2); l->SetLineStyle(7); l->SetLineWidth(2); l->Draw();
        }
        ge[s]->SetMarkerStyle(20); ge[s]->SetMarkerColor(kAzure + 2); ge[s]->SetLineColor(kAzure + 2);
        gp[s]->SetMarkerStyle(21); gp[s]->SetMarkerColor(kOrange + 7); gp[s]->SetLineColor(kOrange + 7);
        gx[s]->SetMarkerStyle(24); gx[s]->SetMarkerColor(kBlack); gx[s]->SetLineColor(kBlack);
        if (ge[s]->GetN()) ge[s]->Draw("P");
        if (gp[s]->GetN()) gp[s]->Draw("P");
        if (gx[s]->GetN()) gx[s]->Draw("P");
        if (s == 3)
        {
            TLegend *lg = new TLegend(0.14, 0.14, 0.90, 0.42);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
            lg->AddEntry(ge[s], "50% efficiency (cut keeps half the fired events)", "p");
            lg->AddEntry(gp[s], "50% purity (half the cut's events fired)", "p");
            lg->AddEntry(gx[s], "x_{50} (bit set in 50% at that factor)", "p");
            lg->AddEntry(thl[0], "period median x_{50}: the period's nominal factor", "l");
            lg->Draw();
        }
        if (s == 0)
        {
            TLatex tn; tn.SetNDC(); tn.SetTextSize(0.040);
            tn.DrawLatex(0.15, 0.16, "one f scales threshold and veto buffer");
        }
    }
    c.SaveAs((outdir + "/trg_time_scan.png").c_str());
    printf("trg_time_scan: %zu groups; %s/trg_time_scan.{txt,png}\n", ngroups, outdir.c_str());
    return 0;
}
