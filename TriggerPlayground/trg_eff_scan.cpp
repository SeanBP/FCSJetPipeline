// Trigger efficiency across Run 22 for the per-jet trigger values (RecoJets
// reco_trg_*, shared/FcsTrigSums.h), at the cut factors 1.15 and 1.3 x the
// nominal threshold (AnConfig.h kTrgThresholdGeV), from unfiltered jet files
// (data_to_jet -t "" -E 1 -F 1, one MuDst file per run).
//
// Per event, the leading jet's trigger value: EM, the largest ECal 4x4 among
// jets passing the HCal-veto margin (em_hmax <= em_et / (4 x 1.15)); JP, the
// largest patch sum. The efficiency of trigger t at factor f is the fraction
// of events with leading value >= f x threshold whose hardware bit is set
// (fcs_tcu_bits, before prescale: EM0-2 bits 3-5; JP2 bit 6; JP1 bits 7-9 and
// JP0 bits 10-12 for regions A, BC, DE, the jet's own region required).
// Reference events: recorded by a trigger whose decision is independent of
// the family measured and cannot imply its bit (whether or not the family
// is also listed; requiring "no EM listed" would bias low, since a fired,
// prescale-accepted EM bit lists the EM trigger). EM: fcsJP*, fcsDiJP*,
// fcsHad*, fcsHHT_*, fcsHTOT_*, fcs_HcalTot (not EHT, ELE, DY, J/psi, EM3:
// built on the same ECal windows). JP: fcsEM0/1 (+ _tpc; thresholds below
// JP0's 4 GeV, so they cannot imply a patch above it), fcsDY, fcsDYAsy,
// fcsJPsi, fcsELE2 (not EM2, 5 GeV, or Had: a window above 4-6 GeV implies
// its patch above JP0). Two versions: "integrated" (every event above f x
// threshold) and "edge" (f to f + 0.15 x threshold, the slice a looser cut
// would add, where a plateau problem shows first).
//
// Turn-on safety (trg_safety_EM.png, trg_safety_JP.png, trg_safety.txt): per
// trigger, consecutive runs are grouped (within a gain period, at least
// [runs per point] runs and kMinTurn reference events with x = value /
// threshold in [0.7, 1.7]); per group the reference-event efficiency vs x is
// fitted (binomial likelihood, 0.3 < x < 2.5) with a + (p - a)/2 (1 +
// erf((x - x50) / (sqrt2 sigma))): a is the floor (the bit set by another
// jet or patch well below threshold, 0-5%), p the plateau (threshold-
// independent losses: veto, hot or dead regions, reference bias); both are
// divided out, leaving the turn-on T(x).
// Rows per trigger:
//   1. the turn-on region over time: box from x10 to x99 (T = 0.1 to 0.99),
//      x50 with its fit error, the cuts 1.15 and 1.3 as lines;
//   2. turn-on loss at each cut: the fraction of jets above the cut that the
//      turn-on removes, 1 - sum S / sum (S / T) over x >= f, S the leading
//      eligible jet's x in events the trigger recorded (listed; EM: veto-safe
//      jets, JP0/JP1: the jet's region listed), which is the true spectrum
//      times T for the triggering jet;
//   3. statistics kept: N(x >= 1.3) / N(x >= 1.15) of that spectrum.
//
// Output: <outdir>/trg_eff_integrated.png, trg_eff_edge.png (per group of
// consecutive runs, both factors, gain-period boundaries 23027047 and
// 23066056 dashed; a group never spans a boundary), trg_eff_periods.txt (per
// period) and trg_eff_runs.txt (per run).
//
// Usage: ./trg_eff_scan <jet file list> <outdir> [runs per point = 5]
// Build: g++ -m64 -O2 -I../shared -I../asymmetry trg_eff_scan.cpp ../shared/VectorDict.cxx -o trg_eff_scan
//            `root-config --cflags --libs` -std=c++11

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <TBinomialEfficiencyFitter.h>
#include <TBox.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "VectorDict.h"
#include "AnConfig.h"

static const int kNT = 6;  // EM0, EM1, EM2, JP0, JP1, JP2
static const int kNF = 2;
static const double kFac[kNF] = {1.15, 1.30};
static const int kPeriodEdge[2] = {23027047, 23066056};

struct Count { long n = 0, k = 0; };
// [trigger][factor][integrated, edge]
struct RunCounts { Count c[kNT][kNF][2]; };
// x = trigger value / threshold in bins of kDX up to kXMax, the last index overflow
static const int kNB = 160;
static const double kDX = 0.025, kXMax = 4.0;
static const long kMinTurn = 200;
// runs whose bit fires on more than kHotFloor of the reference events at x in
// [0.2, 0.55] (below the turn-on in every period, x10 >= 0.63: a hot trigger path) are left out of the
// turn-on fits of that trigger and listed in trg_safety.txt
static const double kHotFloor = 0.2;
static int BinOf(double x) { const int b = (int)(x / kDX); return b < 0 ? 0 : b > kNB ? kNB : b; }
// per run and trigger: reference events (tot, pass) and the recorded-sample spectrum
struct RunHists
{
    std::vector<long> tot, pass, spec;
    RunHists() : tot(kNT * (kNB + 1), 0), pass(kNT * (kNB + 1), 0), spec(kNT * (kNB + 1), 0) {}
};
struct GroupFit
{
    int i0 = 0, i1 = 0;            // run index range [i0, i1]
    int first = 0, last = 0;       // runs
    long nturn = 0, nspec = 0;
    bool ok = false;
    double p = 0, a = 0, mu = 0, sig = 0, emu = 0, x10 = 0, x99 = 0, loss[kNF] = {0, 0}, kept = -1;
};

static int PeriodOf(int run) { return run < kPeriodEdge[0] ? 0 : run < kPeriodEdge[1] ? 1 : 2; }

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <jet file list> <outdir> [runs per point = 5]" << std::endl; return 1; }
    const int group = argc > 3 ? std::atoi(argv[3]) : 5;
    if (group < 1) { std::cerr << "ERROR: bad group size" << std::endl; return 1; }
    const std::string outdir = argv[2];
    gSystem->mkdir(outdir.c_str(), true);
    std::map<int, RunCounts> runs;
    std::map<int, RunHists> rhist;
    RunCounts period[3];

    std::ifstream list(argv[1]);
    std::string fn;
    int nfiles = 0;
    while (list >> fn)
    {
        TFile f(fn.c_str());
        TTree *t = f.IsZombie() ? 0 : (TTree *)f.Get("jetTree");
        if (!t || !t->GetBranch("reco_trg_em_et") || !t->GetBranch("fcs_tcu_bits")) { std::cerr << "skip " << fn << std::endl; continue; }
        // the run: event_run per event (merged files), else from jet_output_<run>_<file>.root
        int file_run = 0, seq = 0;
        const std::string base = fn.substr(fn.find_last_of('/') + 1);
        if (std::sscanf(base.c_str(), "jet_output_%d_%d.root", &file_run, &seq) != 2) file_run = 0;
        Int_t ev_run = 0;
        if (t->GetBranch("event_run")) t->SetBranchAddress("event_run", &ev_run);
        else if (file_run == 0) { std::cerr << "skip (no run) " << fn << std::endl; continue; }
        nfiles++;
        std::vector<float> *em = 0, *hmax = 0, *jp = 0;
        std::vector<int> *patch = 0;
        UInt_t tcu = 0;
        Int_t fem[3] = {0}, fem_tpc[3] = {0}, fjp[7] = {0};
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
        t->SetBranchAddress("reco_trg_jp_patch", &patch);
        t->SetBranchAddress("fcs_tcu_bits", &tcu);
        for (int l = 0; l < 3; l++)
        {
            t->SetBranchAddress(TString::Format("fcsEM%d", l), &fem[l]);
            if (t->GetBranch(TString::Format("fcsEM%d_tpc", l))) t->SetBranchAddress(TString::Format("fcsEM%d_tpc", l), &fem_tpc[l]);
        }
        const char *jn[7] = {"fcsJP2", "fcsJPA1", "fcsJPBC1", "fcsJPDE1", "fcsJPA0", "fcsJPBC0", "fcsJPDE0"};
        for (int i = 0; i < 7; i++) t->SetBranchAddress(jn[i], &fjp[i]);
        for (Long64_t i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);
            const int run = ev_run > 0 ? ev_run : file_run;
            if (run <= 0) continue;
            RunCounts &rc = runs[run];
            RunCounts &pc = period[PeriodOf(run)];
            bool em_ref = false, jp_ref = false;
            for (int k = 0; k < n_ref_em; k++) em_ref |= ref_em[k] != 0;
            for (int k = 0; k < n_ref_jp; k++) jp_ref |= ref_jp[k] != 0;
            double lem = -1, ljp = -1;
            int lpatch = -1;
            for (size_t j = 0; j < em->size(); j++)
            {
                if ((*em)[j] > 0 && (*hmax)[j] >= 0 && (*hmax)[j] * 4 * kEmVetoMargin <= (*em)[j] && (*em)[j] > lem) lem = (*em)[j];
                if ((*jp)[j] > ljp) { ljp = (*jp)[j]; lpatch = (*patch)[j]; }
            }
            const int region = lpatch < 0 ? -1 : (lpatch % 5 == 0 ? 0 : lpatch % 5 <= 2 ? 1 : 2);  // A, BC, DE
            RunHists &rh = rhist[run];
            // recorded-sample spectrum: the leading eligible jet in events the trigger recorded
            for (int tr = 0; tr < kNT; tr++)
            {
                double best = -1;
                for (size_t j = 0; j < em->size(); j++)
                {
                    double v = -1;
                    if (tr < 3)
                    {
                        if (!(fem[tr] || fem_tpc[tr])) break;
                        if ((*em)[j] > 0 && (*hmax)[j] >= 0 && (*hmax)[j] * 4 * kEmVetoMargin <= (*em)[j]) v = (*em)[j];
                    }
                    else if (tr == 5)
                    {
                        if (!fjp[0]) break;
                        v = (*jp)[j];
                    }
                    else
                    {
                        const int pj = (*patch)[j];
                        if (pj < 0) continue;
                        const int rg = pj % 5 == 0 ? 0 : pj % 5 <= 2 ? 1 : 2;
                        if (fjp[(tr == 3 ? 4 : 1) + rg]) v = (*jp)[j];
                    }
                    if (v > best) best = v;
                }
                if (best > 0) rh.spec[tr * (kNB + 1) + BinOf(best / kTrgThresholdGeV[tr])]++;
            }
            for (int tr = 0; tr < kNT; tr++)
            {
                const bool isjp = tr >= 3;
                if (!(isjp ? jp_ref : em_ref)) continue;  // reference events only
                const double v = (isjp ? ljp : lem) / kTrgThresholdGeV[tr];
                if (!(v > 0)) continue;
                bool fired;
                if (tr < 3) fired = (tcu >> (3 + tr)) & 1;
                else if (tr == 5) fired = (tcu >> 6) & 1;
                else
                {
                    if (region < 0) continue;
                    fired = (tcu >> ((tr == 3 ? 10 : 7) + region)) & 1;
                }
                rh.tot[tr * (kNB + 1) + BinOf(v)]++;
                if (fired) rh.pass[tr * (kNB + 1) + BinOf(v)]++;
                for (int fi = 0; fi < kNF; fi++)
                {
                    if (v < kFac[fi]) continue;
                    rc.c[tr][fi][0].n++; rc.c[tr][fi][0].k += fired;
                    pc.c[tr][fi][0].n++; pc.c[tr][fi][0].k += fired;
                    if (v < kFac[fi] + 0.15)
                    {
                        rc.c[tr][fi][1].n++; rc.c[tr][fi][1].k += fired;
                        pc.c[tr][fi][1].n++; pc.c[tr][fi][1].k += fired;
                    }
                }
            }
        }
    }
    if (nfiles == 0) { std::cerr << "ERROR: no usable files" << std::endl; return 1; }
    printf("trg_eff_scan: %d files, %zu runs\n", nfiles, runs.size());

    // tables
    const char *pname[3] = {"before 23027047", "23027047-23066055", "from 23066056"};
    FILE *tp = std::fopen((outdir + "/trg_eff_periods.txt").c_str(), "w");
    std::fprintf(tp, "# per gain period: efficiency (events) of the hardware bit, leading jet trigger value >= f x threshold (integrated) and in [f, f+0.15) (edge)\n");
    std::fprintf(tp, "# %-20s %-4s %22s %22s %22s %22s\n", "period", "trg", "f=1.15 integrated", "f=1.15 edge", "f=1.30 integrated", "f=1.30 edge");
    for (int p = 0; p < 3; p++)
        for (int tr = 0; tr < kNT; tr++)
        {
            std::fprintf(tp, "  %-20s %-4s", pname[p], kTrgLevelName[tr]);
            for (int fi = 0; fi < kNF; fi++)
                for (int e = 0; e < 2; e++)
                {
                    const Count &c = period[p].c[tr][fi][e];
                    std::fprintf(tp, " %12.4f (%7ld)", c.n ? (double)c.k / c.n : -1, c.n);
                }
            std::fprintf(tp, "\n");
        }
    std::fclose(tp);
    FILE *tr_ = std::fopen((outdir + "/trg_eff_runs.txt").c_str(), "w");
    std::fprintf(tr_, "# run trg  [f=1.15 integrated k n] [edge k n] [f=1.30 integrated k n] [edge k n]\n");
    for (std::map<int, RunCounts>::const_iterator r = runs.begin(); r != runs.end(); ++r)
        for (int tr = 0; tr < kNT; tr++)
        {
            std::fprintf(tr_, "%d %s", r->first, kTrgLevelName[tr]);
            for (int fi = 0; fi < kNF; fi++)
                for (int e = 0; e < 2; e++) std::fprintf(tr_, " %ld %ld", r->second.c[tr][fi][e].k, r->second.c[tr][fi][e].n);
            std::fprintf(tr_, "\n");
        }
    std::fclose(tr_);

    // groups of consecutive runs, not spanning a gain-period boundary
    std::vector<std::vector<int> > groups;
    for (std::map<int, RunCounts>::const_iterator r = runs.begin(); r != runs.end(); ++r)
    {
        if (groups.empty() || (int)groups.back().size() >= group || PeriodOf(groups.back().back()) != PeriodOf(r->first)) groups.push_back(std::vector<int>());
        groups.back().push_back(r->first);
    }
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(0.06);
    const int col[kNF] = {kBlue + 1, kRed + 1};
    const int mk[kNF] = {24, 20};
    for (int e = 0; e < 2; e++)
    {
        TCanvas *c = new TCanvas(TString::Format("c%d", e), "", 1800, 1000);
        c->Divide(3, 2);
        for (int tr = 0; tr < kNT; tr++)
        {
            c->cd(tr + 1);
            gPad->SetGridy();
            gPad->SetBottomMargin(0.22); gPad->SetTopMargin(0.1);
            TH1D *fr = new TH1D(TString::Format("fr%d_%d", e, tr), "", groups.size(), 0, groups.size());
            fr->SetTitle(TString::Format("%s (%.1f GeV): %s;first run of each group of %d;efficiency of the hardware bit", kTrgLevelName[tr], kTrgThresholdGeV[tr],
                                         e ? "trigger value in [f, f+0.15) #times thr" : "trigger value #geq f #times thr", group));
            fr->SetMinimum(e ? 0.0 : 0.5); fr->SetMaximum(1.08);
            for (size_t i = 0; i < groups.size(); i++) if (i % 3 == 0) fr->GetXaxis()->SetBinLabel(i + 1, TString::Format("%d", groups[i].front()));
            fr->GetXaxis()->SetLabelSize(0.04);
            fr->GetXaxis()->LabelsOption("v");
            fr->GetXaxis()->SetTitleOffset(2.3);
            fr->Draw("axis");
            TLegend *lg = new TLegend(0.14, 0.25, 0.75, 0.40);
            lg->SetBorderSize(0); lg->SetFillStyle(0);
            for (int fi = 0; fi < kNF; fi++)
            {
                TGraphAsymmErrors *g = new TGraphAsymmErrors();
                long K = 0, N = 0;
                for (size_t i = 0; i < groups.size(); i++)
                {
                    long k = 0, n = 0;
                    for (size_t r = 0; r < groups[i].size(); r++) { k += runs[groups[i][r]].c[tr][fi][e].k; n += runs[groups[i][r]].c[tr][fi][e].n; }
                    K += k; N += n;
                    if (n < 3) continue;
                    const double eff = (double)k / n;
                    const double lo = TEfficiency::ClopperPearson(n, k, 0.683, false), hi = TEfficiency::ClopperPearson(n, k, 0.683, true);
                    const int ip = g->GetN();
                    g->SetPoint(ip, i + 0.5 + (fi ? 0.15 : -0.15), eff);
                    g->SetPointError(ip, 0, 0, eff - lo, hi - eff);
                }
                g->SetMarkerStyle(mk[fi]); g->SetMarkerColor(col[fi]); g->SetLineColor(col[fi]); g->SetMarkerSize(0.8);
                if (g->GetN()) g->Draw("p same");
                lg->AddEntry(g, TString::Format("f = %.2f: all %.3f (%ld events)", kFac[fi], N ? (double)K / N : 0, N), "pe");
            }
            for (int b = 0; b < 2; b++)
            {
                size_t i = 0;
                while (i < groups.size() && groups[i].front() < kPeriodEdge[b]) i++;
                TLine *ln = new TLine(i, e ? 0.0 : 0.5, i, 1.08);
                ln->SetLineStyle(2); ln->SetLineColor(kGray + 2); ln->SetLineWidth(2); ln->Draw();
            }
            lg->Draw();
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.055);
            tl.DrawLatex(0.12, 0.93, TString::Format("%s (%.1f GeV): %s", kTrgLevelName[tr], kTrgThresholdGeV[tr],
                                                    e ? "trigger value in [f, f+0.15) #times thr" : "trigger value #geq f #times thr"));
        }
        c->SaveAs((outdir + (e ? "/trg_eff_edge.png" : "/trg_eff_integrated.png")).c_str());
    }
    // ---------------- turn-on safety ----------------
    std::vector<int> rlist;
    for (std::map<int, RunCounts>::const_iterator r = runs.begin(); r != runs.end(); ++r) rlist.push_back(r->first);
    const int nr = rlist.size();
    std::vector<GroupFit> gf[kNT];
    std::vector<std::string> hot_lines;
    std::map<int, bool> hot[kNT];
    for (int tr = 0; tr < kNT; tr++)
        for (int i = 0; i < nr; i++)
        {
            const RunHists &rh = rhist[rlist[i]];
            double k = 0, n = 0;
            for (int b = BinOf(0.2); b < BinOf(0.55); b++) { k += rh.pass[tr * (kNB + 1) + b]; n += rh.tot[tr * (kNB + 1) + b]; }
            if (n >= 10 && k / n > kHotFloor)
            {
                hot[tr][rlist[i]] = true;
                hot_lines.push_back(TString::Format("# hot: %s run %d fires %.0f of %.0f reference events at x 0.2-0.55 (%.2f); left out of the fits", kTrgLevelName[tr], rlist[i], k, n, k / n).Data());
            }
        }
    for (int tr = 0; tr < kNT; tr++)
    {
        // adaptive groups
        GroupFit cur; cur.i0 = 0;
        for (int i = 0; i < nr; i++)
        {
            if (i > cur.i0 && PeriodOf(rlist[i]) != PeriodOf(rlist[cur.i0]))
            {
                cur.i1 = i - 1; gf[tr].push_back(cur); cur = GroupFit(); cur.i0 = i;
            }
            const RunHists &rh = rhist[rlist[i]];
            if (hot[tr].count(rlist[i])) continue;
            for (int b = BinOf(0.7); b < BinOf(1.7); b++) cur.nturn += rh.tot[tr * (kNB + 1) + b];
            if (cur.nturn >= kMinTurn && i - cur.i0 + 1 >= group)
            {
                cur.i1 = i; gf[tr].push_back(cur); cur = GroupFit(); cur.i0 = i + 1;
            }
        }
        if (cur.i0 < nr) { cur.i1 = nr - 1; gf[tr].push_back(cur); }
        // a small last group of a period joins the previous group of the same period
        for (size_t g = 1; g < gf[tr].size(); g++)
        {
            GroupFit &a = gf[tr][g - 1], &b = gf[tr][g];
            if (b.nturn < kMinTurn / 2 && PeriodOf(rlist[a.i0]) == PeriodOf(rlist[b.i0]))
            {
                a.i1 = b.i1; a.nturn += b.nturn; gf[tr].erase(gf[tr].begin() + g); g--;
            }
        }
        // fit each group
        for (size_t g = 0; g < gf[tr].size(); g++)
        {
            GroupFit &G = gf[tr][g];
            G.first = rlist[G.i0]; G.last = rlist[G.i1];
            TH1D ht(TString::Format("ht_%d_%d", tr, (int)g), "", kNB, 0, kXMax), hp(TString::Format("hp_%d_%d", tr, (int)g), "", kNB, 0, kXMax);
            std::vector<double> S(kNB + 1, 0);
            for (int i = G.i0; i <= G.i1; i++)
            {
                if (hot[tr].count(rlist[i])) continue;
                const RunHists &rh = rhist[rlist[i]];
                for (int b = 0; b < kNB; b++)
                {
                    ht.AddBinContent(b + 1, rh.tot[tr * (kNB + 1) + b]);
                    hp.AddBinContent(b + 1, rh.pass[tr * (kNB + 1) + b]);
                }
                for (int b = 0; b <= kNB; b++) S[b] += rh.spec[tr * (kNB + 1) + b];
            }
            for (int b = 0; b <= kNB; b++) G.nspec += (long)S[b];
            double pk = 0, pn = 0;
            for (int b = BinOf(1.6); b < BinOf(2.5); b++) { pk += hp.GetBinContent(b + 1); pn += ht.GetBinContent(b + 1); }
            TF1 f(TString::Format("f_%d_%d", tr, (int)g), "[3]+([0]-[3])*0.5*(1+TMath::Erf((x-[1])/(sqrt(2.)*[2])))", 0.3, 2.5);
            f.SetParameters(pn > 0 ? std::max(0.5, pk / pn) : 0.95, 1.0, 0.1, 0.01);
            f.SetParLimits(0, 0.2, 1.0); f.SetParLimits(1, 0.3, 2.2); f.SetParLimits(2, 0.01, 0.8); f.SetParLimits(3, 0.0, 0.3);
            if (ht.Integral() < 20) continue;
            TBinomialEfficiencyFitter fitter(&hp, &ht);
            TFitResultPtr r = fitter.Fit(&f, "RSQ");
            if ((int)r != 0) continue;
            G.ok = true;
            G.p = f.GetParameter(0); G.a = f.GetParameter(3); G.mu = f.GetParameter(1); G.sig = f.GetParameter(2); G.emu = f.GetParError(1);
            G.x10 = G.mu - 1.2816 * G.sig; G.x99 = G.mu + 2.3263 * G.sig;
            // turn-on loss above each cut, and the statistics kept
            for (int fi = 0; fi < kNF; fi++)
            {
                double nrec = 0, ntrue = 0;
                for (int b = (int)std::lround(kFac[fi] / kDX); b <= kNB; b++)
                {
                    const double x = b < kNB ? (b + 0.5) * kDX : kXMax;
                    const double T = std::max(0.02, 0.5 * (1 + TMath::Erf((x - G.mu) / (std::sqrt(2.) * G.sig))));
                    nrec += S[b]; ntrue += S[b] / T;
                }
                G.loss[fi] = ntrue > 0 ? 1 - nrec / ntrue : -1;
            }
            double n115 = 0, n130 = 0;
            for (int b = (int)std::lround(1.15 / kDX); b <= kNB; b++) { n115 += S[b]; if (b >= (int)std::lround(1.30 / kDX)) n130 += S[b]; }
            G.kept = n115 > 0 ? n130 / n115 : -1;
        }
    }
    {
        FILE *ts = std::fopen((outdir + "/trg_safety.txt").c_str(), "w");
        std::fprintf(ts, "# turn-on fit per group of runs (reference events, x = trigger value / nominal threshold); T(x) = erf turn-on with the plateau p divided out\n");
        std::fprintf(ts, "# loss_f = fraction of jets above f x thr removed by the turn-on (recorded-sample spectrum, leading eligible jet); kept = N(>=1.3)/N(>=1.15)\n");
        for (size_t h = 0; h < hot_lines.size(); h++) std::fprintf(ts, "%s\n", hot_lines[h].c_str());
        std::fprintf(ts, "# trg first_run last_run n_turn n_spec fit plateau floor x50 +- sigma x10 x99 z115 z130 loss115 loss130 kept\n");
        for (int tr = 0; tr < kNT; tr++)
            for (size_t g = 0; g < gf[tr].size(); g++)
            {
                const GroupFit &G = gf[tr][g];
                std::fprintf(ts, "%-4s %d %d %6ld %7ld %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %6.2f %6.2f %.5f %.5f %.3f\n", kTrgLevelName[tr], G.first, G.last, G.nturn, G.nspec,
                             G.ok, G.p, G.a, G.mu, G.emu, G.sig, G.x10, G.x99, G.ok ? (1.15 - G.mu) / G.sig : 0, G.ok ? (1.30 - G.mu) / G.sig : 0, G.loss[0], G.loss[1], G.kept);
            }
        std::fclose(ts);
    }
    // plots: one canvas per family, columns = levels, rows = turn-on band, loss, kept
    for (int fam = 0; fam < 2; fam++)
    {
        TCanvas *c = new TCanvas(TString::Format("cs%d", fam), "", 1800, 1500);
        c->Divide(3, 3, 0.002, 0.002);
        for (int col3 = 0; col3 < 3; col3++)
        {
            const int tr = fam * 3 + col3;
            for (int row = 0; row < 3; row++)
            {
                c->cd(row * 3 + col3 + 1);
                gPad->SetBottomMargin(0.2); gPad->SetTopMargin(0.1); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
                gPad->SetGridy();
                if (row == 1) gPad->SetLogy();
                TH1D *fr = new TH1D(TString::Format("fs%d_%d", tr, row), "", nr, 0, nr);
                const char *ytit[3] = {"trigger value / nominal threshold", "turn-on loss above the cut", "N(#geq 1.3) / N(#geq 1.15)"};
                fr->GetYaxis()->SetTitle(ytit[row]);
                fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.25);
                if (row == 0) { fr->SetMinimum(0.4); fr->SetMaximum(1.9); }
                if (row == 1) { fr->SetMinimum(1e-4); fr->SetMaximum(1.0); }
                if (row == 2) { fr->SetMinimum(0.0); fr->SetMaximum(1.0); }
                for (int i = 0; i < nr; i += 10) fr->GetXaxis()->SetBinLabel(i + 1, TString::Format("%d", rlist[i]));
                fr->GetXaxis()->SetLabelSize(0.05); fr->GetXaxis()->LabelsOption("v");
                fr->Draw("axis");
                const double y0 = fr->GetMinimum(), y1 = fr->GetMaximum();
                TGraphErrors *gx50 = new TGraphErrors(), *gk = new TGraphErrors(), *gbad = new TGraphErrors();
                TGraphErrors *gl[kNF] = {new TGraphErrors(), new TGraphErrors()};
                for (size_t g = 0; g < gf[tr].size(); g++)
                {
                    const GroupFit &G = gf[tr][g];
                    const double xc = 0.5 * (G.i0 + G.i1 + 1), xw = 0.5 * (G.i1 - G.i0 + 1);
                    if (!G.ok)
                    {
                        if (row == 0) { const int ip = gbad->GetN(); gbad->SetPoint(ip, xc, 0.45); gbad->SetPointError(ip, xw, 0); }
                        continue;
                    }
                    if (row == 0)
                    {
                        TBox *bx = new TBox(G.i0 + 0.1, std::max(y0, G.x10), G.i1 + 0.9, std::min(y1, G.x99));
                        bx->SetFillColor(G.nturn < kMinTurn ? kGray : kAzure - 9);
                        bx->Draw();
                        const int ip = gx50->GetN();
                        gx50->SetPoint(ip, xc, G.mu); gx50->SetPointError(ip, 0, G.emu);
                    }
                    else if (row == 1)
                        for (int fi = 0; fi < kNF; fi++)
                        {
                            const int ip = gl[fi]->GetN();
                            gl[fi]->SetPoint(ip, xc, std::max(1.2e-4, G.loss[fi]));
                            gl[fi]->SetPointError(ip, xw, 0);
                        }
                    else if (G.kept >= 0) { const int ip = gk->GetN(); gk->SetPoint(ip, xc, G.kept); gk->SetPointError(ip, xw, 0); }
                }
                if (row == 0)
                {
                    for (int fi = 0; fi < kNF; fi++)
                    {
                        TLine *ln = new TLine(0, kFac[fi], nr, kFac[fi]);
                        ln->SetLineColor(fi ? kRed + 1 : kBlue + 1); ln->SetLineWidth(2); ln->SetLineStyle(fi ? 1 : 7); ln->Draw();
                    }
                    gx50->SetMarkerStyle(20); gx50->SetMarkerSize(0.7); gx50->SetMarkerColor(kBlack); gx50->SetLineColor(kBlack);
                    if (gx50->GetN()) gx50->Draw("p same");
                    gbad->SetMarkerStyle(5); gbad->SetMarkerColor(kGray + 2); gbad->SetLineColor(kGray + 2);
                    if (gbad->GetN()) gbad->Draw("p same");
                    fr->Draw("axis same");
                }
                if (row == 1)
                    for (int fi = 0; fi < kNF; fi++)
                    {
                        gl[fi]->SetMarkerStyle(fi ? 20 : 24); gl[fi]->SetMarkerSize(0.8);
                        gl[fi]->SetMarkerColor(fi ? kRed + 1 : kBlue + 1); gl[fi]->SetLineColor(fi ? kRed + 1 : kBlue + 1);
                        if (gl[fi]->GetN()) gl[fi]->Draw("p same");
                    }
                if (row == 2)
                {
                    gk->SetMarkerStyle(21); gk->SetMarkerSize(0.8); gk->SetMarkerColor(kGreen + 3); gk->SetLineColor(kGreen + 3);
                    if (gk->GetN()) gk->Draw("p same");
                }
                for (int b = 0; b < 2; b++)
                {
                    int i = 0;
                    while (i < nr && rlist[i] < kPeriodEdge[b]) i++;
                    TLine *ln = new TLine(i, y0, i, y1);
                    ln->SetLineStyle(2); ln->SetLineColor(kGray + 2); ln->SetLineWidth(2); ln->Draw();
                }
                TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
                const char *rt[3] = {"x_{10}-x_{99} box, x_{50}; cut 1.15 - -, 1.3 #minus",
                                     "loss: 1.15 #circ, 1.3 #bullet (floor 10^{-4})",
                                     "kept: N(#geq1.3)/N(#geq1.15)"};
                tl.DrawLatex(0.13, 0.93, TString::Format("#bf{%s} (%g GeV): %s", kTrgLevelName[tr], kTrgThresholdGeV[tr], rt[row]));
            }
        }
        c->SaveAs((outdir + (fam ? "/trg_safety_JP.png" : "/trg_safety_EM.png")).c_str());
    }
    printf("wrote %s/trg_eff_integrated.png, trg_eff_edge.png, trg_eff_periods.txt, trg_eff_runs.txt, trg_safety_EM.png, trg_safety_JP.png, trg_safety.txt\n", outdir.c_str());
    return 0;
}
