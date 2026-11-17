// How the nominal trigger levels are determined, per gain period (Sean
// 2026-11-17): the EM threshold and the EM HCal-veto edge measured
// separately, each with the other held where it cannot matter, and the JP
// thresholds. Reference events only (EM: JP, DiJP, Had, HHT, HTOT, HcalTot;
// JP: EM0/1, DY, J/psi, ELE2); hardware bit: fcs_tcu_bits before prescale,
// JP0/JP1 any region. Per event the leading-E_4x4 jet (with an HCal sum) and
// its veto ratio r = 4 H_max / E_4x4 (hardware setting: veto at r > 1), and
// the largest JP patch sum.
//   EM threshold: bit fraction vs E_4x4 / E_hw (Akio's settings 2 / 3.5 / 5
//                 GeV), events whose leading jet has r <= kVetoHold (far below
//                 every veto edge): the 50% crossing is the nominal threshold;
//   EM veto:      bit fraction vs r, events whose leading jet has E_4x4 >=
//                 kTrgHold x that period's EM0 threshold (on the plateau). One
//                 veto rule for EM0-2 (firmware stage_2_202207), so EM0-2 are
//                 combined with each event once: the bit of the level whose
//                 energy slice the jet falls in ([1.3 T0, 1.3 T1) EM0, [1.3 T1,
//                 1.3 T2) EM1, >= 1.3 T2 EM2). The 50% crossing of the combined
//                 curve is the nominal veto edge; the mid and high slices alone
//                 check that the edge does not depend on energy;
//   JP threshold: bit fraction vs patch / E_hw (4 / 6 / 8 GeV).
// One curve per gain period (23027047, 23066056); the crossings are direct
// (linear interpolation between bins, no fit), errors where the edges of the
// 68% Clopper-Pearson band cross 0.5 (interpolated the same way). Hot JP run 23070019 is left out of JP.
// Output: <outdir>/trg_nominal.png (row 1 EM threshold, row 2 EM veto,
// row 3 JP threshold), trg_nominal.txt (the values per period and level);
// the same nine values over Run 22, per group of consecutive runs (enough
// events near the crossing; never across a gain change): trg_nominal_time.png
// and trg_nominal_time.txt.
// Input: unfiltered jet files (data_to_jet -t "", the 125-run scan) or A_N
// yield files with the trigger reference tallies h_trgref_* (AnHistograms
// with the job's reference file, from 2026-11-17; same binning), mixed freely.
// Usage: ./trg_nominal <jet or yield file list> <outdir>
// Build: g++ -m64 -O2 -I../shared -I../asymmetry trg_nominal.cpp ../shared/VectorDict.cxx -o trg_nominal
//            `root-config --cflags --libs` -std=c++11
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2.h>
#include <TNamed.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <functional>
#include <map>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>

#include "VectorDict.h"
#include "AnConfig.h"

static const double kVetoHold = 0.5;  // threshold curves: leading jet r <= 0.5
static const double kTrgHold = 1.3;   // veto curves: leading jet E_4x4 >= 1.3 x the period's nominal threshold
static const int kPeriodEdge[2] = {23027047, 23066056};
static int PeriodOf(int run) { return run < kPeriodEdge[0] ? 0 : run < kPeriodEdge[1] ? 1 : 2; }
static const char *kPeriodLab[3] = {"< 23027047", "23027047-23066055", "#geq 23066056"};
static const int kPCol[3] = {kAzure + 2, kOrange + 7, kBlack}, kPMk[3] = {20, 21, 22};
// E_4x4 / E_hw (or patch / E_hw) bins and r bins
static const int kNE = 80;   static const double kDE = 0.05;  // 0 - 4
static const int kNR = 40;   static const double kDR = 0.05;  // 0 - 2

struct Curve { std::vector<double> x; std::vector<long> n, k; };

// 0.5 crossing along x (rising or falling), interpolated; err: x range near it whose 68% CP band holds 0.5
static bool Cross(const Curve &c, bool falling, double &x50, double &lo, double &hi, long nmin = 10)
{
    std::vector<double> x, y, ylo, yhi;
    for (size_t i = 0; i < c.x.size(); i++)
    {
        if (c.n[i] < nmin) continue;
        x.push_back(c.x[i]);
        y.push_back(double(c.k[i]) / c.n[i]);
        ylo.push_back(TEfficiency::ClopperPearson(c.n[i], c.k[i], 0.683, false));
        yhi.push_back(TEfficiency::ClopperPearson(c.n[i], c.k[i], 0.683, true));
    }
    // first crossing of 0.5 by a curve in the given direction, interpolated
    auto first = [&](const std::vector<double> &yy, double &out) {
        for (size_t i = 1; i < x.size(); i++)
        {
            const double a = yy[i - 1] - 0.5, b = yy[i] - 0.5;
            if ((falling && a >= 0 && b < 0) || (!falling && a <= 0 && b > 0)) { out = x[i - 1] + (x[i] - x[i - 1]) * (a / (a - b)); return true; }
        }
        return false;
    };
    if (!first(y, x50)) return false;
    // the 68% band's edges cross 0.5 on either side of x50 (upper edge first for a rising curve)
    double a = x50, b = x50;
    if (!first(falling ? ylo : yhi, a)) a = x50;
    if (!first(falling ? yhi : ylo, b)) b = x50;
    lo = std::min(a, b); hi = std::max(a, b);
    lo = std::min(lo, x50); hi = std::max(hi, x50);
    return true;
}

static TGraphAsymmErrors *Graph(const Curve &c, double dx, long nmin = 3)
{
    TGraphAsymmErrors *g = new TGraphAsymmErrors();
    for (size_t i = 0; i < c.x.size(); i++)
    {
        if (c.n[i] < nmin) continue;
        const double e = double(c.k[i]) / c.n[i];
        const int n = g->GetN();
        g->SetPoint(n, c.x[i] + dx, e);
        g->SetPointError(n, 0, 0, e - TEfficiency::ClopperPearson(c.n[i], c.k[i], 0.683, false),
                         TEfficiency::ClopperPearson(c.n[i], c.k[i], 0.683, true) - e);
    }
    return g;
}

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <jet file list> <outdir>" << std::endl; return 1; }
    const std::string outdir = argv[2];
    gSystem->mkdir(outdir.c_str(), true);
    // EM: [period][level][E bin][r bin] reference events and bit set; JP: [period][level][bin]
    std::vector<long> emn(3 * 3 * kNE * kNR, 0), emk(3 * 3 * kNE * kNR, 0), jpn(3 * 3 * kNE, 0), jpk(3 * 3 * kNE, 0);
    auto EI = [](int p, int l, int be, int br) { return ((p * 3 + l) * kNE + be) * kNR + br; };
    // the same per run, for the values over time: EM [level][E bin][r bin], JP [level][bin]
    struct RunH { std::vector<long> en, ek, jn, jk; RunH() : en(3 * kNE * kNR, 0), ek(3 * kNE * kNR, 0), jn(3 * kNE, 0), jk(3 * kNE, 0) {} };
    std::map<int, RunH> rh;
    auto RI = [](int l, int be, int br) { return (l * kNE + be) * kNR + br; };

    std::vector<std::string> filt_names;
    long nfilt_drop = 0;
    if (getenv("TRG_FILTER_EMULATE"))
    {
        std::stringstream ss(getenv("TRG_FILTER_EMULATE"));
        std::string tok;
        while (std::getline(ss, tok, ',')) if (!tok.empty()) filt_names.push_back(tok);
        printf("trg_nominal: emulating the trigger filter (%zu triggers)\n", filt_names.size());
    }
    std::ifstream list(argv[1]);
    std::string fn;
    int nfiles = 0;
    long nev_em = 0, nev_jp = 0;
    while (list >> fn)
    {
        TFile f(fn.c_str());
        // A_N yield files (AnHistograms with a trigger reference list, 2026-11-17): the per-job
        // tallies h_trgref_* (the same binning as here), added to their run
        if (!f.IsZombie() && f.Get("h_trgref_em0_n"))
        {
            TNamed *rn = (TNamed *)f.Get("runs");
            const int run = rn ? std::atoi(rn->GetTitle()) : 0;
            if (run <= 0) { std::cerr << "skip (no run) " << fn << std::endl; continue; }
            nfiles++;
            const int p = PeriodOf(run);
            RunH &h = rh[run];
            for (int l = 0; l < 3; l++)
            {
                TH2 *en = (TH2 *)f.Get(Form("h_trgref_em%d_n", l)), *ek = (TH2 *)f.Get(Form("h_trgref_em%d_k", l));
                TH1 *jn = (TH1 *)f.Get(Form("h_trgref_jp%d_n", l)), *jk = (TH1 *)f.Get(Form("h_trgref_jp%d_k", l));
                if (!en || !ek || !jn || !jk || en->GetNbinsX() != kNE || en->GetNbinsY() != kNR || jn->GetNbinsX() != kNE)
                { std::cerr << "ERROR: " << fn << ": h_trgref_* missing or binned differently from trg_nominal" << std::endl; return 1; }
                for (int be = 0; be < kNE; be++)
                {
                    for (int br = 0; br < kNR; br++)
                    {
                        const long n = (long)std::lround(en->GetBinContent(be + 1, br + 1)), k = (long)std::lround(ek->GetBinContent(be + 1, br + 1));
                        emn[EI(p, l, be, br)] += n; emk[EI(p, l, be, br)] += k;
                        h.en[RI(l, be, br)] += n; h.ek[RI(l, be, br)] += k;
                        if (l == 0) nev_em += n;
                    }
                    const long n = (long)std::lround(jn->GetBinContent(be + 1)), k = (long)std::lround(jk->GetBinContent(be + 1));
                    jpn[(p * 3 + l) * kNE + be] += n; jpk[(p * 3 + l) * kNE + be] += k;
                    h.jn[l * kNE + be] += n; h.jk[l * kNE + be] += k;
                    if (l == 0) nev_jp += n;
                }
            }
            continue;
        }
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
        // env TRG_FILTER_EMULATE="<trigger,...>": keep only events listing one of them
        // (the A_N pass's trigger filter), to check what the filter leaves of the reference sample
        std::vector<Int_t> flt(filt_names.size(), 0);
        for (size_t k = 0; k < filt_names.size(); k++)
            if (t->GetBranch(filt_names[k].c_str())) t->SetBranchAddress(filt_names[k].c_str(), &flt[k]);
        for (Long64_t i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);
            const int run = ev_run > 0 ? ev_run : file_run;
            if (run <= 0) continue;
            if (!filt_names.empty())
            {
                bool keep = false;
                for (size_t k = 0; k < flt.size(); k++) keep |= flt[k] != 0;
                if (!keep) { nfilt_drop++; continue; }
            }
            const int p = PeriodOf(run);
            bool em_ref = false, jp_ref = false;
            for (int k = 0; k < n_ref_em; k++) em_ref |= ref_em[k] != 0;
            for (int k = 0; k < n_ref_jp; k++) jp_ref |= ref_jp[k] != 0;
            double lE = -1, lr = -1, ljp = -1;
            for (size_t j = 0; j < em->size(); j++)
            {
                if ((*em)[j] > lE && (*hmax)[j] >= 0) { lE = (*em)[j]; lr = 4 * (*hmax)[j] / (*em)[j]; }
                ljp = std::max(ljp, (double)(*jp)[j]);
            }
            if (em_ref && lE > 0)
            {
                nev_em++;
                const int br = std::min(kNR - 1, (int)(lr / kDR));
                for (int l = 0; l < 3; l++)
                {
                    const int be = std::min(kNE - 1, (int)(lE / kTrgThresholdGeV[l] / kDE));
                    emn[EI(p, l, be, br)]++;
                    if ((tcu >> (3 + l)) & 1) emk[EI(p, l, be, br)]++;
                    RunH &h = rh[run];
                    h.en[RI(l, be, br)]++;
                    if ((tcu >> (3 + l)) & 1) h.ek[RI(l, be, br)]++;
                }
            }
            if (jp_ref && run != 23070019 && ljp > 0)
            {
                nev_jp++;
                for (int l = 0; l < 3; l++)
                {
                    const int tr = 3 + l;
                    const bool hw = tr == 5 ? ((tcu >> 6) & 1) : (((tcu >> (tr == 3 ? 10 : 7)) & 7) != 0);
                    const int b = std::min(kNE - 1, (int)(ljp / kTrgThresholdGeV[tr] / kDE));
                    jpn[(p * 3 + l) * kNE + b]++;
                    if (hw) jpk[(p * 3 + l) * kNE + b]++;
                    RunH &h = rh[run];
                    h.jn[l * kNE + b]++;
                    if (hw) h.jk[l * kNE + b]++;
                }
            }
        }
    }
    printf("trg_nominal: %d files; reference events EM %ld, JP %ld%s\n", nfiles, nev_em, nev_jp,
           filt_names.empty() ? "" : Form(" (after the filter; %ld events dropped by it)", nfilt_drop));

    // curves and crossings
    // cveto[q][p]: q = 0 combined (EM0-2 by energy slice, the value used), 1 mid slice (EM1 bit),
    // 2 high slice (EM2 bit), 3 low slice (EM0 bit; text only)
    Curve cth[3][3], cveto[4][3], cjp[3][3];  // [level or q][period]
    double xth[3][3], xv[4][3], xjp[3][3], eth[3][3][2], ev[4][3][2], ejp[3][3][2];
    bool okth[3][3], okv[4][3], okjp[3][3];
    const int brHold = (int)std::lround(kVetoHold / kDR);  // r bins [0, brHold)
    for (int l = 0; l < 3; l++)
        for (int p = 0; p < 3; p++)
        {
            Curve &c = cth[l][p];
            for (int be = 0; be < kNE; be++)
            {
                long n = 0, k = 0;
                for (int br = 0; br < brHold; br++) { n += emn[EI(p, l, be, br)]; k += emk[EI(p, l, be, br)]; }
                c.x.push_back((be + 0.5) * kDE); c.n.push_back(n); c.k.push_back(k);
            }
            okth[l][p] = Cross(c, false, xth[l][p], eth[l][p][0], eth[l][p][1]);
            Curve &j = cjp[l][p];
            for (int b = 0; b < kNE; b++) { j.x.push_back((b + 0.5) * kDE); j.n.push_back(jpn[(p * 3 + l) * kNE + b]); j.k.push_back(jpk[(p * 3 + l) * kNE + b]); }
            okjp[l][p] = Cross(j, false, xjp[l][p], ejp[l][p][0], ejp[l][p][1]);
        }

    // ---- veto edge: one veto rule for EM0-2 (firmware stage_2_202207), so each event enters once,
    // with the bit of the level whose energy slice its leading jet falls in (T_m = the period's
    // threshold of level m in GeV): [1.3 T0, 1.3 T1) EM0 bit, [1.3 T1, 1.3 T2) EM1 bit, >= 1.3 T2
    // EM2 bit. The slice edges are rounded to the E bins of each level (0.05 x its hardware
    // setting), so a jet near an edge can fall in a gap or overlap of up to half a bin. ----
    // the E-bin range [lo, hi) of slice m (in level m's bins) for period p
    auto sliceBins = [&](int m, int p, int &lo, int &hi) {
        lo = okth[m][p] ? (int)std::ceil(kTrgHold * xth[m][p] / kDE - 1e-9) : kNE;
        if (m == 2) { hi = kNE; return; }
        const double up = okth[m + 1][p] ? kTrgHold * xth[m + 1][p] * kTrgThresholdGeV[m + 1] / kTrgThresholdGeV[m] : 1e9;
        hi = std::min(kNE, (int)std::ceil(up / kDE - 1e-9));
    };
    static const int kSlices[4][3] = {{0, 1, 2}, {1, -1, -1}, {2, -1, -1}, {0, -1, -1}};  // levels per q
    // veto curve q for period p from n/k accessors (level, E bin, r bin)
    auto vetoCurve = [&](int q, int p, const std::function<long(int, int, int)> &N, const std::function<long(int, int, int)> &K) {
        Curve v;
        for (int br = 0; br < kNR; br++)
        {
            long n = 0, k = 0;
            for (int i = 0; i < 3; i++)
            {
                const int m = kSlices[q][i];
                if (m < 0) continue;
                int lo, hi; sliceBins(m, p, lo, hi);
                for (int be = lo; be < hi; be++) { n += N(m, be, br); k += K(m, be, br); }
            }
            v.x.push_back((br + 0.5) * kDR); v.n.push_back(n); v.k.push_back(k);
        }
        return v;
    };
    for (int q = 0; q < 4; q++)
        for (int p = 0; p < 3; p++)
        {
            cveto[q][p] = vetoCurve(q, p, [&](int m, int be, int br) { return emn[EI(p, m, be, br)]; }, [&](int m, int be, int br) { return emk[EI(p, m, be, br)]; });
            okv[q][p] = Cross(cveto[q][p], true, xv[q][p], ev[q][p][0], ev[q][p][1]);
        }

    // ---------------- values over time ----------------
    // scan s: 0-2 EM threshold (level s), 3-5 EM veto edge (level s-3), 6-8 JP threshold (level s-6).
    // Consecutive runs are grouped (never across a gain change) until the group has kMinWin events
    // in the window around the period's value (thresholds: 0.6-1.4 x; veto: 0.5-1.5 x); a short
    // last group of a period joins the previous one.
    auto curveOf = [&](const std::vector<int> &runs_g, int s, int p) {
        Curve c;
        const int l = s % 3;
        if (s < 3)
            for (int be = 0; be < kNE; be++)
            {
                long n = 0, k = 0;
                for (int r : runs_g) for (int br = 0; br < brHold; br++) { n += rh[r].en[RI(l, be, br)]; k += rh[r].ek[RI(l, be, br)]; }
                c.x.push_back((be + 0.5) * kDE); c.n.push_back(n); c.k.push_back(k);
            }
        else if (s < 6)
            c = vetoCurve(l, p, [&](int m, int be, int br) { long n = 0; for (int r : runs_g) n += rh[r].en[RI(m, be, br)]; return n; },
                          [&](int m, int be, int br) { long k = 0; for (int r : runs_g) k += rh[r].ek[RI(m, be, br)]; return k; });
        else
            for (int b = 0; b < kNE; b++)
            {
                long n = 0, k = 0;
                for (int r : runs_g) { n += rh[r].jn[l * kNE + b]; k += rh[r].jk[l * kNE + b]; }
                c.x.push_back((b + 0.5) * kDE); c.n.push_back(n); c.k.push_back(k);
            }
        return c;
    };
    auto pval = [&](int s, int p) { const int l = s % 3; return s < 3 ? (okth[l][p] ? xth[l][p] : -1) : s < 6 ? (okv[l][p] ? xv[l][p] : -1) : (okjp[l][p] ? xjp[l][p] : -1); };
    // the value the analysis uses: thresholds per level; one veto edge for EM0-2, the combined curve
    auto used = [&](int s, int p) { return (s >= 3 && s < 6) ? pval(3, p) : pval(s, p); };
    auto winCount = [&](const Curve &c, double v, bool veto) {
        long n = 0;
        for (size_t i = 0; i < c.x.size(); i++) if (c.x[i] >= (veto ? 0.5 : 0.6) * v && c.x[i] <= (veto ? 1.5 : 1.4) * v) n += c.n[i];
        return n;
    };
    const long kMinWin[3] = {400, 150, 400};  // threshold, veto, JP
    std::vector<int> runlist;
    for (std::map<int, RunH>::const_iterator it = rh.begin(); it != rh.end(); ++it) runlist.push_back(it->first);
    TGraphAsymmErrors *gt[9];
    FILE *tt = std::fopen((outdir + "/trg_nominal_time.txt").c_str(), "w");
    std::fprintf(tt, "# per group of runs: scan first_run last_run n_runs value lo hi (scan: EMl_thr factor, EMl_veto r, JPl_thr factor)\n");
    const char *sname[9] = {"EM0_thr", "EM1_thr", "EM2_thr", "veto_combined", "veto_midE_EM1bit", "veto_highE_EM2bit", "JP0_thr", "JP1_thr", "JP2_thr"};
    for (int s = 0; s < 9; s++)
    {
        gt[s] = new TGraphAsymmErrors();
        if (s == 4 || s == 5) continue;  // veto energy-slice checks: not shown
        const bool veto = s >= 3 && s < 6;
        const long nmin = kMinWin[s / 3];
        for (int p = 0; p < 3; p++)
        {
            const double v = pval(s, p) > 0 ? pval(s, p) : used(s, p);
            if (v <= 0) continue;
            std::vector<std::vector<int> > gs;
            std::vector<int> cur;
            for (int r : runlist)
            {
                if (PeriodOf(r) != p) continue;
                cur.push_back(r);
                if (winCount(curveOf(cur, s, p), v, veto) >= nmin) { gs.push_back(cur); cur.clear(); }
            }
            if (!cur.empty())
            {
                if (!gs.empty() && winCount(curveOf(cur, s, p), v, veto) < nmin / 2) gs.back().insert(gs.back().end(), cur.begin(), cur.end());
                else gs.push_back(cur);
            }
            for (const std::vector<int> &g : gs)
            {
                double x50, lo, hi;
                if (!Cross(curveOf(g, s, p), veto, x50, lo, hi)) continue;
                const double r0 = g.front(), r1 = g.back(), xc = 0.5 * (r0 + r1);
                const int n = gt[s]->GetN();
                gt[s]->SetPoint(n, xc, x50);
                gt[s]->SetPointError(n, xc - r0 + 0.5, r1 - xc + 0.5, x50 - lo, hi - x50);
                std::fprintf(tt, "%s %d %d %zu %.3f %.3f %.3f\n", sname[s], g.front(), g.back(), g.size(), x50, lo, hi);
            }
        }
    }
    std::fclose(tt);

    FILE *tx = std::fopen((outdir + "/trg_nominal.txt").c_str(), "w");
    std::fprintf(tx, "# nominal trigger levels per gain period (trg_nominal.cpp): 50%% crossing of the hardware-bit fraction, reference events\n");
    std::fprintf(tx, "# EM threshold: factor of the hardware setting (leading jet r <= %.2f); EM veto edge: r (leading jet E_4x4 >= %.2f x nominal threshold);\n", kVetoHold, kTrgHold);
    std::fprintf(tx, "# JP threshold: factor of the hardware setting. Columns: kind level period value lo hi [GeV]\n");
    for (int l = 0; l < 3; l++)
        for (int p = 0; p < 3; p++)
        {
            std::fprintf(tx, "EM_threshold EM%d %d %.3f %.3f %.3f %.3f\n", l, p, okth[l][p] ? xth[l][p] : -1, eth[l][p][0], eth[l][p][1], okth[l][p] ? xth[l][p] * kTrgThresholdGeV[l] : -1);
            const char *qn[4] = {"combined", "check_midE_EM1bit", "check_highE_EM2bit", "check_lowE_EM0bit"};
            if (l == 0) for (int q = 0; q < 4; q++) std::fprintf(tx, "EM_veto_edge %s %d %.3f %.3f %.3f\n", qn[q], p, okv[q][p] ? xv[q][p] : -1, ev[q][p][0], ev[q][p][1]);
            std::fprintf(tx, "JP_threshold JP%d %d %.3f %.3f %.3f %.3f\n", l, p, okjp[l][p] ? xjp[l][p] : -1, ejp[l][p][0], ejp[l][p][1], okjp[l][p] ? xjp[l][p] * kTrgThresholdGeV[3 + l] : -1);
        }
    std::fclose(tx);

    gStyle->SetOptStat(0);
    // one panel: curves per period, the crossings as dashed lines, values in the legend
    auto panel = [&](Curve c[3], double x50[3], bool ok[3], double e[3][2], double xmin, double xmax, const char *xt, const char *title,
                     const char *sub, const char *unit_fmt, bool falling) {
        gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.11); gPad->SetBottomMargin(0.13);
        TH1F *fr = new TH1F(Form("fr%p", (void *)c), Form(";%s;fraction with the hardware bit set", xt), 10, xmin, xmax);
        fr->SetMinimum(0); fr->SetMaximum(1.45);
        fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045);
        fr->GetYaxis()->SetTitleOffset(1.15);
        fr->Draw("AXIS");
        TLine *h = new TLine(xmin, 0.5, xmax, 0.5); h->SetLineColor(kGray + 1); h->SetLineStyle(3); h->Draw();
        TLegend *lg = new TLegend(falling ? 0.28 : 0.14, 0.70, falling ? 0.93 : 0.75, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
        for (int p = 0; p < 3; p++)
        {
            TGraphAsymmErrors *g = Graph(c[p], (p - 1) * 0.006);
            g->SetMarkerStyle(kPMk[p]); g->SetMarkerSize(0.8); g->SetMarkerColor(kPCol[p]); g->SetLineColor(kPCol[p]);
            if (g->GetN()) g->Draw("P");
            if (ok[p])
            {
                TLine *v = new TLine(x50[p], 0, x50[p], 1.05); v->SetLineColor(kPCol[p]); v->SetLineStyle(2); v->SetLineWidth(2); v->Draw();
                lg->AddEntry(g, Form("%s: %s", kPeriodLab[p], Form(unit_fmt, x50[p], 0.5 * (e[p][1] - e[p][0]))), "p");
            }
            else lg->AddEntry(g, Form("%s: no crossing", kPeriodLab[p]), "p");
        }
        lg->Draw();
        TLatex tl; tl.SetNDC(); tl.SetTextFont(62); tl.SetTextSize(0.055); tl.DrawLatex(0.12, 0.935, title);
        TLatex ts; ts.SetNDC(); ts.SetTextSize(0.040); ts.SetTextAlign(31); ts.DrawLatex(0.96, 0.935, sub);
    };
    // layout (both pages): 4 x 2; top EM0-2 threshold + the EM veto edge, bottom JP0-2 + the values used
    const int kPad[9] = {1, 2, 3, 4, -1, -1, 5, 6, 7};  // pad of scan s (veto slices 4, 5 not drawn)
    auto valueTable = [&]() {
        TLatex t; t.SetNDC(); t.SetTextSize(0.050);
        t.SetTextFont(62); t.DrawLatex(0.06, 0.90, "values used per gain period");
        t.SetTextFont(42); t.SetTextSize(0.045);
        const double xc[4] = {0.06, 0.34, 0.56, 0.78};
        t.DrawLatex(xc[1], 0.80, "< 23027047"); t.DrawLatex(xc[2], 0.80, "to 23066055"); t.DrawLatex(xc[3], 0.80, "#geq 23066056");
        const char *rn[7] = {"EM0", "EM1", "EM2", "EM veto r", "JP0", "JP1", "JP2"};
        const int rs[7] = {0, 1, 2, 3, 6, 7, 8};
        for (int i = 0; i < 7; i++)
        {
            const double y = 0.71 - 0.085 * i;
            t.DrawLatex(xc[0], y, rn[i]);
            for (int p = 0; p < 3; p++)
            {
                const double v = used(rs[i], p);
                if (v <= 0) { t.DrawLatex(xc[1 + p], y, "-"); continue; }
                const int lv = rs[i] < 3 ? rs[i] : rs[i] - 3;
                if (rs[i] == 3) t.DrawLatex(xc[1 + p], y, Form("%.3f", v));
                else t.DrawLatex(xc[1 + p], y, Form("%.2f GeV", v * kTrgThresholdGeV[lv]));
            }
        }
        t.SetTextSize(0.038);
        t.DrawLatex(0.06, 0.07, "thresholds: 50% point x hardware setting; veto: r edge");
    };
    {
        TCanvas c("ce", "", 2400, 1000);
        c.Divide(4, 2, 0.003, 0.003);
        double x[3], e[3][2]; bool ok[3];
        for (int l = 0; l < 3; l++)
        {
            c.cd(1 + l);
            for (int p = 0; p < 3; p++) { x[p] = xth[l][p]; ok[p] = okth[l][p]; e[p][0] = eth[l][p][0]; e[p][1] = eth[l][p][1]; }
            panel(cth[l], x, ok, e, 0.3, 2.0, Form("E_{4x4} / %.1f GeV (hardware setting)", kTrgThresholdGeV[l]),
                  Form("EM%d threshold", l), Form("veto held: r #leq %.1f", kVetoHold), "%.3f #pm %.3f", false);
            c.cd(5 + l);
            for (int p = 0; p < 3; p++) { x[p] = xjp[l][p]; ok[p] = okjp[l][p]; e[p][0] = ejp[l][p][0]; e[p][1] = ejp[l][p][1]; }
            panel(cjp[l], x, ok, e, 0.3, 2.0, Form("patch sum / %.1f GeV (hardware setting)", kTrgThresholdGeV[3 + l]),
                  Form("JP%d threshold", l), "", "%.3f #pm %.3f", false);
        }
        c.cd(4);
        for (int p = 0; p < 3; p++) { x[p] = xv[0][p]; ok[p] = okv[0][p]; e[p][0] = ev[0][p][0]; e[p][1] = ev[0][p][1]; }
        panel(cveto[0], x, ok, e, 0.0, 2.0, "r = 4 H_{max} / E_{4x4} (hardware veto at r > 1)", "EM veto edge (EM0-2)",
              Form("E_{4x4} #geq %.1f #times thr", kTrgHold), "%.3f #pm %.3f", true);
        c.cd(8);
        valueTable();
        c.SaveAs((outdir + "/trg_nominal.png").c_str());
    }
    {
        TCanvas c("ct", "", 2400, 1000);
        c.Divide(4, 2, 0.003, 0.003);
        const double rlo = runlist.empty() ? 23000000 : runlist.front() - 2000, rhi = runlist.empty() ? 23110000 : runlist.back() + 2000;
        for (int s = 0; s < 9; s++)
        {
            if (kPad[s] < 0) continue;
            const int l = s % 3;
            const bool veto = s == 3;
            c.cd(kPad[s]);
            gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.11); gPad->SetBottomMargin(0.13);
            const double ylo = veto ? 0.4 : 0.5, yhi = veto ? 1.6 : 1.3;
            TH1F *fr = new TH1F(Form("ft%d", s), Form(";run number;%s", veto ? "veto edge r (hardware: 1)" : "threshold / hardware setting"), 10, rlo, rhi);
            fr->SetMinimum(ylo); fr->SetMaximum(yhi);
            fr->GetXaxis()->SetNoExponent(); fr->GetXaxis()->SetNdivisions(505);
            fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetLabelSize(0.045);
            fr->GetYaxis()->SetTitleOffset(1.15);
            fr->Draw("AXIS");
            for (int k = 0; k < 2; k++)
            {
                TLine *b = new TLine(kPeriodEdge[k], ylo, kPeriodEdge[k], yhi);
                b->SetLineColor(kGray + 2); b->SetLineStyle(7); b->SetLineWidth(2); b->Draw();
            }
            for (int p = 0; p < 3; p++)
            {
                const double v = used(s, p);
                if (v <= 0) continue;
                const double a = p == 0 ? rlo : kPeriodEdge[p - 1], b = p == 2 ? rhi : kPeriodEdge[p];
                TLine *h = new TLine(a, v, b, v); h->SetLineColor(kPCol[p]); h->SetLineWidth(2); h->Draw();
            }
            gt[s]->SetMarkerStyle(20); gt[s]->SetMarkerSize(0.8); gt[s]->SetMarkerColor(kBlack); gt[s]->SetLineColor(kBlack);
            if (gt[s]->GetN()) gt[s]->Draw("P");
            TLatex tl; tl.SetNDC(); tl.SetTextFont(62); tl.SetTextSize(0.055);
            if (s < 3) tl.DrawLatex(0.12, 0.935, Form("EM%d threshold (#times %.1f GeV)", l, kTrgThresholdGeV[l]));
            else if (veto) tl.DrawLatex(0.12, 0.935, "EM veto edge (EM0-2)");
            else tl.DrawLatex(0.12, 0.935, Form("JP%d threshold (#times %.1f GeV)", l, kTrgThresholdGeV[3 + l]));
            TLatex ts; ts.SetNDC(); ts.SetTextSize(0.040); ts.SetTextAlign(31);
            ts.DrawLatex(0.96, 0.935, s < 3 ? Form("veto held: r #leq %.1f", kVetoHold) : veto ? Form("E_{4x4} #geq %.1f #times thr", kTrgHold) : "");
            if (s == 0)
            {
                TLatex tn; tn.SetNDC(); tn.SetTextSize(0.038);
                tn.DrawLatex(0.15, 0.17, "points: groups of runs; lines: period values used");
            }
        }
        c.cd(8);
        valueTable();
        c.SaveAs((outdir + "/trg_nominal_time.png").c_str());
    }
    printf("trg_nominal: wrote %s/trg_nominal{,_time}.png and trg_nominal{,_time}.txt\n", outdir.c_str());
    return 0;
}
