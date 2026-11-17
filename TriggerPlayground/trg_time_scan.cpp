// Trigger thresholds over Run 22 (Sean 2026-11-17): per group of consecutive
// runs and per trigger, the jet trigger value (in units of the NOMINAL
// threshold, AnConfig.h kTrgThresholdGeV, so steps between gain periods show)
// at which
//   efficiency = N(offline and hardware) / N(hardware with a jet)   = 50%
//   purity     = N(offline and hardware) / N(offline)               = 50%
// with the offline decision "some jet has trigger value >= f x nominal", the
// definitions of trg_cut_perf.cpp (reference events only; EM: threshold and
// HCal veto as one cut, the jet's value min(E_4x4, nominal / r), r = 4 H_max /
// E_4x4, so f moves the ECal cut and the veto (nominal edge r = 1) together; hardware bit:
// fcs_tcu_bits before prescale, JP0/JP1 any region). Efficiency falls and
// purity rises with f; each crossing is interpolated on a 0.01 grid, its
// error the range of f where the 68% Clopper-Pearson band contains 0.5.
// Runs are grouped consecutively inside a gain period (23027047, 23066056)
// until a group has [runs per group] runs; the last short group of a period
// joins the previous one.
// Input: jet files with reco_trg_em_et / reco_trg_em_hmax / reco_trg_jp_et,
// fcs_tcu_bits and the reference trigger branches (data_to_jet; the Run 22
// trigger scan trigger_turnon/run22_scan_v1, or a full pass's kept jet trees).
// Output: <outdir>/trg_time_scan.txt (per group and trigger: first/last run,
// counts, f_eff50 +- , f_pur50 +- ) and trg_time_scan.png (one panel per
// trigger, run number on x, both crossings on y, gain changes dashed).
// Hot JP run 23070019 (JP0-2 fire in every reference event far below
// threshold; trg_eff_scan.cpp's hot test) is left out of the JP triggers.
// Optional 4th argument: trg_eff_scan.cpp's trg_safety.txt, whose turn-on
// x50 per group (bit set in 50% of the reference events at that trigger
// value; independent of the jet spectrum) is overlaid as open circles.
// Usage: ./trg_time_scan <jet file list> <outdir> [runs per group = 5] [trg_safety.txt]
// Build: g++ -m64 -O2 -I../shared -I../asymmetry trg_time_scan.cpp ../shared/VectorDict.cxx -o trg_time_scan
//            `root-config --cflags --libs` -std=c++11
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
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
#include <string>
#include <vector>

#include "VectorDict.h"
#include "AnConfig.h"
#include <TLatex.h>
#include <sstream>

static const int kNT = 6;  // EM0, EM1, EM2, JP0, JP1, JP2
static const char *kTrgName[kNT] = {"EM0", "EM1", "EM2", "JP0", "JP1", "JP2"};
static const int kPeriodEdge[2] = {23027047, 23066056};
static int PeriodOf(int run) { return run < kPeriodEdge[0] ? 0 : run < kPeriodEdge[1] ? 1 : 2; }
static const int kNF = 196;
static double FacOf(int i) { return 0.05 + 0.01 * i; }  // 0.05 ... 2.0 x nominal

struct Tally { long hwj = 0, off = 0, both = 0; };
// per run: [trigger][factor]
typedef std::vector<Tally> RunTally;

// the f where the curve (num/den per factor) crosses 0.5, falling (eff) or rising (pur);
// err: the f range where the 68% Clopper-Pearson band contains 0.5
static bool Crossing(const std::vector<Tally> &v, int tr, bool eff, double &f50, double &lo, double &hi)
{
    std::vector<double> x, y, ylo, yhi;
    for (int i = 0; i < kNF; i++)
    {
        const Tally &c = v[tr * kNF + i];
        const long den = eff ? c.hwj : c.off;
        if (den < 5) continue;
        x.push_back(FacOf(i));
        y.push_back(double(c.both) / den);
        ylo.push_back(TEfficiency::ClopperPearson(den, c.both, 0.683, false));
        yhi.push_back(TEfficiency::ClopperPearson(den, c.both, 0.683, true));
    }
    if (x.size() < 3) return false;
    // central crossing: first sign change of (y - 0.5) in the expected direction
    f50 = -1;
    for (size_t i = 1; i < x.size(); i++)
    {
        const double a = y[i - 1] - 0.5, b = y[i] - 0.5;
        if ((eff && a >= 0 && b < 0) || (!eff && a <= 0 && b > 0))
        {
            f50 = x[i - 1] + (x[i] - x[i - 1]) * (a / (a - b));
            break;
        }
    }
    if (f50 < 0) return false;
    lo = hi = f50;
    for (size_t i = 0; i < x.size(); i++)
        if (ylo[i] <= 0.5 && yhi[i] >= 0.5 && std::fabs(x[i] - f50) < 0.3) { lo = std::min(lo, x[i]); hi = std::max(hi, x[i]); }
    return true;
}

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <jet file list> <outdir> [runs per group = 5] [trg_safety.txt]" << std::endl; return 1; }
    const std::string outdir = argv[2];
    const int per_group = argc > 3 ? std::max(1, std::atoi(argv[3])) : 5;
    gSystem->mkdir(outdir.c_str(), true);
    std::map<int, RunTally> runs;

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
            RunTally &rt = runs[run];
            if (rt.empty()) rt.resize(kNT * kNF);
            for (int tr = 0; tr < kNT; tr++)
            {
                const bool isjp = tr >= 3;
                if (!(isjp ? jp_ref : em_ref)) continue;
                if (isjp && run == 23070019) continue;  // hot JP run
                bool hw;
                if (tr < 3) hw = (tcu >> (3 + tr)) & 1;
                else if (tr == 5) hw = (tcu >> 6) & 1;
                else hw = ((tcu >> (tr == 3 ? 10 : 7)) & 7) != 0;
                const double nom = kTrgThresholdGeV[tr];
                double anyv = -1, lead = -1;
                for (size_t j = 0; j < em->size(); j++)
                {
                    anyv = std::max(anyv, (double)(isjp ? (*jp)[j] : (*em)[j]));
                    double v;
                    if (isjp) v = (*jp)[j];
                    // EM threshold and veto as one cut (2026-11-17): the jet's value in nominal-threshold
                    // units is min(E_4x4, nominal x 1 / r), the nominal veto edge being r = 1
                    else if (!((*em)[j] > 0) || (*hmax)[j] < 0) v = -1;
                    else v = (*hmax)[j] > 0 ? std::min((double)(*em)[j], nom * (*em)[j] / (4 * (*hmax)[j])) : (*em)[j];
                    lead = std::max(lead, v);
                }
                // hardware events with a jet near the threshold (trg_cut_perf: >= 0.5 x the
                // measured threshold; here 0.5 x the measured threshold as well)
                const bool hwj = hw && anyv >= 0.5 * TrgThresholdGeV(tr, run);
                for (int fi = 0; fi < kNF; fi++)
                {
                    const bool off = lead >= FacOf(fi) * nom;
                    Tally &c = rt[tr * kNF + fi];
                    c.hwj += hwj; c.off += off; c.both += hw && off;
                }
            }
        }
    }
    printf("trg_time_scan: %d files, %zu runs\n", nfiles, runs.size());

    // groups of consecutive runs within a gain period
    std::vector<std::vector<int> > groups;
    {
        std::vector<int> cur;
        int curp = -1;
        for (std::map<int, RunTally>::const_iterator it = runs.begin(); it != runs.end(); ++it)
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

    std::ofstream txt(outdir + "/trg_time_scan.txt");
    txt << "# trigger value / NOMINAL threshold where efficiency = N(both)/N(hw with a jet) and purity = N(both)/N(offline) cross 50%\n"
        << "# (reference events; EM with the per-period veto margin); err = f range where the 68% Clopper-Pearson band holds 0.5\n"
        << "# trg first_run last_run n_runs n_hw n_offline_at_1 f_eff50 lo hi f_pur50 lo hi pur_at_0.05 pur_at_1\n";
    TGraphAsymmErrors *ge[kNT], *gp[kNT];
    for (int tr = 0; tr < kNT; tr++) { ge[tr] = new TGraphAsymmErrors(); gp[tr] = new TGraphAsymmErrors(); }
    for (size_t g = 0; g < groups.size(); g++)
    {
        std::vector<Tally> sum(kNT * kNF);
        for (size_t r = 0; r < groups[g].size(); r++)
        {
            const RunTally &rt = runs[groups[g][r]];
            for (int k = 0; k < kNT * kNF; k++) { sum[k].hwj += rt[k].hwj; sum[k].off += rt[k].off; sum[k].both += rt[k].both; }
        }
        const int r0 = groups[g].front(), r1 = groups[g].back();
        const double xc = 0.5 * (r0 + r1);
        for (int tr = 0; tr < kNT; tr++)
        {
            double fe = -1, el = 0, eh = 0, fp = -1, pl = 0, ph = 0;
            const bool oke = Crossing(sum, tr, true, fe, el, eh), okp = Crossing(sum, tr, false, fp, pl, ph);
            const int i1 = (int)std::lround((1.0 - 0.05) / 0.01);
            // purity at the lowest f (the bit's overall rate in the reference sample) and at 1 x nominal
            const Tally &c0 = sum[tr * kNF], &c1 = sum[tr * kNF + i1];
            const double pur0 = c0.off ? double(c0.both) / c0.off : -1, pur1 = c1.off ? double(c1.both) / c1.off : -1;
            char line[320];
            std::snprintf(line, sizeof line, "%s %d %d %zu %ld %ld %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", kTrgName[tr], r0, r1, groups[g].size(), sum[tr * kNF].hwj,
                          sum[tr * kNF + i1].off, oke ? fe : -1, oke ? el : -1, oke ? eh : -1, okp ? fp : -1, okp ? pl : -1, okp ? ph : -1, pur0, pur1);
            txt << line;
            if (oke) { const int n = ge[tr]->GetN(); ge[tr]->SetPoint(n, xc, fe); ge[tr]->SetPointError(n, xc - r0, r1 - xc, fe - el, eh - fe); }
            if (okp) { const int n = gp[tr]->GetN(); gp[tr]->SetPoint(n, xc, fp); gp[tr]->SetPointError(n, xc - r0, r1 - xc, fp - pl, ph - fp); }
        }
    }
    txt.close();

    // optional overlay: turn-on x50 per group from trg_safety.txt
    TGraphAsymmErrors *gx[kNT];
    for (int tr = 0; tr < kNT; tr++) gx[tr] = new TGraphAsymmErrors();
    if (argc > 4)
    {
        std::ifstream sf(argv[4]);
        std::string ln;
        while (std::getline(sf, ln))
        {
            if (ln.empty() || ln[0] == '#') continue;
            std::istringstream ss(ln);
            std::string t; long a, b; double nt, ns, fit, pl, fl, x50, ex;
            if (!(ss >> t >> a >> b >> nt >> ns >> fit >> pl >> fl >> x50 >> ex)) continue;
            for (int tr = 0; tr < kNT; tr++)
                if (t == kTrgName[tr]) { const int n = gx[tr]->GetN(); const double xc = 0.5 * (a + b); gx[tr]->SetPoint(n, xc, x50); gx[tr]->SetPointError(n, xc - a, b - xc, ex, ex); }
        }
    }
    gStyle->SetOptStat(0);
    TCanvas c("c", "c", 1800, 1000);
    c.Divide(3, 2, 0.004, 0.004);
    const double rlo = runs.empty() ? 23000000 : runs.begin()->first - 2000, rhi = runs.empty() ? 23110000 : runs.rbegin()->first + 2000;
    for (int tr = 0; tr < kNT; tr++)
    {
        c.cd(tr + 1);
        gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.12);
        TH1F *fr = new TH1F(Form("fr%d", tr), Form("%s (nominal %.1f GeV);run number;jet trigger value / nominal threshold", kTrgName[tr], kTrgThresholdGeV[tr]),
                            10, rlo, rhi);
        fr->SetMinimum(0.0); fr->SetMaximum(1.6);
        fr->GetXaxis()->SetNoExponent(); fr->GetXaxis()->SetNdivisions(505);
        fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetLabelSize(0.045);
        fr->Draw("AXIS");
        TLatex tl; tl.SetNDC(); tl.SetTextFont(62); tl.SetTextSize(0.055);
        tl.DrawLatex(0.13, 0.93, Form("%s  (nominal %.1f GeV)", kTrgName[tr], kTrgThresholdGeV[tr]));
        for (int k = 0; k < 2; k++)
        {
            TLine *l = new TLine(kPeriodEdge[k], 0.0, kPeriodEdge[k], 1.6);
            l->SetLineColor(kGray + 2); l->SetLineStyle(7); l->SetLineWidth(2); l->Draw();
        }
        ge[tr]->SetMarkerStyle(20); ge[tr]->SetMarkerColor(kAzure + 2); ge[tr]->SetLineColor(kAzure + 2);
        gp[tr]->SetMarkerStyle(21); gp[tr]->SetMarkerColor(kOrange + 7); gp[tr]->SetLineColor(kOrange + 7);
        if (ge[tr]->GetN()) ge[tr]->Draw("P");
        if (gp[tr]->GetN()) gp[tr]->Draw("P");
        if (gx[tr]->GetN()) { gx[tr]->SetMarkerStyle(24); gx[tr]->SetMarkerColor(kBlack); gx[tr]->SetLineColor(kBlack); gx[tr]->Draw("P"); }
        if (tr == 1)
        {
            TLegend *lg = new TLegend(0.14, 0.14, 0.80, 0.33);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
            lg->AddEntry(ge[tr], "50% efficiency (cut keeps half the fired events)", "p");
            lg->AddEntry(gp[tr], "50% purity (half the cut's events fired)", "p");
            if (gx[tr]->GetN()) lg->AddEntry(gx[tr], "turn-on x_{50} (bit set at that value 50%)", "p");
            lg->Draw();
        }
    }
    c.SaveAs((outdir + "/trg_time_scan.png").c_str());
    printf("trg_time_scan: %zu groups; %s/trg_time_scan.{txt,png}\n", groups.size(), outdir.c_str());
    return 0;
}
