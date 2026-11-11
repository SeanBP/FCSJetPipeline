// Efficiency and purity of the offline trigger cuts against the hardware
// trigger bits (Sean 2026-11-11), from unfiltered jet files (the Run 22
// trigger scan, data_to_jet -t "" -E 1 -F 1, one MuDst file per run).
//
// Per event, the offline decision for trigger t: some jet has a trigger value
// >= f x the nominal threshold (AnConfig.h kTrgThresholdGeV). EM: the jet's
// largest ECal 4x4 window (reco_trg_em_et, the stage-2 EM sum) and, for the
// HCal veto, the largest neighbouring HCal 4x4 (reco_trg_em_hmax): the jet
// passes the veto with margin m if 4 m hmax <= em (the hardware vetoes at
// 4 hmax > em, i.e. m = 1; m = 0 means no veto cut). JP: the jet's patch sum
// (reco_trg_jp_et, ECal + HCal E_T of the large patch).
// The hardware decision: the trigger's bit in fcs_tcu_bits (stage 3, before
// prescale): EM0-2 bits 3-5, JP2 bit 6, JP1 any of bits 7-9, JP0 any of bits
// 10-12 (event level: any region).
//
//   efficiency = N(offline and hardware) / N(hardware with a jet): how many
//                events the hardware fired the cut keeps, among those with a
//                jet whose trigger value (any veto) is >= 0.5 x threshold
//                (hardware events without such a jet, the firing window
//                outside any reconstructed jet: fiducial cut, hot trigger
//                regions, are the "coverage" loss, printed separately)
//   purity     = N(offline and hardware) / N(offline): how many events the
//                cut keeps the hardware also fired
//
// Only reference events count (as trg_eff_scan.cpp): recorded by a trigger
// that cannot imply the measured family's bit (EM: JP, DiJP, Had, HHT, HTOT,
// HcalTot; JP: EM0/1, DY, J/psi, ELE2), so the bits are not biased by the
// trigger that chose the event.
// Curves per gain period (before 23027047, to 23066056, after: the offline-
// to-online scale differs) and for all runs.
//
// Output (<outdir>): trg_cut_perf_EM.png, trg_cut_perf_JP.png (efficiency and
// purity vs f, the current cut marked), trg_cut_perf_veto.png (EM vs the veto
// margin m at f = kEmCutFactor), trg_cut_perf.txt (the numbers at the current
// cuts).
// Usage: ./trg_cut_perf <jet file list> <outdir>
// Build: g++ -m64 -O2 -I../shared -I../asymmetry trg_cut_perf.cpp ../shared/VectorDict.cxx -o trg_cut_perf
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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "VectorDict.h"
#include "AnConfig.h"

static const int kNT = 6;  // EM0, EM1, EM2, JP0, JP1, JP2
static const int kNP = 4;  // gain periods 0-2, all
static const int kPeriodEdge[2] = {23027047, 23066056};
static int PeriodOf(int run) { return run < kPeriodEdge[0] ? 0 : run < kPeriodEdge[1] ? 1 : 2; }
static const char *kPeriodName[kNP] = {"before 23027047", "23027047-23066055", "from 23066056", "all runs"};
static const int kPeriodCol[kNP] = {kRed + 1, kGreen + 2, kBlue + 1, kBlack};
// cut factor scan f = trigger value / threshold
static const int kNF = 41;
static double FacOf(int i) { return 0.5 + 0.05 * i; }  // 0.5 ... 2.5
// veto margin scan m (0 = no veto cut)
static const int kNM = 31;
static double MarOf(int i) { return 0.1 * i; }          // 0 ... 3.0

struct Tally { long hw = 0, hwj = 0, off = 0, both = 0; };

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <jet file list> <outdir>" << std::endl; return 1; }
    const std::string outdir = argv[2];
    gSystem->mkdir(outdir.c_str(), true);
    // [period][trigger][factor] and, EM only, [period][level][margin]
    std::vector<Tally> tf(kNP * kNT * kNF), tm(kNP * 3 * kNM);
    long nref[kNP][2] = {{0}};

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
            const int per = PeriodOf(run);
            bool em_ref = false, jp_ref = false;
            for (int k = 0; k < n_ref_em; k++) em_ref |= ref_em[k] != 0;
            for (int k = 0; k < n_ref_jp; k++) jp_ref |= ref_jp[k] != 0;
            if (em_ref) { nref[per][0]++; nref[3][0]++; }
            if (jp_ref) { nref[per][1]++; nref[3][1]++; }
            for (int tr = 0; tr < kNT; tr++)
            {
                const bool isjp = tr >= 3;
                if (!(isjp ? jp_ref : em_ref)) continue;
                bool hw;
                if (tr < 3) hw = (tcu >> (3 + tr)) & 1;
                else if (tr == 5) hw = (tcu >> 6) & 1;
                else hw = ((tcu >> (tr == 3 ? 10 : 7)) & 7) != 0;  // any region
                const double thr = kTrgThresholdGeV[tr];
                // the largest trigger value of any jet, no veto: the coverage
                double anyv = -1;
                for (size_t j = 0; j < em->size(); j++) anyv = std::max(anyv, (double)(isjp ? (*jp)[j] : (*em)[j]));
                const bool hwj = hw && anyv >= 0.5 * thr;
                // factor scan, EM at the current veto margin
                double lead = -1;
                for (size_t j = 0; j < em->size(); j++)
                {
                    double v;
                    if (isjp) v = (*jp)[j];
                    else v = ((*em)[j] > 0 && (*hmax)[j] >= 0 && 4 * kEmVetoMargin * (*hmax)[j] <= (*em)[j]) ? (*em)[j] : -1;
                    if (v > lead) lead = v;
                }
                for (int fi = 0; fi < kNF; fi++)
                {
                    const bool off = lead >= FacOf(fi) * thr;
                    for (int p : {per, 3})
                    {
                        Tally &c = tf[(p * kNT + tr) * kNF + fi];
                        c.hw += hw; c.hwj += hwj; c.off += off; c.both += hw && off;
                    }
                }
                // veto-margin scan, EM at the current factor
                if (!isjp)
                    for (int mi = 0; mi < kNM; mi++)
                    {
                        const double m = MarOf(mi);
                        bool off = false;
                        for (size_t j = 0; j < em->size() && !off; j++)
                            off = (*em)[j] >= kEmCutFactor * thr && (*hmax)[j] >= 0 && 4 * m * (*hmax)[j] <= (*em)[j];
                        for (int p : {per, 3})
                        {
                            Tally &c = tm[(p * 3 + tr) * kNM + mi];
                            c.hw += hw; c.hwj += hwj; c.off += off; c.both += hw && off;
                        }
                    }
            }
        }
    }
    printf("trg_cut_perf: %d files; reference events EM %ld, JP %ld\n", nfiles, nref[3][0], nref[3][1]);

    gStyle->SetOptStat(0);
    // efficiency (row 0) or purity (row 1) graph with Clopper-Pearson errors
    auto make = [](const std::vector<Tally> &v, int base, int n, double (*xof)(int), int row) {
        TGraphAsymmErrors *g = new TGraphAsymmErrors();
        for (int i = 0; i < n; i++)
        {
            const Tally &c = v[base + i];
            const long den = row == 0 ? c.hwj : c.off;
            if (den <= 0) continue;
            const double e = double(c.both) / den;
            const int k = g->GetN();
            g->SetPoint(k, xof(i), e);
            g->SetPointError(k, 0, 0, e - TEfficiency::ClopperPearson(den, c.both, 0.683, false),
                             TEfficiency::ClopperPearson(den, c.both, 0.683, true) - e);
        }
        return g;
    };
    std::ofstream txt(outdir + "/trg_cut_perf.txt");
    txt << "# offline trigger cut vs hardware bit, reference events; eff = N(both)/N(hw with a jet >= 0.5 thr), pur = N(both)/N(offline);\n"
        << "# coverage = N(hw with a jet >= 0.5 thr)/N(hw)\n";
    // factor pages: EM (tr 0-2) and JP (tr 3-5)
    for (int fam = 0; fam < 2; fam++)
    {
        TCanvas c("c", "", 1500, 850);
        c.Divide(3, 2, 0.002, 0.002);
        const double fcur = fam ? kJpCutFactor : kEmCutFactor;
        for (int k = 0; k < 3; k++)
        {
            const int tr = 3 * fam + k;
            for (int row = 0; row < 2; row++)
            {
                c.cd(1 + k + 3 * row);
                gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
                TH1D *fr = new TH1D(Form("fr%d%d%d", fam, k, row), Form(";offline %s / threshold (%g GeV);%s", fam ? "JP patch sum" : "ECal 4x4 sum",
                                    kTrgThresholdGeV[tr], row ? "purity" : "efficiency"), 10, 0.5, 2.5);
                fr->SetMinimum(0); fr->SetMaximum(1.05);
                fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleOffset(1.1);
                fr->Draw();
                TLegend *lg = row == 1 && k == 0 ? new TLegend(0.45, 0.15, 0.93, 0.45) : 0;
                for (int p = 0; p < kNP; p++)
                {
                    TGraphAsymmErrors *g = make(tf, (p * kNT + tr) * kNF, kNF, FacOf, row);
                    g->SetLineColor(kPeriodCol[p]); g->SetMarkerColor(kPeriodCol[p]);
                    g->SetMarkerStyle(p == 3 ? 20 : 24); g->SetMarkerSize(0.6); g->SetLineWidth(p == 3 ? 2 : 1);
                    g->Draw("PL SAME");
                    if (lg) lg->AddEntry(g, kPeriodName[p], "pl");
                }
                TLine *l = new TLine(fcur, 0, fcur, 1.05); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
                if (lg) { lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04); lg->Draw(); }
                TLatex tl; tl.SetNDC(); tl.SetTextSize(0.06);
                tl.DrawLatex(0.15, 0.92, Form("#bf{%s} %s", kTrgLevelName[tr], row ? "purity" : "efficiency"));
                tl.SetTextSize(0.045);
                if (row == 0)
                {
                    const Tally &ca = tf[(3 * kNT + tr) * kNF];
                    tl.DrawLatex(0.5, 0.75, Form("coverage (all runs): %.2f", ca.hw ? double(ca.hwj) / ca.hw : 0.));
                }
                if (!fam && row == 0) tl.DrawLatex(0.5, 0.82, Form("HCal veto margin %.2g", kEmVetoMargin));
            }
            for (int p = 0; p < kNP; p++)
            {
                int fi = (int)std::lround((fcur - 0.5) / 0.05);
                const Tally &cc = tf[(p * kNT + tr) * kNF + fi];
                txt << Form("%s f=%.2f %-18s hw %8ld off %8ld both %8ld eff %.3f pur %.3f", kTrgLevelName[tr], fcur, kPeriodName[p], cc.hw, cc.off, cc.both,
                            cc.hwj ? double(cc.both) / cc.hwj : 0., cc.off ? double(cc.both) / cc.off : 0.)
                    << Form("   coverage %.3f\n", cc.hw ? double(cc.hwj) / cc.hw : 0.);
            }
        }
        c.SaveAs((outdir + (fam ? "/trg_cut_perf_JP.png" : "/trg_cut_perf_EM.png")).c_str());
    }
    // veto-margin page
    {
        TCanvas c("cv", "", 1500, 850);
        c.Divide(3, 2, 0.002, 0.002);
        for (int k = 0; k < 3; k++)
            for (int row = 0; row < 2; row++)
            {
                c.cd(1 + k + 3 * row);
                gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
                TH1D *fr = new TH1D(Form("frv%d%d", k, row), Form(";HCal veto margin m (passes: 4m H_{4x4} #leq E_{4x4}, m = 0: none);%s",
                                    row ? "purity" : "efficiency"), 10, 0, 3);
                fr->SetMinimum(0); fr->SetMaximum(1.05);
                fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleOffset(1.1);
                fr->Draw();
                TLegend *lg = row == 1 && k == 0 ? new TLegend(0.45, 0.15, 0.93, 0.45) : 0;
                for (int p = 0; p < kNP; p++)
                {
                    TGraphAsymmErrors *g = make(tm, (p * 3 + k) * kNM, kNM, MarOf, row);
                    g->SetLineColor(kPeriodCol[p]); g->SetMarkerColor(kPeriodCol[p]);
                    g->SetMarkerStyle(p == 3 ? 20 : 24); g->SetMarkerSize(0.6); g->SetLineWidth(p == 3 ? 2 : 1);
                    g->Draw("PL SAME");
                    if (lg) lg->AddEntry(g, kPeriodName[p], "pl");
                }
                TLine *l = new TLine(kEmVetoMargin, 0, kEmVetoMargin, 1.05); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
                TLine *l1 = new TLine(1, 0, 1, 1.05); l1->SetLineStyle(3); l1->SetLineColor(kGray + 2); l1->Draw();
                if (lg) { lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04); lg->Draw(); }
                TLatex tl; tl.SetNDC(); tl.SetTextSize(0.06);
                tl.DrawLatex(0.15, 0.92, Form("#bf{%s} %s vs veto margin", kTrgLevelName[k], row ? "purity" : "efficiency"));
                tl.SetTextSize(0.045);
                if (row == 0) tl.DrawLatex(0.45, 0.82, Form("ECal 4x4 #geq %.2g #times threshold", kEmCutFactor));
            }
        c.SaveAs((outdir + "/trg_cut_perf_veto.png").c_str());
        for (int k = 0; k < 3; k++)
            for (int p = 0; p < kNP; p++)
                for (double m : {0.0, 1.0, kEmVetoMargin, 2.0})
                {
                    const Tally &cc = tm[(p * 3 + k) * kNM + (int)std::lround(m / 0.1)];
                    txt << Form("%s veto m=%.1f (f=%.2f) %-18s hw %8ld off %8ld both %8ld eff %.3f pur %.3f\n", kTrgLevelName[k], m, kEmCutFactor, kPeriodName[p],
                                cc.hw, cc.off, cc.both, cc.hwj ? double(cc.both) / cc.hwj : 0., cc.off ? double(cc.both) / cc.off : 0.);
                }
    }
    printf("trg_cut_perf: wrote %s/trg_cut_perf_{EM,JP,veto}.png and trg_cut_perf.txt\n", outdir.c_str());
    return 0;
}
