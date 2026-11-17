// HCal-veto edge of the FCS EM triggers over Run 22 (Sean 2026-11-17: the EM
// threshold and the veto are one cut, moved together by one factor f; this
// measures the veto's own "threshold"). Reference events for EM (as
// trg_cut_perf.cpp: JP, DiJP, Had, HHT, HTOT, HcalTot) whose jet with the
// largest ECal 4x4 has E_4x4 >= 1.3 x the period's measured EM0 threshold
// (AnConfig.h TrgThresholdGeV, on the turn-on plateau); its veto ratio
// r = 4 H_max / E_4x4 (the hardware vetoes at r > 1 nominally). The EM0 bit
// fraction vs r falls from the plateau to ~0 across the veto; fitted with
//   a + (p - a) / 2 (1 - erf((r - r50) / (sqrt2 s)))
// r50 = the veto edge (in nominal veto-ratio units), per gain period and per
// group of consecutive runs (as trg_time_scan.cpp).
// Output: <outdir>/trg_veto_edge.txt (per period and per group: r50 +- err,
// s, n events) and trg_veto_edge.png (bit fraction vs r per period with fits;
// r50 vs run number).
// Usage: ./trg_veto_edge <jet file list> <outdir> [runs per group = 10]
// Build: g++ -m64 -O2 -I../shared -I../asymmetry trg_veto_edge.cpp ../shared/VectorDict.cxx -o trg_veto_edge
//            `root-config --cflags --libs` -std=c++11
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
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
#include <string>
#include <vector>

#include "VectorDict.h"
#include "AnConfig.h"

static const int kPeriodEdge[2] = {23027047, 23066056};
static int PeriodOf(int run) { return run < kPeriodEdge[0] ? 0 : run < kPeriodEdge[1] ? 1 : 2; }
static const int kNR = 40;                                // r bins 0 ... 2, 0.05 wide
static int RBin(double r) { return r < 0 ? -1 : r >= 2 ? kNR - 1 : (int)(r / 0.05); }
struct Bins { long n[kNR] = {0}, k[kNR] = {0}; };

// fit the falling edge; returns false when the fit fails
static bool FitEdge(const Bins &b, const char *name, double &r50, double &er, double &s, TGraphAsymmErrors *&g, TF1 *&f)
{
    g = new TGraphAsymmErrors();
    long tot = 0;
    for (int i = 0; i < kNR; i++)
    {
        if (b.n[i] < 3) continue;
        tot += b.n[i];
        const double e = double(b.k[i]) / b.n[i], x = 0.05 * (i + 0.5);
        const int m = g->GetN();
        g->SetPoint(m, x, e);
        g->SetPointError(m, 0, 0, e - TEfficiency::ClopperPearson(b.n[i], b.k[i], 0.683, false), TEfficiency::ClopperPearson(b.n[i], b.k[i], 0.683, true) - e);
    }
    if (g->GetN() < 5 || tot < 50) return false;
    f = new TF1(name, "[0] + ([1] - [0]) / 2 * (1 - TMath::Erf((x - [2]) / (sqrt(2) * [3])))", 0, 2);
    f->SetParameters(0.1, 1.0, 0.9, 0.1);
    f->SetParLimits(0, 0, 0.5); f->SetParLimits(1, 0.5, 1.0); f->SetParLimits(2, 0.3, 1.6); f->SetParLimits(3, 0.01, 0.5);
    if (g->Fit(f, "QN0") != 0) return false;
    r50 = f->GetParameter(2); er = f->GetParError(2); s = f->GetParameter(3);
    return er > 0;
}

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <jet file list> <outdir> [runs per group = 10]" << std::endl; return 1; }
    const std::string outdir = argv[2];
    const int per_group = argc > 3 ? std::max(1, std::atoi(argv[3])) : 10;
    gSystem->mkdir(outdir.c_str(), true);
    std::map<int, Bins> runs;
    Bins per[3];

    std::ifstream list(argv[1]);
    std::string fn;
    while (list >> fn)
    {
        TFile f(fn.c_str());
        TTree *t = f.IsZombie() ? 0 : (TTree *)f.Get("jetTree");
        if (!t || !t->GetBranch("reco_trg_em_et") || !t->GetBranch("fcs_tcu_bits")) continue;
        int file_run = 0, seq = 0;
        const std::string base = fn.substr(fn.find_last_of('/') + 1);
        if (std::sscanf(base.c_str(), "jet_output_%d_%d.root", &file_run, &seq) != 2) file_run = 0;
        Int_t ev_run = 0;
        if (t->GetBranch("event_run")) t->SetBranchAddress("event_run", &ev_run);
        std::vector<float> *em = 0, *hmax = 0;
        UInt_t tcu = 0;
        const char *ref_em_names[] = {"fcsJP2", "fcsJPA1", "fcsJPBC1", "fcsJPDE1", "fcsJPA0", "fcsJPBC0", "fcsJPDE0", "fcsDiJP", "fcsDiJPAsy",
                                      "fcsHad0", "fcsHad1", "fcsHad2", "fcsHHT_N", "fcsHHT_S", "fcsHTOT_N", "fcsHTOT_S", "fcs_HcalTot"};
        const int n_ref = sizeof(ref_em_names) / sizeof(ref_em_names[0]);
        std::vector<Int_t> ref(n_ref, 0);
        for (int k = 0; k < n_ref; k++) if (t->GetBranch(ref_em_names[k])) t->SetBranchAddress(ref_em_names[k], &ref[k]);
        t->SetBranchAddress("reco_trg_em_et", &em);
        t->SetBranchAddress("reco_trg_em_hmax", &hmax);
        t->SetBranchAddress("fcs_tcu_bits", &tcu);
        for (Long64_t i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);
            const int run = ev_run > 0 ? ev_run : file_run;
            if (run <= 0) continue;
            bool isref = false;
            for (int k = 0; k < n_ref; k++) isref |= ref[k] != 0;
            if (!isref) continue;
            double aem = -1, ar = -1;
            for (size_t j = 0; j < em->size(); j++)
                if ((*em)[j] > aem && (*hmax)[j] >= 0) { aem = (*em)[j]; ar = 4 * (*hmax)[j] / (*em)[j]; }
            if (aem < 1.3 * TrgThresholdGeV(0, run)) continue;
            const int rb = RBin(ar);
            if (rb < 0) continue;
            const bool bit = (tcu >> 3) & 1;  // EM0
            Bins &rr = runs[run];
            rr.n[rb]++; rr.k[rb] += bit;
            per[PeriodOf(run)].n[rb]++; per[PeriodOf(run)].k[rb] += bit;
        }
    }

    std::ofstream txt(outdir + "/trg_veto_edge.txt");
    txt << "# EM0 HCal-veto edge: bit fraction vs r = 4 H_max / E_4x4 (reference events, leading-ECal jet >= 1.3 x measured EM0 threshold),\n"
        << "# fit a + (p - a)/2 (1 - erf((r - r50)/(sqrt2 s))); r50 in nominal veto-ratio units\n"
        << "# kind first_run last_run r50 err s\n";
    gStyle->SetOptStat(0);
    TCanvas c("c", "c", 1600, 650);
    c.Divide(2, 1, 0.004, 0.004);
    c.cd(1); gPad->SetGridy(); gPad->SetLeftMargin(0.12);
    TH1F *fr = new TH1F("fr", ";veto ratio r = 4 H_{4x4}^{max} / E_{4x4} (nominal veto at 1);fraction with the EM0 bit set", 10, 0, 2);
    fr->SetMinimum(0); fr->SetMaximum(1.2); fr->Draw("AXIS");
    TLegend *lg = new TLegend(0.45, 0.70, 0.95, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0);
    const int col[3] = {kAzure + 2, kOrange + 7, kBlack};
    const char *pname[3] = {"before 23027047", "23027047-23066055", "from 23066056"};
    double pr50[3] = {0, 0, 0};
    for (int p = 0; p < 3; p++)
    {
        double r50, er, s; TGraphAsymmErrors *g; TF1 *f = 0;
        const bool ok = FitEdge(per[p], Form("fp%d", p), r50, er, s, g, f);
        g->SetMarkerStyle(20 + p); g->SetMarkerColor(col[p]); g->SetLineColor(col[p]); g->Draw("P");
        if (ok) { f->SetLineColor(col[p]); f->Draw("same"); pr50[p] = r50; }
        lg->AddEntry(g, ok ? Form("%s: r_{50} = %.3f #pm %.3f", pname[p], r50, er) : pname[p], "p");
        char line[200];
        std::snprintf(line, sizeof line, "period%d %d %d %.4f %.4f %.4f\n", p, p ? kPeriodEdge[p - 1] : 0, p < 2 ? kPeriodEdge[p] - 1 : 99999999, ok ? r50 : -1, ok ? er : -1, ok ? s : -1);
        txt << line;
    }
    lg->Draw();

    // groups of consecutive runs within a period
    std::vector<std::vector<int> > groups;
    {
        std::vector<int> cur; int curp = -1;
        for (std::map<int, Bins>::const_iterator it = runs.begin(); it != runs.end(); ++it)
        {
            const int p = PeriodOf(it->first);
            if (p != curp && !cur.empty()) { groups.push_back(cur); cur.clear(); }
            curp = p; cur.push_back(it->first);
            if ((int)cur.size() == per_group) { groups.push_back(cur); cur.clear(); }
        }
        if (!cur.empty()) groups.push_back(cur);
    }
    TGraphAsymmErrors *gt = new TGraphAsymmErrors();
    for (size_t g = 0; g < groups.size(); g++)
    {
        Bins b;
        for (size_t r = 0; r < groups[g].size(); r++) { const Bins &x = runs[groups[g][r]]; for (int i = 0; i < kNR; i++) { b.n[i] += x.n[i]; b.k[i] += x.k[i]; } }
        double r50, er, s; TGraphAsymmErrors *gg; TF1 *f = 0;
        const bool ok = FitEdge(b, Form("fg%zu", g), r50, er, s, gg, f);
        const int r0 = groups[g].front(), r1 = groups[g].back();
        char line[200];
        std::snprintf(line, sizeof line, "group %d %d %.4f %.4f %.4f\n", r0, r1, ok ? r50 : -1, ok ? er : -1, ok ? s : -1);
        txt << line;
        if (ok) { const int n = gt->GetN(); const double xc = 0.5 * (r0 + r1); gt->SetPoint(n, xc, r50); gt->SetPointError(n, xc - r0, r1 - xc, er, er); }
    }
    txt.close();
    c.cd(2); gPad->SetGridy(); gPad->SetLeftMargin(0.12);
    const double rlo = runs.empty() ? 23000000 : runs.begin()->first - 2000, rhi = runs.empty() ? 23110000 : runs.rbegin()->first + 2000;
    TH1F *f2 = new TH1F("f2", ";run number;veto edge r_{50} (nominal veto ratio units)", 10, rlo, rhi);
    f2->SetMinimum(0.4); f2->SetMaximum(1.4); f2->GetXaxis()->SetNoExponent(); f2->GetXaxis()->SetNdivisions(505); f2->Draw("AXIS");
    for (int k = 0; k < 2; k++) { TLine *l = new TLine(kPeriodEdge[k], 0.4, kPeriodEdge[k], 1.4); l->SetLineStyle(7); l->SetLineColor(kGray + 2); l->Draw(); }
    for (int p = 0; p < 3; p++)
        if (pr50[p] > 0)
        {
            const double a = p ? kPeriodEdge[p - 1] : rlo, b = p < 2 ? kPeriodEdge[p] : rhi;
            TLine *l = new TLine(a, pr50[p], b, pr50[p]); l->SetLineColor(col[p]); l->SetLineWidth(2); l->Draw();
        }
    gt->SetMarkerStyle(20); gt->Draw("P");
    c.cd(0);
    TLatex tl; tl.SetNDC(); tl.SetTextFont(62); tl.SetTextSize(0.03);
    tl.DrawLatex(0.01, 0.965, "EM0 HCal-veto edge (reference events, leading ECal 4x4 #geq 1.3 #times measured threshold)");
    c.SaveAs((outdir + "/trg_veto_edge.png").c_str());
    printf("trg_veto_edge: periods r50 %.3f / %.3f / %.3f; %zu groups; %s/trg_veto_edge.{txt,png}\n", pr50[0], pr50[1], pr50[2], groups.size(), outdir.c_str());
    return 0;
}
