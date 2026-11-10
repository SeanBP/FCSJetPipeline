// Purity / efficiency of the leading-point pi0 probability P_lead as a pi0-jet
// selector, and confusion matrices of the leading point's truth origin, on the
// validation half (odd file_idx) of truth-labelled point-jet dumps
// (pi0_jet_dump), cross-section weighted (Pi0DumpWeights.h).
//
// Truth class of the jet's leading (highest-energy) point:
//   pi0      matched to a photon whose parent is a pi0 (pt_pi0 >= 0)
//   eta      matched to a photon from an eta (pt_par 17)
//   other g  matched to another photon (prompt, omega, ... ; pt_par >= 0)
//   hadron   no truth photon within the match distance, and the nearest
//            particle above 0.3 GeV created before the ECal (pt_ngeid) is a
//            hadron (Geant id >= 8)
//   other    unmatched, nearest particle an e+-, a photon or none
// "pi0 jet" = leading point from a pi0, the class P_lead estimates.
//
// Pages (also PNGs <prefix>_p<N>.png):
//   1  left: purity vs efficiency of the cut P_lead >= t (t scanned) per raw
//      jet-energy bin, with the photon-count cuts n <= 1, 2, 3 as markers;
//      right: purity and efficiency vs the threshold t (all A_N-range jets)
//   3  Collins pair selection: purity / efficiency of P_leadpair >= t for "the
//      chosen pair (leading point + closest-mass partner, Pi0LeadPair) is the
//      leading photon's true pi0 pair", per raw-E bin, with n == 2 as a marker;
//      printed: how often the closest-mass partner is right when the true pair
//      is there
//   2  confusion matrices, A_N-range jets (raw E 28-100 GeV): truth class vs
//      P_lead bin and truth class vs points in the jet; each column normalised
//      (the composition of the selected jets), weighted
//
// Build: g++ -O2 pi0_prob_perf.cpp -o pi0_prob_perf -I../../shared -I../../jet_scale -I. `root-config --cflags --libs` -std=c++11
// Run:   ./pi0_prob_perf "<dump glob>" <sim calib.json> <prefix>
#include "Pi0Prob.h"
#include "Pi0DumpWeights.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace std;

enum { kPi0 = 0, kEta, kOtherG, kHadron, kOther, kNTruth };
static const char *kTruthName[kNTruth] = {"#pi^{0}", "#eta", "other #gamma", "hadron", "other"};
static const char *kTruthText[kNTruth] = {"pi0", "eta", "other g", "hadron", "other"};

struct J { float p, w, e; int n, cls; float pp; int pairOK, avail; };

int main(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "Usage: %s \"<dump glob>\" <sim calib.json> <prefix>\n", argv[0]); return 1; }
    gROOT->SetBatch(1);
    gStyle->SetOptStat(0);
    const string prefix = argv[3];
    Pi0ProbCalib cal;
    string err;
    if (!cal.Load(argv[2], false, &err)) { fprintf(stderr, "%s\n", err.c_str()); return 1; }
    map<int, double> wOf;
    if (!Pi0DumpXSecWeights(argv[1], 0, wOf)) return 1;

    TChain ch("jets");
    ch.Add(argv[1]);
    const int MX = 200;
    Int_t fidx, n, ppi0[MX], ppar[MX], pngeid[MX];
    Float_t jetE, pE[MX], px[MX], py[MX], pz[MX], pdist[MX];
    ch.SetBranchStatus("*", 0);
    for (auto b : {"file_idx", "n", "jet_E", "pt_E", "pt_x", "pt_y", "pt_z", "pt_pi0", "pt_par", "pt_dist", "pt_ngeid"}) ch.SetBranchStatus(b, 1);
    ch.SetBranchAddress("file_idx", &fidx); ch.SetBranchAddress("n", &n); ch.SetBranchAddress("jet_E", &jetE);
    ch.SetBranchAddress("pt_E", pE); ch.SetBranchAddress("pt_x", px); ch.SetBranchAddress("pt_y", py); ch.SetBranchAddress("pt_z", pz);
    ch.SetBranchAddress("pt_pi0", ppi0); ch.SetBranchAddress("pt_par", ppar); ch.SetBranchAddress("pt_dist", pdist);
    ch.SetBranchAddress("pt_ngeid", pngeid);

    vector<J> jets;
    for (Long64_t i = 0; i < ch.GetEntries(); i++)
    {
        ch.GetEntry(i);
        if (fidx % 2 == 0 || n < 1 || n > MX) continue;
        auto it = wOf.find(fidx);
        if (it == wOf.end() || !(it->second > 0)) continue;
        vector<double> e(n), qx(n), qy(n), qz(n);
        int L = 0;
        for (int k = 0; k < n; k++)
        {
            const double r = sqrt(px[k] * px[k] + py[k] * py[k] + pz[k] * pz[k]);
            e[k] = pE[k]; qx[k] = pE[k] * px[k] / r; qy[k] = pE[k] * py[k] / r; qz[k] = pE[k] * pz[k] / r;
            if (pE[k] > pE[L]) L = k;
        }
        const Pi0ProbResult res = Pi0Prob(cal, e, qx, qy, qz, jetE);
        if (res.lead < 0) continue;
        int cls;
        if (ppi0[L] >= 0) cls = kPi0;
        else if (ppar[L] == 17) cls = kEta;
        else if (ppar[L] >= 0) cls = kOtherG;
        else if (pngeid[L] >= 8) cls = kHadron;
        else cls = kOther;
        // the chosen pair is the leading photon's true pi0 pair; and whether that pair is in the jet at all
        int avail = 0;
        if (ppi0[L] >= 0) for (int k = 0; k < n; k++) if (k != L && ppi0[k] == ppi0[L]) avail = 1;
        const int pj = res.lead_partner;
        const int ok = (ppi0[L] >= 0 && pj >= 0 && ppi0[pj] == ppi0[L]) ? 1 : 0;
        jets.push_back({res.lead, (float)it->second, jetE, n, cls, res.lead_pair, ok, avail});
    }
    printf("validation jets (odd file_idx, n >= 1, P_lead computed): %zu\n", jets.size());

    // raw-E bins (about the corrected A_N ranges 35-50, 50-70, 70-100, 100-125)
    const int NE = 5;
    const double elo[NE] = {28, 40, 56, 80, 28}, ehi[NE] = {40, 56, 80, 100, 100};
    const int col[NE] = {kBlue + 1, kGreen + 2, kOrange + 7, kRed + 1, kBlack};
    const int NT = 101;
    TCanvas cv("cv", "", 1300, 600);
    const string pdf = prefix + ".pdf";
    cv.Print((pdf + "[").c_str());

    // ---- page 1: purity vs efficiency ----
    cv.Divide(2, 1);
    vector<TGraph *> roc(NE), pur(NE), eff(NE);
    vector<TGraph *> mk(NE);
    printf("\nraw E bin   pi0 frac  n_eff | cut         eff    purity | P_lead at that eff: purity\n");
    for (int b = 0; b < NE; b++)
    {
        double wA = 0, wA2 = 0, wS = 0;
        for (auto &j : jets) if (j.e >= elo[b] && j.e < ehi[b]) { wA += j.w; wA2 += j.w * j.w; if (j.cls == kPi0) wS += j.w; }
        roc[b] = new TGraph(); pur[b] = new TGraph(); eff[b] = new TGraph(); mk[b] = new TGraph();
        if (!(wA > 0)) continue;
        for (int t = 0; t < NT; t++)
        {
            const double thr = t / (double)(NT - 1);
            double sel = 0, selS = 0;
            for (auto &j : jets) if (j.e >= elo[b] && j.e < ehi[b] && j.p >= thr) { sel += j.w; if (j.cls == kPi0) selS += j.w; }
            if (sel > 0) { roc[b]->SetPoint(roc[b]->GetN(), selS / wS, selS / sel); pur[b]->SetPoint(pur[b]->GetN(), thr, selS / sel); }
            eff[b]->SetPoint(eff[b]->GetN(), thr, selS / wS);
        }
        printf("%3.0f-%-4.0f  %.3f  %6.0f |", elo[b], ehi[b], wS / wA, wA * wA / wA2);
        for (int nc = 1; nc <= 3; nc++)
        {
            double sel = 0, selS = 0;
            for (auto &j : jets) if (j.e >= elo[b] && j.e < ehi[b] && j.n <= nc) { sel += j.w; if (j.cls == kPi0) selS += j.w; }
            if (!(sel > 0)) continue;
            const double ef = selS / wS, pu = selS / sel;
            mk[b]->SetPoint(mk[b]->GetN(), ef, pu);
            // P_lead purity at the same efficiency (interpolated on the scan)
            double pp = -1;
            for (int q = 1; q < roc[b]->GetN(); q++)
            {
                double x0, y0, x1, y1;
                roc[b]->GetPoint(q - 1, x0, y0); roc[b]->GetPoint(q, x1, y1);
                if ((x0 - ef) * (x1 - ef) <= 0 && x0 != x1) { pp = y0 + (y1 - y0) * (ef - x0) / (x1 - x0); break; }
            }
            printf("%s n<=%d %.3f %.3f -> %.3f", nc > 1 ? "\n                           |" : " ", nc, ef, pu, pp);
        }
        printf("\n");
    }
    cv.cd(1);
    gPad->SetGrid();
    TH1F *fr = gPad->DrawFrame(0, 0.5, 1, 1.0, "leading point from a #pi^{0}: cut P_{lead} #geq t;efficiency;purity");
    (void)fr;
    TLegend lg(0.13, 0.13, 0.62, 0.42);
    lg.SetTextSize(0.033);
    for (int b = 0; b < NE; b++)
    {
        roc[b]->SetLineColor(col[b]); roc[b]->SetLineWidth(b == NE - 1 ? 3 : 2); roc[b]->Draw("L");
        mk[b]->SetMarkerColor(col[b]); mk[b]->SetMarkerStyle(b == NE - 1 ? 29 : 24); mk[b]->SetMarkerSize(b == NE - 1 ? 2.0 : 1.3); mk[b]->Draw("P");
        lg.AddEntry(roc[b], Form("raw E %.0f-%.0f GeV%s", elo[b], ehi[b], b == NE - 1 ? " (all)" : ""), "l");
    }
    lg.AddEntry(mk[NE - 1], "photon-count cuts n #leq 1, 2, 3", "p");
    lg.Draw();
    cv.cd(2);
    gPad->SetGrid();
    gPad->DrawFrame(0, 0, 1, 1.0, "raw E 28-100 GeV;threshold t on P_{lead};fraction");
    pur[NE - 1]->SetLineColor(kRed + 1); pur[NE - 1]->SetLineWidth(3); pur[NE - 1]->Draw("L");
    eff[NE - 1]->SetLineColor(kBlue + 1); eff[NE - 1]->SetLineWidth(3); eff[NE - 1]->Draw("L");
    TLegend lg2(0.15, 0.15, 0.5, 0.3);
    lg2.AddEntry(pur[NE - 1], "purity", "l"); lg2.AddEntry(eff[NE - 1], "efficiency", "l");
    lg2.Draw();
    cv.Print(pdf.c_str());
    cv.Print((prefix + "_p1.png").c_str());

    // ---- page 2: confusion matrices (A_N-range jets) ----
    cv.Clear();
    cv.Divide(2, 1);
    const int NP = 4;
    const double pe[NP + 1] = {0, 0.5, 0.7, 0.85, 1.0001};
    const char *pl[NP] = {"<0.5", "0.5-0.7", "0.7-0.85", "#geq0.85"};
    TH2D *cmP = new TH2D("cmP", "truth of the leading point vs P_{lead} (raw E 28-100 GeV);P_{lead};", NP, 0, NP, kNTruth, 0, kNTruth);
    TH2D *cmN = new TH2D("cmN", "truth of the leading point vs points in the jet (raw E 28-100 GeV);points n;", 5, 0, 5, kNTruth, 0, kNTruth);
    const char *nl[5] = {"1", "2", "3", "4", "#geq5"};
    double colP[NP] = {0}, colN[5] = {0}, rowW[kNTruth] = {0};
    for (auto &j : jets)
    {
        if (!(j.e >= 28 && j.e < 100)) continue;
        int bp = 0;
        while (bp < NP - 1 && j.p >= pe[bp + 1]) bp++;
        const int bn = min(j.n, 5) - 1;
        cmP->Fill(bp + 0.5, kNTruth - 1 - j.cls + 0.5, j.w); colP[bp] += j.w;
        cmN->Fill(bn + 0.5, kNTruth - 1 - j.cls + 0.5, j.w); colN[bn] += j.w;
        rowW[j.cls] += j.w;
    }
    double tot = 0;
    for (int c = 0; c < kNTruth; c++) tot += rowW[c];
    printf("\nconfusion (raw E 28-100), column = composition of the selected jets; [share of all jets]\n%-9s", "truth");
    for (int b = 0; b < NP; b++) printf(" %9s", pl[b][0] == '#' ? ">=0.85" : pl[b]);
    printf(" |");
    for (int b = 0; b < 5; b++) printf(" %6s", b == 4 ? "n>=5" : Form("n=%s", nl[b]));
    printf(" | all\n");
    for (int c = 0; c < kNTruth; c++)
    {
        printf("%-9s", kTruthText[c]);
        for (int b = 0; b < NP; b++) printf(" %9.3f", colP[b] > 0 ? cmP->GetBinContent(b + 1, kNTruth - c) / colP[b] : 0);
        printf(" |");
        for (int b = 0; b < 5; b++) printf(" %6.3f", colN[b] > 0 ? cmN->GetBinContent(b + 1, kNTruth - c) / colN[b] : 0);
        printf(" | %.3f\n", rowW[c] / tot);
    }
    printf("%-9s", "[share]");
    for (int b = 0; b < NP; b++) printf(" %9.3f", colP[b] / tot);
    printf(" |");
    for (int b = 0; b < 5; b++) printf(" %6.3f", colN[b] / tot);
    printf("\n");
    for (int b = 0; b < NP; b++) for (int c = 0; c < kNTruth; c++) if (colP[b] > 0) cmP->SetBinContent(b + 1, c + 1, cmP->GetBinContent(b + 1, c + 1) / colP[b]);
    for (int b = 0; b < 5; b++) for (int c = 0; c < kNTruth; c++) if (colN[b] > 0) cmN->SetBinContent(b + 1, c + 1, cmN->GetBinContent(b + 1, c + 1) / colN[b]);
    gStyle->SetPaintTextFormat(".2f");
    for (TH2D *h : {cmP, cmN})
    {
        for (int c = 0; c < kNTruth; c++) h->GetYaxis()->SetBinLabel(kNTruth - c, kTruthName[c]);
        h->SetMinimum(0); h->SetMaximum(1); h->SetMarkerSize(1.6);
        h->GetXaxis()->SetLabelSize(0.05); h->GetYaxis()->SetLabelSize(0.06);
    }
    for (int b = 0; b < NP; b++) cmP->GetXaxis()->SetBinLabel(b + 1, Form("%s [%.0f%%]", pl[b], 100 * colP[b] / tot));
    for (int b = 0; b < 5; b++) cmN->GetXaxis()->SetBinLabel(b + 1, Form("%s [%.0f%%]", nl[b], 100 * colN[b] / tot));
    cv.cd(1); gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.12); cmP->Draw("COLZ TEXT");
    cv.cd(2); gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.12); cmN->Draw("COLZ TEXT");
    cv.Print(pdf.c_str());
    cv.Print((prefix + "_p2.png").c_str());
    // ---- page 3: Collins pair selection with P_leadpair ----
    cv.Clear();
    cv.Divide(2, 1);
    vector<TGraph *> proc(NE), ppur(NE), peff(NE), pmk(NE);
    printf("\nCollins pair (chosen pair = true pi0 pair), cut P_leadpair >= t\n");
    printf("raw E bin  pair right (all n>=2)  pair there  closest-mass right when there | n==2: eff purity\n");
    for (int b = 0; b < NE; b++)
    {
        proc[b] = new TGraph(); ppur[b] = new TGraph(); peff[b] = new TGraph(); pmk[b] = new TGraph();
        double wA = 0, wOK = 0, wAv = 0, wOKav = 0;
        for (auto &j : jets) if (j.e >= elo[b] && j.e < ehi[b] && j.n >= 2) { wA += j.w; wOK += j.pairOK * j.w; wAv += j.avail * j.w; if (j.avail) wOKav += j.pairOK * j.w; }
        if (!(wOK > 0)) continue;
        for (int t = 0; t < NT; t++)
        {
            const double thr = t / (double)(NT - 1);
            double sel = 0, selS = 0;
            for (auto &j : jets) if (j.e >= elo[b] && j.e < ehi[b] && j.n >= 2 && j.pp >= thr) { sel += j.w; selS += j.pairOK * j.w; }
            if (sel > 0) { proc[b]->SetPoint(proc[b]->GetN(), selS / wOK, selS / sel); ppur[b]->SetPoint(ppur[b]->GetN(), thr, selS / sel); }
            peff[b]->SetPoint(peff[b]->GetN(), thr, selS / wOK);
        }
        double s2 = 0, s2S = 0;
        for (auto &j : jets) if (j.e >= elo[b] && j.e < ehi[b] && j.n == 2) { s2 += j.w; s2S += j.pairOK * j.w; }
        if (s2 > 0) pmk[b]->SetPoint(0, s2S / wOK, s2S / s2);
        printf("%3.0f-%-4.0f   %.3f                %.3f       %.3f                      | %.3f %.3f\n", elo[b], ehi[b], wOK / wA, wAv / wA, wAv > 0 ? wOKav / wAv : 0,
               s2 > 0 ? s2S / wOK : 0, s2 > 0 ? s2S / s2 : 0);
    }
    cv.cd(1);
    gPad->SetGrid();
    gPad->DrawFrame(0, 0, 1, 1.0, "chosen pair = true #pi^{0} pair: cut P_{leadpair} #geq t (n #geq 2);efficiency;purity");
    TLegend lg3(0.13, 0.13, 0.62, 0.42);
    lg3.SetTextSize(0.033);
    for (int b = 0; b < NE; b++)
    {
        proc[b]->SetLineColor(col[b]); proc[b]->SetLineWidth(b == NE - 1 ? 3 : 2); proc[b]->Draw("L");
        pmk[b]->SetMarkerColor(col[b]); pmk[b]->SetMarkerStyle(b == NE - 1 ? 29 : 24); pmk[b]->SetMarkerSize(b == NE - 1 ? 2.0 : 1.3); pmk[b]->Draw("P");
        lg3.AddEntry(proc[b], Form("raw E %.0f-%.0f GeV%s", elo[b], ehi[b], b == NE - 1 ? " (all)" : ""), "l");
    }
    lg3.AddEntry(pmk[NE - 1], "two-point jets (n = 2)", "p");
    lg3.Draw();
    cv.cd(2);
    gPad->SetGrid();
    gPad->DrawFrame(0, 0, 1, 1.0, "raw E 28-100 GeV, n #geq 2;threshold t on P_{leadpair};fraction");
    ppur[NE - 1]->SetLineColor(kRed + 1); ppur[NE - 1]->SetLineWidth(3); ppur[NE - 1]->Draw("L");
    peff[NE - 1]->SetLineColor(kBlue + 1); peff[NE - 1]->SetLineWidth(3); peff[NE - 1]->Draw("L");
    TLegend lg4(0.15, 0.15, 0.5, 0.3);
    lg4.AddEntry(ppur[NE - 1], "pair purity", "l"); lg4.AddEntry(peff[NE - 1], "pair efficiency", "l");
    lg4.Draw();
    cv.Print(pdf.c_str());
    cv.Print((prefix + "_p3.png").c_str());
    cv.Print((pdf + "]").c_str());
    printf("wrote %s, %s_p1.png, %s_p2.png, %s_p3.png\n", pdf.c_str(), prefix.c_str(), prefix.c_str(), prefix.c_str());
    return 0;
}
