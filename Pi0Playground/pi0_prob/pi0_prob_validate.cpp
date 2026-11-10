// Validation of the per-jet pi0 probability (shared/Pi0Prob.h) on the
// held-out half of the truth-labelled point jets (pi0_jet_dump, odd
// file_idx), cross-section weighted as pi0_prob_derive:
//   1. calibration: true fraction of jets with >= 1 resolved pi0 against mean
//      predicted P, in P bins, for all jets and n = 2-3, 4-6, 7+
//   2. ROC and AUC of P against two naive scores: "any pair within +-2 sigma
//      of the S peak of its pair energy" (one working point) and the
//      continuous min |m - mu| / sigma over the pairs
//   3. P for jets with and without a resolved pi0
//   4. leading point from a pi0 (jets with n >= 1): calibration of P_lead,
//      ROC against the photon classes (fewer points = more pi0-like) and the
//      leading point being in a pair within 2 sigma, and the purity of n <= 1,
//      2, 3 against a P_lead cut of the same efficiency
// Jets with n >= 2 only (P = 0 for n < 2 by definition); status counts are
// printed. Writes <prefix>.pdf, <prefix>.root (a tree of every jet's P) and a
// table on stdout.
//
// Build:
//   g++ -O2 pi0_prob_validate.cpp -o pi0_prob_validate -I../../shared -I../../jet_scale \
//       `root-config --cflags --libs` -std=c++11
// Run:  ./pi0_prob_validate "<dump glob>" <calib.json> <prefix> [original SimpleTree list; not needed for in-job dumps, which carry sigma_job]
#include "Pi0Prob.h"
#include "Pi0DumpWeights.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;

struct J { float p, base2, minpull; int y, n; double w; float pl, lbase2, negn; int yl; };

// weighted ROC (score high = signal); returns AUC
static double Roc(vector<J> v, double J::*dummy, float J::*score, TGraph *g)
{
    (void)dummy;
    sort(v.begin(), v.end(), [&](const J &a, const J &b) { return a.*score > b.*score; });
    double S = 0, B = 0;
    for (auto &j : v) (j.y ? S : B) += j.w;
    double s = 0, b = 0, auc = 0, ps = 0, pb = 0;
    int ip = 0;
    for (size_t i = 0; i < v.size(); i++)
    {
        (v[i].y ? s : b) += v[i].w;
        const bool last = i + 1 == v.size() || v[i + 1].*score != v[i].*score;
        if (!last) continue;
        auc += (b - pb) / B * 0.5 * (s + ps) / S;
        ps = s; pb = b;
        if (g) g->SetPoint(ip++, b / B, s / S);
    }
    return auc;
}

int main(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "Usage: %s \"<dump glob>\" <calib.json> <prefix> [SimpleTree list]\n", argv[0]); return 1; }
    gROOT->SetBatch(1);
    gStyle->SetOptStat(0);
    Pi0ProbCalib cal;
    string err;
    if (!cal.Load(argv[2], false, &err)) { fprintf(stderr, "%s\n", err.c_str()); return 1; }
    const string prefix = argv[3];

    TChain ch("jets");
    ch.Add(argv[1]);
    map<int, double> wOf;
    if (!Pi0DumpXSecWeights(argv[1], argc > 4 ? argv[4] : 0, wOf)) return 1;

    const int MX = 200;
    Int_t file_idx, n, k_res, ppi0[MX];
    Float_t jet_E, pE[MX], px[MX], py[MX], pz[MX];
    ch.SetBranchStatus("*", 0);
    for (auto b : {"file_idx", "n", "k_res", "jet_E", "pt_E", "pt_x", "pt_y", "pt_z", "pt_pi0"}) ch.SetBranchStatus(b, 1);
    ch.SetBranchAddress("file_idx", &file_idx); ch.SetBranchAddress("n", &n); ch.SetBranchAddress("k_res", &k_res);
    ch.SetBranchAddress("jet_E", &jet_E); ch.SetBranchAddress("pt_E", pE); ch.SetBranchAddress("pt_x", px);
    ch.SetBranchAddress("pt_y", py); ch.SetBranchAddress("pt_z", pz); ch.SetBranchAddress("pt_pi0", ppi0);

    TFile fo((prefix + ".root").c_str(), "RECREATE");
    TTree out("val", "pi0 probability, validation jets");
    Float_t o_p, o_lo, o_E, o_base, o_minpull; Int_t o_y, o_n, o_st; Double_t o_w;
    out.Branch("prob", &o_p); out.Branch("logodds", &o_lo); out.Branch("status", &o_st); out.Branch("n", &o_n);
    out.Branch("jet_E", &o_E); out.Branch("truth", &o_y); out.Branch("w", &o_w);
    out.Branch("base2sig", &o_base); out.Branch("minpull", &o_minpull);

    vector<J> all, lead;
    long st[4] = {0, 0, 0, 0};
    const Long64_t N = ch.GetEntries();
    for (Long64_t i = 0; i < N; i++)
    {
        ch.GetEntry(i);
        if (file_idx % 2 != 1) continue;
        auto it = wOf.find(file_idx);
        if (it == wOf.end() || !(it->second > 0)) continue;
        vector<double> E(n), X(n), Y(n), Z(n);
        for (int a = 0; a < n; a++)
        {
            const double r = sqrt(px[a] * px[a] + py[a] * py[a] + pz[a] * pz[a]);
            E[a] = pE[a]; X[a] = pE[a] * px[a] / r; Y[a] = pE[a] * py[a] / r; Z[a] = pE[a] * pz[a] / r;
        }
        const Pi0ProbResult res = Pi0Prob(cal, E, X, Y, Z, jet_E);
        st[res.status]++;
        // leading point: truth "from a pi0" and the mass-window baseline (in a pair within 2 sigma)
        int L = 0;
        for (int a = 1; a < n; a++) if (pE[a] > pE[L]) L = a;
        if (n >= 1 && res.lead >= 0)
        {
            J jl;
            jl.p = 0; jl.y = 0; jl.base2 = 0; jl.minpull = 0;
            jl.n = n; jl.w = it->second; jl.pl = res.lead; jl.yl = ppi0[L] >= 0; jl.negn = (float)-n;
            jl.lbase2 = 0;
            for (int b = 0; b < n; b++)
            {
                if (b == L) continue;
                const double m = Pi0PairMass(E[L], X[L], Y[L], Z[L], E[b], X[b], Y[b], Z[b]);
                const int e = Pi0ProbCalib::Bin(cal.eS, E[L] + E[b]);
                if (fabs(m - cal.sMean[e]) < 2 * cal.sSigma[e]) jl.lbase2 = 1;
            }
            lead.push_back(jl);
        }
        if (n < 2 || res.status != kPi0ProbOk) continue;
        double minpull = 1e9;
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
            {
                const double m = Pi0PairMass(E[a], X[a], Y[a], Z[a], E[b], X[b], Y[b], Z[b]);
                const int e = Pi0ProbCalib::Bin(cal.eS, E[a] + E[b]);
                minpull = min(minpull, fabs(m - cal.sMean[e]) / cal.sSigma[e]);
            }
        J j;
        j.p = res.prob; j.y = k_res >= 1; j.n = n; j.w = it->second;
        j.base2 = minpull < 2 ? 1.f : 0.f; j.minpull = (float)-minpull;
        all.push_back(j);
        o_p = res.prob; o_lo = res.logodds; o_st = res.status; o_n = n; o_E = jet_E; o_y = j.y; o_w = j.w;
        o_base = j.base2; o_minpull = (float)minpull;
        out.Fill();
    }
    printf("validation jets: status OK %ld, n<2 %ld, too many points %ld, no calib %ld\n", st[0], st[1], st[2], st[3]);

    // ---- calibration ----
    const int NB = 10;
    const char *cname[4] = {"all n#geq2", "n = 2-3", "n = 4-6", "n #geq 7"};
    const int ccol[4] = {kBlack, kBlue + 1, kOrange + 7, kRed + 1};
    TGraphErrors *gc[4];
    double wtot = 0, wsig = 0;
    for (auto &j : all) { wtot += j.w; if (j.y) wsig += j.w; }
    printf("weighted fraction with >= 1 resolved pi0: %.3f\n", wsig / wtot);
    printf("\ncalibration: <P> -> true fraction +- error [effective jets]\n");
    for (int q = 0; q < 4; q++)
    {
        double sw[NB] = {0}, sw2[NB] = {0}, sp[NB] = {0}, sy[NB] = {0};
        for (auto &j : all)
        {
            const int cls = j.n <= 3 ? 1 : (j.n <= 6 ? 2 : 3);
            if (q && cls != q) continue;
            const int b = min(NB - 1, (int)(j.p * NB));
            sw[b] += j.w; sw2[b] += j.w * j.w; sp[b] += j.w * j.p; sy[b] += j.w * j.y;
        }
        gc[q] = new TGraphErrors();
        printf("%s\n", cname[q]);
        for (int b = 0; b < NB; b++)
        {
            if (!(sw[b] > 0)) continue;
            const double neff = sw[b] * sw[b] / sw2[b], mp = sp[b] / sw[b], f = sy[b] / sw[b];
            const double e = sqrt(max(f * (1 - f), 1.0 / neff) / neff);
            printf("   %.3f -> %.3f +- %.3f [%.0f]\n", mp, f, e, neff);
            if (neff < 10) continue;
            const int ip = gc[q]->GetN();
            gc[q]->SetPoint(ip, mp, f); gc[q]->SetPointError(ip, 0, e);
        }
    }
    // ---- ROC ----
    TGraph *gP = new TGraph(), *gM = new TGraph(), *gB = new TGraph();
    const double aucP = Roc(all, 0, &J::p, gP), aucM = Roc(all, 0, &J::minpull, gM), aucB = Roc(all, 0, &J::base2, gB);
    double sS = 0, sB = 0, tS = 0, tB = 0;
    for (auto &j : all) { if (j.y) { sS += j.w; if (j.base2 > 0) tS += j.w; } else { sB += j.w; if (j.base2 > 0) tB += j.w; } }
    printf("\nAUC: P %.3f | min |m-mu|/sigma %.3f | any pair within 2 sigma %.3f (efficiency %.3f at false-positive rate %.3f)\n",
           aucP, aucM, aucB, tS / sS, tB / sB);
    for (int q = 1; q < 4; q++)
    {
        vector<J> sub;
        for (auto &j : all) { const int cls = j.n <= 3 ? 1 : (j.n <= 6 ? 2 : 3); if (cls == q) sub.push_back(j); }
        printf("  %-10s AUC P %.3f | min pull %.3f | 2 sigma %.3f\n", cname[q], Roc(sub, 0, &J::p, 0), Roc(sub, 0, &J::minpull, 0), Roc(sub, 0, &J::base2, 0));
    }

    // ---- leading point from a pi0 ----
    TGraphErrors *gl[5];
    const char *lname[5] = {"all n#geq1", "n = 1", "n = 2-3", "n = 4-6", "n #geq 7"};
    {
        double a = 0, y = 0;
        for (auto &j : lead) { a += j.w; if (j.yl) y += j.w; }
        printf("\nleading point from a pi0 (weighted): %.3f of %zu jets\ncalibration of P_lead: <P_lead> -> true fraction +- error [effective jets]\n", y / a, lead.size());
        for (int q = 0; q < 5; q++)
        {
            double sw[NB] = {0}, sw2[NB] = {0}, sp[NB] = {0}, sy[NB] = {0};
            for (auto &j : lead)
            {
                const int cls = j.n == 1 ? 1 : (j.n <= 3 ? 2 : (j.n <= 6 ? 3 : 4));
                if (q && cls != q) continue;
                const int b = min(NB - 1, (int)(j.pl * NB));
                sw[b] += j.w; sw2[b] += j.w * j.w; sp[b] += j.w * j.pl; sy[b] += j.w * j.yl;
            }
            gl[q] = new TGraphErrors();
            printf("%s\n", lname[q]);
            for (int b = 0; b < NB; b++)
            {
                if (!(sw[b] > 0)) continue;
                const double neff = sw[b] * sw[b] / sw2[b], mp = sp[b] / sw[b], f = sy[b] / sw[b];
                const double e = sqrt(max(f * (1 - f), 1.0 / neff) / neff);
                printf("   %.3f -> %.3f +- %.3f [%.0f]\n", mp, f, e, neff);
                if (neff < 10) continue;
                const int ip = gl[q]->GetN();
                gl[q]->SetPoint(ip, mp, f); gl[q]->SetPointError(ip, 0, e);
            }
        }
    }
    // swap the lead fields into p / y for the ROC helper
    vector<J> LJ = lead;
    for (auto &j : LJ) { j.y = j.yl; j.p = j.pl; j.base2 = j.lbase2; j.minpull = j.negn; }
    TGraph *glP = new TGraph(), *glN = new TGraph();
    const double aucL = Roc(LJ, 0, &J::p, glP), aucN = Roc(LJ, 0, &J::minpull, glN), aucLB = Roc(LJ, 0, &J::base2, 0);
    printf("\nleading point AUC: P_lead %.3f | fewer points (-n, the photon classes) %.3f | leading point in a pair within 2 sigma %.3f\n", aucL, aucN, aucLB);
    // purity of the photon classes against P_lead cuts at the same efficiency
    {
        double S = 0, Bk = 0;
        for (auto &j : LJ) (j.y ? S : Bk) += j.w;
        for (int ncut : {1, 2, 3})
        {
            double s = 0, b = 0;
            for (auto &j : LJ) if (j.n <= ncut) (j.y ? s : b) += j.w;
            // the P_lead cut with the same signal efficiency
            vector<J> v = LJ;
            sort(v.begin(), v.end(), [](const J &a, const J &c) { return a.p > c.p; });
            double ss = 0, bb = 0;
            for (auto &j : v) { if (ss >= s) break; (j.y ? ss : bb) += j.w; }
            printf("  n <= %d: efficiency %.3f, purity %.3f | P_lead at the same efficiency: purity %.3f\n", ncut, s / S, s / (s + b), ss / (ss + bb));
        }
    }

    // ---- plots ----
    TCanvas cv("cv", "", 1500, 500);
    cv.Divide(3, 1);
    cv.cd(1);
    TH1D fr("fr", "calibration (validation half);mean predicted P;true fraction with #geq1 resolved #pi^{0}", 10, 0, 1);
    fr.SetMaximum(1); fr.SetMinimum(0); fr.Draw();
    TLine dg(0, 0, 1, 1); dg.SetLineStyle(2); dg.Draw();
    TLegend L1(0.15, 0.62, 0.5, 0.88); L1.SetBorderSize(0); L1.SetFillStyle(0);
    for (int q = 0; q < 4; q++)
    {
        gc[q]->SetMarkerStyle(q ? 24 + q : 20); gc[q]->SetMarkerColor(ccol[q]); gc[q]->SetLineColor(ccol[q]);
        gc[q]->Draw("P same"); L1.AddEntry(gc[q], cname[q], "p");
    }
    L1.Draw();
    cv.cd(2);
    TH1D fr2("fr2", "ROC (validation half);false-positive rate (no resolved #pi^{0});efficiency (#geq1 resolved #pi^{0})", 10, 0, 1);
    fr2.SetMaximum(1); fr2.SetMinimum(0); fr2.Draw();
    gP->SetLineColor(kBlack); gP->SetLineWidth(2); gP->Draw("L same");
    gM->SetLineColor(kBlue + 1); gM->SetLineWidth(2); gM->SetLineStyle(2); gM->Draw("L same");
    TGraph wp(1); wp.SetPoint(0, tB / sB, tS / sS); wp.SetMarkerStyle(29); wp.SetMarkerSize(2); wp.SetMarkerColor(kRed + 1); wp.Draw("P same");
    TLegend L2(0.35, 0.15, 0.89, 0.4); L2.SetBorderSize(0); L2.SetFillStyle(0);
    L2.AddEntry(gP, Form("P (AUC %.3f)", aucP), "l");
    L2.AddEntry(gM, Form("min |m-#mu|/#sigma (AUC %.3f)", aucM), "l");
    L2.AddEntry(&wp, "any pair within 2#sigma", "p");
    L2.Draw();
    cv.cd(3);
    TH1D hs("hs", "P (validation half);P;fraction of jets", 25, 0, 1), hb("hb", "", 25, 0, 1);
    for (auto &j : all) (j.y ? hs : hb).Fill(j.p, j.w);
    hs.Scale(1.0 / hs.Integral()); hb.Scale(1.0 / hb.Integral());
    hs.SetLineColor(kRed + 1); hb.SetLineColor(kBlue + 1); hs.SetLineWidth(2); hb.SetLineWidth(2);
    hs.SetMaximum(1.2 * max(hs.GetMaximum(), hb.GetMaximum()));
    hs.Draw("hist"); hb.Draw("hist same");
    TLegend L3(0.3, 0.7, 0.75, 0.88); L3.SetBorderSize(0); L3.SetFillStyle(0);
    L3.AddEntry(&hs, "#geq1 resolved #pi^{0}", "l"); L3.AddEntry(&hb, "no resolved #pi^{0}", "l"); L3.Draw();
    cv.Print((prefix + ".pdf(").c_str());
    cv.Clear(); cv.Divide(2, 1);
    cv.cd(1);
    TH1D frl("frl", "leading point: calibration;mean predicted P_{lead};true fraction, leading point from a #pi^{0}", 10, 0, 1);
    frl.SetMaximum(1); frl.SetMinimum(0); frl.Draw();
    dg.Draw();
    TLegend Ll(0.15, 0.6, 0.5, 0.88); Ll.SetBorderSize(0); Ll.SetFillStyle(0);
    const int lcol[5] = {kBlack, kGreen + 2, kBlue + 1, kOrange + 7, kRed + 1};
    for (int q = 0; q < 5; q++)
    {
        gl[q]->SetMarkerStyle(q ? 23 + q : 20); gl[q]->SetMarkerColor(lcol[q]); gl[q]->SetLineColor(lcol[q]);
        gl[q]->Draw("P same"); Ll.AddEntry(gl[q], lname[q], "p");
    }
    Ll.Draw();
    cv.cd(2);
    TH1D frl2("frl2", "leading point: ROC;false-positive rate;efficiency (leading point from a #pi^{0})", 10, 0, 1);
    frl2.SetMaximum(1); frl2.SetMinimum(0); frl2.Draw();
    glP->SetLineColor(kBlack); glP->SetLineWidth(2); glP->Draw("L same");
    glN->SetLineColor(kBlue + 1); glN->SetLineWidth(2); glN->SetLineStyle(2); glN->SetMarkerStyle(20); glN->SetMarkerColor(kBlue + 1); glN->Draw("LP same");
    TLegend Ll2(0.35, 0.15, 0.89, 0.35); Ll2.SetBorderSize(0); Ll2.SetFillStyle(0);
    Ll2.AddEntry(glP, Form("P_{lead} (AUC %.3f)", aucL), "l");
    Ll2.AddEntry(glN, Form("photon classes, n #leq cut (AUC %.3f)", aucN), "lp");
    Ll2.Draw();
    cv.Print((prefix + ".pdf)").c_str());
    fo.cd(); out.Write(); hs.Write(); hb.Write();
    for (int q = 0; q < 4; q++) gc[q]->Write(Form("calib_%d", q));
    gP->Write("roc_P"); gM->Write("roc_minpull");
    fo.Close();
    printf("wrote %s.pdf, %s.root\n", prefix.c_str(), prefix.c_str());
    return 0;
}
