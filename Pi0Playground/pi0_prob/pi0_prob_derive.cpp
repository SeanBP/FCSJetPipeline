// Derives the inputs of the per-jet pi0 probability (shared/Pi0Prob.h) from
// truth-labelled point jets (pi0_jet_dump), training half only (even
// file_idx; the odd half is kept for pi0_prob_validate):
//   S(m | E_pair)  Crystal Ball fit to the true pi0 pairs of two-point jets,
//                  per pair-energy bin (Sean 2026-11-09: S from two-photon
//                  jets); the true pairs of jets with >= 3 points are fitted
//                  too, for comparison only
//   B(m | E_pair)  every other pair of every jet with >= 2 points, smoothed
//                  histogram density per pair-energy bin; drawn against its
//                  parts: B0 pairs from jets with no pi0 photon at all, B1 one
//                  point from a pi0, B2 points from two different pi0
//   pi_k(n, E_jet) fraction of jets with k resolved pi0 among n points, per
//                  raw jet energy bin (an n row merged over energy where a
//                  cell has fewer than kMinNeff effective jets, the previous
//                  n row's where the whole row has; pi_0 floored at 1%)
//   q(n, E_jet)    leading point: fraction from a pi0 among jets whose leading
//                  point is in no resolved pi0 pair (merged pi0, partner lost),
//                  n from 1, same fallbacks
// The prior is weighted by the MC cross section of each jet's file
// (jet_scale/JetXSecWeights.h on the SimpleTrees named in the dumps); the S
// and B shapes per pair-energy bin are unweighted (detector-driven, and the
// weights leave too few effective pairs), weighted if PI0PROB_SHAPE_WEIGHTED is set.
//
// Output: <prefix>.json (the calibration Pi0ProbCalib::Load reads),
// <prefix>.root (histograms), <prefix>.pdf (fits and overlays), and a summary
// on stdout.
//
// Build:
//   g++ -O2 pi0_prob_derive.cpp -o pi0_prob_derive -I../../shared -I../../jet_scale \
//       `root-config --cflags --libs` -std=c++11
// Run:  ./pi0_prob_derive "<dump glob>" <prefix> [version] [original SimpleTree list; not needed for in-job dumps, which carry sigma_job]
#include "Pi0Prob.h"
#include "Pi0DumpWeights.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;

static const double kEdgesE[] = {2, 6, 10, 15, 20, 30, 40, 60, 100, 300};
static const int kNE = sizeof(kEdgesE) / sizeof(double) - 1;
static const double kNLo[] = {2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
static const int kNN = sizeof(kNLo) / sizeof(double);
static const double kEdgesJ[] = {0, 15, 25, 35, 50, 70, 100, 10000};
static const int kNJ = sizeof(kEdgesJ) / sizeof(double) - 1;
static const int kKmax = 6;
// q(n, E_jet) n bins for the leading point, from one-point jets
static const double kQNLo[] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
static const int kQN = sizeof(kQNLo) / sizeof(double);
static const double kMinNeff = 30;
static const double kPi0Floor = 0.01;
static bool kShapeWeighted = false;
static const double kMmax = 1.2;
static const int kNmS = 600, kNmB = 120; // 2 MeV, 10 MeV

static double CB(double *x, double *p)
{
    const double t = (x[0] - p[1]) / p[2], a = fabs(p[3]), n = p[4];
    if (t > -a) return p[0] * exp(-0.5 * t * t);
    const double A = pow(n / a, n) * exp(-0.5 * a * a), B = n / a - a;
    return p[0] * A * pow(B - t, -n);
}

struct SFit { double mu, sig, alpha, n, emu, esig, chi2ndf, neff; bool ok; };

// fixAl/fixN > 0: the CB tail (alpha, n) is fixed to these values, and only
// the peak (mu, sigma) is fitted
static SFit FitS(TH1D *h, TF1 *keep, double fixAl = -1, double fixN = -1)
{
    SFit r = {0.135, 0.02, 1.5, 3, 0, 0, 0, 0, false};
    const double sw = h->GetSumOfWeights();
    double sw2 = 0;
    for (int i = 1; i <= h->GetNbinsX(); i++) sw2 += pow(h->GetBinError(i), 2);
    r.neff = sw2 > 0 ? sw * sw / sw2 : 0;
    if (r.neff < 50) return r;
    const int b0 = h->FindBin(0.06), b1 = h->FindBin(0.25);
    double best = 0, at = 0.13;
    for (int i = b0; i <= b1; i++) if (h->GetBinContent(i) > best) { best = h->GetBinContent(i); at = h->GetBinCenter(i); }
    TF1 f("cb", CB, 0.0, 0.40, 5);
    f.SetParameters(best, at, 0.015, 1.5, 3);
    f.SetParLimits(1, 0.08, 0.22); f.SetParLimits(2, 0.004, 0.08);
    f.SetParLimits(3, 0.3, 5); f.SetParLimits(4, 1.01, 30);
    if (fixAl > 0) { f.FixParameter(3, fixAl); f.FixParameter(4, fixN); }
    h->Fit(&f, "QRN0", "", 0.02, 0.40);
    r.mu = f.GetParameter(1); r.sig = f.GetParameter(2); r.alpha = fabs(f.GetParameter(3)); r.n = f.GetParameter(4);
    r.emu = f.GetParError(1); r.esig = f.GetParError(2);
    r.chi2ndf = f.GetNDF() > 0 ? f.GetChisquare() / f.GetNDF() : 0;
    r.ok = true;
    if (keep) for (int i = 0; i < 5; i++) keep->SetParameter(i, f.GetParameter(i));
    return r;
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "Usage: %s \"<dump glob>\" <prefix> [version] [SimpleTree list]\n", argv[0]); return 1; }
    const string prefix = argv[2], version = argc > 3 ? argv[3] : "pi0prob_v1";
    if (getenv("PI0PROB_SHAPE_WEIGHTED")) kShapeWeighted = true;
    gROOT->SetBatch(1);
    gStyle->SetOptStat(0);

    TChain ch("jets");
    ch.Add(argv[1]);
    // cross-section weight per file_idx (Pi0DumpWeights.h)
    map<int, double> wOf;
    if (!Pi0DumpXSecWeights(argv[1], argc > 4 ? argv[4] : 0, wOf)) return 1;

    const int MX = 200;
    Int_t file_idx, n, k_res, n_merged, n_lost, ppi0[MX];
    Float_t jet_E, pE[MX], px[MX], py[MX], pz[MX];
    ch.SetBranchStatus("*", 0);
    for (auto b : {"file_idx", "n", "k_res", "n_merged", "n_lost", "jet_E", "pt_E", "pt_x", "pt_y", "pt_z", "pt_pi0"}) ch.SetBranchStatus(b, 1);
    ch.SetBranchAddress("file_idx", &file_idx); ch.SetBranchAddress("n", &n); ch.SetBranchAddress("k_res", &k_res);
    ch.SetBranchAddress("n_merged", &n_merged); ch.SetBranchAddress("n_lost", &n_lost);
    ch.SetBranchAddress("jet_E", &jet_E); ch.SetBranchAddress("pt_E", pE); ch.SetBranchAddress("pt_x", px);
    ch.SetBranchAddress("pt_y", py); ch.SetBranchAddress("pt_z", pz); ch.SetBranchAddress("pt_pi0", ppi0);

    TFile fo((prefix + ".root").c_str(), "RECREATE");
    TH1D *hS2[kNE], *hS3[kNE], *hB[kNE], *hB0[kNE], *hB1[kNE], *hB2[kNE];
    for (int e = 0; e < kNE; e++)
    {
        TString t = TString::Format("%g-%g GeV", kEdgesE[e], kEdgesE[e + 1]);
        hS2[e] = new TH1D(Form("hS2_%d", e), "true pi0 pairs, 2-point jets, " + t + ";m_{#gamma#gamma} (GeV)", kNmS, 0, kMmax);
        hS3[e] = new TH1D(Form("hS3_%d", e), "true pi0 pairs, #geq3-point jets, " + t + ";m_{#gamma#gamma} (GeV)", kNmS, 0, kMmax);
        hB[e] = new TH1D(Form("hB_%d", e), "background pairs, " + t + ";m_{#gamma#gamma} (GeV)", kNmB, 0, kMmax);
        hB0[e] = new TH1D(Form("hB0_%d", e), "B0, " + t, kNmB, 0, kMmax);
        hB1[e] = new TH1D(Form("hB1_%d", e), "B1, " + t, kNmB, 0, kMmax);
        hB2[e] = new TH1D(Form("hB2_%d", e), "B2, " + t, kNmB, 0, kMmax);
        for (TH1D *h : {hS2[e], hS3[e], hB[e], hB0[e], hB1[e], hB2[e]}) h->Sumw2();
    }
    // prior sums [n][E][k] and sum w^2 per cell
    vector<double> pw(kNN * kNJ * (kKmax + 1), 0), cw(kNN * kNJ, 0), cw2(kNN * kNJ, 0);
    double wJets = 0, wRes = 0, wMergedOnly = 0, wLostOnly = 0;
    vector<double> qw(kQN * kNJ, 0), qw2(kQN * kNJ, 0), qpi(kQN * kNJ, 0), leadAll(kQN * kNJ, 0), leadRes(kQN * kNJ, 0), leadPi0(kQN * kNJ, 0);
    long nTrain = 0, nPairs = 0;
    TH1D hn("h_n", "points per jet (training, weighted);n;jets", 30, 0, 30);

    const Long64_t N = ch.GetEntries();
    for (Long64_t i = 0; i < N; i++)
    {
        ch.GetEntry(i);
        if (file_idx % 2 != 0) continue;
        auto it = wOf.find(file_idx);
        if (it == wOf.end() || !(it->second > 0)) continue;
        const double w = it->second;
        // mass shapes per pair-energy bin: unweighted unless -w (the shape is
        // detector-driven; the weights leave ~70 effective pairs per bin)
        const double ws = kShapeWeighted ? w : 1.0;
        nTrain++;
        hn.Fill(n, w);
        // leading point: from a pi0, and in a resolved pi0 pair of the jet?
        if (n >= 1)
        {
            int L = 0;
            for (int a = 1; a < n; a++) if (pE[a] > pE[L]) L = a;
            bool inRes = false;
            for (int a = 0; a < n; a++) if (a != L && ppi0[L] >= 0 && ppi0[a] == ppi0[L]) inRes = true;
            const int iq = Pi0ProbCalib::LoBin(vector<double>(kQNLo, kQNLo + kQN), n);
            const int ie = Pi0ProbCalib::Bin(vector<double>(kEdgesJ, kEdgesJ + kNJ + 1), jet_E);
            leadAll[iq * kNJ + ie] += w;
            if (inRes) leadRes[iq * kNJ + ie] += w;
            else
            {
                qw[iq * kNJ + ie] += w; qw2[iq * kNJ + ie] += w * w;
                if (ppi0[L] >= 0) qpi[iq * kNJ + ie] += w;
            }
            if (ppi0[L] >= 0) leadPi0[iq * kNJ + ie] += w;
        }
        if (n < 2) continue;
        wJets += w;
        if (k_res > 0) wRes += w;
        else if (n_merged > 0) wMergedOnly += w;
        else if (n_lost > 0) wLostOnly += w;
        bool anyPi0 = false;
        for (int a = 0; a < n; a++) if (ppi0[a] >= 0) anyPi0 = true;
        // 4-vectors along the direction from the origin, as the jets
        vector<double> E(n), X(n), Y(n), Z(n);
        for (int a = 0; a < n; a++)
        {
            const double r = sqrt(px[a] * px[a] + py[a] * py[a] + pz[a] * pz[a]);
            E[a] = pE[a]; X[a] = pE[a] * px[a] / r; Y[a] = pE[a] * py[a] / r; Z[a] = pE[a] * pz[a] / r;
        }
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
            {
                const double m = Pi0PairMass(E[a], X[a], Y[a], Z[a], E[b], X[b], Y[b], Z[b]);
                const int e = Pi0ProbCalib::Bin(vector<double>(kEdgesE, kEdgesE + kNE + 1), E[a] + E[b]);
                nPairs++;
                if (ppi0[a] >= 0 && ppi0[a] == ppi0[b]) { (n == 2 ? hS2[e] : hS3[e])->Fill(m, ws); continue; }
                hB[e]->Fill(m, ws);
                const int np = (ppi0[a] >= 0) + (ppi0[b] >= 0);
                if (np == 0 && !anyPi0) hB0[e]->Fill(m, ws);
                else if (np == 1) hB1[e]->Fill(m, ws);
                else if (np == 2) hB2[e]->Fill(m, ws);
            }
        const int in = Pi0ProbCalib::LoBin(vector<double>(kNLo, kNLo + kNN), n);
        const int ie = Pi0ProbCalib::Bin(vector<double>(kEdgesJ, kEdgesJ + kNJ + 1), jet_E);
        const int k = k_res > kKmax ? kKmax : k_res;
        pw[(in * kNJ + ie) * (kKmax + 1) + k] += w;
        cw[in * kNJ + ie] += w; cw2[in * kNJ + ie] += w * w;
    }
    printf("training jets %ld (all n), pairs %ld\n", nTrain, nPairs);
    printf("jets with n>=2 (weighted): >=1 resolved pi0 %.3f; no resolved, a merged pi0 %.3f; only a lost pi0 photon %.3f\n",
           wRes / wJets, wMergedOnly / wJets, wLostOnly / wJets);

    // ---- S ----
    vector<double> sMu(kNE), sSig(kNE), sAl(kNE), sN(kNE);
    TCanvas cv("cv", "", 1200, 900);
    const string pdf = prefix + ".pdf";
    cv.Print((pdf + "[").c_str());
    cv.Divide(3, 3);
    printf("\nS: Crystal Ball fits, true pi0 pairs (2-point jets | >=3-point jets)\n");
    printf("%-10s %8s %8s %8s %6s %6s %7s | %8s %8s %7s\n", "E_pair", "mu", "sigma", "sig/mu", "alpha", "n", "neff", "mu3", "sigma3", "neff3");
    // S per pair-energy bin: the 2-point-jet fit where it is good (n_eff >=
    // kMinNeffS and sigma off its 0.08 limit). Otherwise all true pairs (2-point
    // and >= 3-point jets): CB tail free where the nearest good 2-point bin is
    // not available (lowest bins), else the tail (alpha, n) fixed to that bin's
    // and only mu, sigma fitted; with too few pairs even then, the previous
    // bin's S. (pi0prob_v1 copied the previous bin or a 0.135/0.020 default.)
    const double kMinNeffS = 300;
    SFit last = {0.135, 0.02, 1.5, 3, 0, 0, 0, 0, false};
    bool haveGood = false;
    for (int e = 0; e < kNE; e++)
    {
        cv.cd(e + 1);
        TF1 *f2 = new TF1(Form("fS2_%d", e), CB, 0, 0.4, 5), *f3 = new TF1(Form("fS3_%d", e), CB, 0, 0.4, 5);
        SFit s2 = FitS(hS2[e], f2), s3 = FitS(hS3[e], f3);
        const bool good2 = s2.ok && s2.neff >= kMinNeffS && s2.sig < 0.075;
        if (!good2)
        {
            TH1D *hall = (TH1D *)hS2[e]->Clone(Form("hSall_%d", e));
            hall->Add(hS3[e]);
            SFit sa = haveGood ? FitS(hall, f2, last.alpha, last.n) : FitS(hall, f2);
            if (sa.ok && sa.neff >= kMinNeffS && sa.sig < 0.075)
            {
                printf("  (bin %d: 2-point fit %s; using all true pairs, n_eff %.0f%s)\n", e, s2.ok ? "poor" : "too few pairs", sa.neff,
                       haveGood ? ", CB tail fixed to the previous good bin" : "");
                s2 = sa;
            }
            else { s2 = last; printf("  (bin %d: too few pairs even with >= 3-point jets, using the previous bin's S)\n", e); }
        }
        else haveGood = true;
        last = s2;
        sMu[e] = s2.mu; sSig[e] = s2.sig; sAl[e] = s2.alpha; sN[e] = s2.n;
        printf("%4g-%-5g %8.4f %8.4f %8.3f %6.2f %6.2f %7.0f | %8.4f %8.4f %7.0f\n", kEdgesE[e], kEdgesE[e + 1], s2.mu, s2.sig, s2.sig / s2.mu,
               s2.alpha, s2.n, s2.neff, s3.ok ? s3.mu : 0, s3.ok ? s3.sig : 0, s3.neff);
        TH1D *a = (TH1D *)hS2[e]->Clone(), *b = (TH1D *)hS3[e]->Clone();
        a->Rebin(4); b->Rebin(4);
        a->GetXaxis()->SetRangeUser(0, 0.4);
        if (a->Integral() > 0) a->Scale(1.0 / a->Integral());
        if (b->Integral() > 0) b->Scale(1.0 / b->Integral());
        a->SetLineColor(kBlue + 1); b->SetLineColor(kOrange + 7);
        a->SetMaximum(1.3 * max(a->GetMaximum(), b->GetMaximum()));
        a->Draw("hist"); b->Draw("hist same");
        if (s2.ok)
        {
            TF1 *g = new TF1(Form("fd_%d", e), CB, 0, 0.4, 5);
            for (int p = 0; p < 5; p++) g->SetParameter(p, f2->GetParameter(p));
            // scale to the rebinned, normalised histogram
            g->SetParameter(0, f2->GetParameter(0) * 4 / hS2[e]->Integral());
            g->SetLineColor(kBlue + 1); g->SetLineStyle(2); g->Draw("same");
        }
        TLegend *L = new TLegend(0.45, 0.65, 0.89, 0.89);
        L->SetBorderSize(0); L->SetFillStyle(0);
        L->AddEntry(a, Form("2-pt jets: #sigma/#mu %.2f", s2.sig / s2.mu), "l");
        L->AddEntry(b, s3.ok ? Form("#geq3-pt jets: #sigma/#mu %.2f", s3.sig / s3.mu) : "#geq3-pt jets", "l");
        L->Draw();
    }
    cv.Print(pdf.c_str());

    // ---- B ----
    vector<double> mB(kNmB + 1);
    for (int i = 0; i <= kNmB; i++) mB[i] = kMmax * i / kNmB;
    vector<vector<double> > bpdf(kNE, vector<double>(kNmB, 0));
    cv.Clear(); cv.Divide(3, 3);
    printf("\nB: background pair shapes; mean mass and fraction below 0.25 GeV (all | B0 | B1 | B2), weighted\n");
    for (int e = 0; e < kNE; e++)
    {
        TH1D *s = (TH1D *)hB[e]->Clone(Form("hBs_%d", e));
        if (s->GetSumOfWeights() > 0)
        {
            s->Smooth(1);
            const double I = s->Integral() * (kMmax / kNmB);
            for (int i = 0; i < kNmB; i++) bpdf[e][i] = I > 0 ? s->GetBinContent(i + 1) / I : 0;
        }
        auto frac = [](TH1D *h) { double t = h->Integral(); return t > 0 ? h->Integral(1, h->FindBin(0.2499)) / t : 0.0; };
        printf("%4g-%-5g mean %.3f %.3f %.3f %.3f   <0.25 %.2f %.2f %.2f %.2f   (eff. pairs %0.f %0.f %0.f %0.f)\n", kEdgesE[e], kEdgesE[e + 1],
               hB[e]->GetMean(), hB0[e]->GetMean(), hB1[e]->GetMean(), hB2[e]->GetMean(), frac(hB[e]), frac(hB0[e]), frac(hB1[e]), frac(hB2[e]),
               hB[e]->GetEffectiveEntries(), hB0[e]->GetEffectiveEntries(), hB1[e]->GetEffectiveEntries(), hB2[e]->GetEffectiveEntries());
        cv.cd(e + 1);
        TH1D *hs[4] = {(TH1D *)hB[e]->Clone(), (TH1D *)hB0[e]->Clone(), (TH1D *)hB1[e]->Clone(), (TH1D *)hB2[e]->Clone()};
        const int col[4] = {kBlack, kGreen + 2, kRed + 1, kViolet + 1};
        const char *lab[4] = {"all background", "B0: jets with no #pi^{0} #gamma", "B1: one #gamma from a #pi^{0}", "B2: two different #pi^{0}"};
        double mx = 0;
        for (int q = 0; q < 4; q++)
        {
            hs[q]->Rebin(2);
            if (hs[q]->Integral() > 0) hs[q]->Scale(1.0 / hs[q]->Integral());
            hs[q]->SetLineColor(col[q]); hs[q]->SetLineWidth(q == 0 ? 3 : 2);
            mx = max(mx, hs[q]->GetMaximum());
        }
        hs[0]->SetMaximum(1.3 * mx); hs[0]->GetXaxis()->SetRangeUser(0, 1.0);
        hs[0]->SetTitle(TString::Format("background pairs, E_{pair} %g-%g GeV;m_{#gamma#gamma} (GeV);fraction", kEdgesE[e], kEdgesE[e + 1]));
        hs[0]->Draw("hist");
        for (int q = 1; q < 4; q++) hs[q]->Draw("hist same");
        if (e == 0)
        {
            TLegend *L = new TLegend(0.35, 0.55, 0.89, 0.89);
            L->SetBorderSize(0); L->SetFillStyle(0);
            for (int q = 0; q < 4; q++) L->AddEntry(hs[q], lab[q], "l");
            L->Draw();
        }
    }
    cv.Print(pdf.c_str());

    // ---- prior ----
    vector<double> prior(kNN * kNJ * (kKmax + 1), 0);
    printf("\nprior pi_k(n, E_jet): pi_0 per cell [effective jets]; '*' = merged over E (n row), '^' = the previous n row's\n%-6s", "n\\E");
    for (int j = 0; j < kNJ; j++) printf("  %6g-%-6g", kEdgesJ[j], kEdgesJ[j + 1]);
    printf("\n");
    for (int in = 0; in < kNN; in++)
    {
        // the n row summed over E, the fallback
        vector<double> row(kKmax + 1, 0);
        double rw = 0;
        for (int ie = 0; ie < kNJ; ie++)
        {
            for (int k = 0; k <= kKmax; k++) row[k] += pw[(in * kNJ + ie) * (kKmax + 1) + k];
            rw += cw[in * kNJ + ie];
        }
        double rw2 = 0;
        for (int ie = 0; ie < kNJ; ie++) rw2 += cw2[in * kNJ + ie];
        const bool rowOk = rw2 > 0 && rw * rw / rw2 >= kMinNeff;
        printf("%-6s", in + 1 < kNN ? Form("%g-%g", kNLo[in], kNLo[in + 1] - 1) : Form("%g+", kNLo[in]));
        for (int ie = 0; ie < kNJ; ie++)
        {
            const int c = in * kNJ + ie;
            const double neff = cw2[c] > 0 ? cw[c] * cw[c] / cw2[c] : 0;
            const bool own = neff >= kMinNeff;
            double *pk = &prior[c * (kKmax + 1)];
            if (own) for (int k = 0; k <= kKmax; k++) pk[k] = pw[c * (kKmax + 1) + k] / cw[c];
            else if (rowOk) for (int k = 0; k <= kKmax; k++) pk[k] = row[k] / rw;
            else if (in > 0) for (int k = 0; k <= kKmax; k++) pk[k] = prior[((in - 1) * kNJ + ie) * (kKmax + 1) + k];
            else { pk[0] = 1; for (int k = 1; k <= kKmax; k++) pk[k] = 0; }
            // floor pi_0 at kPi0Floor: pi_0 = 0 would make P = 1 whatever the masses
            if (pk[0] < kPi0Floor)
            {
                double t = 0;
                for (int k = 1; k <= kKmax; k++) t += pk[k];
                for (int k = 1; k <= kKmax; k++) pk[k] *= (1 - kPi0Floor) / (t > 0 ? t : 1);
                pk[0] = kPi0Floor;
            }
            printf("  %5.3f%s[%5.0f]", pk[0], own ? " " : (rowOk ? "*" : "^"), neff);
        }
        printf("\n");
    }
    printf("mean k per n row:");
    for (int in = 0; in < kNN; in++)
    {
        double s = 0, t = 0;
        for (int ie = 0; ie < kNJ; ie++) for (int k = 0; k <= kKmax; k++) { s += k * pw[(in * kNJ + ie) * (kKmax + 1) + k]; t += pw[(in * kNJ + ie) * (kKmax + 1) + k]; }
        printf(" n%g: %.2f", kNLo[in], t > 0 ? s / t : 0);
    }
    printf("\n");

    // ---- q(n, E) for the leading point ----
    vector<double> qv(kQN * kNJ, 0);
    printf("\nleading point: q = P(from a pi0 | not in a resolved pair) per cell [effective jets]; '*' n row, '^' previous row\n%-6s", "n\\E");
    for (int j = 0; j < kNJ; j++) printf("  %6g-%-6g", kEdgesJ[j], kEdgesJ[j + 1]);
    printf("\n");
    for (int iq = 0; iq < kQN; iq++)
    {
        double rw = 0, rw2 = 0, rp = 0;
        for (int ie = 0; ie < kNJ; ie++) { rw += qw[iq * kNJ + ie]; rw2 += qw2[iq * kNJ + ie]; rp += qpi[iq * kNJ + ie]; }
        const bool rowOk = rw2 > 0 && rw * rw / rw2 >= kMinNeff;
        printf("%-6s", iq + 1 < kQN ? Form("%g-%g", kQNLo[iq], kQNLo[iq + 1] - 1) : Form("%g+", kQNLo[iq]));
        for (int ie = 0; ie < kNJ; ie++)
        {
            const int c = iq * kNJ + ie;
            const double neff = qw2[c] > 0 ? qw[c] * qw[c] / qw2[c] : 0;
            const bool own = neff >= kMinNeff;
            qv[c] = own ? qpi[c] / qw[c] : (rowOk ? rp / rw : (iq > 0 ? qv[(iq - 1) * kNJ + ie] : 0.5));
            printf("  %5.3f%s[%5.0f]", qv[c], own ? " " : (rowOk ? "*" : "^"), neff);
        }
        printf("\n");
    }
    {
        double a = 0, r = 0, p = 0;
        for (int c = 0; c < kQN * kNJ; c++) { a += leadAll[c]; r += leadRes[c]; p += leadPi0[c]; }
        printf("leading point (weighted, all jets): from a pi0 %.3f, of which in a resolved pair %.3f\n", p / a, r / p);
    }

    cv.Clear(); cv.cd(); hn.Draw("hist"); cv.SetLogy(1);
    cv.Print(pdf.c_str());
    cv.SetLogy(0);
    cv.Print((pdf + "]").c_str());
    fo.Write();
    fo.Close();

    // ---- JSON ----
    FILE *js = fopen((prefix + ".json").c_str(), "w");
    auto arr = [&](const char *key, const vector<double> &v, bool comma = true) {
        fprintf(js, "  \"%s\": [", key);
        for (size_t i = 0; i < v.size(); i++) fprintf(js, "%s%.6g", i ? ", " : "", v[i]);
        fprintf(js, "]%s\n", comma ? "," : "");
    };
    fprintf(js, "{\n  \"version\": \"%s\",\n", version.c_str());
    fprintf(js, "  \"source\": \"jes_pixring_v1 SimpleTrees (k = 0.947), even files; pi0_jet_dump + pi0_prob_derive\",\n");
    fprintf(js, "  \"mass\": \"point 4-vectors along the direction from the origin, as the jets\",\n");
    fprintf(js, "  \"eps\": 1e-6,\n  \"max_points\": 20,\n  \"kmax\": %d,\n", kKmax);
    arr("S_pairE_edges", vector<double>(kEdgesE, kEdgesE + kNE + 1));
    arr("S_mean", sMu); arr("S_sigma", sSig); arr("S_alpha", sAl); arr("S_n", sN);
    arr("B_pairE_edges", vector<double>(kEdgesE, kEdgesE + kNE + 1));
    arr("B_mass_edges", mB);
    fprintf(js, "  \"B_pdf\": [\n");
    for (int e = 0; e < kNE; e++)
    {
        fprintf(js, "    [");
        for (int i = 0; i < kNmB; i++) fprintf(js, "%s%.6g", i ? ", " : "", bpdf[e][i]);
        fprintf(js, "]%s\n", e + 1 < kNE ? "," : "");
    }
    fprintf(js, "  ],\n");
    arr("prior_n_lo", vector<double>(kNLo, kNLo + kNN));
    arr("prior_jetE_edges", vector<double>(kEdgesJ, kEdgesJ + kNJ + 1));
    arr("prior", prior);
    arr("lead_q_n_lo", vector<double>(kQNLo, kQNLo + kQN));
    arr("lead_q", qv, false);
    fprintf(js, "}\n");
    fclose(js);
    printf("\nwrote %s.json, %s.root, %s.pdf\n", prefix.c_str(), prefix.c_str(), prefix.c_str());
    return 0;
}
