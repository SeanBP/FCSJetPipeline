// Sign check of the A_N chain on pi0s: EM jets with exactly two constituent
// points, their invariant mass (point momenta from the origin), and the blue
// and yellow left-right A_N in mass windows: pi0 (0.10-0.20 GeV), sideband
// (0.30-0.40), eta (0.45-0.65). Forward pi0 A_N for the polarized beam at
// x_F > 0 is positive at RHIC (FPD, FMS, PHENIX MPC, BRAHMS, at 62-510 GeV),
// so a negative pi0-window A_N means the sign convention of the chain
// (spin state -> beam spin, side -> left/right) is inverted.
//
// Same selection as AnHistograms: no LED, a listed EM trigger, Spin_config in
// {5,6,9,10}, reco_E_corr > 0, physical cut, corrected pT > 2 GeV and
// uncorrected pT > 1.15 x the threshold of a listed level. Same conventions
// as AnAsymmetry: blue left = x > 0, yellow left = x < 0, the other beam's
// spin held fixed (one cross ratio per state, inverse-variance mean). P is a
// constant (the sign is the question), --pol to change it.
//
// Usage: ./pi0_sign_check <jet file list> [P = 0.52]
// Build (in asymmetry/sivers): g++ -m64 -O2 pi0_sign_check.cpp ../../shared/VectorDict.cxx -o pi0_sign_check -I.. -I../../shared `root-config --cflags --libs` -std=c++11

#include <TChain.h>
#include <TH1D.h>
#include <TString.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "VectorDict.h"
#include "AnConfig.h"

static const int kUp[2][2] = {{0, 1}, {0, 2}}, kDown[2][2] = {{2, 3}, {1, 3}};
static const int kLeftSide[2] = {1, 0};

int main(int argc, char **argv)
{
    if (argc < 2) { std::cerr << "Usage: " << argv[0] << " <jet file list> [P]" << std::endl; return 1; }
    const double P = argc > 2 ? std::atof(argv[2]) : 0.52;
    TChain t("jetTree");
    { std::ifstream in(argv[1]); std::string f; while (in >> f) t.Add(f.c_str()); }
    Int_t spin_config = -1, led = 0, em[kNEmLevels] = {0}, em_tpc[kNEmLevels] = {0};
    Double_t sqrt_s = 0;
    std::vector<float> *eta = 0, *x = 0, *y = 0, *pt_raw = 0, *e_corr = 0, *xf_corr = 0;
    std::vector<float> *cpx = 0, *cpy = 0, *cpz = 0, *cE = 0;
    std::vector<int> *cj = 0;
    t.SetBranchAddress("Spin_config", &spin_config);
    t.SetBranchAddress("fcs_led", &led);
    t.SetBranchAddress("event_sqrt_s_gev", &sqrt_s);
    for (int l = 0; l < kNEmLevels; l++)
    {
        t.SetBranchAddress(TString::Format("fcsEM%d", l), &em[l]);
        if (t.GetBranch(TString::Format("fcsEM%d_tpc", l))) t.SetBranchAddress(TString::Format("fcsEM%d_tpc", l), &em_tpc[l]);
    }
    t.SetBranchAddress("reco_eta", &eta);
    t.SetBranchAddress("reco_x", &x);
    t.SetBranchAddress("reco_y", &y);
    t.SetBranchAddress("reco_pt", &pt_raw);
    t.SetBranchAddress("reco_E_corr", &e_corr);
    t.SetBranchAddress("reco_x_F_corr", &xf_corr);
    t.SetBranchAddress("reco_const_px", &cpx);
    t.SetBranchAddress("reco_const_py", &cpy);
    t.SetBranchAddress("reco_const_pz", &cpz);
    t.SetBranchAddress("reco_const_E", &cE);
    t.SetBranchAddress("reco_const_jetIndex", &cj);

    // mass windows x x_F bins x side x spin state
    const int nW = 3, nX = 4;
    const char *wname[nW] = {"pi0 0.10-0.20", "side 0.30-0.40", "eta 0.45-0.65"};
    const double wlo[nW] = {0.10, 0.30, 0.45}, whi[nW] = {0.20, 0.40, 0.65};
    const double xe[nX + 1] = {0.10, 0.15, 0.20, 0.30, 0.45};
    double N[nW][nX][2][4] = {};
    TH1D hm("hm", ";M_{#gamma#gamma} (GeV);jets", 50, 0, 1);
    long n2 = 0;
    const Long64_t ne = t.GetEntries();
    for (Long64_t i = 0; i < ne; i++)
    {
        t.GetEntry(i);
        if (led) continue;
        bool listed[kNEmLevels], any = false;
        for (int l = 0; l < kNEmLevels; l++) { listed[l] = em[l] || em_tpc[l]; any = any || listed[l]; }
        const int spin = SpinState(spin_config);
        if (!any || spin < 0) continue;
        for (size_t j = 0; j < eta->size(); j++)
        {
            const double E = (*e_corr)[j];
            if (!(E > 0) || !pass_physical_jet_cut(E, (*eta)[j], sqrt_s)) continue;
            if (E / std::cosh((*eta)[j]) <= kJetPtCorrMin) continue;
            bool trig = false;
            for (int l = 0; l < kNEmLevels; l++) trig = trig || (listed[l] && (*pt_raw)[j] > kTriggerPtFactor * kEmThresholdGeV[l]);
            if (!trig) continue;
            std::vector<int> c;
            for (size_t k = 0; k < cj->size(); k++) if ((*cj)[k] == (int)j) c.push_back(k);
            if (c.size() != 2) continue;
            const int a = c[0], b = c[1];
            const double sE = (*cE)[a] + (*cE)[b], spx = (*cpx)[a] + (*cpx)[b], spy = (*cpy)[a] + (*cpy)[b], spz = (*cpz)[a] + (*cpz)[b];
            const double m = std::sqrt(std::max(0.0, sE * sE - spx * spx - spy * spy - spz * spz));
            const double xf = (*xf_corr)[j];
            n2++;
            if (xf > 0.1) hm.Fill(m);
            int ix = -1;
            for (int k = 0; k < nX; k++) if (xf >= xe[k] && xf < xe[k + 1]) ix = k;
            if (ix < 0) continue;
            const int side = computeSide((*x)[j]) > 0 ? 1 : 0;
            for (int w = 0; w < nW; w++) if (m >= wlo[w] && m < whi[w]) N[w][ix][side][spin] += 1;
        }
    }
    printf("pi0_sign_check: %lld events, %ld two-point jets; P = %.2f\n", ne, n2, P);
    printf("two-point jets with x_F > 0.1 by mass (GeV):\n");
    for (int b = 1; b <= hm.GetNbinsX(); b += 2)
        printf("  %.2f-%.2f %8.0f\n", hm.GetBinLowEdge(b), hm.GetBinLowEdge(b) + 2 * hm.GetBinWidth(b), hm.GetBinContent(b) + hm.GetBinContent(b + 1));
    for (int bm = 0; bm < 2; bm++)
    {
        printf("%s beam (left = %s), A_N = eps / P (no <|cos phi|>):\n", bm ? "yellow" : "blue", bm ? "x < 0" : "x > 0");
        for (int w = 0; w < nW; w++)
            for (int k = 0; k < nX; k++)
            {
                const int sl = kLeftSide[bm], sr = 1 - sl;
                double sw = 0, swe = 0, nj = 0;
                for (int q = 0; q < 2; q++)
                {
                    const double Lu = N[w][k][sl][kUp[bm][q]], Ld = N[w][k][sl][kDown[bm][q]], Ru = N[w][k][sr][kUp[bm][q]], Rd = N[w][k][sr][kDown[bm][q]];
                    nj += Lu + Ld + Ru + Rd;
                    if (!(Lu > 0 && Ld > 0 && Ru > 0 && Rd > 0)) continue;
                    const double p = std::sqrt(Lu * Rd), mm = std::sqrt(Ld * Ru);
                    const double e = (p - mm) / (p + mm), de = p * mm / ((p + mm) * (p + mm)) * std::sqrt(1 / Lu + 1 / Ld + 1 / Ru + 1 / Rd);
                    sw += 1 / (de * de); swe += e / (de * de);
                }
                if (sw > 0)
                    printf("  %-15s x_F %.2f-%.2f  A_N %+.4f +- %.4f  (%+.1f sigma, %.0f jets)\n", wname[w], xe[k], xe[k + 1], swe / sw / P,
                           1 / std::sqrt(sw) / P, swe / sw * std::sqrt(sw), nj);
            }
    }
    return 0;
}
