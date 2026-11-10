#ifndef COLLINS_HISTOGRAMS_H
#define COLLINS_HISTOGRAMS_H

// Collins yields filled by ../AnHistograms.cpp in the same pass as the A_N
// yields (one event loop, one output file). Jets with a leading-point pi0
// candidate pair (reco_pi0_lead_gup computed; shared/Pi0Prob.h Pi0LeadPair),
// E <= sqrt(s)/2, at the default trigger factor (sample 0, no suffix) or JP
// not EM at it (sample 1, _jponly), fill
//   h_xFg_<level>_<side>_<spin state>_<lp0|lp1|lp2>[_jponly]
// (TH3D: x_F 0-0.6 in 0.05 bins, g = sin(phi_S - phi_H) for spin up in 20 bins
// over [-1, 1], z = E_pair / jet E in 0-0.3-0.5-0.7-1; lp0/1/2 = P_leadpair
// (reco_pi0_lead_pair) < 0.5, 0.5-0.8, >= 0.8), and
//   h_absg_<level>_<lp>[_jponly], h_xfsumg_<level>_<lp>[_jponly]
// (TH2D x_F x z weighted with |g| / x_F, summed over sides and spin states).
// Needs AnConfig.h (levels, sides, spin states) included first.

#include "CollinsConfig.h"

#include <TH2D.h>
#include <TH3D.h>
#include <TString.h>

#include <cmath>

struct CollinsHistograms
{
    TH3D *g[2][kNTrgLevels][kNSides][kNSpinStates][kNPleadClsG] = {};
    TH2D *gabs[2][kNTrgLevels][kNPleadClsG] = {}, *gxfs[2][kNTrgLevels][kNPleadClsG] = {};

    // sample mg (0 default, 1 _jponly) of level l; sfx = the sample's suffix
    void Book(int mg, int l, const TString &sfx)
    {
        double xe[kNXfBinsG + 1], ge[kNGBins + 1];
        for (int i = 0; i <= kNXfBinsG; i++) xe[i] = kXfMinG + i * (kXfMaxG - kXfMinG) / kNXfBinsG;
        for (int i = 0; i <= kNGBins; i++) ge[i] = -1.0 + i * 2.0 / kNGBins;
        for (int c = 0; c < kNPleadClsG; c++)
        {
            for (int s = 0; s < kNSides; s++)
                for (int p = 0; p < kNSpinStates; p++)
                {
                    TString gn = TString::Format("h_xFg_%s_%s_%s_%s", kTrgLevelName[l], kSideName[s], kSpinStateName[p], kPleadClsNameG[c]) + sfx;
                    g[mg][l][s][p][c] = new TH3D(gn, gn + ";x_{F};g = sin(#phi_{S}-#phi_{H});z", kNXfBinsG, xe, kNGBins, ge, kNZBinsG, kZEdgesG);
                }
            TString gn = TString::Format("h_absg_%s_%s", kTrgLevelName[l], kPleadClsNameG[c]) + sfx;
            gabs[mg][l][c] = new TH2D(gn, gn + ";x_{F};z", kNXfBinsG, xe, kNZBinsG, kZEdgesG);
            gn = TString::Format("h_xfsumg_%s_%s", kTrgLevelName[l], kPleadClsNameG[c]) + sfx;
            gxfs[mg][l][c] = new TH2D(gn, gn + ";x_{F};z", kNXfBinsG, xe, kNZBinsG, kZEdgesG);
        }
    }
    // one jet of sample mg, level l: gg = g, zz = z_pair, plp = P_leadpair
    void Fill(int mg, int l, int side, int spin, double xf, double gg, double zz, double plp)
    {
        int pc = 0;
        while (pc < kNPleadClsG - 1 && plp >= kPleadEdgesG[pc + 1]) pc++;
        g[mg][l][side][spin][pc]->Fill(xf, gg, zz);
        gabs[mg][l][pc]->Fill(xf, zz, std::fabs(gg));
        gxfs[mg][l][pc]->Fill(xf, zz, xf);
    }
    // a jet has a Collins candidate: g computed and a P_leadpair
    static bool Has(double gg, double plp) { return gg >= -1.0001 && plp >= 0; }
};

#endif
