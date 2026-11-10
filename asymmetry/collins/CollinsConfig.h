#ifndef COLLINS_CONFIG_H
#define COLLINS_CONFIG_H

// Collins-analysis binning (asymmetry/collins/, see README.md there). The A_N
// (Sivers-type) settings are in ../AnConfig.h, which this follows.

// h_xFg_* yields: x_F (coarse) x g = sin(phi_S - phi_H) (reco_pi0_lead_gup,
// spin-up phi_S; shared/Pi0Prob.h Pi0LeadPair) x z = E_pair / jet E
// (reco_pi0_lead_zpair), in classes of P_leadpair (reco_pi0_lead_pair: the
// leading point is in a resolved pi0 pair; Collins needs the pair, so merged
// and partner-lost pi0, which P_lead also counts, are not wanted here). g in 20 bins keeps both the
// sign split (cross ratio of g > 0 / g < 0 within a side, divided by P <|g|>)
// and a binned fit open.
static const int kNXfBinsG = 12;
static const double kXfMinG = 0.0, kXfMaxG = 0.6;
static const int kNGBins = 20;
static const int kNZBinsG = 4;
static const double kZEdgesG[kNZBinsG + 1] = {0.0, 0.3, 0.5, 0.7, 1.0001};
static const int kNPleadClsG = 3;
static const double kPleadEdgesG[kNPleadClsG + 1] = {0.0, 0.5, 0.8, 1.0001};
static const char *const kPleadClsNameG[kNPleadClsG] = {"lp0", "lp1", "lp2"};  // P_leadpair < 0.5, 0.5-0.8, >= 0.8

#endif
