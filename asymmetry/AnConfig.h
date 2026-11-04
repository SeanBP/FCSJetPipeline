#ifndef AN_CONFIG_H
#define AN_CONFIG_H

// Selection and binning of the A_N yield histograms (AnHistograms.cpp).

#include <cmath>

// ------------------ EM trigger levels ------------------
// A jet counts if its uncorrected pT (reco_pt: no UE subtraction, no JES, the
// energy the hardware trigger sees) passes kTriggerPtFactor x the threshold
// of any EM trigger listed in its event (the FMS EM-jet rule); its level is
// the lowest such. The 2 GeV minimum (kJetPtCorrMin) is on corrected pT. Listed = recorded the event (data: nominal trigger ids, so
// after prescale); fcsEMk_tpc counts as fcsEMk (same condition, own id; run
// 23007048: 187 events list fcsEM2_tpc without fcsEM2, of 155506). fcsEM3
// (runs 23101059-23108024 only) is not used. Trigger emulation cannot stand
// in for the list: at run 23007048's settings it fires for only 86/71/47% of
// the events listing EM0/1/2.
static const int kNEmLevels = 3;
static const char *const kEmLevelName[kNEmLevels] = {"EM0", "EM1", "EM2"};

// EM trigger thresholds, GeV, compared with uncorrected jet pT (reco_pt: the
// hardware trigger sees raw energy). The nominal online values, taken as
// given for now (Sean, 2026-11-03); to be checked against measured turn-ons
// (TriggerPlayground/turn_on.cpp). Online (A. Ogawa's Run 22 FCS trigger
// pages): each EM bit is a stage-2 ECal 4x4 + nearest HCal 4x4 E_T sum above
// EMTHR with ECal 4x4 > 0.25 x HCal 4x4, ORed north/south at stage 3; 0.0316
// GeV E_T per count, EM0/1/2 = 64/112/160 counts = 2.0/3.5/5.0 GeV from run
// 22359013 (63/111/158 before; stage_params.txt still has those). The online
// E_T scale moved with the ET gain files (new at 23027047, removing an old
// 0.7 factor) and the ECal electronics gain (halved at 23066056), so the
// offline turn-ons can differ between those periods. Earlier rough estimate
// on run 23007048: uncorrected point-jet pT turns on near 2.6 (EM1) and 3.7
// GeV (EM2), about 0.74 of the nominal 3.5 and 5.0.
static const double kEmThresholdGeV[kNEmLevels] = {2.0, 3.5, 5.0};

// A jet counts only if its uncorrected pT is above this factor times the
// threshold of an EM trigger listed in its event (above the turn-on); 1.15 as
// the FMS EM-jet analysis (UnderlyingEventPlayground/fms_emjet_reference.txt)
static const double kTriggerPtFactor = 1.15;

// Corrected jet pT minimum (GeV), after UE subtraction and calibration, as
// the FMS analysis; the same value as the reco cut at jet finding
static const double kJetPtCorrMin = 2.0;

// ------------------ Point multiplicity ------------------
// reco_nphotons (FCS points above point_energy_min in the jet) classes:
// 0, 1, 2, 3, 4, >= 5
static const int kNPhotonClasses = 6;
inline int PhotonClass(int nphotons) { return nphotons < 0 ? 0 : (nphotons >= kNPhotonClasses - 1 ? kNPhotonClasses - 1 : nphotons); }
inline const char *PhotonClassName(int k)
{
    static const char *const n[kNPhotonClasses] = {"np0", "np1", "np2", "np3", "np4", "np5plus"};
    return n[k];
}

// ------------------ Spin states ------------------
// Spin_config is StSpinDbMaker::spin4usingBX7(). Its bits carry CAD's "+" and
// "-", and CAD "+" is spin DOWN at STAR, for both beams (STAR spin-pattern
// note, drupal.star.bnl.gov/STAR/blog/oleg/spin-patterns-and-polarization-direction;
// used by the Run 15/17 FMS EM-jet A_N code, BrightSTAR, and D. Kapukchyan's
// Run 22 FCS pi0 A_N, StFwdDataEvent::BlueSpin/YellowSpin):
//   10 = blue up, yellow up     9 = blue up, yellow down
//    6 = blue down, yellow up   5 = blue down, yellow down
// Only these four enter; -1 (no spin DB entry) and the rest are dropped.
// Until 2026-11-04 this mapping took StSpinDbMaker's labels at face value
// (5 = blue up, yellow up), which flips both beams and so the sign of both
// A_N; yield files written then have no "spin_convention" entry and
// AnAsymmetry relabels them (kSpinConventionTag).
static const int kNSpinStates = 4;
static const char *const kSpinStateName[kNSpinStates] = {"BuYu", "BuYd", "BdYu", "BdYd"};
inline int SpinState(int spin_config)
{
    switch (spin_config)
    {
        case 10: return 0;
        case 9: return 1;
        case 6: return 2;
        case 5: return 3;
        default: return -1;
    }
}
// Written by AnHistograms as TNamed "spin_convention"; a yield file without
// it used the old mapping, whose state p is the physical state 3 - p
static const char *const kSpinConventionTag = "physical: Spin_config 10 = BuYu, 9 = BuYd, 6 = BdYu, 5 = BdYd (CAD + = down)";

// ------------------ Detector side ------------------
// computeSide() (JetParameters.h) on reco_x: -1 for x <= 0, +1 for x > 0.
// Left/right relative to each beam is assigned when A_N is formed
// (AnAsymmetry.cpp): blue left = xpos, yellow left = xneg.
static const int kNSides = 2;
static const char *const kSideName[kNSides] = {"xneg", "xpos"};

// ------------------ x_F binning ------------------
// Fine bins, to be rebinned when A_N is formed
static const int kNXfBins = 40;
static const double kXfMin = 0.0, kXfMax = 1.0;

#endif
