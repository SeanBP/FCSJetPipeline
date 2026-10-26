#ifndef AN_CONFIG_H
#define AN_CONFIG_H

// Selection and binning of the A_N yield histograms (AnHistograms.cpp).

#include <cmath>

// ------------------ EM trigger levels ------------------
// Highest EM trigger fired, as the underlying-event step (ProxyEmLevel in
// underlying_event/ProxySelection.h): fcsEM2 > fcsEM1 > fcsEM0. With
// prescales a lower flag can be missing, so the highest flag present counts.
// fcsEM3 (runs 23101059-23108024 only) and the fcsEM*_tpc copies are not
// used; no event of run 23007048 fired only a _tpc copy.
static const int kNEmLevels = 3;
static const char *const kEmLevelName[kNEmLevels] = {"EM0", "EM1", "EM2"};

// EM trigger thresholds, GeV of jet pT. PLACEHOLDER, to be confirmed: the
// Run 22 stage-2 thresholds are 63/111/158 counts (TriggerPlayground/
// stage_params.txt), and the leading uncalibrated point-jet pT of run
// 23007048 turns on near 2.6 GeV for EM1 and 3.7 GeV for EM2, about 0.023
// GeV per count; EM0 is scaled from those.
static const double kEmThresholdGeV[kNEmLevels] = {1.5, 2.6, 3.7};

// A jet counts only if its corrected pT is above this fraction of its
// event's highest EM trigger threshold (above the trigger turn-on)
static const double kTriggerPtFactor = 1.1;

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
// Spin_config is StSpinDbMaker::spin4usingBX7(): blue up = {4,5,6}, blue
// down = {8,9,10}, yellow up = {1,5,9}, yellow down = {2,6,10}. Only the four
// both-polarized states enter; -1 (no spin DB entry) and the rest are dropped.
static const int kNSpinStates = 4;
static const char *const kSpinStateName[kNSpinStates] = {"BuYu", "BuYd", "BdYu", "BdYd"};
inline int SpinState(int spin_config)
{
    switch (spin_config)
    {
        case 5: return 0;
        case 6: return 1;
        case 9: return 2;
        case 10: return 3;
        default: return -1;
    }
}

// ------------------ Detector side ------------------
// computeSide() (JetParameters.h) on reco_x: -1 for x <= 0, +1 for x > 0.
// Left/right relative to each beam is assigned when A_N is formed.
static const int kNSides = 2;
static const char *const kSideName[kNSides] = {"xneg", "xpos"};

// ------------------ x_F binning ------------------
// Fine bins, to be rebinned when A_N is formed
static const int kNXfBins = 40;
static const double kXfMin = 0.0, kXfMax = 1.0;

#endif
