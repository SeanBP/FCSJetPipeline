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
// EMTHR with ECal 4x4 > 0.25 x HCal 4x4, ORed north/south at stage 3. The
// emulator code differs by version (StRoot/RTS/src/TRG_FCS): stage_2_202203
// compares ECal + HCal with EMTHR when HCal <= ECal/4; stage_2_202207 (EM0-2
// naming, as the Run 22 trigger ids) compares the ECal 4x4 alone with EMTHR,
// with the largest neighbouring HCal 4x4 <= ECal/4 as a veto. Which one the
// Run 22 firmware ran is not in the code; 0.0316
// GeV E_T per count, EM0/1/2 = 64/112/160 counts = 2.0/3.5/5.0 GeV from run
// 22359013 (63/111/158 before; stage_params.txt still has those). The online
// E_T scale moved with the ET gain files (new at 23027047, removing an old
// 0.7 factor) and the ECal electronics gain (halved at 23066056), so the
// offline turn-ons can differ between those periods. Earlier rough estimate
// on run 23007048: uncorrected point-jet pT turns on near 2.6 (EM1) and 3.7
// GeV (EM2), about 0.74 of the nominal 3.5 and 5.0.
static const double kEmThresholdGeV[kNEmLevels] = {2.0, 3.5, 5.0};

// ------------------ Trigger levels: EM and JP (from 2026-11-05) ------------------
// A jet "fired" a trigger level when the quantity that trigger compared with
// its threshold, rebuilt for the jet from offline hit E_T (RecoJets
// reco_trg_*: shared/FcsTrigSums.h), is above CutFactor() x the nominal
// threshold (A. Ogawa's Run 22 trigger page, the same value for every run and
// in simulation), for a level listed in the event:
//   EM0-2: reco_trg_em_et, the jet's largest ECal 4x4 window, and the stage-2
//          HCal veto with a margin: reco_trg_em_hmax <= reco_trg_em_et /
//          (4 x kEmVetoMargin). The Run 22 firmware ran the 202207 logic
//          (ECal-only window, veto when the largest neighbouring HCal 4x4 >
//          ECal / 4): jets failing the veto fire 0-57%, veto-safe jets above
//          1.3 x threshold 97-100% (one file per gain period, 2026-11-05).
//   JP0-2: reco_trg_jp_et, the largest JP patch (ECal + HCal) containing the
//          jet axis; JP0 and JP1 are per region (fcsJPA*, fcsJPBC*, fcsJPDE*:
//          patch A, B or C, D or E, either side) and need a listed region
//          among the jet's patches (reco_trg_jp_mask); fcsJP2 is any patch.
// The two families are judged independently: a jet has an EM level (the
// lowest listed EM level it fires) and, separately, a JP level (the lowest
// listed JP level it fires), and fills both, so the A_N of one level (e.g.
// JP0) has every jet that fired it, whatever the jet did in the other
// family. Within a family a jet is in one level only, so the EM levels add
// up, and so do the JP levels. For an EM + JP union without double counting,
// the JP histograms are split by the same jet's EM result (kEmTagSuffix:
// _noem, no EM level fired, veto failures included; _emt115 ... _emt160, the
// EM slice it reached): the union at factor f is the EM levels plus the JP
// jets whose EM slice is below f (AnAsymmetry -J). Jets above an EM
// threshold that fail only the veto go to their own histograms (kVetoSuffix). Default
// factor 1.3 for both families; see kEmCutFactor / kJpCutFactor. The yields are
// kept per slice of trigger value / threshold (kTrgSliceEdge, suffixes
// kTrgSliceSuffix), so AnAsymmetry -f can choose another factor afterwards.
static const int kNTrgLevels = 6;
static const char *const kTrgLevelName[kNTrgLevels] = {"EM0", "EM1", "EM2", "JP0", "JP1", "JP2"};
static const double kTrgThresholdGeV[kNTrgLevels] = {2.0, 3.5, 5.0, 4.0, 6.0, 8.0};
inline bool IsJpLevel(int l) { return l >= kNEmLevels; }
// Default cut factor per family: 1.3 for both (Sean 2026-11-05, for the 20%
// pass; EM was 1.15 for a few hours). At 1.3 the turn-on loss
// (trg_safety.txt), jet-weighted over each gain period, is <= 0.4% for every
// trigger; single run groups reach 1.5% (EM1, 23028001-23034017) and 1.1%
// (EM0, 23035007-23037065), both wide turn-on fits from ~200 reference
// events; JP <= 0.2% in every group. A sharp cut in simulation is then good
// to about 1% for EM; 1.15 stays available (AnAsymmetry -fe 1.15). Run 22 scan (125 runs, TriggerPlayground/trg_eff_scan.cpp):
// hardware-bit efficiency above 1.15 x threshold EM0 0.958, EM1 0.904, EM2
// 0.98; JP0 0.976 with the slice 1.15-1.3 only 0.80-0.95 between 23027047
// and 23066055, while JP >= 1.3 is 0.999 in every period.
// Runs whose JP triggers are not used (EM unaffected): the JP levels count as
// not listed there. 23070019: every JP bit (JP0, JP1, JP2) set on every event,
// even far below threshold (TriggerPlayground/trg_eff_scan.cpp hot-run test;
// excluded by Sean 2026-11-05).
static const int kNJpBadRuns = 1;
static const int kJpBadRuns[kNJpBadRuns] = {23070019};
inline bool IsJpBadRun(int run)
{
    for (int i = 0; i < kNJpBadRuns; i++) if (kJpBadRuns[i] == run) return true;
    return false;
}
static const double kEmCutFactor = 1.3;
static const double kJpCutFactor = 1.3;
inline double CutFactor(int level) { return IsJpLevel(level) ? kJpCutFactor : kEmCutFactor; }
// EM HCal-veto margin: a jet counts only if its largest neighbouring HCal 4x4
// <= ECal 4x4 / (4 x margin). 1.3 since 2026-11-05 (Sean: same number as the
// cut factor, uniform): EM above 1.3 x threshold then fires 96.7% (EM0) /
// 97.4% (EM1) before 23027047 and >= 99.6% later; 1.15 gave 90% early (the
// online HCal/ECal ratio was ~1.4-1.5 x the offline one before 23027047;
// 1.5 would give >= 99.7% everywhere). TriggerPlayground/trg_eff_scan.cpp
// trg_veto_periods.png. Not stored in slices: changing it reruns AnHistograms
// on the kept jet trees (submit_rehist.sh).
static const double kEmVetoMargin = 1.3;
static const int kNTrgSlices = 4;
static const double kTrgSliceEdge[kNTrgSlices] = {1.15, 1.30, 1.45, 1.60};  // lower edges; the last runs to infinity
static const char *const kTrgSliceSuffix[kNTrgSlices] = {"_t115", "_t130", "_t145", "_t160"};
static const char *const kVetoSuffix = "_veto";
static const int kNEmTags = kNTrgSlices + 1;  // 0 no EM, 1 + g: EM slice g
static const char *const kEmTagSuffix[kNEmTags] = {"_noem", "_emt115", "_emt130", "_emt145", "_emt160"};

// A jet counts only if its uncorrected pT is above this factor times the
// threshold of an EM trigger listed in its event (above the turn-on); 1.15 as
// the FMS EM-jet analysis (UnderlyingEventPlayground/fms_emjet_reference.txt)
static const double kTriggerPtFactor = 1.15;

// Corrected jet pT minimum (GeV), after UE subtraction and calibration, as
// the FMS analysis; the same value as the reco cut at jet finding
static const double kJetPtCorrMin = 2.0;

// ------------------ Point multiplicity ------------------
// Photons in the jet = its point constituents (reco_nconst of point jets, FCS
// points above point_energy_min), the FMS EM-jet definition (BrightSTAR:
// jet->GetNumberOfTowers(), whose "towers" are FMS points); classes:
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
// Fine bins of 0.025, to be rebinned when A_N is formed. The range reaches
// past 1 because jets with corrected E above sqrt(s)/2 are kept (in their
// own _hiE histograms, see below), and x_F <= 2E/sqrt(s), so every jet with
// x_F > 1 is one of them. Yield files made before 2026-11-04 have 40 bins
// from 0 to 1 and cannot be summed with these.
static const int kNXfBins = 60;
static const double kXfMin = 0.0, kXfMax = 1.5;

// ------------------ Jets above sqrt(s)/2 ------------------
// The FMS EM-jet analysis dropped jets with E > sqrt(s)/2 (odd behaviour
// there). AnHistograms keeps them apart instead, in h_xF_*_hiE, h_abscos_*_hiE
// and h_xfsum_*_hiE (the same binning), and the default yields hold only
// E <= sqrt(s)/2. AnAsymmetry -E adds them in; the default leaves them out,
// as the FMS. Decide from the QA (AnPlots: qa_jets, yields_hiE) whether they
// behave.
static const char *const kHiESuffix = "_hiE";
// ------------------ Companion jets ------------------
// Jets passing every cut but the trigger threshold (kTriggerPtFactor) in an
// event where another jet passes it. They did not fire the trigger, so they
// are kept apart: the same histograms with this suffix (before _hiE), at the
// highest level among the event's passing jets. AnAsymmetry -C forms their
// A_N alone, -A adds them to the triggering jets; the default leaves them out.
static const char *const kCompSuffix = "_comp";
// ------------------ Two-photon jets ------------------
// Invariant mass of the two point constituents of a two-photon jet: QA (qa_mgg_*) and x_F vs mass yields (h_xFm_*) for a pi0
// window chosen later. 25 MeV bins; the data pi0 peak sits 8-17% below sim
// (point energy scale).
static const int kNMassBins = 40;
static const double kMassMin = 0.0, kMassMax = 1.0;

#endif
