#ifndef JETPARAMETERS_H
#define JETPARAMETERS_H

#include <vector>
#include <cmath>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

// ------------------ Jet clustering ------------------
static const double R = 0.7;

// ------------------ Projection ------------------
static const double z_proj = 721.5402;

// ------------------ MIP values ------------------
static const float ecal_mip = 3.1;
static const float hcal_mip = 15.0;

// ------------------ Threshold ------------------
// Tower hits (reco jets built from hits): minimum energy as a fraction of
// the MIP values above.
static const float mip_threshold = 0.1;

// FCS points (reco jets built from points, -F 1): minimum point energy,
// GeV. Same value as the STAR FMS EM-jet analysis's photon-candidate cut
// (kFmsEngMin, see UnderlyingEventPlayground/fms_emjet_reference.txt),
// which that analysis chose to reject low-energy photons / hadronic
// contributions and to match data better. The MIP-based tower cut above
// does not apply to points. Also the energy floor of reco_nphotons.
// Applied through pass_point_energy_cut() below, together with the upper
// limit.
static const float point_energy_min = 1.0f;

// ------------------ Geometry cut (symmetric two-window detector) ------------------
// HCAL's own projected footprint (north/south combined into one shared
// rectangle -- see below), inner column (nearest the beam pipe) bounding
// the y-range, at the ECAL front plane. Reproducible directly from the
// STAR conditions DB via GeometryPlayground/run_number_to_fiducial.sh
// <run_number> (no simulation needed): stage 1 (get_fcs_geometry.C) pulls
// StFcsDb's detector-position geometry for that run's era, stage 2
// (fiducial_from_corners.cpp, plain C++) reproduces this exact rectangle.
//
// Investigated (Aug 2026) whether north and south actually need separate
// values: they don't, currently -- StFcsDb's detector-position DB table is
// exactly mirror-symmetric (bit-identical at full double precision, both
// before and after the one real position update in Run 22's history), and
// a direct check against real reconstructed hits in a production muDst
// (run 22359013, ~1.9M ECAL hits) found identical extremal hit positions
// AND identical observed column/row extent (both sides fully saturate the
// max 22x34 tower grid) on both arms. So a single shared rectangle is
// correct, not an approximation -- there's no known source of real
// geometric north/south asymmetry to split this on. (Real per-side hit
// *counts* do differ by a few percent in that same check, but that's an
// occupancy/rate effect, not a geometric one, and doesn't belong in a
// position-based fiducial cut.)
//
// Detector position itself is constant across essentially all of Run 22:
// checked at DB level, hour-granularity, from a pre-survey placeholder
// value through 2021-12-20 18:00 UTC, to the real surveyed position from
// 2021-12-21 onward through the y2023 nominal BFC date -- one single
// transition, no drift after. get_fcs_geometry.C warns (but does not
// block) if pointed at a run before that transition; these values are
// only valid from 2021-12-21 onward.
static const float cut_x_inner = 24.15591212f;
static const float cut_x_outer = 131.82035143f;

static const float cut_y_min = -83.19687629f;
static const float cut_y_max =  86.41376856f;

// ------------------ ECAL-only geometry cut (EM-only reco jets) ------------------
// Same derivation as cut_x_inner/etc above (GeometryPlayground/
// run_number_to_fiducial.sh, ECAL front plane, run 22359013's geometry --
// see fcs_corners_run22.txt), but using ECAL's OWN footprint instead of
// being clipped to HCAL's smaller projected footprint: when reco jets are
// built from ECAL hits only (no HCAL), there is no reason to cut them down
// to HCAL's acceptance. ECAL's front face is tilted, not perpendicular to
// z (confirmed from the real corner data -- the inner-column corners sit
// at larger z than the outer-column ones), so these values are ECAL's own
// corners projected onto the shared z_ecal plane like everything else
// here, not their raw (x,y) -- only the outer-column corners happen to
// already sit at z_ecal itself. North/south again confirmed bit-identical
// (same DB-level symmetry as the default rectangle above), so one shared
// rectangle applies here too. Regenerate via
// GeometryPlayground/fiducial_from_corners.cpp's "ECAL-only rectangle"
// output if the detector position ever changes.
static const float cut_x_inner_ecal = 20.53604584f;
static const float cut_x_outer_ecal = 137.59500000f;

static const float cut_y_min_ecal = -96.72439178f;
static const float cut_y_max_ecal =  86.25565167f;

// ------------------ Buffer ------------------
static const float cut_buffer = 0.0f;
static const float R_buffer = 0.0f * R / 4.0f;

// margin for the simple rectangle cut (pass_fiducial_cut): reco jets use a
// margin since the jet energy response is poorly calibrated near the
// fiducial edge. Truth jets use pass_jet_scale_cut instead, extending this
// same boundary outward in eta-phi (see below): by half a jet radius for
// ECAL+HCAL jets, by a full jet radius (truth_fiducial_margin_em) for
// EM-only jets.
static const float reco_fiducial_buffer = 10.0f;

// Same role as reco_fiducial_buffer, but for EM-only reco jets against the
// ECAL-only rectangle above. Characterized via the EM-only JES calibration
// diagnostics fork (JetScalePlayground/): a survey of high-statistics,
// high-chi2 position bins found the poor-fit failure modes concentrated
// within ~20cm of the fiducial edge, so bins inside that margin were
// excluded from the calibration survey entirely -- 20.0f bakes that same
// exclusion into jet-finding itself, rather than only filtering for it
// downstream in the calibration/analysis stage.
static const float reco_fiducial_buffer_ecal = 20.0f;

// ------------------ EM-only reco jets ------------------
// Minimum pT (GeV) of an EM-only reco jet at jet finding. Same value as the
// STAR FMS EM-jet analysis (kPtMin, see
// UnderlyingEventPlayground/fms_emjet_reference.txt). That analysis applies
// the same 2 GeV again after its UE subtraction and energy/pT correction;
// that second, analysis-level cut is not applied here.
static const double reco_jet_pt_min_em = 2.0;

// ------------------ EM-only truth jets ------------------
// Minimum pT (GeV) of an EM-only truth jet at jet finding. Same value as
// the STAR FMS EM-jet analysis's particle-level jets (AntiKtR070Particle,
// see UnderlyingEventPlayground/fms_emjet_reference.txt).
static const double truth_jet_pt_min_em = 1.0;

// How far outside the reco fiducial region (kFiducialRectEcal inset by
// reco_fiducial_buffer_ecal) an EM-only truth jet's axis may sit, in
// eta-phi: one full jet radius. A truth jet further out than R has no
// overlap with any jet whose axis is inside the fiducial region, so it
// cannot affect what is observed there and is dropped.
static const float truth_fiducial_margin_em = (float)R;

// ------------------ Feynman x / detector side ------------------
// Real Run 22 pp sqrt(s), GeV -- confirmed from RunLog.beamInfo (the
// authoritative real-time beam-energy readback DB), investigated
// 2026-10. Run 22's forward-production period actually had TWO distinct
// beam energies, not one constant: blueEnergy=yellowEnergy=254.867 GeV
// for runs <=23002045 (~2021-12-15 to 2022-01-02, 321/3368 forward
// runs / ~8.3% of events) and =254.213 GeV for runs >=23002052
// (2022-01-02 onward, 3047/3368 runs / ~91.7% of events) -- cross-checked
// externally against 2022 IPAC accelerator-physics proceedings, which
// describe RHIC changing the Run 22 collision energy from ~255 to 254.2
// GeV mid-run (snake-magnet damage mitigation). No production runs fall
// between those two run numbers, so there's no real ambiguity right at
// the threshold. kSqrtS_MainRun22_GeV is also what sim_to_jet's
// generator-level SetFrame("CMS",508.4) already uses (fixed previously;
// see starsim_pythia8_filter.C/starsim_pythia6_filter.C) -- sim has no
// real per-event run number, so it always uses this one fixed value
// rather than the per-run lookup below.
//
// The OLD SQRT_S=500.0 constant below is kept only so any code that
// hasn't been updated to pass a real sqrt(s) into computeFeynmanX() yet
// fails loudly (no default parameter) rather than silently reusing a
// value that was wrong by ~1.7-1.9% for every x_F this pipeline has ever
// written (sim and real data both) -- see PRODUCTION_TRACKER.md / the
// 2026-10 x_F investigation for the full writeup. Do not reintroduce a
// default here.
static const double kSqrtS_EarlyRun22_GeV = 509.734; // runs <= 23002045
static const double kSqrtS_MainRun22_GeV  = 508.426; // runs >= 23002052 (and sim's fixed value)
static const int    kSqrtSRunThreshold    = 23002045;

inline double computeSqrtSForRun(int runNumber)
{
    return (runNumber <= kSqrtSRunThreshold) ? kSqrtS_EarlyRun22_GeV : kSqrtS_MainRun22_GeV;
}

// ------------------ Shared functions ------------------

// x_F = 2*pz/sqrt(s). Uses jet energy (not pT) for pz, so that an energy
// correction applied downstream (see ApplyCorrections.cpp) propagates into
// a genuinely different corrected x_F rather than reproducing the
// uncorrected value -- pT and eta are unaffected by the energy-only
// correction, but E*tanh(eta) is not. sqrtS is a required argument (no
// default) -- call with computeSqrtSForRun(runNumber) for real data, or
// kSqrtS_MainRun22_GeV for simulation; see the comment above.
inline double computeFeynmanX(double E, double eta, double sqrtS)
{
    double pz = E * std::tanh(eta);
    return 2.0 * pz / sqrtS;
}

// Energy window for an FCS point to be used (reco jets built from points,
// and reco_nphotons): above point_energy_min and no more than one beam's
// energy, sqrt(s)/2. A single photon cannot carry more than the beam
// energy, so anything above it is unphysical (hot/miscalibrated channels,
// beam background) -- the same role as the STAR FMS EM-jet analysis's
// upper photon-energy cut (kFmsEngMax = 100 GeV at sqrt(s) = 200 GeV).
// Shared by sim_to_jet and data_to_jet so both apply exactly the same
// code; it removes nothing physical in simulation. Pass
// computeSqrtSForRun(runNumber) for real data, kSqrtS_MainRun22_GeV for
// simulation.
inline bool pass_point_energy_cut(float e, double sqrtS)
{
    return e > point_energy_min && e <= 0.5 * sqrtS;
}

// Rejects unphysical reco jets, as the STAR FMS EM-jet analysis does: a
// jet cannot carry more than one beam's energy (E <= sqrt(s)/2), and its
// x_F must lie in (0, 1]. ANALYSIS-LEVEL ONLY: these limits are only
// meaningful for the final jet energy, i.e. after underlying-event
// subtraction and energy calibration, so nothing in sim_to_jet or
// data_to_jet calls this -- it lives here so the analysis has one
// definition to use. Pass the corrected E.
inline bool pass_physical_jet_cut(double E, double eta, double sqrtS)
{
    double xF = computeFeynmanX(E, eta, sqrtS);
    return E <= 0.5 * sqrtS && xF > 0.0 && xF <= 1.0;
}

// Detector side from the ECAL-front-plane x position (same north/south
// split convention as FcsJetFilter): +1 for x>0, -1 for x<=0.
inline int computeSide(double x)
{
    return (x > 0.0) ? 1 : -1;
}

// Projects a detector-plane (x, y) position at z_proj into (eta, phi)
inline void computeEtaPhi(double x, double y, double &eta, double &phi)
{
    double r = std::sqrt(x*x + y*y + z_proj*z_proj);

    double px = x / r;
    double py = y / r;
    double pz = z_proj / r;

    eta = std::asinh(pz / std::sqrt(px*px + py*py));
    phi = std::atan2(py, px);
}

inline double computeTauN(const std::vector<fastjet::PseudoJet> &constituents, int N, double R)
{
    if (constituents.empty()) return -1;
    if (constituents.size() < (size_t)N) return 0.;

    fastjet::ClusterSequence cs(
        constituents,
        fastjet::JetDefinition(fastjet::kt_algorithm, R)
    );

    std::vector<fastjet::PseudoJet> subjets = cs.exclusive_jets(N);

    double num = 0., den = 0.;

    for (const auto &c : constituents) {
        double minDR = 1e6;

        for (const auto &sj : subjets) {
            double dphi = std::fabs(c.phi() - sj.phi());
            if (dphi > M_PI) dphi = 2*M_PI - dphi;

            double deta = c.eta() - sj.eta();
            double dr = std::sqrt(deta*deta + dphi*dphi);

            if (dr < minDR) minDR = dr;
        }

        num += c.perp() * minDR;
        den += c.perp();
    }

    return (den > 0) ? num / (den * R) : -1;
}

// minimum eta-phi distance from (jetXE, jetYE) to the boundary of the
// rectangle [x_inner, x_outer] x [y_min, y_max] (sampled along the edges),
// signed positive if the point is inside the rectangle and negative if it
// is outside
inline double etaPhiSignedDistance(float jetXE, float jetYE,
                                    double x_inner, double x_outer,
                                    double y_min, double y_max)
{
    // Reflect to the right detector
    double x = std::fabs(jetXE);
    double y = jetYE;

    bool inside = (x >= x_inner && x <= x_outer && y >= y_min && y <= y_max);

    double r = std::sqrt(x*x + y*y + z_proj*z_proj);

    double px = x / r;
    double py = y / r;
    double pz = z_proj / r;

    double jetEta = std::asinh(pz / std::sqrt(px*px + py*py));
    double jetPhi = std::atan2(py, px);

    auto deltaPhi = [](double a, double b)
    {
        double d = a - b;
        while (d >  M_PI) d -= 2.0*M_PI;
        while (d < -M_PI) d += 2.0*M_PI;
        return d;
    };

    auto updateDistance = [&](double bx, double by, double &best)
    {
        double br = std::sqrt(bx*bx + by*by + z_proj*z_proj);

        double bpx = bx / br;
        double bpy = by / br;
        double bpz = z_proj / br;

        double eta = std::asinh(bpz / std::sqrt(bpx*bpx + bpy*bpy));
        double phi = std::atan2(bpy, bpx);

        double deta = jetEta - eta;
        double dphi = deltaPhi(jetPhi, phi);

        double dr = std::sqrt(deta*deta + dphi*dphi);

        if (dr < best)
            best = dr;
    };

    const int nSamples = 200;
    double best = 1e30;

    // Bottom and top
    for (int i = 0; i <= nSamples; ++i)
    {
        double t = double(i) / nSamples;
        double bx = x_inner + t * (x_outer - x_inner);

        updateDistance(bx, y_min, best);
        updateDistance(bx, y_max, best);
    }

    // Left and right
    for (int i = 0; i <= nSamples; ++i)
    {
        double t = double(i) / nSamples;
        double by = y_min + t * (y_max - y_min);

        updateDistance(x_inner, by, best);
        updateDistance(x_outer, by, best);
    }

    return inside ? best : -best;
}

// A fiducial rectangle's four bounds, so pass_jet_scale_cut/pass_fiducial_cut
// below can be pointed at either the default (ECAL+HCAL-limited) rectangle
// or the ECAL-only one (EM-only reco jets), instead of always hardcoding
// cut_x_inner/etc.
struct FiducialRect { float x_inner, x_outer, y_min, y_max; };

static const FiducialRect kFiducialRect     = { cut_x_inner,      cut_x_outer,      cut_y_min,      cut_y_max };
static const FiducialRect kFiducialRectEcal = { cut_x_inner_ecal, cut_x_outer_ecal, cut_y_min_ecal, cut_y_max_ecal };

// symmetric two-window acceptance, with an eta-phi margin relative to a
// rectangle boundary (rect inset by cmBuffer). A positive etaPhiMargin
// erodes the acceptance region: the jet must clear the boundary by that
// much in eta-phi (used to define where the jet-energy-scale calibration
// map is valid). A negative etaPhiMargin dilates it: the jet may fall
// outside the boundary by up to that much in eta-phi (used to accept truth
// jets whose reco counterpart can still land inside the reco fiducial cut).
// rect defaults to the standard ECAL+HCAL rectangle; pass kFiducialRectEcal
// for EM-only reco jets.
inline bool pass_jet_scale_cut(float jetXE, float jetYE, float cmBuffer, float etaPhiMargin,
                                const FiducialRect &rect = kFiducialRect)
{
    double x_inner = rect.x_inner + cmBuffer;
    double x_outer = rect.x_outer - cmBuffer;
    double y_min   = rect.y_min   + cmBuffer;
    double y_max   = rect.y_max   - cmBuffer;

    return etaPhiSignedDistance(jetXE, jetYE, x_inner, x_outer, y_min, y_max) >= etaPhiMargin;
}


// symmetric two-window acceptance at ECAL front plane, used for reco jet
// selection (with an adjustable margin). rect defaults to the standard
// ECAL+HCAL rectangle; pass kFiducialRectEcal for EM-only reco jets.
inline bool pass_fiducial_cut(float jetXE, float jetYE, float buffer,
                               const FiducialRect &rect = kFiducialRect)
{
    const float ax = std::fabs(jetXE);
    const bool in_x = ax >= (rect.x_inner + buffer) && ax <= (rect.x_outer - buffer);
    const bool in_y = jetYE >= (rect.y_min + buffer) && jetYE <= (rect.y_max - buffer);
    return in_x && in_y;
}

#endif