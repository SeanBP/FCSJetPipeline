#ifndef JETPARAMETERS_H
#define JETPARAMETERS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <string>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"

// ------------------ Jet clustering ------------------
static const double R = 0.7;

// Jet areas, as the FMS EM-jet analysis (STAR StJetFinder defaults with ghost
// rapidity 5): active area, ghost area 0.04, |ghost rapidity| < 5, one
// repetition, grid scatter 1, pT scatter 0.1, mean ghost pT 1e-100. The UE
// subtracted from a jet is the UE density (mean cone energy / pi R^2) times
// the jet's area (underlying_event/UnderlyingEvent.cpp).
inline fastjet::AreaDefinition JetAreaDefinition()
{
    return fastjet::AreaDefinition(fastjet::active_area, fastjet::GhostedAreaSpec(5.0, 1, 0.04, 1.0, 0.1, 1e-100));
}

// ------------------ Projection ------------------
static const double z_proj = 721.5402;

// ------------------ MIP values ------------------
static const float ecal_mip = 3.1;
static const float hcal_mip = 15.0;

// ------------------ Thresholds ------------------
// Tower hits: minimum energy as a fraction of the MIP values above
static const float mip_threshold = 0.1;

// FCS points (-F 1): minimum point energy, GeV, as the FMS EM-jet analysis
// (kFmsEngMin, UnderlyingEventPlayground/fms_emjet_reference.txt). Also the
// energy floor of the points counted as photons. Applied via
// pass_point_energy_cut() (no upper cut on points; see there).
static const float point_energy_min = 1.0f;

// ------------------ Geometry cut (symmetric two-window detector) ------------------
// HCAL's projected footprint at the ECAL front plane, one rectangle shared by
// north and south (the DB geometry is mirror-symmetric). Valid for runs from
// 2021-12-21 onward (surveyed detector position). Regenerate with
// GeometryPlayground/run_number_to_fiducial.sh <run_number>.
static const float cut_x_inner = 24.15591212f;
static const float cut_x_outer = 131.82035143f;

static const float cut_y_min = -83.19687629f;
static const float cut_y_max =  86.41376856f;

// ------------------ ECAL-only geometry cut (EM-only reco jets) ------------------
// ECAL's own footprint, its tilted front-face corners projected onto the
// z_proj plane. Regenerate with GeometryPlayground/fiducial_from_corners.cpp
// ("ECAL-only rectangle").
static const float cut_x_inner_ecal = 20.53604584f;
static const float cut_x_outer_ecal = 137.59500000f;

static const float cut_y_min_ecal = -96.72439178f;
static const float cut_y_max_ecal =  86.25565167f;

// ------------------ Buffer ------------------
static const float cut_buffer = 0.0f;
static const float R_buffer = 0.0f * R / 4.0f;

// Reco jet margin inside the rectangle (pass_fiducial_cut), cm: the energy
// response is poorly calibrated near the edge. Truth jets use
// pass_jet_scale_cut, extending the boundary outward in eta-phi instead.
static const float reco_fiducial_buffer = 10.0f;

// Same, for EM-only reco jets against the ECAL-only rectangle. JES fits
// degrade within ~20 cm of the edge.
static const float reco_fiducial_buffer_ecal = 20.0f;

// ------------------ EM-only reco jets ------------------
// Minimum pT (GeV) of an EM-only reco jet at jet finding, as the FMS EM-jet
// analysis (kPtMin). Its second, post-correction 2 GeV cut is analysis level.
static const double reco_jet_pt_min_em = 2.0;

// ------------------ EM-only truth jets ------------------
// Minimum pT (GeV) of an EM-only truth jet, as the FMS particle-level jets
static const double truth_jet_pt_min_em = 1.0;

// EM-only truth acceptance. Truth jets: only those that could be matched to a
// fiducial reco jet (dR < R/2), i.e. axis within R/2 in eta-phi of the reco
// fiducial region (pass_jet_scale_cut with -truth_jet_margin_em). Truth
// particles: a full-azimuth eta band around the eta range of that region
// (EmFiducialEtaRange() below), 1.5R wide, so every particle of such a jet,
// and of its off-axis UE cones at the same eta, is kept.
static const float truth_jet_margin_em = 0.5f * (float)R;
static const double truth_particle_eta_margin_em = 1.5 * R;

// EM-only truth particles: photons created within this distance (cm) of the
// primary vertex, i.e. prompt and pi0/eta-decay photons but not GEANT shower
// secondaries. Shared by the truth jets (sim_to_jet/JetMatcher.cpp) and the
// truth-level UE (underlying_event/UnderlyingEvent.cpp).
static const double em_truth_displacement_max = 0.01;

// ------------------ Full (ECAL+HCAL) jets ------------------
// Analogous to the EM-only cuts above (Sean 2026-11-16), on the region both
// calorimeters cover (kFiducialRectFull below: the HCAL footprint at the ECAL
// front plane, cut_* above, clipped to the ECAL rectangle; only y_max moves,
// by 0.16 cm).
// Reco jet margin inside that rectangle, cm: as EM-only's 20 cm until the
// full-jet JES test shows where the fits degrade (was 10 cm on the HCAL
// rectangle).
static const float reco_fiducial_buffer_full = 20.0f;
// Truth jets: axis within R/2 in eta-phi of the reco fiducial region, as EM-only.
static const float truth_jet_margin_full = 0.5f * (float)R;
// Truth particles: a full-azimuth eta band 1.5R wide around that region's eta
// range, as EM-only (FiducialEtaRange() below).
static const double truth_particle_eta_margin_full = 1.5 * R;

// Geant3 particle code of the photon
static const int kGeant3Photon = 1;

// The primary vertex: the creation point of the first particle made at
// StMcVertex 1 (idVtx == 1). False if there is none.
inline bool FindEmTruthPrimaryVertex(int n, const int *idVtx, const float *vx, const float *vy, const float *vz,
                                     float &x, float &y, float &z)
{
    for (int i = 0; i < n; i++)
        if (idVtx[i] == 1) { x = vx[i]; y = vy[i]; z = vz[i]; return true; }
    return false;
}

// Whether a particle enters the EM-only truth: a photon created within
// em_truth_displacement_max of the primary vertex (px, py, pz)
inline bool IsEmTruthParticle(int geid, float vx, float vy, float vz, float px, float py, float pz)
{
    if (geid != kGeant3Photon) return false;
    const double dx = vx - px, dy = vy - py, dz = vz - pz;
    return std::sqrt(dx * dx + dy * dy + dz * dz) <= em_truth_displacement_max;
}

// ------------------ Feynman x / detector side ------------------
// Run 22 pp sqrt(s), GeV, from RunLog.beamInfo: RHIC lowered the beam energy
// from 254.867 to 254.213 GeV between runs 23002045 and 23002052.
// Simulation uses kSqrtS_MainRun22_GeV (as the generator frame does).
static const double kSqrtS_EarlyRun22_GeV = 509.734; // runs <= 23002045
static const double kSqrtS_MainRun22_GeV  = 508.426; // runs >= 23002052 (and sim's fixed value)
static const int    kSqrtSRunThreshold    = 23002045;

// ------------------ JES evaluation range ------------------
// The EM-only log-cubic JES (ApplyCorrections_em_only, jes_closure) is
// evaluated at the jet's UE-subtracted energy held at kJesEvalMinGeV from
// below and at sqrt(s)/2 from above (Sean 2026-11-06). 2 GeV is about the
// curve's minimum on jesmap_pixring_v1 (median scale 0.797; 0.801 at the
// map's lowest fitted response, 1.29 GeV), so the hold changes nothing
// measurable; analysis jets are above ~12 GeV.
static const double kJesEvalMinGeV = 2.0;

inline double computeSqrtSForRun(int runNumber)
{
    return (runNumber <= kSqrtSRunThreshold) ? kSqrtS_EarlyRun22_GeV : kSqrtS_MainRun22_GeV;
}

// ------------------ Shared functions ------------------

// x_F = 2 pz / sqrt(s), with pz = E tanh(eta) so an energy correction
// propagates into x_F. Pass computeSqrtSForRun(run) for data,
// kSqrtS_MainRun22_GeV for simulation.
inline double computeFeynmanX(double E, double eta, double sqrtS)
{
    double pz = E * std::tanh(eta);
    return 2.0 * pz / sqrtS;
}

// FCS point energy cut: above point_energy_min. Shared by sim_to_jet,
// data_to_jet and underlying_event. Until 2026-11-06 points above sqrt(s)/2
// were dropped here as well; that silently rebuilt a jet with a hot channel
// from its remaining points. Now such a jet keeps its energy and is handled as
// a whole by the jet-level cut (pass_physical_jet_cut; AnHistograms keeps jets
// above sqrt(s)/2 apart, AnAsymmetry -E adds them) (Sean 2026-11-06).
inline bool pass_point_energy_cut(float e)
{
    return e > point_energy_min;
}

// Analysis-level cut on the corrected jet (after UE subtraction and energy
// calibration), as the FMS analysis: E <= sqrt(s)/2 and 0 < x_F <= 1.
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

// Run 22 beamline (Sean 2026-11-16): the TPC primary-vertex beamline,
// event-weighted over 40 runs (P25ib st_physics; vertex_tpc_v1/code/
// beamline_runs.C): x = kBeamX0 + kBeamDxDz z, y = kBeamY0 + kBeamDyDz z (cm).
// The EPD / BBC TAC give only z; the event vertex is (beamX(z), beamY(z), z),
// and with no vertex the nominal IP (beamX, beamY at z = -7.68 cm, the TPC
// beam-spot mean; StSimpleReaderMaker::mNominalVz). The simulation generates
// its vertices on the same beamline (starsim SetVertex / SetSlope).
static const double kBeamX0 = -0.061, kBeamY0 = -0.203;
static const double kBeamDxDz = -1.04e-3, kBeamDyDz = 1.95e-4;
inline double beamX(double vz) { return kBeamX0 + kBeamDxDz * vz; }
inline double beamY(double vz) { return kBeamY0 + kBeamDyDz * vz; }

// The same seen from a vertex (vx, vy, vz) (cm): the direction of the
// detector point (x, y, z_proj) from the vertex
inline void computeEtaPhi(double x, double y, double vx, double vy, double vz, double &eta, double &phi)
{
    const double dx = x - vx, dy = y - vy, dz = z_proj - vz;
    eta = std::asinh(dz / std::sqrt(dx*dx + dy*dy));
    phi = std::atan2(dy, dx);
}
// ... from a vertex at z = vz on the z axis (x = y = 0)
inline void computeEtaPhi(double x, double y, double vz, double &eta, double &phi)
{
    computeEtaPhi(x, y, 0.0, 0.0, vz, eta, phi);
}

// Event vertex in reco jet kinematics (2026-11-16): on by default (Sean,
// 2026-11-16); FCS_JET_VERTEX=0 in the environment turns it off (origin-based
// kinematics, as before). RecoJets / JetMatcher take reco point directions from the
// event vertex (beamX(z), beamY(z), z) with z = event_vertex_z (the nominal z
// when no method has one) and project jets back to the detector plane from
// it, so reco_x / reco_y stay detector positions; truth uses the generated
// vertex (mc_vx, mc_vy, mc_vz). Note: the JES map in use (jesmap_pixring_v1)
// was built with origin-based reco jets; the next JES production uses this.
inline bool useJetVertex()
{
    const char* v = std::getenv("FCS_JET_VERTEX");
    return !(v && std::string(v) == "0");
}

// Event vertex cut (Sean 2026-11-16): |event_vertex_z| < kVertexZCut on the
// chosen vertex as measured (StSimpleReaderMaker: EPD, else BBC TAC, else
// the nominal z, which passes; no rescaling to the TPC scale, the methods
// read ~0.8 z_TPC, so 100 cm is ~125 cm true). 100 cm keeps 91.2% of st_fwd
// FCS-triggered events (runs from 23050000), close to the share FMS's
// |vz| < 80 cm kept in Run 15 (93.4%: 76% of their events had no vertex and
// sat at z = 0). The simulation emulates the same response, so the cut acts
// alike in both.
// Trees from before the 2026-11-16 vertex code (no epd_source / mc_vz branch:
// event_vertex_source meant TPC / VPD / old BBC then) are not cut; callers
// apply it only to trees with an epd_source (data) or mc_vz (simulation) branch.
static const double kVertexZCut = 100.0;
inline bool pass_vertex_cut(float event_vertex_z, int event_vertex_source)
{
    if (event_vertex_source <= 0) return true; // not set: simulation without the emulation
    return std::fabs(event_vertex_z) < kVertexZCut;
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

// A fiducial rectangle's bounds (|x| for both detector halves)
struct FiducialRect { float x_inner, x_outer, y_min, y_max; };

static const FiducialRect kFiducialRect     = { cut_x_inner,      cut_x_outer,      cut_y_min,      cut_y_max };
static const FiducialRect kFiducialRectEcal = { cut_x_inner_ecal, cut_x_outer_ecal, cut_y_min_ecal, cut_y_max_ecal };
// Full (ECAL+HCAL) jets: the HCAL rectangle clipped to the ECAL rectangle
static const FiducialRect kFiducialRectFull = {
    cut_x_inner > cut_x_inner_ecal ? cut_x_inner : cut_x_inner_ecal,
    cut_x_outer < cut_x_outer_ecal ? cut_x_outer : cut_x_outer_ecal,
    cut_y_min > cut_y_min_ecal ? cut_y_min : cut_y_min_ecal,
    cut_y_max < cut_y_max_ecal ? cut_y_max : cut_y_max_ecal};

// The reco fiducial rectangle, inset buffer (cm), truth-jet eta-phi margin and
// truth-particle eta margin of a jet type (em_only: ECAL-only; else ECAL+HCAL)
inline const FiducialRect &JetFiducialRect(bool em_only) { return em_only ? kFiducialRectEcal : kFiducialRectFull; }
inline float JetFiducialBuffer(bool em_only) { return em_only ? reco_fiducial_buffer_ecal : reco_fiducial_buffer_full; }
inline float TruthJetMargin(bool em_only) { return em_only ? truth_jet_margin_em : truth_jet_margin_full; }
inline double TruthParticleEtaMargin(bool em_only) { return em_only ? truth_particle_eta_margin_em : truth_particle_eta_margin_full; }

// Acceptance with an eta-phi margin from the boundary of rect inset by
// cmBuffer. A positive margin requires the jet to clear the boundary by that
// much (JES map validity); a negative one accepts jets up to that far outside
// (truth jets whose reco partner can still land inside).
inline bool pass_jet_scale_cut(float jetXE, float jetYE, float cmBuffer, float etaPhiMargin,
                                const FiducialRect &rect = kFiducialRect)
{
    double x_inner = rect.x_inner + cmBuffer;
    double x_outer = rect.x_outer - cmBuffer;
    double y_min   = rect.y_min   + cmBuffer;
    double y_max   = rect.y_max   - cmBuffer;

    return etaPhiSignedDistance(jetXE, jetYE, x_inner, x_outer, y_min, y_max) >= etaPhiMargin;
}


// Reco jet acceptance: inside rect inset by buffer (cm) at the ECAL front plane
inline bool pass_fiducial_cut(float jetXE, float jetYE, float buffer,
                               const FiducialRect &rect = kFiducialRect)
{
    const float ax = std::fabs(jetXE);
    const bool in_x = ax >= (rect.x_inner + buffer) && ax <= (rect.x_outer - buffer);
    const bool in_y = jetYE >= (rect.y_min + buffer) && jetYE <= (rect.y_max - buffer);
    return in_x && in_y;
}

// Eta range of a jet type's reco fiducial region (its rectangle minus its
// buffer), from its nearest and farthest points to the beam
inline void FiducialEtaRange(bool em_only, double &eta_lo, double &eta_hi)
{
    const FiducialRect &r = JetFiducialRect(em_only);
    const double b = JetFiducialBuffer(em_only);
    const double x_lo = r.x_inner + b, x_hi = r.x_outer - b, y_lo = r.y_min + b, y_hi = r.y_max - b;
    const double y_near = std::min(std::max(0.0, y_lo), y_hi);
    const double r_min = std::sqrt(x_lo * x_lo + y_near * y_near);
    const double r_max = std::sqrt(x_hi * x_hi + std::max(y_lo * y_lo, y_hi * y_hi));
    eta_lo = std::asinh(z_proj / r_max);
    eta_hi = std::asinh(z_proj / r_min);
}

inline void EmFiducialEtaRange(double &eta_lo, double &eta_hi) { FiducialEtaRange(true, eta_lo, eta_hi); }

// Truth particles kept for a jet type's truth jets (and truth UE)
inline bool InTruthParticleEtaRange(double eta, bool em_only = true)
{
    double lo, hi;
    FiducialEtaRange(em_only, lo, hi);
    const double m = TruthParticleEtaMargin(em_only);
    return eta >= lo - m && eta <= hi + m;
}

#endif