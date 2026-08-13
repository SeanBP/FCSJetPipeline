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
static const float mip_threshold = 0.1;

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

// ------------------ Buffer ------------------
static const float cut_buffer = 0.0f;
static const float R_buffer = 0.0f * R / 4.0f;

// margin for the simple rectangle cut (pass_fiducial_cut): reco jets use a
// margin since the jet energy response is poorly calibrated near the
// fiducial edge. Truth jets use pass_jet_scale_cut instead, extending this
// same boundary outward by half a jet radius in eta-phi (see below).
static const float reco_fiducial_buffer = 10.0f;

// ------------------ Feynman x / detector side ------------------
static const double SQRT_S = 500.0; // GeV, pp sqrt(s) for this dataset

// ------------------ Shared functions ------------------

// x_F = 2*pz/sqrt(s). Uses jet energy (not pT) for pz, so that an energy
// correction applied downstream (see ApplyCorrections.cpp) propagates into
// a genuinely different corrected x_F rather than reproducing the
// uncorrected value -- pT and eta are unaffected by the energy-only
// correction, but E*tanh(eta) is not.
inline double computeFeynmanX(double E, double eta)
{
    double pz = E * std::tanh(eta);
    return 2.0 * pz / SQRT_S;
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

// symmetric two-window acceptance, with an eta-phi margin relative to a
// rectangle boundary (the nominal fiducial rectangle inset by cmBuffer). A
// positive etaPhiMargin erodes the acceptance region: the jet must clear the
// boundary by that much in eta-phi (used to define where the jet-energy-scale
// calibration map is valid). A negative etaPhiMargin dilates it: the jet may
// fall outside the boundary by up to that much in eta-phi (used to accept
// truth jets whose reco counterpart can still land inside the reco fiducial
// cut).
inline bool pass_jet_scale_cut(float jetXE, float jetYE, float cmBuffer, float etaPhiMargin)
{
    double x_inner = cut_x_inner + cmBuffer;
    double x_outer = cut_x_outer - cmBuffer;
    double y_min   = cut_y_min   + cmBuffer;
    double y_max   = cut_y_max   - cmBuffer;

    return etaPhiSignedDistance(jetXE, jetYE, x_inner, x_outer, y_min, y_max) >= etaPhiMargin;
}


// symmetric two-window acceptance at ECAL front plane, used for reco jet
// selection (with an adjustable margin)
inline bool pass_fiducial_cut(float jetXE, float jetYE, float buffer)
{
    const float ax = std::fabs(jetXE);
    const bool in_x = ax >= (cut_x_inner + buffer) && ax <= (cut_x_outer - buffer);
    const bool in_y = jetYE >= (cut_y_min + buffer) && jetYE <= (cut_y_max - buffer);
    return in_x && in_y;
}

#endif