#ifndef UE_SELECTION_H
#define UE_SELECTION_H

// Geometry of the underlying-event (UE) event selection, shared by
// max_dr_jet_placement.cpp (which draws it) and ue_map.cpp (which applies
// it), so both use exactly the same definition.
//
// For a point P on the FCS surface where the UE is assessed:
//   P' line   the line at 180 deg in azimuth from P, at any eta: a jet
//             back-to-back with one at P is opposite in azimuth but need not
//             be at the same eta, so it can be anywhere along this line
//   allowed   a jet axis J with zero overlap with a jet at P, i.e.
//             dR(J, P) >= ue_min_dr_to_point = 2R
//   distance  of a jet axis J:  |dphi(J, P' line)|, the eta-phi distance from
//             J to the nearest point of the line (the one at J's own eta).
//             The larger it is, the less a jet at J can overlap a jet
//             back-to-back with one at P.
//   best      the largest distance any allowed J inside the EM-only reco
//             fiducial region of the half OPPOSITE P can reach. Restricting
//             the reference to the opposite half keeps it continuous in P:
//             taken over both halves it jumps where axes on P's own side,
//             near P in azimuth, first clear the 2R requirement.
//   tag jet   an allowed jet (on either half) whose distance is within
//             ue_tag_tolerance_in_R * R of the best one, or at least
//             ue_tag_distance_cap = 2R: beyond 2R a jet has zero overlap with
//             any jet on the line, so more distance gains nothing
// By symmetry |dphi(J, P' line)| = |dphi(P, J' line)|, so the distance is
// equally how far P is in azimuth from the line back-to-back with the jet.

#include <vector>
#include <cmath>
#include <algorithm>

#include "JetParameters.h"

// A jet axis closer than this to P would overlap a jet at P
static const double ue_min_dr_to_point = 2.0 * R;

// Allowed shortfall, in units of R, of a tag jet's distance relative to the
// best achievable one
static const double ue_tag_tolerance_in_R = 0.3;

// A jet axis at least this far from the P' line has zero overlap with a jet
// anywhere on it, so no larger distance is ever required of a tag jet
static const double ue_tag_distance_cap = 2.0 * R;

struct SurfacePoint
{
    double x, y, eta, phi;
};

inline double DeltaPhi(double a, double b)
{
    double d = a - b;
    while (d >  M_PI) d -= 2.0 * M_PI;
    while (d < -M_PI) d += 2.0 * M_PI;
    return d;
}

inline double DeltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = DeltaPhi(phi1, phi2);
    return std::sqrt(deta * deta + dphi * dphi);
}

// Distance of a jet axis from the P' line (180 deg in azimuth from P): |dphi|
inline double JetDistance(double jet_phi, double p_phi)
{
    return std::fabs(DeltaPhi(jet_phi, p_phi + M_PI));
}

// Candidate jet-axis points: the EM-only reco fiducial region on both
// halves, sampled every `spacing` cm. The last point of each axis is pinned
// to the far edge so the corners, where the extremes are, are sampled exactly.
inline std::vector<SurfacePoint> FiducialJetAxisCandidates(double spacing)
{
    const FiducialRect &rect = kFiducialRectEcal;
    const float buf = reco_fiducial_buffer_ecal;

    const double x_lo = rect.x_inner + buf;
    const double x_hi = rect.x_outer - buf;
    const double y_lo = rect.y_min + buf;
    const double y_hi = rect.y_max - buf;

    std::vector<double> xs, ys;
    for (double x = x_lo; x < x_hi; x += spacing) xs.push_back(x);
    xs.push_back(x_hi);
    for (double y = y_lo; y < y_hi; y += spacing) ys.push_back(y);
    ys.push_back(y_hi);

    std::vector<SurfacePoint> pts;
    for (int side = -1; side <= 1; side += 2)
    for (size_t ix = 0; ix < xs.size(); ix++)
    for (size_t iy = 0; iy < ys.size(); iy++)
    {
        SurfacePoint j = {side * xs[ix], ys[iy], 0, 0};
        computeEtaPhi(j.x, j.y, j.eta, j.phi);
        pts.push_back(j);
    }
    return pts;
}

// Best achievable distance for the point (p_eta, p_phi): the maximum of
// JetDistance over the candidates on the half opposite P that are at least
// ue_min_dr_to_point from P. Returns -1 if no candidate qualifies; best_index (if given) receives the
// index of the candidate that achieves it.
inline double BestJetDistance(double p_eta, double p_phi,
                              const std::vector<SurfacePoint> &candidates,
                              size_t *best_index = 0)
{
    double best = -1.0;
    const double p_side = std::cos(p_phi);

    for (size_t j = 0; j < candidates.size(); j++)
    {
        if (candidates[j].x * p_side >= 0.0)
            continue;

        double d_p = DeltaR(candidates[j].eta, candidates[j].phi, p_eta, p_phi);
        if (d_p < ue_min_dr_to_point)
            continue;

        double d = JetDistance(candidates[j].phi, p_phi);
        if (d > best)
        {
            best = d;
            if (best_index)
                *best_index = j;
        }
    }
    return best;
}

// Smallest distance a tag jet may have, given the best achievable one
inline double MinTagDistance(double best_distance)
{
    return std::min(best_distance - ue_tag_tolerance_in_R * R, ue_tag_distance_cap);
}

#endif
