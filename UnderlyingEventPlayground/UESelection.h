#ifndef UE_SELECTION_H
#define UE_SELECTION_H

// Geometry of the underlying-event (UE) event selection, shared by
// max_dr_jet_placement.cpp (which draws it) and ue_map.cpp (which applies
// it), so both use exactly the same definition.
//
// For a point P on the FCS surface where the UE is assessed:
//   P'        the point 180 deg in azimuth from P at the same eta, where a
//             jet back-to-back with one at P would be
//   distance  of a jet axis J:  min( dR(J, P), dR(J, P') )
//   best      the largest distance any J inside the EM-only reco fiducial
//             region (either half) can reach while having zero overlap with
//             a jet at P, i.e. dR(J, P) >= ue_min_dr_to_point = 2R
//   tag jet   a jet whose distance is within ue_tag_tolerance_in_R * R of
//             the best one (and at least 2R from P)
// By symmetry dR(J, P') = dR(P, J'), so the distance is equally how far P is
// from the nearer of the jet and the point back-to-back with the jet.

#include <vector>
#include <cmath>
#include <algorithm>

#include "JetParameters.h"

// A jet axis closer than this to P would overlap a jet at P
static const double ue_min_dr_to_point = 2.0 * R;

// Allowed shortfall, in units of R, of a tag jet's distance relative to the
// best achievable one
static const double ue_tag_tolerance_in_R = 0.3;

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

// min( dR(J, P), dR(J, P') )
inline double JetDistance(double jet_eta, double jet_phi, double p_eta, double p_phi)
{
    return std::min(DeltaR(jet_eta, jet_phi, p_eta, p_phi),
                    DeltaR(jet_eta, jet_phi, p_eta, p_phi + M_PI));
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
// JetDistance over the candidates that are at least ue_min_dr_to_point from
// P. Returns -1 if no candidate qualifies; best_index (if given) receives the
// index of the candidate that achieves it.
inline double BestJetDistance(double p_eta, double p_phi,
                              const std::vector<SurfacePoint> &candidates,
                              size_t *best_index = 0)
{
    double best = -1.0;

    for (size_t j = 0; j < candidates.size(); j++)
    {
        double d_p = DeltaR(candidates[j].eta, candidates[j].phi, p_eta, p_phi);
        if (d_p < ue_min_dr_to_point || d_p <= best)
            continue;

        double d = std::min(d_p, DeltaR(candidates[j].eta, candidates[j].phi, p_eta, p_phi + M_PI));
        if (d > best)
        {
            best = d;
            if (best_index)
                *best_index = j;
        }
    }
    return best;
}

#endif
