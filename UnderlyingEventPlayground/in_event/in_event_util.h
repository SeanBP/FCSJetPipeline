// Small helpers of the in-event UE studies, from the retired proxy UE's
// underlying_event/ProxySelection.h (git history): the reco jet fiducial
// test, Delta R and the highest EM trigger level.

#ifndef IN_EVENT_UTIL_H
#define IN_EVENT_UTIL_H

#include <cmath>

#include "JetParameters.h"

// Whether a detector-plane position is inside the EM-only reco fiducial region
inline bool ProxyInFiducial(double x, double y)
{
    const FiducialRect &r = kFiducialRectEcal;
    const double b = reco_fiducial_buffer_ecal;
    return std::fabs(x) >= r.x_inner + b && std::fabs(x) <= r.x_outer - b && y >= r.y_min + b && y <= r.y_max - b;
}

inline double ProxyDeltaR(double eta1, double phi1, double eta2, double phi2)
{
    const double dphi = std::remainder(phi1 - phi2, 2.0 * M_PI);
    return std::sqrt((eta1 - eta2) * (eta1 - eta2) + dphi * dphi);
}

// Highest EM trigger level fired (0, 1, 2), -1 if none. With prescales a lower
// flag can be missing, so the highest flag present is used.
inline int ProxyEmLevel(int em0, int em1, int em2)
{
    if (em2) return 2;
    if (em1) return 1;
    if (em0) return 0;
    return -1;
}

#endif
