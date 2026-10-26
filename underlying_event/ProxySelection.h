#ifndef PROXY_SELECTION_H
#define PROXY_SELECTION_H

// Selection rules of the near off-axis proxy-event underlying-event (UE)
// estimate, applied by ProxySearch.h. Studies behind the numbers are in
// UnderlyingEventPlayground/.
//
// The UE under a jet is sampled in other events (proxy events), in a cone at
// the jet's position (the UE cone). A proxy event needs its own reco jet (the
// proxy jet) placed so that it stays out of the UE cone:
//   separation   dR(proxy jet, UE cone) >= proxy_min_dr = 2R + R/2: no cone
//                overlap, plus the truth-reco matching radius
//   near off-axis  | |dphi| - 90 deg | <= half-width, at any eta
//   fiducial     proxy jet inside the EM-only reco fiducial region
//   pT           proxy jet pT within proxy_pt_window_fraction of the jet's
//   trigger      same highest EM trigger level (fcsEM0/1/2), which sums ECal
//                ET only; the jet-patch triggers include uncalibrated HCal
//   group        same run in data (same detector state; luminosity changes
//                1-4% within a run), same ptHat bin in simulation
//
// The half-width depends on the UE position: the smallest that gives the
// allowed region the fixed eta-phi area ProxyTargetArea(), capped at
// proxy_half_width_max_deg to stay away from the recoil direction. Eta-phi
// area is used because jet axes are close to uniform in eta-phi.

#include <vector>
#include <cmath>
#include <algorithm>

#include "JetParameters.h"

// Minimum eta-phi distance between the proxy jet axis and the UE cone axis,
// 2.5R; building with -DPROXY_MIN_DR_R=<x> sets it to x R instead (0 = no
// separation: the band is then sized without the exclusion and the jet veto,
// which uses the same distance, does nothing; UnderlyingEventPlayground's
// separation test)
#ifndef PROXY_MIN_DR_R
#define PROXY_MIN_DR_R 2.5
#endif
static const double proxy_min_dr = PROXY_MIN_DR_R * R;

// Veto: reject a proxy event if any of its reco jets (pT >= reco_jet_pt_min_em)
// lies within proxy_min_dr of the UE cone. Off by default: the cone keeps
// whatever it contains, as the truth UE's same-event cones do. In simulation
// (UnderlyingEventPlayground/truth_decomp, unfiltered sample) separation
// without the veto matches the truth UE within a few percent, while the veto
// removes real UE and comes out ~10% low. Building with -DPROXY_JET_VETO=1
// turns it on (the separation test's sep_veto variant).
#ifndef PROXY_JET_VETO
#define PROXY_JET_VETO 0
#endif
static const bool proxy_jet_veto = PROXY_JET_VETO;

// Target allowed area, as a fraction of one detector half's reco fiducial
// area in eta-phi (1.70), i.e. 0.17
static const double proxy_target_area_fraction = 0.10;

// Largest half-width of the band about 90 deg, degrees (the target area
// needs at most 46)
static const double proxy_half_width_max_deg = 50.0;

// Half-width of the pT window, as a fraction of the pT of the jet being
// corrected: the proxy jet must have pT within jet pT * (1 +- this).
static const double proxy_pt_window_fraction = 0.05;

// Every accepted proxy event is used (proxy_n_max = 0; a positive value stops
// the search there, nearest events first). A jet with fewer than proxy_n_min
// has its pT window widened in steps of proxy_pt_window_step, up to
// proxy_pt_window_max, until it reaches proxy_n_min.
static const int proxy_n_max = 0;
static const int proxy_n_min = 4;
static const double proxy_pt_window_step = 0.05;
static const double proxy_pt_window_max = 0.50;

// Side of the square cells the fiducial region is divided into to measure area, cm
static const double proxy_cell_cm = 1.0;

struct ProxyCell
{
    double x, y, eta, phi;
    double area;                                  // eta-phi area of the cell
};

// The EM-only reco fiducial region as cells, built once
inline const std::vector<ProxyCell> &ProxyFiducialCells()
{
    static std::vector<ProxyCell> cells;
    if (!cells.empty())
        return cells;

    const FiducialRect &rect = kFiducialRectEcal;
    const double buf = reco_fiducial_buffer_ecal;
    const double x_lo = rect.x_inner + buf, x_hi = rect.x_outer - buf;
    const double y_lo = rect.y_min + buf, y_hi = rect.y_max - buf;
    const int nx = (int)std::ceil((x_hi - x_lo) / proxy_cell_cm), ny = (int)std::ceil((y_hi - y_lo) / proxy_cell_cm);
    const double dx = (x_hi - x_lo) / nx, dy = (y_hi - y_lo) / ny;

    for (int side = -1; side <= 1; side += 2)
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
            {
                ProxyCell c;
                c.x = side * (x_lo + (i + 0.5) * dx);
                c.y = y_lo + (j + 0.5) * dy;
                computeEtaPhi(c.x, c.y, c.eta, c.phi);
                // r = z / sinh(eta), so dx dy = r^2 coth(eta) deta dphi
                const double r2 = c.x * c.x + c.y * c.y;
                c.area = dx * dy * z_proj / (r2 * std::sqrt(r2 + z_proj * z_proj));
                cells.push_back(c);
            }
    return cells;
}

// The target area in eta-phi: proxy_target_area_fraction of one half's reco fiducial area
inline double ProxyTargetArea()
{
    static double target = -1.0;
    if (target < 0.0)
    {
        const std::vector<ProxyCell> &cells = ProxyFiducialCells();
        double total = 0.0;
        for (size_t i = 0; i < cells.size(); i++) total += cells[i].area;
        target = proxy_target_area_fraction * 0.5 * total;
    }
    return target;
}

// Whether a detector-plane position is inside the EM-only reco fiducial region
inline bool ProxyInFiducial(double x, double y)
{
    const FiducialRect &r = kFiducialRectEcal;
    const double b = reco_fiducial_buffer_ecal;
    return std::fabs(x) >= r.x_inner + b && std::fabs(x) <= r.x_outer - b && y >= r.y_min + b && y <= r.y_max - b;
}

// How far, in degrees, an azimuth is from 90 deg away from the UE cone: | |dphi| - 90 |
inline double ProxyOffAxisDeviationDeg(double ue_phi, double phi)
{
    const double dphi = std::fabs(std::remainder(phi - ue_phi, 2.0 * M_PI));
    return std::fabs(dphi * 180.0 / M_PI - 90.0);
}

inline double ProxyDeltaR(double eta1, double phi1, double eta2, double phi2)
{
    const double dphi = std::remainder(phi1 - phi2, 2.0 * M_PI);
    return std::sqrt((eta1 - eta2) * (eta1 - eta2) + dphi * dphi);
}

// Area (eta-phi) of the allowed region for a UE cone at (ue_eta, ue_phi) with
// the given half-width
inline double ProxyAllowedArea(double ue_eta, double ue_phi, double half_width_deg)
{
    const std::vector<ProxyCell> &cells = ProxyFiducialCells();
    double area = 0.0;
    for (size_t i = 0; i < cells.size(); i++)
    {
        if (ProxyOffAxisDeviationDeg(ue_phi, cells[i].phi) > half_width_deg) continue;
        if (ProxyDeltaR(cells[i].eta, cells[i].phi, ue_eta, ue_phi) < proxy_min_dr) continue;
        area += cells[i].area;
    }
    return area;
}

// Half-width (deg) giving the allowed region the area target_area
// (ProxyTargetArea() if negative) for a UE cone at (ue_eta, ue_phi), capped at
// max_deg; area_reached receives the area obtained. Sorts every cell on each
// call, so per-jet use goes through ProxySearch.h's ProxyWidthTable.
inline double ProxyHalfWidthDeg(double ue_eta, double ue_phi, double *area_reached = 0,
                                double target_area = -1.0, double max_deg = proxy_half_width_max_deg)
{
    if (target_area < 0.0) target_area = ProxyTargetArea();
    const std::vector<ProxyCell> &cells = ProxyFiducialCells();
    std::vector<std::pair<double, double> > dev_area;
    for (size_t i = 0; i < cells.size(); i++)
    {
        const double dev = ProxyOffAxisDeviationDeg(ue_phi, cells[i].phi);
        if (dev > max_deg) continue;
        if (ProxyDeltaR(cells[i].eta, cells[i].phi, ue_eta, ue_phi) < proxy_min_dr) continue;
        dev_area.push_back(std::make_pair(dev, cells[i].area));
    }
    std::sort(dev_area.begin(), dev_area.end());

    double area = 0.0;
    for (size_t i = 0; i < dev_area.size(); i++)
    {
        area += dev_area[i].second;
        if (area >= target_area)
        {
            if (area_reached) *area_reached = area;
            return dev_area[i].first;
        }
    }
    if (area_reached) *area_reached = area;
    return max_deg;
}

// Whether a proxy jet axis is allowed for a UE cone at (ue_eta, ue_phi), given
// that UE position's half-width. The caller checks that the axis is inside
// the reco fiducial region.
inline bool ProxyJetAllowed(double ue_eta, double ue_phi, double half_width_deg, double jet_eta, double jet_phi)
{
    return ProxyOffAxisDeviationDeg(ue_phi, jet_phi) <= half_width_deg &&
           ProxyDeltaR(jet_eta, jet_phi, ue_eta, ue_phi) >= proxy_min_dr;
}

// Whether a proxy jet's pT is close enough to that of the jet being corrected
inline bool ProxyPtAllowed(double jet_pt, double proxy_jet_pt, double window_fraction = proxy_pt_window_fraction)
{
    return std::fabs(proxy_jet_pt - jet_pt) <= window_fraction * jet_pt;
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

// Whether a proxy event's trigger matches that of the event being corrected:
// the same highest EM level, and an EM trigger fired at all
inline bool ProxyTriggerAllowed(int em_level, int proxy_em_level)
{
    return em_level >= 0 && proxy_em_level == em_level;
}

// Whether a proxy event is from the same run as the event being corrected
inline bool ProxyRunAllowed(int run, int proxy_run)
{
    return proxy_run == run;
}

#endif
