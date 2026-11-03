// The in-event underlying-event (UE) estimate of a jet: the pixel ring.
//
// Ring: the jet's eta band |eta - eta_jet| < R, at least 2R from the jet axis
// (kUeJetExcl), within kUeDphiMax of the jet in azimuth (the 30 deg around
// 180 are left out: the recoil side), other jets not excluded. The ring is
// weighted so that its area-per-eta matches the jet cone's: with h_jet(b) and
// h_ring(b) the (live) eta-phi areas of the jet cone and the ring in eta bins
// b of kUeBin, a point in bin b weighs r(b) = h_jet(b) / h_ring(b). Bins the
// ring covers under kUeMinRatio of the jet's area carry no weight. The UE
// density is
//   rho = sum_points E r(b(point)) / sum_b r(b) h_ring(b),
// whose denominator is the jet's area in the bins the ring reaches, so empty
// parts of the ring count as area with no energy; the density of the
// reached bins stands in for the others. The jet's UE is rho times the jet
// cone's area: its live area on the ECal at reco level, pi R^2 at truth level.
//
// Reco (RecoPixelRing): the ECal (both halves of kFiducialRectEcal) as
// kUeCell pixels, each split into kUeSub x kUeSub sub-pixels for the areas;
// ring membership per pixel (its centre), weights per point (its own eta bin).
// Truth (TruthPixelRing): truth photons with full azimuth; the areas are
// exact, so the weights depend only on eta - eta_jet and are computed once.
//
// Chosen on FMS Run 15 and FCS simulation, 2026-10 (UnderlyingEventPlayground/
// in_event): the azimuthal scans put the jet's edge at about 2R and the
// recoil rise from about 120-130 deg; on the FMS the FCS-like ring sits within
// about 10% of the FMS two-cone (+-90 deg) estimate.

#ifndef PIXEL_RING_UE_H
#define PIXEL_RING_UE_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "JetParameters.h"

static const double kUeJetExcl = 2.0 * R;           // ring starts this far from the jet axis
static const double kUeDphiMax = 150.0 * M_PI / 180; // ring ends here in |dphi|
static const double kUeBin = 0.05;                  // eta bins of the area matching
static const double kUeMinRatio = 0.1;              // a bin is reached if h_ring >= kUeMinRatio h_jet
static const double kUeCell = 3.0;                  // cm, reco pixel size
static const int kUeSub = 3;                        // sub-pixels per pixel side
static const int kUeNB = 32;                        // reco eta bins: 2R + 0.2 (half a pixel each side)

inline double UeDeltaR(double eta1, double phi1, double eta2, double phi2)
{
    const double dphi = std::remainder(phi1 - phi2, 2.0 * M_PI);
    return std::sqrt((eta1 - eta2) * (eta1 - eta2) + dphi * dphi);
}

// whether a direction is in the ring of a jet at (jeta, jphi)
inline bool InUeRing(double eta, double phi, double jeta, double jphi)
{
    return std::fabs(eta - jeta) < R && UeDeltaR(eta, phi, jeta, jphi) >= kUeJetExcl &&
           std::fabs(std::remainder(phi - jphi, 2.0 * M_PI)) <= kUeDphiMax;
}

struct UeResult
{
    double rho, rho_pt;  // UE density (energy, pT) per unit eta-phi area; -1: no ring at the jet's eta
    double area;         // the jet cone's area the UE is scaled to
    double cov;          // share of the jet's area in eta bins the ring reaches
    double E() const { return rho >= 0 ? rho * area : -1; }
    double pt() const { return rho_pt >= 0 ? rho_pt * area : -1; }
};

// The ECal as pixels: kUeCell squares over both halves, with sub-pixel areas
struct UeGeometry
{
    struct Pixel { float eta, phi; };
    struct Sub { float eta, phi, area; };
    std::vector<Pixel> pix;
    std::vector<Sub> sub;  // kUeSub^2 per pixel
    int nx = 0, ny = 0;
    double dx = 0, dy = 0;
    void Build()
    {
        const FiducialRect &r = kFiducialRectEcal;
        nx = (int)std::ceil((r.x_outer - r.x_inner) / kUeCell); ny = (int)std::ceil((r.y_max - r.y_min) / kUeCell);
        dx = (r.x_outer - r.x_inner) / nx; dy = (r.y_max - r.y_min) / ny;
        pix.clear(); sub.clear();
        for (int side = 0; side < 2; side++)
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    const double s = side ? 1 : -1;
                    double eta, phi;
                    computeEtaPhi(s * (r.x_inner + (i + 0.5) * dx), r.y_min + (j + 0.5) * dy, eta, phi);
                    pix.push_back(Pixel{(float)eta, (float)phi});
                    for (int a = 0; a < kUeSub; a++)
                        for (int b = 0; b < kUeSub; b++)
                        {
                            const double px = s * (r.x_inner + (i + (a + 0.5) / kUeSub) * dx), py = r.y_min + (j + (b + 0.5) / kUeSub) * dy;
                            double pe, pp;
                            computeEtaPhi(px, py, pe, pp);
                            const double r2 = px * px + py * py;
                            sub.push_back(Sub{(float)pe, (float)pp, (float)(dx * dy / (kUeSub * kUeSub) * z_proj / (r2 * std::sqrt(r2 + z_proj * z_proj)))});
                        }
                }
    }
    // pixel of an ECal position (cm at z_proj), -1 outside
    int PixelOf(double x, double y) const
    {
        const FiducialRect &r = kFiducialRectEcal;
        const int side = x > 0 ? 1 : 0;
        const int i = (int)std::floor((std::fabs(x) - r.x_inner) / dx), j = (int)std::floor((y - r.y_min) / dy);
        if (i < 0 || i >= nx || j < 0 || j >= ny) return -1;
        return (side * nx + i) * ny + j;
    }
};

// Reco: points (eta, phi, E, pixel) of the event's ECal; the UE of a jet at (jeta, jphi)
inline UeResult RecoPixelRing(const UeGeometry &g, double jeta, double jphi, const std::vector<float> &p_eta,
                              const std::vector<float> &p_phi, const std::vector<float> &p_E, const std::vector<int> &p_pix)
{
    UeResult res = {-1, -1, 0, 0};
    const double eta0 = jeta - R - 0.1;
    double hj[kUeNB] = {0}, hr[kUeNB] = {0};
    std::vector<char> in(g.pix.size(), 0);
    const int ns = kUeSub * kUeSub;
    for (size_t c = 0; c < g.pix.size(); c++)
    {
        if (std::fabs(g.pix[c].eta - jeta) > R + 0.2) continue;
        in[c] = InUeRing(g.pix[c].eta, g.pix[c].phi, jeta, jphi);
        for (int k = 0; k < ns; k++)
        {
            const UeGeometry::Sub &q = g.sub[c * ns + k];
            const int b = (int)std::floor((q.eta - eta0) / kUeBin);
            if (b < 0 || b >= kUeNB) continue;
            if (UeDeltaR(q.eta, q.phi, jeta, jphi) < R) { hj[b] += q.area; res.area += q.area; }
            if (in[c]) hr[b] += q.area;
        }
    }
    double ratio[kUeNB], den = 0;
    for (int b = 0; b < kUeNB; b++)
    {
        ratio[b] = hj[b] > 0 && hr[b] >= kUeMinRatio * hj[b] ? hj[b] / hr[b] : 0;
        den += ratio[b] * hr[b];
    }
    res.cov = res.area > 0 ? den / res.area : 0;
    if (!(den > 0)) return res;
    double sE = 0, sP = 0;
    for (size_t i = 0; i < p_E.size(); i++)
    {
        if (p_pix[i] < 0 || !in[p_pix[i]]) continue;
        const int b = (int)std::floor((p_eta[i] - eta0) / kUeBin);
        if (b < 0 || b >= kUeNB) continue;
        sE += p_E[i] * ratio[b];
        sP += p_E[i] / std::cosh(p_eta[i]) * ratio[b];
    }
    res.rho = sE / den; res.rho_pt = sP / den;
    return res;
}

// Truth: full azimuth, so the weights are a function of eta - eta_jet only
struct TruthPixelRing
{
    static const int kNB = 28;  // 2R / kUeBin
    double ratio[kNB], den;
    TruthPixelRing()
    {
        // exact areas from fine eta steps: at offset d the jet cone spans |dphi| < sqrt(R^2 - d^2), the
        // ring kUeJetExcl-distance-limited from sqrt(kUeJetExcl^2 - d^2) to kUeDphiMax, on both sides
        den = 0;
        const int nsub = 50;
        for (int b = 0; b < kNB; b++)
        {
            double hj = 0, hr = 0;
            for (int q = 0; q < nsub; q++)
            {
                const double d = -R + (b + (q + 0.5) / nsub) * kUeBin, de = kUeBin / nsub;
                hj += 2 * std::sqrt(std::max(0.0, R * R - d * d)) * de;
                hr += 2 * std::max(0.0, kUeDphiMax - std::sqrt(std::max(0.0, kUeJetExcl * kUeJetExcl - d * d))) * de;
            }
            ratio[b] = hj > 0 && hr >= kUeMinRatio * hj ? hj / hr : 0;
            den += ratio[b] * hr;
        }
    }
    // photons (eta, phi, E) of the event; the UE of a truth jet at (jeta, jphi), scaled to pi R^2
    UeResult Estimate(double jeta, double jphi, const std::vector<float> &p_eta, const std::vector<float> &p_phi,
                      const std::vector<float> &p_E) const
    {
        UeResult res = {-1, -1, M_PI * R * R, den / (M_PI * R * R)};
        if (!(den > 0)) return res;
        double sE = 0, sP = 0;
        for (size_t i = 0; i < p_E.size(); i++)
        {
            if (!InUeRing(p_eta[i], p_phi[i], jeta, jphi)) continue;
            const int b = (int)std::floor((p_eta[i] - jeta + R) / kUeBin);
            if (b < 0 || b >= kNB) continue;
            sE += p_E[i] * ratio[b];
            sP += p_E[i] / std::cosh(p_eta[i]) * ratio[b];
        }
        res.rho = sE / den; res.rho_pt = sP / den;
        return res;
    }
};

#endif
