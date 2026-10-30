// ECal cell geometry and the eta-weighted annulus estimate of the in-event UE
// (in_event_ue.cpp, annulus_diagram.cpp).
//
// Annulus: the cells of the jet's eta band |eta - eta_jet| < R, at least 2R
// from the jet and R from every other reco jet, in a phi range (whole ring,
// transverse 60-120 deg on both sides, or away side > 120 deg). The band is cut
// in eta strips 0.1 wide; a point in strip s is weighted by
//   w_s = A(jet cone, live; s) / A(annulus, live; s)
// so the weighted annulus has the eta distribution of the jet cone's live
// area. Strips with less live annulus area than kStripMinLive of a full cone
// slice are uncovered; the covered strips' share of the jet cone's template
// sum scales the estimate up to the whole live jet cone.

#ifndef IN_EVENT_GEO_H
#define IN_EVENT_GEO_H

#include <cmath>
#include <map>
#include <vector>

#include "JetParameters.h"
#include "in_event_util.h"

static const double kCell = 3.0;      // cm
static const double kBucket = 0.1;    // eta-phi bucket size of the cell lookup
static const double kStrip = 0.1;     // eta strips of the strip estimates
static const int kNS = 14;            // 2 R / kStrip
static const double kStripMinLive = 0.2;

// eta-phi area of the slice of a radius-R disk between eta offsets a and b from its centre
static double DiskSlice(double a, double b)
{
    a = std::max(a, -R); b = std::min(b, R);
    if (b <= a) return 0;
    struct F { static double I(double u) { return u * std::sqrt(R * R - u * u) + R * R * std::asin(u / R); } };
    return F::I(b) - F::I(a);
}
static int StripOf(double deta) { const int s = (int)std::floor((deta + R) / kStrip); return s >= 0 && s < kNS ? s : -1; }

// ECal detector cells (no fiducial buffer), with their eta-phi area
struct Cell { float x, y, eta, phi, area; };
// sub-pixels of the cells (kSub x kSub per cell) for the pixel ring estimate
static const int kSub = 3;
struct Pix { float eta, phi, area; };
struct Geometry
{
    std::vector<Cell> cells;
    int nx, ny;
    double x_lo, y_lo;
    std::map<std::pair<int, int>, std::vector<int> > buckets;
    std::vector<Pix> pix;
    void Build()
    {
        const FiducialRect &r = kFiducialRectEcal;
        x_lo = r.x_inner; y_lo = r.y_min;
        nx = (int)std::ceil((r.x_outer - r.x_inner) / kCell); ny = (int)std::ceil((r.y_max - r.y_min) / kCell);
        const double dx = (r.x_outer - r.x_inner) / nx, dy = (r.y_max - r.y_min) / ny;
        for (int side = -1; side <= 1; side += 2)
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    Cell c;
                    c.x = side * (x_lo + (i + 0.5) * dx); c.y = y_lo + (j + 0.5) * dy;
                    double eta, phi; computeEtaPhi(c.x, c.y, eta, phi);
                    c.eta = eta; c.phi = phi;
                    const double r2 = (double)c.x * c.x + (double)c.y * c.y;
                    c.area = dx * dy * z_proj / (r2 * std::sqrt(r2 + z_proj * z_proj));
                    buckets[Key(eta, phi)].push_back((int)cells.size());
                    cells.push_back(c);
                    for (int a = 0; a < kSub; a++)
                        for (int b = 0; b < kSub; b++)
                        {
                            const double px = side * (x_lo + (i + (a + 0.5) / kSub) * dx), py = y_lo + (j + (b + 0.5) / kSub) * dy;
                            double pe, pp; computeEtaPhi(px, py, pe, pp);
                            const double p2 = px * px + py * py;
                            Pix q = {(float)pe, (float)pp, (float)(dx * dy / (kSub * kSub) * z_proj / (p2 * std::sqrt(p2 + z_proj * z_proj)))};
                            pix.push_back(q);
                        }
                }
    }
    static std::pair<int, int> Key(double eta, double phi)
    {
        return std::make_pair((int)std::floor(eta / kBucket), (int)std::floor((phi + M_PI) / kBucket));
    }
    // cell of a detector position, -1 outside
    int CellOf(double x, double y) const
    {
        const int side = x > 0 ? 1 : 0;
        const int i = (int)std::floor((std::fabs(x) - x_lo) / ((kFiducialRectEcal.x_outer - x_lo) / nx));
        const int j = (int)std::floor((y - y_lo) / ((kFiducialRectEcal.y_max - y_lo) / ny));
        if (i < 0 || i >= nx || j < 0 || j >= ny) return -1;
        return side * nx * ny + i * ny + j;
    }
    // cells with centres within dr of (eta, phi)
    void Near(double eta, double phi, double dr, std::vector<int> &out) const
    {
        out.clear();
        const int nb = (int)std::ceil(dr / kBucket) + 1, nphi = (int)std::ceil(2 * M_PI / kBucket);
        const std::pair<int, int> k0 = Key(eta, phi);
        for (int a = -nb; a <= nb; a++)
            for (int b = -nb; b <= nb; b++)
            {
                const int kb = ((k0.second + b) % nphi + nphi) % nphi;
                std::map<std::pair<int, int>, std::vector<int> >::const_iterator it = buckets.find(std::make_pair(k0.first + a, kb));
                if (it == buckets.end()) continue;
                for (size_t q = 0; q < it->second.size(); q++)
                {
                    const Cell &c = cells[it->second[q]];
                    if (ProxyDeltaR(c.eta, c.phi, eta, phi) < dr) out.push_back(it->second[q]);
                }
            }
    }
};

// a projected (eta, phi) direction: the live ECal cell it hits, -1 if none
static int CellOfDirection(const Geometry &g, double eta, double phi)
{
    if (!(eta > 0)) return -1;
    const double r = z_proj / std::sinh(eta);
    return g.CellOf(r * std::cos(phi), r * std::sin(phi));
}

static const int kNA = 3;  // annulus phi ranges: 0 whole ring, 1 transverse, 2 away side
static const double kAnnJetExcl = 2 * R, kAnnOtherExcl = R;
static const double kAnnPhiLo[kNA] = {0.0, M_PI / 3, 2 * M_PI / 3}, kAnnPhiHi[kNA] = {M_PI, 2 * M_PI / 3, M_PI};

struct Annulus
{
    std::vector<float> w;  // per cell: weight of a point in it (0 outside the annulus)
    double fill;           // template scale-up for uncovered strips (0 if nothing covered)
    double area_live;      // live annulus area (eta-phi)
    bool cov[kNS];         // strip covered
};

// A region of the jet's eta band: d_lo <= dR(jet) < d_hi, phi_lo <= |dphi| <= phi_hi,
// >= kAnnOtherExcl from others; strips covered if their live area is at least
// min_live of a full cone slice. Per-cell weights as for the annulus. Ajs[s]
// is the jet cone's live area in strip s and jet_T, Tjs[s] its template sums;
// others: (eta, phi) of the event's other reco jets.
static void BuildRegion(const Geometry &g, const std::vector<double> &rho, double jeta, double jphi,
                        const std::vector<std::pair<float, float> > &others, const double *Ajs, double jet_T,
                        const double *Tjs, double phi_lo, double phi_hi, double d_lo, double d_hi, double min_live, Annulus &an)
{
    an.w.assign(g.cells.size(), 0.f);
    std::vector<int> in;
    double As[kNS] = {0};
    for (size_t c = 0; c < g.cells.size(); c++)
    {
        const Cell &cl = g.cells[c];
        const int s = StripOf(cl.eta - jeta);
        if (s < 0) continue;
        const double dphi = std::fabs(std::remainder(cl.phi - jphi, 2 * M_PI));
        if (dphi < phi_lo || dphi > phi_hi) continue;
        const double d = ProxyDeltaR(cl.eta, cl.phi, jeta, jphi);
        if (d < d_lo || d >= d_hi) continue;
        bool near = false;
        for (size_t o = 0; o < others.size() && !near; o++)
            near = ProxyDeltaR(cl.eta, cl.phi, others[o].first, others[o].second) < kAnnOtherExcl;
        if (near) continue;
        in.push_back((int)c);
        As[s] += cl.area;
    }
    double ws[kNS], tcov = 0;
    an.area_live = 0;
    for (int s = 0; s < kNS; s++)
    {
        const double full = DiskSlice(-R + s * kStrip, -R + (s + 1) * kStrip);
        const bool cov = full > 0 && As[s] >= min_live * full;
        ws[s] = cov ? Ajs[s] / As[s] : 0;
        an.cov[s] = cov;
        if (cov) tcov += Tjs[s];
        an.area_live += As[s];
    }
    an.fill = tcov > 0 ? jet_T / tcov : 0;
    for (size_t q = 0; q < in.size(); q++)
        an.w[in[q]] = ws[StripOf(g.cells[in[q]].eta - jeta)];
}

// The annulus of range v for a jet (>= 2R from it)
static void BuildAnnulus(const Geometry &g, const std::vector<double> &rho, double jeta, double jphi,
                         const std::vector<std::pair<float, float> > &others, const double *Ajs, double jet_T,
                         const double *Tjs, int v, Annulus &an)
{
    BuildRegion(g, rho, jeta, jphi, others, Ajs, jet_T, Tjs, kAnnPhiLo[v], kAnnPhiHi[v], kAnnJetExcl, 1e9, kStripMinLive, an);
}

// Pixel ring: the cells (pixels) of the whole ring (cell centre |deta| < R,
// >= d_lo from the jet, >= kAnnOtherExcl from others) each get a weight w_i
// such that the weighted area-per-eta of the ring matches the jet cone's live
// area-per-eta: with h_jet(b), h_ring(b) the areas in eta bins b (from
// sub-pixels, so a cell spreads over the bins it overlaps),
//   w_i = sum_b A_i(b) h_jet(b) / h_ring(b) / A_i
// The UE density is the weighted average over the ring, each point weighted
// by the ratio of its own eta bin (Weight: exact where a cell spans bins of
// different ratio, as at the beam hole; ring membership stays per cell),
//   rho = sum_points E r(b(point)) / sum_b r(b) h_ring(b),  r = h_jet / h_ring,
// with sum_b r(b) h_ring(b) = sum_i w_i A_i; the cell weights w_i remain for
// pictures. The jet's UE is rho * A_jet (live). Bins of the jet the ring does not
// reach (h_ring < kPixMinRatio h_jet) carry no weight; the normalisation
// extrapolates the matched density to the whole jet area. No jet is dropped.
static const double kPixBin = 0.05, kPixMinRatio = 0.1;
static const int kPixNB = 32;  // (2 R + 0.2) / kPixBin: the band plus half a cell on each side
struct PixelRing
{
    double eta0;            // first bin's lower edge (eta_jet - R - 0.1)
    std::vector<float> w;   // per cell weight (0 outside the ring)
    double den;             // sum_i w_i A_i: the jet area in matched bins
    double area_jet, area_ring, w_max;
    double ratio[kPixNB];   // h_jet / h_ring per eta bin (0: not reached)
    std::vector<char> in;   // per cell: in the ring
    int Bin(double eta) const { const int b = (int)std::floor((eta - eta0) / kPixBin); return b >= 0 && b < kPixNB ? b : -1; }
    // a point's weight: the ratio of its own eta bin if its cell is in the ring
    double Weight(int cell, double eta) const { if (cell < 0 || !in[cell]) return 0; const int b = Bin(eta); return b >= 0 ? ratio[b] : 0; }
};
// dphi_max: cells with |dphi| to the jet above it are left out (a cut around the recoil side)
static bool InPixelRing(double eta, double phi, double jeta, double jphi, double d_lo, const std::vector<std::pair<float, float> > &others, double dphi_max = M_PI)
{
    if (std::fabs(eta - jeta) >= R || ProxyDeltaR(eta, phi, jeta, jphi) < d_lo) return false;
    if (std::fabs(std::remainder(phi - jphi, 2 * M_PI)) > dphi_max) return false;
    for (size_t o = 0; o < others.size(); o++)
        if (ProxyDeltaR(eta, phi, others[o].first, others[o].second) < kAnnOtherExcl) return false;
    return true;
}
static void BuildPixelRing(const Geometry &g, double jeta, double jphi, const std::vector<std::pair<float, float> > &others, double d_lo, PixelRing &pr,
                           double dphi_max = M_PI)
{
    pr.eta0 = jeta - R - 0.1;
    pr.w.assign(g.cells.size(), 0.f);
    pr.in.assign(g.cells.size(), 0);
    double hj[kPixNB] = {0}, hr[kPixNB] = {0};
    std::vector<int> ring;
    pr.area_jet = pr.area_ring = 0;
    const int ns = kSub * kSub;
    for (size_t c = 0; c < g.cells.size(); c++)
    {
        const Cell &cl = g.cells[c];
        if (std::fabs(cl.eta - jeta) > R + 0.2) continue;
        const bool in_ring = InPixelRing(cl.eta, cl.phi, jeta, jphi, d_lo, others, dphi_max);
        if (in_ring) { ring.push_back((int)c); pr.in[c] = 1; }
        for (int k = 0; k < ns; k++)
        {
            const Pix &p = g.pix[c * ns + k];
            const int b = pr.Bin(p.eta);
            if (b < 0) continue;
            if (ProxyDeltaR(p.eta, p.phi, jeta, jphi) < R) { hj[b] += p.area; pr.area_jet += p.area; }
            if (in_ring) { hr[b] += p.area; pr.area_ring += p.area; }
        }
    }
    double *ratio = pr.ratio;
    for (int b = 0; b < kPixNB; b++) ratio[b] = hj[b] > 0 && hr[b] >= kPixMinRatio * hj[b] ? hj[b] / hr[b] : 0;
    pr.den = 0; pr.w_max = 0;
    for (size_t q = 0; q < ring.size(); q++)
    {
        const int c = ring[q];
        double sw = 0;
        for (int k = 0; k < ns; k++)
        {
            const Pix &p = g.pix[c * ns + k];
            const int b = pr.Bin(p.eta);
            if (b >= 0) sw += p.area * ratio[b];
        }
        pr.w[c] = sw / g.cells[c].area;
        pr.den += sw;
        pr.w_max = std::max(pr.w_max, (double)pr.w[c]);
    }
}

#endif
