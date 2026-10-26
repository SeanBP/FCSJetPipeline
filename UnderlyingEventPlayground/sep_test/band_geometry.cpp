// Proxy-band geometry of the separation test (ProxySelection.h), per UE cone
// position over one half of the reco fiducial region:
//   hw_sep    the pipeline half-width: the band widened until its area
//             outside 2.5R of the cone reaches the target (5% of a half)
//   hw_nosep  the half-width with no exclusion (PROXY_MIN_DR_R=0)
//   left      the fraction of the target that stays allowed inside the
//             hw_nosep band once proxy jets within 2.5R of the cone are
//             excluded (no separation with a 2.5R veto on all jets)
// Printed per jet-eta slice: the means and ranges over positions, and the
// share of positions where nothing or under 25% of the target is left.
//
// Build and run from this folder:
//   g++ -m64 -O2 band_geometry.cpp -o band_geometry -I../../shared -I../../underlying_event `root-config --cflags --libs` `fastjet-config --cxxflags` -std=c++11
//   ./band_geometry

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "ProxySelection.h"

// Half-width giving the target area with no exclusion
static double HalfWidthNoExclusion(double ue_phi, double target)
{
    const std::vector<ProxyCell> &cells = ProxyFiducialCells();
    std::vector<std::pair<double, double> > dev_area;
    for (size_t i = 0; i < cells.size(); i++)
    {
        const double dev = ProxyOffAxisDeviationDeg(ue_phi, cells[i].phi);
        if (dev <= proxy_half_width_max_deg) dev_area.push_back(std::make_pair(dev, cells[i].area));
    }
    std::sort(dev_area.begin(), dev_area.end());
    double area = 0.0;
    for (size_t i = 0; i < dev_area.size(); i++)
        if ((area += dev_area[i].second) >= target) return dev_area[i].first;
    return proxy_half_width_max_deg;
}

int main()
{
    const double target = ProxyTargetArea();
    const FiducialRect &r = kFiducialRectEcal;
    const double b = reco_fiducial_buffer_ecal;
    const double edges[] = {2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.6};
    const int ns = 6;
    struct S { long n = 0, zero = 0, low = 0; double sep = 0, nosep = 0, left = 0, sep_min = 1e9, sep_max = 0, ns_min = 1e9, ns_max = 0; } s[ns];
    printf("proxy_min_dr %.3f (%.1fR), target area %.4f\n", proxy_min_dr, proxy_min_dr / R, target);
    for (double x = r.x_inner + b; x <= r.x_outer - b; x += 4.0)
        for (double y = r.y_min + b; y <= r.y_max - b; y += 4.0)
        {
            double eta, phi;
            computeEtaPhi(x, y, eta, phi);
            int k = -1;
            for (int i = 0; i < ns; i++) if (eta >= edges[i] && eta < edges[i + 1]) k = i;
            if (k < 0) continue;
            const double hs = ProxyHalfWidthDeg(eta, phi);
            const double hn = HalfWidthNoExclusion(phi, target);
            const double left = ProxyAllowedArea(eta, phi, hn) / target;
            S &a = s[k];
            a.n++; a.sep += hs; a.nosep += hn; a.left += left;
            a.sep_min = std::min(a.sep_min, hs); a.sep_max = std::max(a.sep_max, hs);
            a.ns_min = std::min(a.ns_min, hn); a.ns_max = std::max(a.ns_max, hn);
            if (left <= 0) a.zero++;
            if (left < 0.25) a.low++;
        }
    printf("  jet eta    positions  hw_sep mean [min-max]   hw_nosep mean [min-max]   left mean   left=0   left<25%%\n");
    for (int i = 0; i < ns; i++)
    {
        const S &a = s[i];
        if (!a.n) continue;
        printf("  %.1f-%.1f   %5ld      %5.1f [%4.1f-%4.1f]        %5.1f [%4.1f-%4.1f]          %4.2f      %4.0f%%     %4.0f%%\n",
               edges[i], edges[i + 1], a.n, a.sep / a.n, a.sep_min, a.sep_max, a.nosep / a.n, a.ns_min, a.ns_max,
               a.left / a.n, 100.0 * a.zero / a.n, 100.0 * a.low / a.n);
    }
    return 0;
}
