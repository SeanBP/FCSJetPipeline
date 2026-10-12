// For every point on the FCS surface inside the EM-only reco fiducial region,
// finds how far from back-to-back a point on the OPPOSITE half can be: the
// minimum azimuthal separation |dphi| between the point and any point of the
// opposite half's fiducial region. 180 deg is exactly back-to-back (where the
// recoil of a jet at this point goes); the smaller the value, the further
// from the recoil direction an underlying-event cone on the other half can
// be placed.
//
// Two maps are made:
//   1. unconstrained -- any opposite-half fiducial point is allowed
//   2. no overlap    -- the opposite point must also be at least min_dR away
//                       in eta-phi (default 2R, so two cones of radius R
//                       around the two points do not overlap)
//
// The fiducial region is the one reco EM-only jets are selected with
// (pass_fiducial_cut with kFiducialRectEcal / reco_fiducial_buffer_ecal),
// taken straight from shared/JetParameters.h.
//
// Build and run from this folder:
//   g++ -m64 min_opposite_angle.cpp -o min_opposite_angle -I../shared \
//       `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//   ./min_opposite_angle [min_dR] [grid_spacing_cm]
//
// Writes min_opposite_angle.pdf/.png/.root in the current directory.

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TLatex.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "JetParameters.h"

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

// Outlines of the ECAL-only rectangle (dashed) and the reco fiducial region
// inside it (solid), on both halves.
void DrawOutlines()
{
    const FiducialRect &rect = kFiducialRectEcal;
    const float buf = reco_fiducial_buffer_ecal;

    for (int side = -1; side <= 1; side += 2)
    {
        TBox *ecal = new TBox(side * rect.x_inner, rect.y_min, side * rect.x_outer, rect.y_max);
        ecal->SetFillStyle(0);
        ecal->SetLineStyle(2);
        ecal->SetLineColor(kGray + 2);
        ecal->Draw("same");

        TBox *fid = new TBox(side * (rect.x_inner + buf), rect.y_min + buf,
                             side * (rect.x_outer - buf), rect.y_max - buf);
        fid->SetFillStyle(0);
        fid->SetLineColor(kBlack);
        fid->Draw("same");
    }
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);

    double min_dR = 2.0 * R;
    double grid_spacing = 1.0;

    if (argc > 1 && argv[1][0] != '\0')
        min_dR = std::atof(argv[1]);
    if (argc > 2 && argv[2][0] != '\0')
        grid_spacing = std::atof(argv[2]);

    const FiducialRect &rect = kFiducialRectEcal;
    const float buf = reco_fiducial_buffer_ecal;

    const double x_lo = rect.x_inner + buf;
    const double x_hi = rect.x_outer - buf;
    const double y_lo = rect.y_min + buf;
    const double y_hi = rect.y_max - buf;

    // Candidate points on the opposite (x < 0) half. The last point of each
    // axis is pinned to the far edge so the corners, where the extreme
    // angles are, are sampled exactly.
    std::vector<double> xs, ys;
    for (double x = x_lo; x < x_hi; x += grid_spacing) xs.push_back(x);
    xs.push_back(x_hi);
    for (double y = y_lo; y < y_hi; y += grid_spacing) ys.push_back(y);
    ys.push_back(y_hi);

    std::vector<SurfacePoint> far_pts;
    for (size_t ix = 0; ix < xs.size(); ix++)
    for (size_t iy = 0; iy < ys.size(); iy++)
    {
        SurfacePoint q = {-xs[ix], ys[iy], 0, 0};
        computeEtaPhi(q.x, q.y, q.eta, q.phi);
        far_pts.push_back(q);
    }

    // Histogram spans both halves. The map is evaluated at the centre of
    // every bin that lies inside the fiducial region.
    const double hx = std::ceil(rect.x_outer / 10.0) * 10.0 + 10.0;
    const double hy_lo = std::floor(rect.y_min / 10.0) * 10.0 - 10.0;
    const double hy_hi = std::ceil(rect.y_max / 10.0) * 10.0 + 10.0;
    const int nbx = (int)std::round(2.0 * hx / grid_spacing);
    const int nby = (int)std::round((hy_hi - hy_lo) / grid_spacing);

    TH2D h_free("h_min_dphi", "Any opposite-half point;x (cm);y (cm);min |#Delta#phi| (deg)",
                nbx, -hx, hx, nby, hy_lo, hy_hi);
    TH2D h_sep("h_min_dphi_no_overlap",
               Form("Opposite point at #DeltaR #geq %.2f;x (cm);y (cm);min |#Delta#phi| (deg)", min_dR),
               nbx, -hx, hx, nby, hy_lo, hy_hi);

    double free_min = 1e9, free_max = -1e9, sep_min = 1e9, sep_max = -1e9;
    int n_points = 0;
    int n_no_solution = 0;

    // x > 0 half only; the mirrored bin gets the same value since the
    // geometry is mirror-symmetric in x.
    for (int bx = 1; bx <= nbx; bx++)
    for (int by = 1; by <= nby; by++)
    {
        SurfacePoint p = {h_free.GetXaxis()->GetBinCenter(bx), h_free.GetYaxis()->GetBinCenter(by), 0, 0};
        if (p.x < 0 || !pass_fiducial_cut((float)p.x, (float)p.y, buf, rect))
            continue;

        computeEtaPhi(p.x, p.y, p.eta, p.phi);
        n_points++;

        const int bx_mirror = nbx + 1 - bx;

        double best_free = 1e9;
        double best_sep = 1e9;

        for (size_t j = 0; j < far_pts.size(); j++)
        {
            const SurfacePoint &q = far_pts[j];

            double dphi = std::fabs(DeltaPhi(p.phi, q.phi));
            if (dphi < best_free)
                best_free = dphi;

            double deta = p.eta - q.eta;
            if (dphi < best_sep && std::sqrt(deta * deta + dphi * dphi) >= min_dR)
                best_sep = dphi;
        }

        best_free *= 180.0 / M_PI;
        free_min = std::min(free_min, best_free);
        free_max = std::max(free_max, best_free);

        h_free.SetBinContent(bx, by, best_free);
        h_free.SetBinContent(bx_mirror, by, best_free);

        if (best_sep > 1e8)
        {
            n_no_solution++;
            continue;
        }

        best_sep *= 180.0 / M_PI;
        sep_min = std::min(sep_min, best_sep);
        sep_max = std::max(sep_max, best_sep);

        h_sep.SetBinContent(bx, by, best_sep);
        h_sep.SetBinContent(bx_mirror, by, best_sep);
    }

    std::cout << "Fiducial region (one half): x in [" << x_lo << ", " << x_hi << "] cm, y in ["
              << y_lo << ", " << y_hi << "] cm, " << n_points << " map points at "
              << grid_spacing << " cm spacing" << std::endl;
    std::cout << "min |dphi|, any opposite point:      " << free_min << " to " << free_max
              << " deg" << std::endl;
    std::cout << "min |dphi|, opposite point dR >= " << min_dR << ": " << sep_min << " to "
              << sep_max << " deg (" << n_no_solution << " points with no allowed opposite point)"
              << std::endl;

    // Central point of the half, for a quick sanity number.
    {
        SurfacePoint c = {0.5 * (x_lo + x_hi), 0.0, 0, 0};
        std::cout << "At (x, y) = (" << c.x << ", 0): any = "
                  << h_free.GetBinContent(h_free.FindBin(c.x, c.y)) << " deg, no overlap = "
                  << h_sep.GetBinContent(h_sep.FindBin(c.x, c.y)) << " deg" << std::endl;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(53);
    gStyle->SetNumberContours(100);

    // Shared colour scale so the two panels can be compared directly; the
    // lower bound sits just below the smallest value so no filled cell is
    // drawn as empty.
    double z_lo = std::floor(std::min(free_min, sep_min)) - 1.0;
    double z_hi = std::ceil(std::max(free_max, sep_max)) + 1.0;
    h_free.GetZaxis()->SetRangeUser(z_lo, z_hi);
    h_sep.GetZaxis()->SetRangeUser(z_lo, z_hi);

    TCanvas c("c", "", 1500, 650);
    c.Divide(2, 1);

    TH2D *hists[2] = {&h_free, &h_sep};
    for (int i = 0; i < 2; i++)
    {
        c.cd(i + 1);
        gPad->SetRightMargin(0.17);
        gPad->SetLeftMargin(0.11);
        hists[i]->GetZaxis()->SetTitleOffset(1.3);
        hists[i]->GetYaxis()->SetTitleOffset(1.2);
        hists[i]->Draw("colz");
        DrawOutlines();
    }

    c.SaveAs("min_opposite_angle.pdf");
    c.SaveAs("min_opposite_angle.png");

    TFile fout("min_opposite_angle.root", "RECREATE");
    h_free.Write();
    h_sep.Write();
    fout.Close();

    return 0;
}
