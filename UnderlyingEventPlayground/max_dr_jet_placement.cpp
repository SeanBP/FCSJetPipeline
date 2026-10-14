// For a point P on the FCS surface where the underlying event is to be
// assessed, finds where a jet axis J can be placed inside the EM-only reco
// fiducial region of the half opposite P so that
//   - a jet at J has ZERO overlap with a jet at P: dR(J, P) >= 2R (a hard
//     requirement -- a jet there would overlap the UE cone directly), and
//   - subject to that, J is as far as possible, in eta-phi, from the P'
//     line, the line 180 deg in azimuth from P at any eta (a jet back-to-back
//     with one at P can be anywhere along it), so that a jet at J overlaps
//     such a jet as little as possible:
//         d(P) = max over allowed J of |dphi(J, P' line)|,
//     the distance from J to the nearest point of the line, the one at J's
//     own eta.
// J is taken from the opposite half only so that d(P) varies continuously
// with P (see UESelection.h); a tag jet may be on either half.
//
// Two panels are drawn:
//   left  -- d(P) for every P in the fiducial region
//   right -- the worst case (the P with the smallest d(P)): R-radius circles
//            in eta-phi around P and the best J, projected onto the surface,
//            with the P' line and the band within R of it in azimuth.
//            Shaded: where else the jet axis may sit (still with zero overlap
//            with P) if its distance to the P' line is allowed to
//            fall short of the best value by up to 0.1 R, 0.2 R or 0.3 R (see
//            tolerances_in_R below).
//
// The fiducial region is the one reco EM-only jets are selected with
// (pass_fiducial_cut with kFiducialRectEcal / reco_fiducial_buffer_ecal),
// taken straight from shared/JetParameters.h.
//
// Build and run from this folder:
//   g++ -m64 max_dr_jet_placement.cpp -o max_dr_jet_placement -I. -I../shared \
//       `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//   ./max_dr_jet_placement [grid_spacing_cm]
//
// Writes max_dr_jet_placement.pdf/.png/.root in the current directory.

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TLegend.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "JetParameters.h"
#include "UESelection.h"

// Allowed shortfall, in units of R, of a jet axis's distance (to the P'
// line) relative to the best achievable one; one shaded region each.
// Ascending, so the regions are nested.
static const int n_tolerances = 3;
static const double tolerances_in_R[n_tolerances] = {0.1, 0.2, 0.3};

// Inverse of computeEtaPhi(): (eta, phi) -> position on the z_proj plane
inline void etaPhiToSurface(double eta, double phi, double &x, double &y)
{
    double r = z_proj / std::sinh(eta);
    x = r * std::cos(phi);
    y = r * std::sin(phi);
}

// Outlines of the ECAL-only rectangle (grey) and the reco fiducial region
// inside it (black), on both halves.
void DrawOutlines()
{
    const FiducialRect &rect = kFiducialRectEcal;
    const float buf = reco_fiducial_buffer_ecal;

    for (int side = -1; side <= 1; side += 2)
    {
        TBox *ecal = new TBox(side * rect.x_inner, rect.y_min, side * rect.x_outer, rect.y_max);
        ecal->SetFillStyle(0);
        ecal->SetLineColor(kGray + 1);
        ecal->Draw("same");

        TBox *fid = new TBox(side * (rect.x_inner + buf), rect.y_min + buf,
                             side * (rect.x_outer - buf), rect.y_max - buf);
        fid->SetFillStyle(0);
        fid->SetLineColor(kBlack);
        fid->Draw("same");
    }
}

// Circle of the given radius in eta-phi around (eta, phi), projected onto the
// surface. Returns the polyline (not yet drawn) and grows the given extent.
TPolyLine *ProjectedCircle(double eta, double phi, double radius, int color, double &extent)
{
    const int n = 360;
    std::vector<double> px(n + 1), py(n + 1);

    for (int i = 0; i <= n; i++)
    {
        double t = 2.0 * M_PI * i / n;
        etaPhiToSurface(eta + radius * std::cos(t), phi + radius * std::sin(t), px[i], py[i]);
        extent = std::max(extent, std::max(std::fabs(px[i]), std::fabs(py[i])));
    }

    TPolyLine *line = new TPolyLine(n + 1, &px[0], &py[0]);
    line->SetLineColor(color);
    line->SetLineWidth(2);
    return line;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);

    double grid_spacing = 1.0;
    if (argc > 1 && argv[1][0] != '\0')
        grid_spacing = std::atof(argv[1]);

    const FiducialRect &rect = kFiducialRectEcal;
    const float buf = reco_fiducial_buffer_ecal;

    const double x_lo = rect.x_inner + buf;
    const double x_hi = rect.x_outer - buf;
    const double y_lo = rect.y_min + buf;
    const double y_hi = rect.y_max - buf;

    // Candidate jet-axis points J on both halves (see UESelection.h)
    std::vector<SurfacePoint> jet_pts = FiducialJetAxisCandidates(grid_spacing);

    // Histogram spans both halves. The map is evaluated at the centre of
    // every bin that lies inside the fiducial region.
    const double hx = std::ceil(rect.x_outer / 10.0) * 10.0 + 10.0;
    const double hy_lo = std::floor(rect.y_min / 10.0) * 10.0 - 10.0;
    const double hy_hi = std::ceil(rect.y_max / 10.0) * 10.0 + 10.0;
    const int nbx = (int)std::round(2.0 * hx / grid_spacing);
    const int nby = (int)std::round((hy_hi - hy_lo) / grid_spacing);

    // A jet axis closer than this to P would overlap a jet at P
    const double min_dr_to_p = ue_min_dr_to_point;

    TH2D h_map("h_max_dr", Form("Best achievable jet distance, J on the half opposite P (R = %.1f);x (cm);y (cm);max |#Delta#phi(J,P' line)| / R   [#DeltaR(J,P) #geq 2R]", R),
               nbx, -hx, hx, nby, hy_lo, hy_hi);

    double map_min = 1e9, map_max = -1e9;
    int n_no_solution = 0;
    SurfacePoint worst_p = {0, 0, 0, 0}, worst_j = {0, 0, 0, 0};

    // x > 0 half only; the mirrored bin gets the same value since the
    // geometry is mirror-symmetric in x.
    for (int bx = 1; bx <= nbx; bx++)
    for (int by = 1; by <= nby; by++)
    {
        SurfacePoint p = {h_map.GetXaxis()->GetBinCenter(bx), h_map.GetYaxis()->GetBinCenter(by), 0, 0};
        if (p.x < 0 || !pass_fiducial_cut((float)p.x, (float)p.y, buf, rect))
            continue;

        computeEtaPhi(p.x, p.y, p.eta, p.phi);

        size_t best_j = 0;
        double best = BestJetDistance(p.eta, p.phi, jet_pts, &best_j);

        if (best < 0)
        {
            n_no_solution++;
            continue;
        }

        h_map.SetBinContent(bx, by, best);
        h_map.SetBinContent(nbx + 1 - bx, by, best);

        map_max = std::max(map_max, best);
        if (best < map_min)
        {
            map_min = best;
            worst_p = p;
            worst_j = jet_pts[best_j];
        }
    }

    // Largest step of the map between neighbouring bins (x > 0 half), as a
    // check that the best distance varies continuously with P
    double max_step = 0.0;
    for (int bx = 1; bx <= nbx; bx++)
    for (int by = 1; by <= nby; by++)
    {
        if (h_map.GetXaxis()->GetBinCenter(bx) < 0 || h_map.GetBinContent(bx, by) <= 0)
            continue;
        if (bx < nbx && h_map.GetBinContent(bx + 1, by) > 0)
            max_step = std::max(max_step, std::fabs(h_map.GetBinContent(bx + 1, by) - h_map.GetBinContent(bx, by)));
        if (by < nby && h_map.GetBinContent(bx, by + 1) > 0)
            max_step = std::max(max_step, std::fabs(h_map.GetBinContent(bx, by + 1) - h_map.GetBinContent(bx, by)));
    }
    std::cout << "Largest step between neighbouring points: " << max_step / R << " R" << std::endl;

    // P' line: 180 deg in azimuth from the worst-case P
    const double back_phi = worst_p.phi + M_PI;

    double dr_to_p = DeltaR(worst_j.eta, worst_j.phi, worst_p.eta, worst_p.phi);
    double dr_to_back = JetDistance(worst_j.phi, worst_p.phi);

    std::cout << "Fiducial region (one half): x in [" << x_lo << ", " << x_hi << "] cm, y in ["
              << y_lo << ", " << y_hi << "] cm, " << grid_spacing << " cm spacing" << std::endl;
    std::cout << "Best achievable |dphi(J, P' line)| over the surface, with dR(J, P) >= " << min_dr_to_p << ": "
              << map_min << " to " << map_max << " (" << n_no_solution
              << " points with no allowed jet axis)" << std::endl;
    std::cout << "Worst case:" << std::endl;
    std::cout << "  UE point P       (x, y) = (" << worst_p.x << ", " << worst_p.y << ") cm, eta = "
              << worst_p.eta << ", phi = " << worst_p.phi * 180.0 / M_PI << " deg" << std::endl;
    std::cout << "  P' line at phi = " << DeltaPhi(back_phi, 0.0) * 180.0 / M_PI << " deg" << std::endl;
    std::cout << "  best jet axis J  (x, y) = (" << worst_j.x << ", " << worst_j.y << ") cm, eta = "
              << worst_j.eta << ", phi = " << worst_j.phi * 180.0 / M_PI << " deg" << std::endl;
    std::cout << "  dR(J, P) = " << dr_to_p << " (" << dr_to_p / R << " R), |dphi(J, P' line)| = " << dr_to_back
              << " (" << dr_to_back / R << " R)" << std::endl;

    // Tolerance regions for the worst case: every fiducial bin (both halves)
    // whose centre, taken as the jet axis, has zero overlap with P and whose
    // distance to the P' line is within the allowed shortfall of
    // the best one. Filled with 1 so the "box" draw option gives full cells.
    TH2D *h_tol[n_tolerances];
    int tol_bins[n_tolerances];
    int fid_bins = 0;

    for (int k = 0; k < n_tolerances; k++)
    {
        h_tol[k] = new TH2D(Form("h_tolerance_%d", k),
                            Form("Jet axis within %.1f R of the best distance;x (cm);y (cm)", tolerances_in_R[k]),
                            nbx, -hx, hx, nby, hy_lo, hy_hi);
        tol_bins[k] = 0;
    }

    for (int bx = 1; bx <= nbx; bx++)
    for (int by = 1; by <= nby; by++)
    {
        SurfacePoint j = {h_map.GetXaxis()->GetBinCenter(bx), h_map.GetYaxis()->GetBinCenter(by), 0, 0};
        if (!pass_fiducial_cut((float)j.x, (float)j.y, buf, rect))
            continue;

        computeEtaPhi(j.x, j.y, j.eta, j.phi);
        fid_bins++;

        double d_p = DeltaR(j.eta, j.phi, worst_p.eta, worst_p.phi);
        if (d_p < min_dr_to_p)
            continue;

        double d = JetDistance(j.phi, worst_p.phi);

        for (int k = 0; k < n_tolerances; k++)
        {
            if (d < std::min(map_min - tolerances_in_R[k] * R, ue_tag_distance_cap))
                continue;

            h_tol[k]->SetBinContent(bx, by, 1.0);
            tol_bins[k]++;
        }
    }

    std::cout << "  Jet-axis regions allowed by a tolerance on |dphi(J, P' line)|:" << std::endl;
    for (int k = 0; k < n_tolerances; k++)
    {
        double d_min = std::min(map_min - tolerances_in_R[k] * R, ue_tag_distance_cap);

        // Overlap of two R-radius circles whose centres are d_min apart, as a
        // fraction of one circle's area -- the largest overlap this
        // tolerance admits.
        double overlap = 0.0;
        if (d_min < 2.0 * R)
            overlap = (2.0 * R * R * std::acos(d_min / (2.0 * R)) -
                       0.5 * d_min * std::sqrt(4.0 * R * R - d_min * d_min)) / (M_PI * R * R);

        std::cout << "    " << tolerances_in_R[k] << " R: distance >= " << d_min << ", area = "
                  << tol_bins[k] * grid_spacing * grid_spacing << " cm^2 ("
                  << 100.0 * tol_bins[k] / fid_bins << "% of the fiducial region), max cone overlap with a jet on the P' line = "
                  << 100.0 * overlap << "%" << std::endl;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(53);
    gStyle->SetNumberContours(100);

    TCanvas c("c", "", 1500, 850);
    c.Divide(2, 1);

    // ---- left: map of the best achievable distance, in units of R ----
    // (everything drawn is in units of R; h_map is rescaled here, after all
    // the calculations above that use it in absolute dR)
    h_map.Scale(1.0 / R);

    c.cd(1);
    gPad->SetRightMargin(0.17);
    gPad->SetLeftMargin(0.11);
    h_map.GetZaxis()->SetRangeUser(std::floor(map_min / R * 10.0) / 10.0 - 0.05, std::ceil(map_max / R * 10.0) / 10.0);
    h_map.GetZaxis()->SetTitleOffset(1.3);
    h_map.GetYaxis()->SetTitleOffset(1.2);
    h_map.Draw("colz");
    DrawOutlines();

    TMarker worst_marker(worst_p.x, worst_p.y, 29);
    worst_marker.SetMarkerColor(kAzure + 1);
    worst_marker.SetMarkerSize(2.0);
    worst_marker.Draw("same");

    // ---- right: the worst case ----
    c.cd(2);
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.04);

    const int col_p = kRed + 1, col_back = kAzure + 1, col_j = kGreen + 2;

    double extent = hx;
    TPolyLine *circle_p = ProjectedCircle(worst_p.eta, worst_p.phi, R, col_p, extent);
    TPolyLine *circle_j = ProjectedCircle(worst_j.eta, worst_j.phi, R, col_j, extent);
    extent = std::ceil(extent / 10.0) * 10.0 + 10.0;

    // Equal cm-per-pixel on both axes, so the projected shapes are not
    // stretched: the y span follows from the x span and the frame's pixel
    // aspect ratio. The frame starts at -extent in y, leaving the spare
    // room at the top for the legend.
    const double frame_w = gPad->GetWw() * gPad->GetAbsWNDC() * (1.0 - gPad->GetLeftMargin() - gPad->GetRightMargin());
    const double frame_h = gPad->GetWh() * gPad->GetAbsHNDC() * (1.0 - gPad->GetTopMargin() - gPad->GetBottomMargin());
    const double y_span = 2.0 * extent * frame_h / frame_w;

    TH1 *frame = gPad->DrawFrame(-extent, -extent, extent, -extent + y_span);
    frame->SetTitle(Form("Worst case: R = %.1f jets projected on the surface;x (cm);y (cm)", R));
    frame->GetYaxis()->SetTitleOffset(1.5);

    // Widest tolerance first (lightest), so the tighter ones sit on top.
    const int tol_colors[n_tolerances] = {kGreen + 2, kGreen - 6, kGreen - 10};
    for (int k = n_tolerances - 1; k >= 0; k--)
    {
        h_tol[k]->SetFillColor(tol_colors[k]);
        h_tol[k]->SetLineColor(tol_colors[k]);
        h_tol[k]->Draw("box same");
    }

    // P' line (solid) and the edges of the band within R of it in azimuth
    // (dashed), the envelope of R-radius jets anywhere along the line. All
    // three are rays from the beam axis, clipped by the frame.
    TLine *line_back = 0, *line_band = 0;
    for (int k = -1; k <= 1; k++)
    {
        double phi = back_phi + k * R;
        const double cx = std::cos(phi), cy = std::sin(phi);

        // Length at which the ray reaches the frame edge
        double ray_len = 1e9;
        if (cx != 0.0) ray_len = std::min(ray_len, extent / std::fabs(cx));
        if (cy > 0.0)  ray_len = std::min(ray_len, (-extent + y_span) / cy);
        if (cy < 0.0)  ray_len = std::min(ray_len, extent / -cy);

        TLine *ray = new TLine(0.0, 0.0, ray_len * cx, ray_len * cy);
        ray->SetLineColor(col_back);
        ray->SetLineWidth(2);
        ray->SetLineStyle(k == 0 ? 1 : 2);
        ray->Draw("same");
        (k == 0 ? line_back : line_band) = ray;
    }

    DrawOutlines();
    circle_p->Draw("same");
    circle_j->Draw("same");

    TMarker m_p(worst_p.x, worst_p.y, 20);
    TMarker m_j(worst_j.x, worst_j.y, 20);
    m_p.SetMarkerColor(col_p);
    m_j.SetMarkerColor(col_j);
    m_p.Draw("same");
    m_j.Draw("same");

    TLegend leg(0.13, 0.67, 0.94, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.030);
    leg.AddEntry(circle_p, Form("UE point P (%.0f, %.0f) cm", worst_p.x, worst_p.y), "l");
    leg.AddEntry(line_back, "P' line: 180#circ in #phi from P, any #eta", "l");
    leg.AddEntry(line_band, "within R of the P' line", "l");
    leg.AddEntry(circle_j, Form("Best jet axis J (%.0f, %.0f) cm", worst_j.x, worst_j.y), "l");
    leg.AddEntry((TObject *)0, Form("#DeltaR(J,P) = %.2f R  (required #geq 2R)", dr_to_p / R), "");
    leg.AddEntry((TObject *)0, Form("|#Delta#phi(J,P' line)| = %.2f R", dr_to_back / R), "");
    leg.Draw();

    TLegend leg_tol(0.56, 0.12, 0.955, 0.26);
    leg_tol.SetBorderSize(0);
    leg_tol.SetFillStyle(0);
    leg_tol.SetTextSize(0.027);
    leg_tol.SetHeader("Jet axis allowed, shortfall up to");
    for (int k = 0; k < n_tolerances; k++)
        leg_tol.AddEntry(h_tol[k], Form("%.1f R  (distance #geq %.2f R)", tolerances_in_R[k],
                                        map_min / R - tolerances_in_R[k]), "f");
    leg_tol.Draw();

    c.SaveAs("max_dr_jet_placement.pdf");
    c.SaveAs("max_dr_jet_placement.png");

    TFile fout("max_dr_jet_placement.root", "RECREATE");
    h_map.Write();
    for (int k = 0; k < n_tolerances; k++)
        h_tol[k]->Write();
    fout.Close();

    return 0;
}
