// Maps the allowed proxy-jet region of the near off-axis proxy-event method
// with separation (../underlying_event/ProxySelection.h) over UE cone positions.
//
// For every UE cone position on a grid inside the EM-only reco fiducial
// region it finds the half-width of the azimuthal band about 90 deg that
// gives the allowed region the fixed area ProxyTargetArea(), and, for
// comparison, the area a fixed band of 90 +- fixed_deg gives.
//
// Build and run from this folder:
//   g++ -m64 -O2 proxy_region_area.cpp -o proxy_region_area -I. -I../shared -I../underlying_event `root-config --cflags --libs` -std=c++11
//   ./proxy_region_area [output_base] [grid_spacing_cm] [fixed_deg]
//
// Writes <output_base>.pdf/.png/.root (default proxy_region_area). The .root
// file holds half_width (deg) and area_reached against the UE position in cm,
// which serves as the lookup table for applying the selection.

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TLatex.h>
#include <TColor.h>

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>

#include "ProxySelection.h"

static double Percentile(std::vector<double> v, double p)
{
    std::sort(v.begin(), v.end());
    return v[(size_t)(p * (v.size() - 1))];
}

static void DrawDetector()
{
    const FiducialRect &rect = kFiducialRectEcal;
    const double buf = reco_fiducial_buffer_ecal;
    for (int side = -1; side <= 1; side += 2)
    {
        TBox *det = new TBox(side * rect.x_inner, rect.y_min, side * rect.x_outer, rect.y_max);
        det->SetFillStyle(0); det->SetLineWidth(2); det->Draw("l");
        TBox *fid = new TBox(side * (rect.x_inner + buf), rect.y_min + buf, side * (rect.x_outer - buf), rect.y_max - buf);
        fid->SetFillStyle(0); fid->SetLineWidth(1); fid->Draw("l");
    }
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    const std::string out_base = argc > 1 ? argv[1] : "proxy_region_area";
    const double spacing = argc > 2 ? std::atof(argv[2]) : 2.0;
    const double fixed_deg = argc > 3 ? std::atof(argv[3]) : 30.0;
    const double proxy_target_area = ProxyTargetArea();

    const FiducialRect &rect = kFiducialRectEcal;
    const double buf = reco_fiducial_buffer_ecal;
    const double x_lo = rect.x_inner + buf, x_hi = rect.x_outer - buf, y_lo = rect.y_min + buf, y_hi = rect.y_max - buf;
    const int nx = (int)((x_hi - x_lo) / spacing), ny = (int)((y_hi - y_lo) / spacing);
    const double dx = (x_hi - x_lo) / nx, dy = (y_hi - y_lo) / ny;

    // One histogram spanning both halves; bins in the gap between them stay empty
    const int nx_all = (int)(2.0 * x_hi / dx + 0.5);
    TFile fout((out_base + ".root").c_str(), "RECREATE");
    TH2D h_width("half_width", ";UE cone x (cm);UE cone y (cm);half-width about 90#circ (deg)", nx_all, -x_hi, x_hi, ny, y_lo, y_hi);
    TH2D h_area("area_reached", ";UE cone x (cm);UE cone y (cm);allowed area (#eta-#phi)", nx_all, -x_hi, x_hi, ny, y_lo, y_hi);
    TH2D h_fixed("area_fixed", ";UE cone x (cm);UE cone y (cm);allowed area (#eta-#phi)", nx_all, -x_hi, x_hi, ny, y_lo, y_hi);

    std::vector<double> widths, fixed_areas;
    int n_capped = 0, n_fixed_empty = 0;
    for (int side = -1; side <= 1; side += 2)
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
            {
                const double x = side * (x_lo + (i + 0.5) * dx), y = y_lo + (j + 0.5) * dy;
                double eta, phi;
                computeEtaPhi(x, y, eta, phi);

                double area = 0.0;
                const double width = ProxyHalfWidthDeg(eta, phi, &area);
                const double fixed_area = ProxyAllowedArea(eta, phi, fixed_deg);
                if (area < proxy_target_area) n_capped++;
                if (fixed_area <= 0.0) n_fixed_empty++;
                widths.push_back(width); fixed_areas.push_back(fixed_area);

                const int bin = h_width.FindBin(x, y);
                h_width.SetBinContent(bin, width);
                h_area.SetBinContent(bin, area);
                h_fixed.SetBinContent(bin, std::max(fixed_area, 1e-6)); // so a zero is drawn, not left blank
            }

    double total = 0.0;
    const std::vector<ProxyCell> &cells = ProxyFiducialCells();
    for (size_t i = 0; i < cells.size(); i++) total += cells[i].area;

    std::cout << "Reco fiducial region: eta-phi area " << total << "; target allowed area " << proxy_target_area
              << " (" << 100.0 * proxy_target_area_fraction << " % of one half); separation " << proxy_min_dr / R << " R" << std::endl;
    std::cout << widths.size() << " UE positions, " << spacing << " cm grid" << std::endl;
    std::cout << "Half-width for the target area (deg): min " << Percentile(widths, 0.0) << ", 5% " << Percentile(widths, 0.05)
              << ", median " << Percentile(widths, 0.5) << ", 95% " << Percentile(widths, 0.95) << ", max " << Percentile(widths, 1.0) << std::endl;
    std::cout << "UE positions short of the target at the " << proxy_half_width_max_deg << " deg cap: "
              << 100.0 * n_capped / widths.size() << " %" << std::endl;
    std::cout << "Fixed 90 +- " << fixed_deg << " deg band, allowed area: min " << Percentile(fixed_areas, 0.0) << ", 5% " << Percentile(fixed_areas, 0.05)
              << ", median " << Percentile(fixed_areas, 0.5) << ", 95% " << Percentile(fixed_areas, 0.95) << ", max " << Percentile(fixed_areas, 1.0)
              << "; no allowed region at all: " << 100.0 * n_fixed_empty / widths.size() << " %" << std::endl;

    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    gStyle->SetNumberContours(99);
    TCanvas c("c", "", 2000, 900);
    c.Divide(2, 1);
    TH2D *maps[2] = {&h_fixed, &h_width};
    const std::string titles[2] = {std::string(Form("Fixed band 90#circ #pm %.0f#circ: area varies", fixed_deg)),
                                   std::string(Form("Fixed area (%.0f%% of fiducial #eta#minus#phi area): band varies", 100.0 * proxy_target_area_fraction / 2.0))};
    for (int p = 0; p < 2; p++)
    {
        c.cd(p + 1);
        gPad->SetRightMargin(0.19); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.09);
        TH2D *frame = maps[p];
        frame->GetXaxis()->SetRangeUser(-rect.x_outer - 5.0, rect.x_outer + 5.0);
        frame->GetXaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetTitleSize(0.045); frame->GetZaxis()->SetTitleSize(0.045);
        frame->GetXaxis()->SetLabelSize(0.040); frame->GetYaxis()->SetLabelSize(0.040); frame->GetZaxis()->SetLabelSize(0.040);
        frame->GetZaxis()->SetTitleOffset(1.25); frame->GetYaxis()->SetTitleOffset(1.2);
        if (p == 0) frame->SetMinimum(0.0);
        else { frame->SetMinimum(10.0); frame->SetMaximum(proxy_half_width_max_deg); }
        frame->Draw("COLZ");
        DrawDetector();
        TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.042);
        t.DrawLatex(0.12, 0.935, titles[p].c_str());
    }
    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());

    h_width.Write(); h_area.Write(); h_fixed.Write();
    fout.Close();
    return 0;
}
