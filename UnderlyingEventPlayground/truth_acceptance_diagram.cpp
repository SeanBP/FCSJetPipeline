// Diagram of the EM-only truth acceptance on the FCS face (ECal front plane),
// each boundary lettered: A the ECal, B the reco fiducial region, C B dilated
// by R/2 in eta-phi, D/E the innermost/outermost rings touching C (the truth
// jet eta band), F/G the rings R/2 beyond D/E, and J1/J2 the R cones of truth
// jets with axes on D/E. Constant-eta rings are circles on the face; the
// cones are circles in eta-phi, so distorted. Left: full view; right: zoom.
//
// Build and run from this folder:
//   g++ -m64 -O2 truth_acceptance_diagram.cpp -o truth_acceptance_diagram -I../shared \
//       `root-config --cflags --libs` `fastjet-config --cxxflags` -std=c++11
//   ./truth_acceptance_diagram [output_base]

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TBox.h>
#include <TEllipse.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMarker.h>

#include <cmath>
#include <string>
#include <vector>

#include "JetParameters.h"

// Radius on the face of a given eta (direction from the origin)
double RadiusOfEta(double eta) { return z_proj / std::sinh(eta); }

// Outline on the face of a circle of radius r in eta-phi around (eta0, phi0)
TGraph *EtaPhiCircle(double eta0, double phi0, double r)
{
    TGraph *g = new TGraph();
    for (int i = 0; i <= 360; i++)
    {
        const double t = 2.0 * M_PI * i / 360.0;
        const double eta = eta0 + r * std::cos(t), phi = phi0 + r * std::sin(t);
        const double rho = RadiusOfEta(eta);
        g->SetPoint(i, rho * std::cos(phi), rho * std::sin(phi));
    }
    return g;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    const std::string out = argc > 1 ? argv[1] : "diagrams/truth_acceptance";

    double eta_lo, eta_hi;
    EmFiducialEtaRange(eta_lo, eta_hi);
    const double half = 0.5 * R;
    const double eta_jet_in = eta_hi + half, eta_jet_out = eta_lo - half;   // rings touching the dilated boundary
    const double eta_far_in = eta_hi + R, eta_far_out = eta_lo - R;         // R/2 beyond those

    const FiducialRect &e = kFiducialRectEcal;
    const double b = reco_fiducial_buffer_ecal;

    // Signed eta-phi distance to the fiducial region on a grid; its -R/2
    // contour is the dilated boundary (both detector halves)
    const double ext = 320.0, step = 1.0;
    const int n = (int)(2 * ext / step);
    TH2D *hd = new TH2D("hd", "", n, -ext, ext, n, -ext, ext);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
        {
            const double x = hd->GetXaxis()->GetBinCenter(i), y = hd->GetYaxis()->GetBinCenter(j);
            hd->SetBinContent(i, j, etaPhiSignedDistance(x, y, e.x_inner + b, e.x_outer - b, e.y_min + b, e.y_max - b));
        }
    double level[1] = {-half};
    hd->SetContour(1, level);
    hd->SetLineStyle(2);
    hd->SetLineWidth(3);

    // A truth jet on each jet ring: inner at phi = 0, outer towards the
    // farthest fiducial corner
    const double phi_out = std::atan2(e.y_min + b, e.x_outer - b);
    TGraph *cone_in = EtaPhiCircle(eta_jet_in, 0.0, R);
    TGraph *cone_out = EtaPhiCircle(eta_jet_out, phi_out, R);
    for (TGraph *g : {cone_in, cone_out}) g->SetLineWidth(3);

    // Labelled boundaries: A ECal, B fiducial region, C its R/2 dilation,
    // D/E inner/outer jet rings, F/G inner/outer particle rings, J1/J2 the
    // jet cones (circles in eta-phi, so distorted on the face)
    const int colA = kGray + 2, colB = kAzure + 1, colC = kViolet + 1, colD = kGreen + 2, colE = kTeal + 3,
              colF = kOrange + 7, colG = kMagenta + 2, colJ1 = kRed + 1, colJ2 = kPink + 9;
    hd->SetLineColor(colC);
    cone_in->SetLineColor(colJ1);
    cone_out->SetLineColor(colJ2);
    const double etas[4] = {eta_jet_in, eta_jet_out, eta_far_in, eta_far_out};
    const int ring_cols[4] = {colD, colE, colF, colG};
    const char *ring_names[4] = {"D", "E", "F", "G"};

    TCanvas c("c", "", 2000, 1000);
    c.Divide(2, 1, 0.002, 0.002);
    const double views[2] = {450.0, 160.0};
    std::vector<TObject *> keep;
    for (int v = 0; v < 2; v++)
    {
        c.cd(v + 1);
        gPad->SetLeftMargin(0.1); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.06); gPad->SetBottomMargin(0.1);
        const double w = views[v];
        TH1F *frame = gPad->DrawFrame(-w, -w, w, w);
        frame->GetXaxis()->SetTitle("x (cm)");
        frame->GetYaxis()->SetTitle("y (cm)");
        frame->GetYaxis()->SetTitleOffset(1.2);

        TBox *ecal_b = 0, *fid_b = 0;
        for (int s = -1; s <= 1; s += 2)
        {
            const double x1 = s * e.x_inner, x2 = s * e.x_outer;
            TBox *ecal = new TBox(std::min(x1, x2), e.y_min, std::max(x1, x2), e.y_max);
            ecal->SetFillColor(kGray); ecal->SetLineColor(colA); ecal->SetLineWidth(2); ecal->Draw("l");
            const double f1 = s * (e.x_inner + b), f2 = s * (e.x_outer - b);
            TBox *fid = new TBox(std::min(f1, f2), e.y_min + b, std::max(f1, f2), e.y_max - b);
            fid->SetFillColor(kAzure - 9); fid->SetLineColor(colB); fid->SetLineWidth(2); fid->Draw("l");
            keep.push_back(ecal); keep.push_back(fid);
            ecal_b = ecal; fid_b = fid;
        }
        hd->Draw("cont3 same");

        TEllipse *rings[4];
        for (int k = 0; k < 4; k++)
        {
            const double r = RadiusOfEta(etas[k]);
            rings[k] = new TEllipse(0, 0, r, r);
            rings[k]->SetFillStyle(0); rings[k]->SetLineColor(ring_cols[k]); rings[k]->SetLineWidth(3);
            rings[k]->Draw();
            keep.push_back(rings[k]);
        }
        cone_in->Draw("l");
        cone_out->Draw("l");

        // Jet axes
        const double ax[2] = {RadiusOfEta(eta_jet_in), RadiusOfEta(eta_jet_out) * std::cos(phi_out)};
        const double ay[2] = {0.0, RadiusOfEta(eta_jet_out) * std::sin(phi_out)};
        for (int k = 0; k < 2; k++)
        {
            TMarker *m = new TMarker(ax[k], ay[k], 20);
            m->SetMarkerColor(k == 0 ? colJ1 : colJ2); m->SetMarkerSize(1.4); m->Draw();
            keep.push_back(m);
        }
        TMarker *beam = new TMarker(0, 0, 5);
        beam->SetMarkerSize(2); beam->Draw();
        keep.push_back(beam);

        // Letter tags next to each boundary
        TLatex *tag = new TLatex();
        tag->SetTextFont(62); tag->SetTextSize(0.04);
        auto Tag = [&](const char *txt, int col, double x, double y) {
            if (std::fabs(x) > 0.97 * w || std::fabs(y) > 0.97 * w) return;
            tag->SetTextColor(col); tag->DrawLatex(x, y, txt);
        };
        Tag("A", colA, e.x_outer + 2, e.y_max - 12);
        Tag("B", colB, e.x_outer - b - 14, e.y_max - b - 12);
        Tag("C", colC, 0.0, 0.0 + RadiusOfEta(eta_lo - half) * 0.0 + (v == 0 ? 125 : 108));
        const double tag_phi[4] = {2.2, 2.2, 2.6, 2.6};
        for (int k = 0; k < 4; k++)
        {
            const double r = RadiusOfEta(etas[k]) + (k == 0 || k == 2 ? -0.0 : 0.0);
            const double ph = (k == 0 || k == 2) ? (k == 0 ? 1.2 : 2.0) : tag_phi[k];
            Tag(ring_names[k], ring_cols[k], r * std::cos(ph) + 3, r * std::sin(ph) + 3);
        }
        Tag("J1", colJ1, ax[0] + 8, ay[0] + 25);
        Tag("J2", colJ2, ax[1] + 20, ay[1] + 30);
        keep.push_back(tag);

        TLatex *t = new TLatex();
        t->SetNDC(); t->SetTextFont(62); t->SetTextSize(0.035);
        t->DrawLatex(0.11, 0.955, v == 0 ? "FCS face (z = 721.5 cm), full view" : "Zoom on the detector");
        keep.push_back(t);

        if (v == 0)
        {
            TLegend *leg = new TLegend(0.11, 0.11, 0.80, 0.40);
            leg->SetBorderSize(1); leg->SetFillColor(kWhite); leg->SetTextSize(0.022);
            leg->AddEntry(ecal_b, "A  ECal", "f");
            leg->AddEntry(fid_b, Form("B  reco fiducial region (#eta %.2f-%.2f)", eta_lo, eta_hi), "f");
            leg->AddEntry(hd, "C  B dilated by R/2 in #eta-#phi", "l");
            leg->AddEntry(rings[0], Form("D  inner ring touching C (#eta %.2f)", eta_jet_in), "l");
            leg->AddEntry(rings[1], Form("E  outer ring touching C (#eta %.2f)", eta_jet_out), "l");
            leg->AddEntry(rings[2], Form("F  R/2 inside D (#eta %.2f)", eta_far_in), "l");
            leg->AddEntry(rings[3], Form("G  R/2 outside E (#eta %.2f)", eta_far_out), "l");
            leg->AddEntry(cone_in, Form("J1  R cone of a jet with its axis on D (reaches #eta %.2f)", eta_jet_in + R), "l");
            leg->AddEntry(cone_out, Form("J2  R cone of a jet with its axis on E (reaches #eta %.2f)", eta_jet_out - R), "l");
            leg->Draw();
            keep.push_back(leg);
        }
    }
    c.SaveAs((out + ".png").c_str());
    c.SaveAs((out + ".pdf").c_str());
    return 0;
}
