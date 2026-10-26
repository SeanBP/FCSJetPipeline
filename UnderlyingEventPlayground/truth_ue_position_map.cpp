// Truth-level counterpart of ue_position_map.cpp: the mean truth UE estimate
// (truth_ue_E, same-event off-axis cones, written by UnderlyingEvent) against
// the truth jet's position on the detector face, in ranges of truth jet pT.
//
// Same selection and binning as the reco map: cross-section-weighted jets in
// EM-triggered events (fcsEM0/1/2; "all" as the third argument drops this),
// 6 cm bins, at least min_jets effective jets per bin, one colour scale. The
// map covers the truth jet acceptance (the fiducial region dilated by R/2 in
// eta-phi), with the fiducial boxes and that boundary drawn.
//
// Build and run from this folder:
//   g++ -m64 -O2 truth_ue_position_map.cpp ../shared/VectorDict.cxx -o truth_ue_position_map -I../shared -I../jet_scale \
//       `root-config --cflags --libs` `fastjet-config --cxxflags` -std=c++11
//   ./truth_ue_position_map <jet file list> [output_base] [em|all] [bin_cm] [min_jets]
//
// Writes <output_base>.pdf/.png/.root (default truth_ue_position_map).

#include <TROOT.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TLatex.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2)
    {
        std::cerr << "usage: ./truth_ue_position_map <jet file list> [output_base] [em|all] [bin_cm] [min_jets]" << std::endl;
        return 1;
    }
    const std::string out_base = argc > 2 ? argv[2] : "truth_ue_position_map";
    const bool em_only_events = !(argc > 3 && std::string(argv[3]) == "all");
    const double bin_cm = argc > 4 ? std::atof(argv[4]) : 6.0;
    const int min_jets = argc > 5 ? std::atoi(argv[5]) : 5;

    std::vector<std::string> files;
    {
        std::ifstream in(argv[1]);
        std::string s;
        while (in >> s) files.push_back(s);
    }
    std::vector<float> file_weight;
    if (!ComputeXSecWeights(files, "jetTree", file_weight)) return 1;

    const int kNPt = 6;
    const double kLo[kNPt] = {2.0, 3.0, 4.0, 5.0, 6.0, 2.0}, kHi[kNPt] = {3.0, 4.0, 5.0, 6.0, 12.0, 12.0}; // the last is all of them

    // Fiducial region and the truth jet acceptance (signed eta-phi distance >= -R/2)
    const FiducialRect &rect = kFiducialRectEcal;
    const double buf = reco_fiducial_buffer_ecal;
    const double x_in = rect.x_inner + buf, x_hi = rect.x_outer - buf, y_lo = rect.y_min + buf, y_hi = rect.y_max - buf;
    double ext_x = 0, ext_y_lo = 0, ext_y_hi = 0;
    for (double x = 0; x <= 260; x += 1.0)
        for (double y = -260; y <= 260; y += 1.0)
            if (etaPhiSignedDistance(x, y, x_in, x_hi, y_lo, y_hi) >= -truth_jet_margin_em)
            {
                ext_x = std::max(ext_x, x); ext_y_lo = std::min(ext_y_lo, y); ext_y_hi = std::max(ext_y_hi, y);
            }
    const int nx = 2 * (int)std::ceil(ext_x / bin_cm), ny = (int)std::ceil((ext_y_hi - ext_y_lo) / bin_cm);
    const double xr = 0.5 * nx * bin_cm, y0 = ext_y_lo;

    TFile fout((out_base + ".root").c_str(), "RECREATE");
    TH2D *h_sum[kNPt], *h_n[kNPt], *h_mean[kNPt];
    for (int k = 0; k < kNPt; k++)
    {
        h_sum[k] = new TH2D(Form("ue_sum_%d", k), "", nx, -xr, xr, ny, y0, y0 + ny * bin_cm);
        h_n[k] = new TH2D(Form("n_jets_%d", k), ";truth jet x (cm);truth jet y (cm);jets", nx, -xr, xr, ny, y0, y0 + ny * bin_cm);
        h_n[k]->Sumw2();
        h_mean[k] = new TH2D(Form("ue_mean_%d", k), ";Truth jet x (cm);Truth jet y (cm);Mean truth UE energy per cone (GeV)",
                             nx, -xr, xr, ny, y0, y0 + ny * bin_cm);
    }
    double tot_w[kNPt] = {0}, tot_wx[kNPt] = {0};
    long tot_n[kNPt] = {0};

    for (size_t f = 0; f < files.size(); f++)
    {
        if (file_weight[f] < 0) continue;
        TFile *fin = TFile::Open(files[f].c_str());
        TTree *t = (fin && !fin->IsZombie()) ? (TTree *)fin->Get("jetTree") : 0;
        if (!t || !t->GetBranch("truth_ue_E") || !t->GetBranch("truth_ue_ncone"))
        {
            std::cerr << "skipping " << files[f] << " (no jetTree or no truth_ue_E/truth_ue_ncone)" << std::endl;
            if (fin) { fin->Close(); delete fin; }
            continue;
        }
        std::vector<float> *tx = 0, *ty = 0, *tpt = 0, *tue = 0;
        std::vector<int> *tnc = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        t->SetBranchStatus("*", 0);
        const char *br[] = {"truth_x", "truth_y", "truth_pt", "truth_ue_E", "truth_ue_ncone", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int i = 0; i < 8; i++) t->SetBranchStatus(br[i], 1);
        t->SetBranchAddress("truth_x", &tx); t->SetBranchAddress("truth_y", &ty); t->SetBranchAddress("truth_pt", &tpt);
        t->SetBranchAddress("truth_ue_E", &tue); t->SetBranchAddress("truth_ue_ncone", &tnc);
        t->SetBranchAddress("fcsEM0", &em0); t->SetBranchAddress("fcsEM1", &em1); t->SetBranchAddress("fcsEM2", &em2);
        const double w = file_weight[f];
        for (Long64_t ev = 0; ev < t->GetEntries(); ev++)
        {
            t->GetEntry(ev);
            if (em_only_events && !(em0 || em1 || em2)) continue;
            for (size_t j = 0; j < tpt->size(); j++)
            {
                if (tnc->at(j) <= 0) continue;
                const double pt = tpt->at(j);
                for (int k = 0; k < kNPt; k++)
                {
                    if (pt < kLo[k] || pt >= kHi[k]) continue;
                    h_sum[k]->Fill(tx->at(j), ty->at(j), w * tue->at(j));
                    h_n[k]->Fill(tx->at(j), ty->at(j), w);
                    tot_w[k] += w; tot_wx[k] += w * tue->at(j); tot_n[k]++;
                }
            }
        }
        fin->Close();
        delete fin;
    }

    for (int k = 0; k < kNPt; k++)
    {
        int n_filled = 0;
        for (int bx = 1; bx <= nx; bx++)
            for (int by = 1; by <= ny; by++)
            {
                const double n = h_n[k]->GetBinContent(bx, by), e = h_n[k]->GetBinError(bx, by);
                if (n <= 0 || e <= 0 || n * n / (e * e) < min_jets) continue; // effective number of jets
                h_mean[k]->SetBinContent(bx, by, std::max(h_sum[k]->GetBinContent(bx, by) / n, 1e-4));
                n_filled++;
            }
        printf("truth jet pT %4.0f-%-4.0f GeV: %7ld jets, mean truth UE %.2f GeV, %d position bins filled\n", kLo[k], kHi[k], tot_n[k],
               tot_w[k] > 0 ? tot_wx[k] / tot_w[k] : 0.0, n_filled);
    }

    // The truth jet acceptance boundary, as a contour of the signed distance
    TH2D *hd = new TH2D("hd", "", 260, -xr, xr, 260, y0, y0 + ny * bin_cm);
    for (int i = 1; i <= hd->GetNbinsX(); i++)
        for (int j = 1; j <= hd->GetNbinsY(); j++)
            hd->SetBinContent(i, j, etaPhiSignedDistance(hd->GetXaxis()->GetBinCenter(i), hd->GetYaxis()->GetBinCenter(j),
                                                         x_in, x_hi, y_lo, y_hi));
    double level[1] = {-truth_jet_margin_em};
    hd->SetContour(1, level);
    hd->SetLineColor(kViolet + 1); hd->SetLineStyle(2); hd->SetLineWidth(3);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    gStyle->SetNumberContours(99);
    TCanvas c("c", "", 2700, 1700);
    c.Divide(3, 2);
    const double z_top = 6.0; // the reco map's colour scale
    for (int k = 0; k < kNPt; k++)
    {
        c.cd(k + 1);
        gPad->SetRightMargin(0.17); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.13);
        TH2D *h = h_mean[k];
        h->SetMinimum(0.0); h->SetMaximum(z_top);
        h->GetXaxis()->SetTitleSize(0.048); h->GetYaxis()->SetTitleSize(0.048); h->GetZaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.042); h->GetYaxis()->SetLabelSize(0.042); h->GetZaxis()->SetLabelSize(0.042);
        h->GetZaxis()->SetTitleOffset(1.15); h->GetYaxis()->SetTitleOffset(1.15);
        h->Draw("COLZ");
        for (int side = -1; side <= 1; side += 2)
        {
            TBox *fid = new TBox(side * x_in, y_lo, side * x_hi, y_hi);
            fid->SetFillStyle(0); fid->SetLineWidth(2); fid->Draw("l");
        }
        hd->Draw("cont3 same");
        TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.046);
        t.DrawLatex(0.12, 0.935, Form("Truth jet p_{T} %.0f#minus%.0f GeV", kLo[k], kHi[k]));
        t.SetTextFont(42); t.SetTextSize(0.036);
        t.DrawLatex(0.12, 0.885, Form("%ld jets, mean %.2f GeV", tot_n[k], tot_w[k] > 0 ? tot_wx[k] / tot_w[k] : 0.0));
    }
    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());

    fout.cd();
    for (int k = 0; k < kNPt; k++) { h_mean[k]->Write(); h_n[k]->Write(); }
    fout.Close();
    return 0;
}
