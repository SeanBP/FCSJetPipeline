// Mean UE estimate against jet position for all jets of 2-12 GeV pT, three
// panels side by side: simulation truth jets (truth_ue_E, same-event cones),
// simulation reco jets and data reco jets (reco_ue_E, proxy events), as
// stored by UnderlyingEvent. Jets in EM-triggered events (fcsEM0/1/2) with a
// UE estimate (but see truth_all_events); simulation jets cross-section
// weighted. Same bins, colour scale and minimum effective jets per bin as
// ue_position_map.cpp; every panel
// covers the truth jet acceptance (the fiducial region dilated by R/2 in
// eta-phi, drawn on the truth panel). With trig_reweight = 1, simulation
// events are also weighted to the data's EM trigger mix (TriggerReweight.h).
// With truth_all_events = 1, the truth panel takes every event, triggered or
// not (cross-section weights only): all stored truth jets.
//
// Build and run from this folder:
//   g++ -m64 -O2 ue_position_compare.cpp ../shared/VectorDict.cxx -o ue_position_compare -I. -I../shared -I../jet_scale \
//       `root-config --cflags --libs` `fastjet-config --cxxflags` -std=c++11
//   ./ue_position_compare <sim jet file list> <data jet file list> [output_base] [bin_cm] [min_jets] [trig_reweight=0] [truth_all_events=0]
//
// Writes <output_base>.pdf/.png/.root (default ue_position_compare).

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
#include "TriggerReweight.h"

static const double kPtLo = 2.0, kPtHi = 12.0;

static std::vector<std::string> ReadList(const char *path)
{
    std::vector<std::string> files;
    std::ifstream in(path);
    std::string s;
    while (in >> s) files.push_back(s);
    return files;
}

struct Panel
{
    TH2D *sum, *n, *mean;
    double tot_w, tot_wx;
    long tot_n;
};

// Fills a panel from jets of one level: prefix "truth" (needs truth_ue_ncone > 0)
// or "reco" (needs reco_ue_nproxy > 0)
// trig_w: per-trigger-level weight factors, empty for none; all_events: no
// EM trigger requirement
static void FillPanel(const std::vector<std::string> &files, const std::vector<float> &weight, bool truth, Panel &p,
                      const std::vector<double> &trig_w, bool all_events)
{
    const std::string pre = truth ? "truth" : "reco";
    const std::string count = truth ? "truth_ue_ncone" : "reco_ue_nproxy";
    for (size_t f = 0; f < files.size(); f++)
    {
        if (weight[f] < 0) continue;
        TFile *fin = TFile::Open(files[f].c_str());
        TTree *t = (fin && !fin->IsZombie()) ? (TTree *)fin->Get("jetTree") : 0;
        if (!t || !t->GetBranch((pre + "_ue_E").c_str()) || !t->GetBranch(count.c_str()))
        {
            std::cerr << "skipping " << files[f] << " (no jetTree or no " << pre << "_ue_E/" << count << ")" << std::endl;
            if (fin) { fin->Close(); delete fin; }
            continue;
        }
        std::vector<float> *x = 0, *y = 0, *pt = 0, *ue = 0;
        std::vector<int> *nc = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        const std::string br[] = {pre + "_x", pre + "_y", pre + "_pt", pre + "_ue_E", count, "fcsEM0", "fcsEM1", "fcsEM2"};
        t->SetBranchStatus("*", 0);
        for (int i = 0; i < 8; i++) t->SetBranchStatus(br[i].c_str(), 1);
        t->SetBranchAddress(br[0].c_str(), &x); t->SetBranchAddress(br[1].c_str(), &y); t->SetBranchAddress(br[2].c_str(), &pt);
        t->SetBranchAddress(br[3].c_str(), &ue); t->SetBranchAddress(br[4].c_str(), &nc);
        t->SetBranchAddress("fcsEM0", &em0); t->SetBranchAddress("fcsEM1", &em1); t->SetBranchAddress("fcsEM2", &em2);
        const double w = weight[f];
        for (Long64_t ev = 0; ev < t->GetEntries(); ev++)
        {
            t->GetEntry(ev);
            if (!all_events && !(em0 || em1 || em2)) continue;
            const double we = trig_w.empty() ? w : w * trig_w[TriggerLevel(em0, em1, em2)];
            for (size_t j = 0; j < pt->size(); j++)
            {
                if (nc->at(j) <= 0 || pt->at(j) < kPtLo || pt->at(j) >= kPtHi) continue;
                p.sum->Fill(x->at(j), y->at(j), we * ue->at(j));
                p.n->Fill(x->at(j), y->at(j), we);
                p.tot_w += we; p.tot_wx += we * ue->at(j); p.tot_n++;
            }
        }
        fin->Close();
        delete fin;
    }
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3)
    {
        std::cerr << "usage: ./ue_position_compare <sim jet file list> <data jet file list> [output_base] [bin_cm] [min_jets] [trig_reweight] [truth_all_events]" << std::endl;
        return 1;
    }
    const std::vector<std::string> sim_files = ReadList(argv[1]), data_files = ReadList(argv[2]);
    const std::string out_base = argc > 3 ? argv[3] : "ue_position_compare";
    const double bin_cm = argc > 4 ? std::atof(argv[4]) : 6.0;
    const int min_jets = argc > 5 ? std::atoi(argv[5]) : 5;
    const bool trig_reweight = argc > 6 && std::atoi(argv[6]) == 1;
    const bool truth_all_events = argc > 7 && std::atoi(argv[7]) == 1;

    std::vector<float> sim_weight, data_weight(data_files.size(), 1.0f);
    if (!ComputeXSecWeights(sim_files, "jetTree", sim_weight)) return 1;
    std::vector<double> trig_w;
    if (trig_reweight) trig_w = TriggerReweightFactors(TriggerFractions(sim_files, sim_weight), TriggerFractions(data_files, data_weight));

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
    const double xr = 0.5 * nx * bin_cm, y0 = ext_y_lo, y1 = y0 + ny * bin_cm;

    const int kNPanel = 3;
    const char *name[kNPanel] = {"sim_truth", "sim_reco", "data_reco"};
    const char *title[kNPanel] = {"Simulation, truth jets", "Simulation, reco jets", "Data, reco jets"};
    const char *axis[kNPanel] = {"Truth jet", "Reco jet", "Reco jet"};
    TFile fout((out_base + ".root").c_str(), "RECREATE");
    Panel p[kNPanel];
    for (int k = 0; k < kNPanel; k++)
    {
        p[k].sum = new TH2D(Form("ue_sum_%s", name[k]), "", nx, -xr, xr, ny, y0, y1);
        p[k].n = new TH2D(Form("n_jets_%s", name[k]), Form(";%s x (cm);%s y (cm);jets", axis[k], axis[k]), nx, -xr, xr, ny, y0, y1);
        p[k].n->Sumw2();
        p[k].mean = new TH2D(Form("ue_mean_%s", name[k]), Form(";%s x (cm);%s y (cm);Mean UE energy per cone (GeV)", axis[k], axis[k]),
                             nx, -xr, xr, ny, y0, y1);
        p[k].tot_w = p[k].tot_wx = 0; p[k].tot_n = 0;
    }
    FillPanel(sim_files, sim_weight, true, p[0], truth_all_events ? std::vector<double>() : trig_w, truth_all_events);
    FillPanel(sim_files, sim_weight, false, p[1], trig_w, false);
    FillPanel(data_files, data_weight, false, p[2], std::vector<double>(), false);

    for (int k = 0; k < kNPanel; k++)
    {
        int n_filled = 0;
        for (int bx = 1; bx <= nx; bx++)
            for (int by = 1; by <= ny; by++)
            {
                const double n = p[k].n->GetBinContent(bx, by), e = p[k].n->GetBinError(bx, by);
                if (n <= 0 || e <= 0 || n * n / (e * e) < min_jets) continue; // effective number of jets
                p[k].mean->SetBinContent(bx, by, std::max(p[k].sum->GetBinContent(bx, by) / n, 1e-4));
                n_filled++;
            }
        printf("%-24s pT %.0f-%.0f GeV: %8ld jets, mean UE %.2f GeV, %d position bins filled\n", title[k], kPtLo, kPtHi, p[k].tot_n,
               p[k].tot_w > 0 ? p[k].tot_wx / p[k].tot_w : 0.0, n_filled);
    }

    // The truth jet acceptance boundary, as a contour of the signed distance
    TH2D *hd = new TH2D("hd", "", 260, -xr, xr, 260, y0, y1);
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
    TCanvas c("c", "", 2700, 900);
    c.Divide(3, 1);
    const double z_top = 6.0; // as ue_position_map.cpp
    for (int k = 0; k < kNPanel; k++)
    {
        c.cd(k + 1);
        gPad->SetRightMargin(0.17); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.13);
        TH2D *h = p[k].mean;
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
        if (k == 0) hd->Draw("cont3 same");
        TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.046);
        t.DrawLatex(0.12, 0.935, Form("%s, p_{T} %.0f#minus%.0f GeV", title[k], kPtLo, kPtHi));
        t.SetTextFont(42); t.SetTextSize(0.036);
        const char *note = (k == 0 && truth_all_events) ? ", all events" : (trig_reweight && k < 2) ? ", trigger mix reweighted to data" : "";
        t.DrawLatex(0.12, 0.885, Form("%ld jets, mean %.2f GeV%s", p[k].tot_n, p[k].tot_w > 0 ? p[k].tot_wx / p[k].tot_w : 0.0, note));
    }
    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());

    fout.cd();
    for (int k = 0; k < kNPanel; k++) { p[k].mean->Write(); p[k].n->Write(); }
    fout.Close();
    return 0;
}
