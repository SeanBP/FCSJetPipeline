// Distribution of the per-jet reco UE estimate (reco_ue_E, as stored by
// UnderlyingEvent) for reco jets near the beam (eta above a cut), simulation
// (cross-section weighted) against data, each normalized to unit area; linear
// and log y. Reco jets of 2-12 GeV pT in EM-triggered events with a UE
// estimate, as ue_position_compare.cpp. With trig_reweight = 1, simulation
// events are also weighted to the data's EM trigger mix (TriggerReweight.h).
//
// Build and run from this folder:
//   g++ -m64 -O2 ue_inner_distribution.cpp ../shared/VectorDict.cxx -o ue_inner_distribution -I. -I../shared -I../jet_scale \
//       `root-config --cflags --libs` `fastjet-config --cxxflags` -std=c++11
//   ./ue_inner_distribution <sim jet file list> <data jet file list> [output_base] [eta_min=3.2] [trig_reweight=0]
//
// Writes <output_base>.pdf/.png/.root (default ue_inner_distribution).

#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
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

// trig_w: per-trigger-pattern weight factors, empty for none
static long FillHist(const std::vector<std::string> &files, const std::vector<float> &weight, double eta_min, TH1D *h,
                     const std::vector<double> &trig_w)
{
    long n_jets = 0;
    for (size_t f = 0; f < files.size(); f++)
    {
        if (weight[f] < 0) continue;
        TFile *fin = TFile::Open(files[f].c_str());
        TTree *t = (fin && !fin->IsZombie()) ? (TTree *)fin->Get("jetTree") : 0;
        if (!t || !t->GetBranch("reco_ue_E") || !t->GetBranch("reco_ue_nproxy"))
        {
            std::cerr << "skipping " << files[f] << " (no jetTree or no reco_ue_E/reco_ue_nproxy)" << std::endl;
            if (fin) { fin->Close(); delete fin; }
            continue;
        }
        std::vector<float> *eta = 0, *pt = 0, *ue = 0;
        std::vector<int> *np = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        const char *br[] = {"reco_eta", "reco_pt", "reco_ue_E", "reco_ue_nproxy", "fcsEM0", "fcsEM1", "fcsEM2"};
        t->SetBranchStatus("*", 0);
        for (int i = 0; i < 7; i++) t->SetBranchStatus(br[i], 1);
        t->SetBranchAddress("reco_eta", &eta); t->SetBranchAddress("reco_pt", &pt);
        t->SetBranchAddress("reco_ue_E", &ue); t->SetBranchAddress("reco_ue_nproxy", &np);
        t->SetBranchAddress("fcsEM0", &em0); t->SetBranchAddress("fcsEM1", &em1); t->SetBranchAddress("fcsEM2", &em2);
        for (Long64_t ev = 0; ev < t->GetEntries(); ev++)
        {
            t->GetEntry(ev);
            if (!(em0 || em1 || em2)) continue;
            const double we = trig_w.empty() ? weight[f] : weight[f] * trig_w[TriggerLevel(em0, em1, em2)];
            for (size_t j = 0; j < pt->size(); j++)
            {
                if (np->at(j) <= 0 || pt->at(j) < kPtLo || pt->at(j) >= kPtHi || eta->at(j) < eta_min) continue;
                h->Fill(ue->at(j), we);
                n_jets++;
            }
        }
        fin->Close();
        delete fin;
    }
    return n_jets;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3)
    {
        std::cerr << "usage: ./ue_inner_distribution <sim jet file list> <data jet file list> [output_base] [eta_min] [trig_reweight]" << std::endl;
        return 1;
    }
    const std::vector<std::string> sim_files = ReadList(argv[1]), data_files = ReadList(argv[2]);
    const std::string out_base = argc > 3 ? argv[3] : "ue_inner_distribution";
    const double eta_min = argc > 4 ? std::atof(argv[4]) : 3.2;
    const bool trig_reweight = argc > 5 && std::atoi(argv[5]) == 1;

    std::vector<float> sim_weight, data_weight(data_files.size(), 1.0f);
    if (!ComputeXSecWeights(sim_files, "jetTree", sim_weight)) return 1;
    std::vector<double> trig_w;
    if (trig_reweight) trig_w = TriggerReweightFactors(TriggerFractions(sim_files, sim_weight), TriggerFractions(data_files, data_weight));

    TFile fout((out_base + ".root").c_str(), "RECREATE");
    const char *xt = ";Reco jet UE estimate per cone (GeV);Fraction of jets per 0.2 GeV";
    TH1D *h_sim = new TH1D("ue_sim", xt, 60, 0.0, 12.0), *h_data = new TH1D("ue_data", xt, 60, 0.0, 12.0);
    h_sim->Sumw2(); h_data->Sumw2();
    const long n_sim = FillHist(sim_files, sim_weight, eta_min, h_sim, trig_w);
    const long n_data = FillHist(data_files, data_weight, eta_min, h_data, std::vector<double>());
    TH1D *h[2] = {h_sim, h_data};
    const char *lab[2] = {trig_reweight ? "Sim. (trigger reweighted)" : "Simulation", "Data"};
    const long n[2] = {n_sim, n_data};
    for (int k = 0; k < 2; k++)
    {
        // Overflow included in the normalization and the moments
        const double tot = h[k]->Integral(0, h[k]->GetNbinsX() + 1);
        printf("%-36s eta > %.2f: %8ld jets, mean %.2f GeV, RMS %.2f GeV, skewness %.2f, overflow %.2f%%\n", lab[k], eta_min, n[k],
               h[k]->GetMean(), h[k]->GetRMS(), h[k]->GetSkewness(), tot > 0 ? 100.0 * h[k]->GetBinContent(h[k]->GetNbinsX() + 1) / tot : 0.0);
        if (tot > 0) h[k]->Scale(1.0 / tot);
    }
    h_sim->SetLineColor(kAzure + 2); h_sim->SetMarkerColor(kAzure + 2); h_sim->SetMarkerStyle(21);
    h_data->SetLineColor(kBlack); h_data->SetMarkerColor(kBlack); h_data->SetMarkerStyle(20);

    gStyle->SetOptStat(0);
    TCanvas c("c", "", 2000, 800);
    c.Divide(2, 1);
    for (int p = 0; p < 2; p++)
    {
        c.cd(p + 1);
        gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.1); gPad->SetBottomMargin(0.12);
        gPad->SetLogy(p == 1);
        const double ymax = 1.3 * std::max(h_sim->GetMaximum(), h_data->GetMaximum());
        TH1D *hs = (TH1D *)h_sim->Clone(Form("ue_sim_pad%d", p)); // own copy per pad, so the y ranges stay separate
        hs->SetMaximum(p == 1 ? 3.0 * ymax : ymax);
        hs->SetMinimum(p == 1 ? 1e-5 : 0.0);
        hs->GetXaxis()->SetTitleSize(0.048); hs->GetYaxis()->SetTitleSize(0.048);
        hs->GetXaxis()->SetLabelSize(0.042); hs->GetYaxis()->SetLabelSize(0.042);
        hs->GetYaxis()->SetTitleOffset(1.3);
        hs->Draw("E1");
        h_data->Draw("E1 same");
        TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.045);
        t.DrawLatex(0.13, 0.93, Form("Reco jets, #eta > %.2f, p_{T} %.0f#minus%.0f GeV", eta_min, kPtLo, kPtHi));
        if (p == 0)
        {
            TLegend *leg = new TLegend(0.42, 0.72, 0.95, 0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
            for (int k = 0; k < 2; k++)
                leg->AddEntry(h[k], Form("%s: %ld jets, mean %.2f GeV", lab[k], n[k], h[k]->GetMean()), "lep");
            leg->Draw();
        }
    }
    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());
    fout.cd();
    h_sim->Write(); h_data->Write();
    fout.Close();
    return 0;
}
