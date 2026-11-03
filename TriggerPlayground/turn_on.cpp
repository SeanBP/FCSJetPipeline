// FCS EM trigger turn-ons in data against the uncorrected jet pT, from the
// hardware trigger bits (fcs_tcu_bits, StSimpleReaderMaker: the stage-3
// output before prescale; bits 3-5 = EM0-2).
//
// Jets: re-clustered here from the SimpleTree's ECal points as the pipeline
// does (anti-kT R, points with pass_point_energy_cut, pT from the point
// position), but with no pT minimum, so the turn-on is seen below the 2 GeV
// jet-finding cut. Uncorrected pT: no UE subtraction, no JES (the trigger
// sees raw energy). Per event: the leading jet anywhere on the ECal, and the
// leading jet inside the reco fiducial region (pass_fiducial_cut, 20 cm).
//
// Efficiency of EMk = events with bit k / events, per leading-jet pT bin, in
// a reference sample whose selection does not depend on condition k:
//   EM1, EM2: events recorded by fcsEM0 (listed, _tpc copy included); EM2
//             also in events recorded by fcsEM1
//   EM0, EM1, EM2: events recorded only by HCal-only FCS triggers (names
//             with HHT, HTOT or HcalTot, e.g. FCSS7_HTOT_S) and no other
//             FCS trigger
// Each turn-on is fitted with plateau x (1 + erf((pT - mu) / (sqrt2 sigma))) / 2
// from 0.5 to 10 GeV; printed: plateau, mu (50% point), sigma, 1.15 mu (the
// FMS rule: 15% above nominal), and the pT where the fit reaches 95% and 99%
// of its plateau.
//
// Plot: top row the leading jet, bottom row the leading fiducial jet, one
// column per reference sample; dashed lines at each fit's 99% point.
// Usage: ./turn_on <SimpleTree list> <out base>   (writes <out base>.root, .txt, .png)
// Build: g++ -m64 -O2 turn_on.cpp ../shared/VectorDict.cxx -o turn_on -I../shared -I../data_to_jet
//            `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include "JetParameters.h"
#include "FcsTriggerDefs.h"

using namespace fastjet;

static const int kNRef = 3;
static const char *const kRefName[kNRef] = {"fcsEM0 recorded", "fcsEM1 recorded", "HCal-only recorded"};
static const char *const kRefTag[kNRef] = {"refEM0", "refEM1", "refHcal"};
static const int kNJet = 2;
static const char *const kJetTag[kNJet] = {"lead", "leadfid"};
static const char *const kJetName[kNJet] = {"leading jet", "leading fiducial jet"};
static const int kNPt = 100;
static const double kPtMax = 10.0;

static int FlagIndex(const char *name)
{
    for (int i = 0; i < kNTrigFlags; i++)
        if (std::string(kTrigFlagName[i]) == name) return i;
    return -1;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::fprintf(stderr, "usage: ./turn_on <SimpleTree list> <out base>\n"); return 1; }
    const std::string base = argv[2];

    // flag indices: EMk and their _tpc copies; the HCal-only and all other FCS triggers
    int em[3], em_tpc[3];
    for (int k = 0; k < 3; k++)
    {
        em[k] = FlagIndex(Form("fcsEM%d", k));
        em_tpc[k] = FlagIndex(Form("fcsEM%d_tpc", k));
    }
    std::vector<int> is_hcal(kNTrigFlags, 0);
    for (int i = 0; i < kNTrigFlags; i++)
    {
        const std::string n = kTrigFlagName[i];
        is_hcal[i] = n.find("HHT") != std::string::npos || n.find("HTOT") != std::string::npos || n.find("HcalTot") != std::string::npos;
    }
    const int led = FlagIndex("fcs_led");
    is_hcal[led] = 0;

    TFile out((base + ".root").c_str(), "RECREATE");
    TH1D *hd[kNRef][kNJet], *hn[kNRef][kNJet][3];
    for (int r = 0; r < kNRef; r++)
        for (int j = 0; j < kNJet; j++)
        {
            hd[r][j] = new TH1D(Form("den_%s_%s", kRefTag[r], kJetTag[j]), Form("%s;uncorrected %s p_{T} (GeV);events", kRefName[r], kJetName[j]), kNPt, 0, kPtMax);
            for (int k = 0; k < 3; k++)
                hn[r][j][k] = new TH1D(Form("num_EM%d_%s_%s", k, kRefTag[r], kJetTag[j]), "", kNPt, 0, kPtMax);
        }

    std::ifstream in(argv[1]);
    std::string f;
    long n_ev = 0, n_hw = 0, n_ref[kNRef] = {0};
    while (in >> f)
    {
        TFile tf(f.c_str());
        TTree *t = (TTree *)tf.Get("data");
        if (!t || !t->GetBranch("fcs_tcu_bits")) { std::fprintf(stderr, "skipping %s (no fcs_tcu_bits)\n", f.c_str()); continue; }
        const int MAX_PTS = 200;
        Int_t np = 0, det[MAX_PTS], src = 0, flag[kNTrigFlags] = {0};
        Float_t pE[MAX_PTS], px[MAX_PTS], py[MAX_PTS], pz[MAX_PTS];
        UInt_t bits = 0;
        Double_t sqrt_s = kSqrtS_MainRun22_GeV;
        t->SetBranchStatus("*", 0);
        const char *b[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid", "fcs_tcu_bits", "fcs_tcu_source", "event_sqrt_s_gev"};
        for (int i = 0; i < 9; i++) t->SetBranchStatus(b[i], 1);
        t->SetBranchAddress("Cal_npts", &np); t->SetBranchAddress("Cal_pt_energy", pE); t->SetBranchAddress("Cal_pt_x", px);
        t->SetBranchAddress("Cal_pt_y", py); t->SetBranchAddress("Cal_pt_z", pz); t->SetBranchAddress("Cal_pt_detid", det);
        t->SetBranchAddress("fcs_tcu_bits", &bits); t->SetBranchAddress("fcs_tcu_source", &src);
        t->SetBranchAddress("event_sqrt_s_gev", &sqrt_s);
        for (int i = 0; i < kNTrigFlags; i++)
            if (t->GetBranch(kTrigFlagName[i])) { t->SetBranchStatus(kTrigFlagName[i], 1); t->SetBranchAddress(kTrigFlagName[i], &flag[i]); }

        for (Long64_t ev = 0; ev < t->GetEntries(); ev++)
        {
            t->GetEntry(ev);
            n_ev++;
            if (src != 1 || (led >= 0 && flag[led])) continue;
            n_hw++;
            bool ref[kNRef];
            ref[0] = flag[em[0]] || (em_tpc[0] >= 0 && flag[em_tpc[0]]);
            ref[1] = flag[em[1]] || (em_tpc[1] >= 0 && flag[em_tpc[1]]);
            bool any_hcal = false, any_other = false;
            for (int i = 0; i < kNTrigFlags; i++)
                if (flag[i]) { if (is_hcal[i]) any_hcal = true; else any_other = true; }
            ref[2] = any_hcal && !any_other;
            if (!ref[0] && !ref[1] && !ref[2]) continue;

            std::vector<PseudoJet> pts;
            for (int i = 0; i < std::min((int)np, MAX_PTS); i++)
            {
                if ((det[i] != 0 && det[i] != 1) || !pass_point_energy_cut(pE[i], sqrt_s)) continue;
                const double norm = std::sqrt((double)px[i] * px[i] + (double)py[i] * py[i] + (double)pz[i] * pz[i]);
                pts.emplace_back(pE[i] * px[i] / norm, pE[i] * py[i] / norm, pE[i] * pz[i] / norm, pE[i]);
            }
            double lead[kNJet] = {0, 0};
            if (!pts.empty())
            {
                ClusterSequence cs(pts, JetDefinition(antikt_algorithm, R));
                for (const PseudoJet &jet : cs.inclusive_jets(0.0))
                {
                    lead[0] = std::max(lead[0], jet.pt());
                    const float xe = z_proj * jet.px() / jet.pz(), ye = z_proj * jet.py() / jet.pz();
                    if (pass_fiducial_cut(xe, ye, reco_fiducial_buffer_ecal, kFiducialRectEcal)) lead[1] = std::max(lead[1], jet.pt());
                }
            }
            for (int r = 0; r < kNRef; r++)
            {
                if (!ref[r]) continue;
                n_ref[r]++;
                for (int j = 0; j < kNJet; j++)
                {
                    if (!(lead[j] > 0)) continue;
                    hd[r][j]->Fill(lead[j]);
                    for (int k = 0; k < 3; k++)
                        if ((bits >> (3 + k)) & 1) hn[r][j][k]->Fill(lead[j]);
                }
            }
        }
    }

    FILE *txt = std::fopen((base + ".txt").c_str(), "w");
    FILE *outs[2] = {stdout, txt};
    for (int o = 0; o < 2; o++)
        std::fprintf(outs[o], "turn_on: %ld events, %ld with hardware bits (no LED); reference events: %s %ld, %s %ld, %s %ld\n",
                     n_ev, n_hw, kRefName[0], n_ref[0], kRefName[1], n_ref[1], kRefName[2], n_ref[2]);
    gStyle->SetOptStat(0);
    TCanvas c("c", "c", 1500, 1000);
    c.Divide(3, 2);
    const int col[3] = {kBlue + 1, kRed + 1, kGreen + 2};
    for (int j = 0; j < kNJet; j++)
        for (int r = 0; r < kNRef; r++)
        {
            c.cd(j * kNRef + r + 1);
            gPad->SetGridy();
            TH1D *frame = new TH1D(Form("frame_%d_%d", j, r), Form("%s;uncorrected %s p_{T} (GeV);efficiency", kRefName[r], kJetName[j]), kNPt, 0, kPtMax);
            frame->SetMinimum(0); frame->SetMaximum(1.1);
            frame->Draw();
            TLegend *lg = new TLegend(0.40, 0.13, 0.89, 0.36);
            lg->SetHeader("dashed: 99% of plateau");
            for (int k = 0; k < 3; k++)
            {
                if ((r == 0 && k == 0) || (r == 1 && k < 2)) continue;  // the reference's own condition
                if (hd[r][j]->GetEntries() < 10) continue;
                TGraphAsymmErrors *g = new TGraphAsymmErrors(hn[r][j][k], hd[r][j], "cp");
                g->SetName(Form("eff_EM%d_%s_%s", k, kRefTag[r], kJetTag[j]));
                g->SetMarkerStyle(20); g->SetMarkerSize(0.6); g->SetMarkerColor(col[k]); g->SetLineColor(col[k]);
                TF1 *fit = new TF1(Form("fit_EM%d_%s_%s", k, kRefTag[r], kJetTag[j]), "[0]*0.5*(1+TMath::Erf((x-[1])/(sqrt(2.)*[2])))", 0.5, kPtMax);
                fit->SetParameters(1.0, 1.5 + k, 0.5);
                fit->SetParLimits(0, 0.0, 1.05); fit->SetParLimits(1, 0.0, 10.0); fit->SetParLimits(2, 0.05, 5.0);
                fit->SetLineColor(col[k]);
                g->Fit(fit, "QR");
                g->Draw("P SAME");
                out.cd();
                g->Write();
                c.cd(j * kNRef + r + 1);
                const double p0 = fit->GetParameter(0), mu = fit->GetParameter(1), sg = fit->GetParameter(2);
                const double x95 = mu + sg * std::sqrt(2.0) * TMath::ErfInverse(2 * 0.95 - 1);
                const double x99 = mu + sg * std::sqrt(2.0) * TMath::ErfInverse(2 * 0.99 - 1);
                lg->AddEntry(g, Form("EM%d: 50%% %.2f, 99%% %.2f GeV", k, mu, x99), "p");
                // the 99% point, a candidate cut
                TLine *l99 = new TLine(x99, 0, x99, 1.1);
                l99->SetLineColor(col[k]); l99->SetLineStyle(2);
                if (x99 < kPtMax) l99->Draw();
                for (int o = 0; o < 2; o++)
                    std::fprintf(outs[o], "  EM%d, %-18s %-20s plateau %.3f  mu %.2f +- %.2f  sigma %.2f  1.15 mu %.2f  95%% at %.2f  99%% at %.2f GeV  (chi2/ndf %.1f/%d)\n",
                                 k, kRefName[r], kJetName[j], p0, mu, fit->GetParError(1), sg, 1.15 * mu, x95, x99, fit->GetChisquare(), fit->GetNDF());
            }
            lg->Draw();
        }
    c.SaveAs((base + ".png").c_str());
    std::fclose(txt);
    out.Write();
    out.Close();
    return 0;
}
