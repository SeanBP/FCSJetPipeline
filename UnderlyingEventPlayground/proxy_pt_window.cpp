// How wide a jet-pT window can a proxy event be allowed, and what does it cost?
//
// Input: a folder of sim_to_jet -E 1 -s 1 output, i.e. SimpleTree_<s>.root
// (calorimeter points) next to jet_output_<s>.root (reco jets and the
// cross-section weights, jet_scale/JetXSecWeights.h).
//
// Two things are measured, both with the cross-section weights so they
// describe the physical (steeply falling) spectrum and not the flat-ish
// sampling in ptHat bins:
//   spectrum   the reco-jet pT spectrum. For a jet of pT p and a window
//              around p, the fraction of all reco jets inside the window is
//              how often a candidate proxy jet passes the pT requirement, and
//              the mean pT inside the window minus p is how far the proxy
//              jets sit, on average, from the jet they stand in for.
//   UE vs pT   the ECal point energy away from the jet (at least proxy_min_dr
//              from its axis and within 90 +- 45 deg of it in azimuth),
//              scaled from the eta-phi area the detector offers there to the
//              area of one R cone. How fast this changes with jet pT is how
//              much a pT mismatch between jet and proxy jet matters.
// Multiplying the UE slope by the mean pT shift of a window gives the bias in
// the UE estimate that the window causes.
//
// Build and run from this folder:
//   g++ -m64 -O2 proxy_pt_window.cpp ../shared/VectorDict.cxx -o proxy_pt_window -I. -I../shared -I../jet_scale -I../underlying_event \
//       `root-config --cflags --libs` -std=c++11
//   ./proxy_pt_window <input_dir/> [output_base]
//
// Writes <output_base>.pdf/.png/.root (default proxy_pt_window).

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"
#include "ProxySelection.h"

static const double kBandHalfWidthDeg = 45.0;  // azimuthal band about 90 deg where the UE is sampled
static const double kMinArea = 0.3;            // least eta-phi area the detector must offer there
static const double kCellCm = 4.0;             // cells used to measure that area

static const int kNSpec = 300;                 // spectrum bins of 0.1 GeV from 0 to 30 GeV
static const double kSpecMax = 30.0;

static const int kNPt = 8;                     // jet pT bins for the UE measurement
static const double kPt[kNPt + 1] = {2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 12.0};

struct Window { const char *name; double rel, abs; }; // half-width: rel * pT, or abs GeV
static const int kNWin = 6;
static const Window kWin[kNWin] = {{"#pm 5%", 0.05, 0}, {"#pm 10%", 0.10, 0}, {"#pm 20%", 0.20, 0},
                                   {"#pm 0.25 GeV", 0, 0.25}, {"#pm 0.5 GeV", 0, 0.5}, {"#pm 1 GeV", 0, 1.0}};
static const char *kWinText[kNWin] = {"+-5%", "+-10%", "+-20%", "+-0.25 GeV", "+-0.5 GeV", "+-1 GeV"};

struct Cell { double eta, phi, area; };

// Weighted content of the spectrum between lo and hi, and the pT-weighted sum
static void Integrate(const TH1D &h, double lo, double hi, double &sum, double &sum_pt)
{
    sum = 0.0; sum_pt = 0.0;
    for (int b = 1; b <= h.GetNbinsX(); b++)
    {
        const double b_lo = h.GetBinLowEdge(b), b_hi = b_lo + h.GetBinWidth(b);
        const double a = std::max(lo, b_lo), z = std::min(hi, b_hi);
        if (z <= a) continue;
        const double part = h.GetBinContent(b) * (z - a) / (b_hi - b_lo);
        sum += part; sum_pt += part * 0.5 * (a + z);
    }
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2)
    {
        std::cerr << "usage: ./proxy_pt_window <input_dir/> [output_base]" << std::endl;
        return 1;
    }
    std::string input_dir = argv[1];
    if (input_dir[input_dir.size() - 1] != '/') input_dir += '/';
    const std::string out_base = argc > 2 ? argv[2] : "proxy_pt_window";

    std::vector<std::string> jet_files, simple_files;
    {
        TSystemDirectory dir("input_dir", input_dir.c_str());
        TList *files = dir.GetListOfFiles();
        if (!files) { std::cerr << "ERROR: cannot list " << input_dir << std::endl; return 1; }
        TIter next(files);
        TSystemFile *file;
        std::vector<std::string> names;
        while ((file = (TSystemFile *)next()))
        {
            const std::string n = file->GetName();
            if (n.compare(0, 10, "jet_output") == 0 && n.size() > 5 && n.compare(n.size() - 5, 5, ".root") == 0)
                names.push_back(n);
        }
        std::sort(names.begin(), names.end());
        for (size_t i = 0; i < names.size(); i++)
        {
            jet_files.push_back(input_dir + names[i]);
            simple_files.push_back(input_dir + "SimpleTree" + names[i].substr(10));
        }
    }
    std::cout << "Found " << jet_files.size() << " file pairs in " << input_dir << std::endl;

    std::vector<float> file_weight;
    if (!ComputeXSecWeights(jet_files, "jetTree", file_weight))
        return 1;

    // The whole ECal face as cells, to measure the area available away from a jet
    std::vector<Cell> cells;
    {
        const FiducialRect &rect = kFiducialRectEcal;
        for (int side = -1; side <= 1; side += 2)
            for (double x = rect.x_inner; x < rect.x_outer; x += kCellCm)
                for (double y = rect.y_min; y < rect.y_max; y += kCellCm)
                {
                    const double dx = std::min(kCellCm, rect.x_outer - x), dy = std::min(kCellCm, rect.y_max - y);
                    const double cx = side * (x + dx / 2), cy = y + dy / 2, r2 = cx * cx + cy * cy;
                    Cell c;
                    computeEtaPhi(cx, cy, c.eta, c.phi);
                    c.area = dx * dy * z_proj / (r2 * std::sqrt(r2 + z_proj * z_proj));
                    cells.push_back(c);
                }
    }

    TH1D h_spec("spectrum", ";Reco jet p_{T} (GeV);Jets per 0.1 GeV (weighted, arbitrary units)", kNSpec, 0.0, kSpecMax);
    TH1D h_raw("spectrum_raw", ";Reco jet p_{T} (GeV);Jets per 0.1 GeV (unweighted)", kNSpec, 0.0, kSpecMax);
    h_spec.Sumw2();

    double ue_w[kNPt] = {0}, ue_wx[kNPt] = {0}, ue_wxx[kNPt] = {0}, ue_w2[kNPt] = {0}, ue_wpt[kNPt] = {0}, ue_wzero[kNPt] = {0};
    long ue_n[kNPt] = {0};

    const int MAX_PTS = 200; // matches StSimpleReaderMaker.h's Cal_pt_* array size
    Long64_t n_events = 0, n_jets = 0;

    for (size_t f = 0; f < jet_files.size(); f++)
    {
        if (file_weight[f] < 0) continue;
        const double w = file_weight[f];

        TFile *fj = TFile::Open(jet_files[f].c_str());
        TFile *fs = TFile::Open(simple_files[f].c_str());
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        if (!jt || !st || jt->GetEntries() != st->GetEntries() || !st->GetBranch("Cal_npts"))
        {
            std::cerr << "WARNING: skipping " << jet_files[f] << std::endl;
            if (fj) { fj->Close(); delete fj; }
            if (fs) { fs->Close(); delete fs; }
            continue;
        }

        std::vector<float> *reco_pt = 0, *reco_eta = 0, *reco_phi = 0;
        jt->SetBranchStatus("*", 0);
        jt->SetBranchStatus("reco_pt", 1); jt->SetBranchStatus("reco_eta", 1); jt->SetBranchStatus("reco_phi", 1);
        jt->SetBranchAddress("reco_pt", &reco_pt);
        jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi);

        Int_t Cal_npts = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        Int_t Cal_pt_detid[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pt_branches[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pt_branches[i], 1);
        st->SetBranchAddress("Cal_npts", &Cal_npts);
        st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x);
        st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z);
        st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);

        const Long64_t n = jt->GetEntries();
        for (Long64_t ev = 0; ev < n; ev++)
        {
            if (jt->GetEntry(ev) <= 0 || st->GetEntry(ev) <= 0) continue;
            n_events++;

            for (size_t j = 0; j < reco_pt->size(); j++)
            {
                const double pt = reco_pt->at(j), eta = reco_eta->at(j), phi = reco_phi->at(j);
                n_jets++;
                h_spec.Fill(pt, w);
                h_raw.Fill(pt);

                int pb = -1;
                for (int k = 0; k < kNPt; k++) if (pt >= kPt[k] && pt < kPt[k + 1]) pb = k;
                if (pb < 0) continue;

                // No other reco jet may reach into the sampled region
                bool other_jet = false;
                for (size_t o = 0; o < reco_pt->size(); o++)
                    if (o != j) other_jet = true;
                if (other_jet) continue;

                double area = 0.0;
                for (size_t c = 0; c < cells.size(); c++)
                    if (ProxyJetAllowed(eta, phi, kBandHalfWidthDeg, cells[c].eta, cells[c].phi))
                        area += cells[c].area;
                if (area < kMinArea) continue;

                // Points: same selection and direction convention as the reco jets (use_points)
                double e_sum = 0.0;
                for (int i = 0; i < Cal_npts; i++)
                {
                    if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                    if (!pass_point_energy_cut(Cal_pt_energy[i], kSqrtS_MainRun22_GeV)) continue;
                    const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                    const double p_eta = std::asinh(Cal_pt_z[i] / rho), p_phi = std::atan2(Cal_pt_y[i], Cal_pt_x[i]);
                    if (ProxyJetAllowed(eta, phi, kBandHalfWidthDeg, p_eta, p_phi))
                        e_sum += Cal_pt_energy[i];
                }
                const double ue = e_sum * M_PI * R * R / area; // scaled to one cone
                ue_w[pb] += w; ue_wx[pb] += w * ue; ue_wxx[pb] += w * ue * ue; ue_w2[pb] += w * w; ue_wpt[pb] += w * pt; ue_n[pb]++;
                if (e_sum <= 0.0) ue_wzero[pb] += w;
            }
        }
        fj->Close(); delete fj;
        fs->Close(); delete fs;
    }
    std::cout << "Read " << n_events << " events, " << n_jets << " reco jets" << std::endl;

    // ---------------- UE against jet pT ----------------
    TGraphErrors g_ue;
    std::cout << "\nUE energy scaled to one R = " << R << " cone, against reco jet pT (single-jet events):" << std::endl;
    for (int k = 0; k < kNPt; k++)
    {
        if (ue_w[k] <= 0) continue;
        const double mean = ue_wx[k] / ue_w[k], var = std::max(0.0, ue_wxx[k] / ue_w[k] - mean * mean);
        const double n_eff = ue_w[k] * ue_w[k] / ue_w2[k], err = std::sqrt(var / n_eff), pt_mean = ue_wpt[k] / ue_w[k];
        const int n = g_ue.GetN();
        g_ue.SetPoint(n, pt_mean, mean);
        g_ue.SetPointError(n, 0.0, err);
        printf("  pT %4.1f-%4.1f GeV: UE %.3f +- %.3f GeV  (spread %.2f, zero in %.0f%%, %ld jets, %.0f effective)\n",
               kPt[k], kPt[k + 1], mean, err, std::sqrt(var), 100.0 * ue_wzero[k] / ue_w[k], ue_n[k], n_eff);
    }
    TF1 line("line", "pol1", kPt[0], kPt[kNPt]);
    g_ue.Fit(&line, "Q");
    const double slope = line.GetParameter(1), slope_err = line.GetParError(1);
    printf("  linear fit: UE = %.3f + (%.4f +- %.4f) * pT, chi2/ndf %.1f/%d\n", line.GetParameter(0), slope, slope_err,
           line.GetChisquare(), line.GetNDF());

    // ---------------- windows ----------------
    double all = 0.0, dummy = 0.0;
    Integrate(h_spec, reco_jet_pt_min_em, kSpecMax, all, dummy);
    const double centres[] = {2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
    std::cout << "\nWindows about a jet of pT p. Each entry: fraction of all reco jets inside the window; mean pT inside minus p; "
                 "resulting UE bias from the fitted slope" << std::endl;
    printf("  %-11s", "p (GeV)");
    for (int i = 0; i < kNWin; i++) printf(" | %-26s", kWinText[i]);
    printf("\n");
    for (size_t c = 0; c < sizeof(centres) / sizeof(centres[0]); c++)
    {
        const double p = centres[c];
        printf("  %-11.1f", p);
        for (int i = 0; i < kNWin; i++)
        {
            const double hw = kWin[i].rel > 0 ? kWin[i].rel * p : kWin[i].abs;
            double s, spt;
            Integrate(h_spec, std::max(p - hw, (double)reco_jet_pt_min_em), p + hw, s, spt);
            if (s <= 0) { printf(" | %-26s", "-"); continue; }
            const double shift = spt / s - p;
            printf(" | %7.3f%% %+6.2f GeV %+6.3f", 100.0 * s / all, shift, slope * shift);
        }
        printf("\n");
    }

    // ---------------- plot ----------------
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TCanvas c("c", "", 2400, 800);
    c.Divide(3, 1);
    TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.045);

    c.cd(1);
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13); gPad->SetLogy();
    TH1D *spec = (TH1D *)h_spec.Clone("spec_draw");
    spec->Rebin(5);
    spec->Scale(1.0 / spec->GetMaximum());
    spec->SetTitle(";Reco jet p_{T} (GeV);Jets per 0.5 GeV (peak = 1)");
    spec->GetXaxis()->SetRangeUser(reco_jet_pt_min_em, 16.0);
    spec->GetXaxis()->SetTitleSize(0.050); spec->GetYaxis()->SetTitleSize(0.050);
    spec->GetXaxis()->SetLabelSize(0.045); spec->GetYaxis()->SetLabelSize(0.045); spec->GetYaxis()->SetTitleOffset(1.55);
    spec->SetLineColor(kAzure + 1); spec->SetLineWidth(2); spec->SetMarkerColor(kAzure + 1); spec->SetMarkerStyle(20);
    spec->Draw("E1");
    t.DrawLatex(0.16, 0.92, "Reco jet spectrum, weighted");

    c.cd(2);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13); gPad->SetLogy();
    const int cols[kNWin] = {kAzure + 1, kOrange + 7, kGreen + 2, kAzure + 1, kOrange + 7, kGreen + 2};
    TH1D frame("frame_win", ";Jet p_{T} (GeV);Reco jets searched per jet in the window", 10, 2.0, 12.0);
    frame.GetYaxis()->SetRangeUser(1.0, 3e5);
    frame.GetXaxis()->SetTitleSize(0.050); frame.GetYaxis()->SetTitleSize(0.044);
    frame.GetXaxis()->SetLabelSize(0.045); frame.GetYaxis()->SetLabelSize(0.045); frame.GetYaxis()->SetTitleOffset(1.3);
    frame.Draw("AXIS");
    TLegend leg(0.17, 0.62, 0.62, 0.89);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.038); leg.SetNColumns(2);
    for (int i = 0; i < kNWin; i++)
    {
        TGraph *g = new TGraph();
        for (double p = 2.5; p <= 12.0; p += 0.25)
        {
            const double hw = kWin[i].rel > 0 ? kWin[i].rel * p : kWin[i].abs;
            double s, spt;
            Integrate(h_spec, std::max(p - hw, (double)reco_jet_pt_min_em), p + hw, s, spt);
            if (s > 0) g->SetPoint(g->GetN(), p, all / s);
        }
        g->SetLineColor(cols[i]); g->SetLineWidth(3); g->SetLineStyle(kWin[i].rel > 0 ? 1 : 2);
        g->Draw("L SAME");
        leg.AddEntry(g, kWin[i].name, "l");
    }
    leg.Draw();
    t.DrawLatex(0.14, 0.92, "Cost of the p_{T} window");

    c.cd(3);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13);
    TH1D frame_ue("frame_ue", ";Reco jet p_{T} (GeV);UE energy per cone (GeV)", 10, 2.0, 12.0);
    double ue_max = 0.0;
    for (int i = 0; i < g_ue.GetN(); i++) ue_max = std::max(ue_max, g_ue.GetY()[i] + g_ue.GetEY()[i]);
    frame_ue.GetYaxis()->SetRangeUser(0.0, 1.4 * ue_max);
    frame_ue.GetXaxis()->SetTitleSize(0.050); frame_ue.GetYaxis()->SetTitleSize(0.050);
    frame_ue.GetXaxis()->SetLabelSize(0.045); frame_ue.GetYaxis()->SetLabelSize(0.045); frame_ue.GetYaxis()->SetTitleOffset(1.3);
    frame_ue.Draw("AXIS");
    g_ue.SetMarkerStyle(20); g_ue.SetMarkerSize(1.4); g_ue.SetLineWidth(2);
    line.SetLineColor(kRed + 1); line.SetLineWidth(2);
    g_ue.Draw("P SAME");
    t.DrawLatex(0.14, 0.92, "UE energy away from the jet");
    t.SetTextFont(42); t.SetTextSize(0.040);
    t.DrawLatex(0.18, 0.84, Form("slope %.3f #pm %.3f GeV per GeV of jet p_{T}", slope, slope_err));

    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());

    TFile fout((out_base + ".root").c_str(), "RECREATE");
    h_spec.Write(); h_raw.Write(); g_ue.Write("ue_vs_pt");
    fout.Close();
    return 0;
}
