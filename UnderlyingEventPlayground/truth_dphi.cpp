// Azimuthal correlation between truth jets in the FCS pseudorapidity range:
// is a second truth jet more likely to be back-to-back with the first?
//
// Input: a folder of sim_to_jet -E 1 -s 1 output, i.e. SimpleTree_<s>.root
// (truth particles) next to jet_output_<s>.root (used only for the
// cross-section weights, jet_scale/JetXSecWeights.h).
//
// Truth jets are rebuilt here from the SimpleTree so that no fiducial cut is
// applied: the same inputs and algorithm as sim_to_jet's EM-only truth jets
// (photons created within 100 um of the primary vertex, anti-kT R = 0.7,
// pT >= truth_jet_pt_min_em), kept if eta_min < eta < eta_max at any azimuth.
//
// For every ordered pair of truth jets in an event, |dphi| is histogrammed
// with the event weight and divided by the weighted number of truth jets,
// giving the mean number of other truth jets per truth jet in each |dphi|
// bin. The acceptance is the full azimuth, so with no correlation this is
// flat; a rise toward 180 deg is recoil.
//
// Build and run from this folder:
//   g++ -m64 truth_dphi.cpp ../shared/VectorDict.cxx -o truth_dphi -I. -I../shared -I../jet_scale \
//       `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//   ./truth_dphi <input_dir/> [output_base] [eta_min] [eta_max]
//   ./truth_dphi replot <output_base> [eta_min] [eta_max]     (redraw from the saved histograms)
//
// Writes <output_base>.pdf/.png/.root (default truth_dphi) and the cone-probability
// plots <output_base>_cone_prob.pdf/.png (cone at the jet's eta) and
// <output_base>_band_prob.pdf/.png (any eta in the range), and the increase in
// probability of the opposite band over the sideways band, <output_base>_band_increase.pdf/.png
// (in pT) and <output_base>_band_increase_energy.pdf/.png (in energy). Default eta range 2.5 to 4.0.

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "fastjet/ClusterSequence.hh"

#include "JetParameters.h"
#include "JetXSecWeights.h"

static const int kNDphi = 36;                    // 5 deg bins
static const double kDecayDisplacementMax = 0.01; // cm, as in JetMatcher.cpp
static const int kGeantPhoton = 1;

struct Sel
{
    const char *name;
    double pt_min;                               // both jets of a pair
    double sw, sw2;                              // weighted truth jets
    long n_jets;
    std::vector<double> pair, pair_w2;
    std::vector<long> pair_n;
    double sw_events_ge2;                        // weighted events with >= 2 such jets
    Sel(const char *n, double p) : name(n), pt_min(p), sw(0), sw2(0), n_jets(0),
        pair(kNDphi, 0.0), pair_w2(kNDphi, 0.0), pair_n(kNDphi, 0), sw_events_ge2(0) {}
};

// Bin widths used when drawing, in degrees, one per pT threshold: wider for
// the harder, sparser selections
static const int kNSel = 3;
static const double kPtMin[kNSel] = {1.0, 2.0, 3.0};
static const int kDrawWidth[kNSel] = {5, 15, 30};

// Draws the saved 5 deg histograms (content = other truth jets per truth jet
// in that 5 deg bin). Each is merged into wider bins and divided by the number
// of 5 deg bins merged, so every curve stays "per 5 deg of |dphi|" and the
// curves are directly comparable whatever their bin width.
static void Plot(const std::vector<TH1D *> &h5, const std::string &out_base, double eta_min, double eta_max)
{
    gStyle->SetOptStat(0);
    TCanvas c("c", "", 1200, 900);
    c.SetLeftMargin(0.13); c.SetRightMargin(0.04); c.SetBottomMargin(0.13); c.SetLogy();

    const int cols[kNSel] = {kAzure + 1, kOrange + 7, kGreen + 2};
    std::vector<TH1D *> h;
    double ymin = 1e9, ymax = 0;
    for (int s = 0; s < kNSel && s < (int)h5.size(); s++)
    {
        const int merge = kDrawWidth[s] / 5, nb = 180 / kDrawWidth[s];
        TH1D *hh = new TH1D(Form("draw_%d", s), "", nb, 0, 180);
        for (int b = 0; b < nb; b++)
        {
            double sum = 0, err2 = 0;
            for (int k = 0; k < merge; k++)
            {
                sum += h5[s]->GetBinContent(b * merge + k + 1);
                err2 += std::pow(h5[s]->GetBinError(b * merge + k + 1), 2);
            }
            hh->SetBinContent(b + 1, sum / merge);
            hh->SetBinError(b + 1, std::sqrt(err2) / merge);
            if (sum > 0) ymin = std::min(ymin, sum / merge);
            ymax = std::max(ymax, sum / merge);
        }
        hh->SetLineColor(cols[s]); hh->SetMarkerColor(cols[s]); hh->SetMarkerStyle(20); hh->SetLineWidth(2);
        h.push_back(hh);
    }
    h[0]->SetTitle(";|#Delta#phi| between two truth jets (deg);Other truth jets per truth jet, per 5#circ of |#Delta#phi|");
    h[0]->GetYaxis()->SetRangeUser(0.3 * ymin, 60.0 * ymax);
    h[0]->GetXaxis()->SetTitleSize(0.050); h[0]->GetXaxis()->SetLabelSize(0.045);
    h[0]->GetYaxis()->SetTitleSize(0.043); h[0]->GetYaxis()->SetLabelSize(0.045); h[0]->GetYaxis()->SetTitleOffset(1.4);
    h[0]->Draw("E1");
    for (size_t s = 1; s < h.size(); s++) h[s]->Draw("E1 SAME");

    TLegend leg(0.17, 0.72, 0.72, 0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.038);
    for (size_t s = 0; s < h.size(); s++)
        leg.AddEntry(h[s], Form("both jets p_{T} #geq %.0f GeV  (%d#circ bins)", kPtMin[s], kDrawWidth[s]), "lp");
    leg.Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.040); t.SetTextFont(62);
    t.DrawLatex(0.13, 0.94, Form("Truth-jet pairs with %.1f < #eta < %.1f, simulation", eta_min, eta_max));

    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());
}

// Bins of trigger truth-jet pT and thresholds on the jet seen, for the cone-probability plot
static const int kNPB = 5, kNSeen = 3;
static const double kPB[kNPB + 1] = {2.0, 3.0, 4.0, 5.0, 7.0, 12.0};
static const double kSeen[kNSeen] = {1.0, 2.0, 3.0};
// The same in energy: bins of trigger truth-jet energy and thresholds on the energy of the jet seen
static const double kEB[kNPB + 1] = {20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
static const double kSeenE[kNSeen] = {10.0, 20.0, 30.0};
// Truth jets are clustered down to this pT so the energy thresholds are complete
// over the eta range (E = 10 GeV at eta = 4 is pT = 0.37 GeV); every pT-based
// count applies its own threshold of at least truth_jet_pt_min_em
static const double kClusterPtMin = 0.3;

// Draws the probability (%) of another truth-jet axis inside an R cone opposite
// and sideways to a trigger truth jet, against the trigger pT. hc holds, per
// seen-jet threshold, the opposite then the sideways histogram. With band set
// the histograms are instead for the band |dphi| < R about those azimuths, at
// any eta in the range.
static void PlotCone(const std::vector<TH1D *> &hc, const std::string &out_base, double eta_min, double eta_max, bool band)
{
    gStyle->SetOptStat(0);
    TCanvas c(band ? "cb" : "cc", "", 1600, 800);
    c.Divide(2, 1);
    for (int q = 0; q < 2; q++) // the two lowest seen-jet thresholds
    {
        c.cd(q + 1);
        gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.09);
        TH1D *ho = hc[2 * q], *hs = hc[2 * q + 1];
        ho->SetLineColor(kRed + 1); ho->SetMarkerColor(kRed + 1); ho->SetMarkerStyle(20); ho->SetLineWidth(2); ho->SetMarkerSize(1.4);
        hs->SetLineColor(kAzure + 1); hs->SetMarkerColor(kAzure + 1); hs->SetMarkerStyle(24); hs->SetLineWidth(2); hs->SetMarkerSize(1.4);
        double ymax = 0;
        for (int k = 1; k <= kNPB; k++) ymax = std::max(ymax, ho->GetBinContent(k) + ho->GetBinError(k));
        ho->SetTitle(band ? ";Truth jet p_{T} (GeV);Probability of another truth jet in the band (%)"
                          : ";Truth jet p_{T} (GeV);Probability of another truth jet in the cone (%)");
        ho->GetYaxis()->SetRangeUser(0.0, 1.45 * ymax);
        ho->GetXaxis()->SetTitleSize(0.050); ho->GetXaxis()->SetLabelSize(0.045);
        ho->GetYaxis()->SetTitleSize(0.046); ho->GetYaxis()->SetLabelSize(0.045); ho->GetYaxis()->SetTitleOffset(1.35);
        ho->Draw("E1");
        hs->Draw("E1 SAME");
        TLegend *leg = new TLegend(0.42, 0.75, 0.95, 0.89);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
        leg->AddEntry(ho, band ? "opposite band (180#circ)" : "opposite cone (180#circ)", "lp");
        leg->AddEntry(hs, band ? "sideways band (#pm 90#circ)" : "sideways cone (#pm 90#circ)", "lp");
        leg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.042); t.SetTextFont(62);
        if (band) t.DrawLatex(0.14, 0.935, Form("Other truth jet p_{T} #geq %.0f GeV, any #eta, within #pm%.0f#circ", kSeen[q], R * 180.0 / M_PI));
        else t.DrawLatex(0.14, 0.935, Form("Other truth jet p_{T} #geq %.0f GeV, R = %.1f cone", kSeen[q], R));
        t.SetTextFont(42); t.SetTextSize(0.036);
        t.DrawLatex(0.18, 0.84, Form("%.1f < #eta < %.1f", eta_min, eta_max));
        t.DrawLatex(0.18, 0.79, "simulation");
    }
    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());
}

// Draws the increase in probability of the opposite band over the sideways
// band, opposite - sideways in percentage points, against the trigger pT, one
// series per seen-jet threshold. hb is as for PlotCone. With energy set the
// trigger bins and seen-jet thresholds are in energy (kEB, kSeenE), not pT.
static void PlotIncrease(const std::vector<TH1D *> &hb, const std::string &out_base, double eta_min, double eta_max, bool energy)
{
    const double *kPB = energy ? kEB : ::kPB, *kSeen = energy ? kSeenE : ::kSeen;
    const char *var = energy ? "E" : "pT";
    gStyle->SetOptStat(0);
    TCanvas c(energy ? "cie" : "ci", "", 1200, 900);
    c.SetLeftMargin(0.13); c.SetRightMargin(0.04); c.SetBottomMargin(0.13); c.SetTopMargin(0.08);

    const int cols[kNSeen] = {kAzure + 1, kOrange + 7, kGreen + 2};
    const int marks[kNSeen] = {20, 21, 22};
    std::vector<TGraphErrors *> g;
    double ymax = 0, ymin = 0;
    for (int q = 0; q < kNSeen; q++)
    {
        TGraphErrors *gr = new TGraphErrors();
        for (int k = 1; k <= kNPB; k++)
        {
            const double o = hb[2 * q]->GetBinContent(k), eo = hb[2 * q]->GetBinError(k);
            const double sd = hb[2 * q + 1]->GetBinContent(k), es = hb[2 * q + 1]->GetBinError(k);
            const double d = o - sd, ed = std::sqrt(eo * eo + es * es);
            const double width = kPB[k] - kPB[k - 1];
            const int n = gr->GetN();
            gr->SetPoint(n, 0.5 * (kPB[k - 1] + kPB[k]) + (q - 1) * 0.08 * width, d); // offset so the series do not overlap
            gr->SetPointError(n, 0.0, ed);
            std::cout << "  increase, other jet " << var << " >= " << kSeen[q] << " GeV, trigger " << var << " " << kPB[k - 1] << "-" << kPB[k] << " GeV: "
                      << d << " +- " << ed << " percentage points  (opposite " << o << " %, sideways " << sd << " %)" << std::endl;
            ymax = std::max(ymax, d + ed); ymin = std::min(ymin, d - ed);
        }
        gr->SetLineColor(cols[q]); gr->SetMarkerColor(cols[q]); gr->SetMarkerStyle(marks[q]); gr->SetMarkerSize(1.6); gr->SetLineWidth(2);
        g.push_back(gr);
    }
    TH1D frame(energy ? "frame_inc_e" : "frame_inc", energy ? ";Truth jet energy (GeV);Opposite #minus sideways band probability (percentage points)"
                                                              : ";Truth jet p_{T} (GeV);Opposite #minus sideways band probability (percentage points)", kNPB, kPB);
    const double ylo = ymin - 0.1 * ymax, yhi = 1.45 * ymax;
    frame.GetYaxis()->SetRangeUser(ylo, yhi);
    frame.GetXaxis()->SetTitleSize(0.050); frame.GetXaxis()->SetLabelSize(0.045);
    frame.GetYaxis()->SetTitleSize(0.038); frame.GetYaxis()->SetLabelSize(0.045); frame.GetYaxis()->SetTitleOffset(1.55);
    frame.Draw("AXIS");
    TLine zero(kPB[0], 0.0, kPB[kNPB], 0.0); zero.SetLineStyle(2); zero.Draw();
    for (int k = 1; k < kNPB; k++)
    {
        TLine *edge = new TLine(kPB[k], ylo, kPB[k], yhi);
        edge->SetLineStyle(3); edge->SetLineColor(kGray + 1); edge->Draw();
    }
    for (size_t q = 0; q < g.size(); q++) if (g[q]->GetN() > 0) g[q]->Draw("P SAME");

    TLegend leg(0.52, 0.70, 0.95, 0.89);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.038);
    for (size_t q = 0; q < g.size(); q++) leg.AddEntry(g[q], Form(energy ? "other truth jet E #geq %.0f GeV" : "other truth jet p_{T} #geq %.0f GeV", kSeen[q]), "lp");
    leg.Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.038); t.SetTextFont(62);
    t.DrawLatex(0.13, 0.94, Form("Truth jets, %.1f < #eta < %.1f, bands within #pm%.0f#circ of 180#circ and #pm90#circ, simulation", eta_min, eta_max, R * 180.0 / M_PI));

    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2)
    {
        std::cerr << "usage: ./truth_dphi <input_dir/> [output_base] [eta_min] [eta_max]" << std::endl;
        return 1;
    }
    // "replot <output_base> [eta_min] [eta_max]": redraw from an existing <output_base>.root
    if (std::string(argv[1]) == "replot" && argc > 2)
    {
        TFile fin((std::string(argv[2]) + ".root").c_str());
        std::vector<TH1D *> h;
        for (int s = 0; s < kNSel; s++)
        {
            TH1D *hh = (TH1D *)fin.Get(Form("dphi_%d", s));
            if (!hh) { std::cerr << "ERROR: dphi_" << s << " not in " << argv[2] << ".root" << std::endl; return 1; }
            h.push_back(hh);
        }
        Plot(h, argv[2], argc > 3 ? std::atof(argv[3]) : 2.5, argc > 4 ? std::atof(argv[4]) : 4.0);
        std::vector<TH1D *> hc;
        for (int q = 0; q < kNSeen; q++)
            for (int place = 0; place < 2; place++)
                hc.push_back((TH1D *)fin.Get(Form("cone_%s_%d", place == 0 ? "opp" : "side", q)));
        if (hc[0])
            PlotCone(hc, std::string(argv[2]) + "_cone_prob", argc > 3 ? std::atof(argv[3]) : 2.5, argc > 4 ? std::atof(argv[4]) : 4.0, false);
        std::vector<TH1D *> hb;
        for (int q = 0; q < kNSeen; q++)
            for (int place = 0; place < 2; place++)
                hb.push_back((TH1D *)fin.Get(Form("band_%s_%d", place == 0 ? "opp" : "side", q)));
        if (hb[0])
        {
            PlotCone(hb, std::string(argv[2]) + "_band_prob", argc > 3 ? std::atof(argv[3]) : 2.5, argc > 4 ? std::atof(argv[4]) : 4.0, true);
            if (hb[2 * kNSeen - 1])
                PlotIncrease(hb, std::string(argv[2]) + "_band_increase", argc > 3 ? std::atof(argv[3]) : 2.5, argc > 4 ? std::atof(argv[4]) : 4.0, false);
            std::vector<TH1D *> he;
            for (int q = 0; q < kNSeen; q++)
                for (int place = 0; place < 2; place++)
                    he.push_back((TH1D *)fin.Get(Form("eband_%s_%d", place == 0 ? "opp" : "side", q)));
            if (he[0])
                PlotIncrease(he, std::string(argv[2]) + "_band_increase_energy", argc > 3 ? std::atof(argv[3]) : 2.5, argc > 4 ? std::atof(argv[4]) : 4.0, true);
        }
        return 0;
    }

    std::string input_dir = argv[1];
    if (input_dir[input_dir.size() - 1] != '/') input_dir += '/';
    const std::string out_base = argc > 2 ? argv[2] : "truth_dphi";
    const double eta_min = argc > 3 ? std::atof(argv[3]) : 2.5, eta_max = argc > 4 ? std::atof(argv[4]) : 4.0;

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

    std::vector<Sel> sels;
    sels.push_back(Sel("truth pT >= 1 GeV", 1.0));
    sels.push_back(Sel("truth pT >= 2 GeV", 2.0));
    sels.push_back(Sel("truth pT >= 3 GeV", 3.0));

    const int MAX = 10000;
    static Float_t mcpart_px[MAX], mcpart_py[MAX], mcpart_pz[MAX], mcpart_E[MAX];
    static Float_t mcpart_Vtx_x[MAX], mcpart_Vtx_y[MAX], mcpart_Vtx_z[MAX];
    static Int_t mcpart_geid[MAX], mcpart_idVtx[MAX];
    Int_t mcpart_num = 0;

    // ---- cone probabilities ----
    // For a "trigger" truth jet (pT above kTrig) at (eta_t, phi_t): is there
    // another truth jet (pT above kSeen, FCS eta range) whose axis lies
    // within R of the point opposite it, (eta_t, phi_t + 180 deg)? The same
    // is counted for the two sideways points (eta_t, phi_t +- 90 deg), and
    // for an unrelated event the chance of a truth-jet axis in a cone at
    // eta_t follows from the eta distribution of truth jets in all events.
    const int kNTrig = 3;
    const double kTrig[kNTrig] = {2.0, 3.0, 5.0};
    const int kNEta = 400; const double kEtaLo = 1.0, kEtaHi = 5.0; // 0.01 bins
    double c_w[kNTrig] = {0}, c_w2[kNTrig] = {0}; long c_n[kNTrig] = {0};
    double c_opp[kNTrig][kNSeen] = {{0}}, c_opp2[kNTrig][kNSeen] = {{0}}, c_side[kNTrig][kNSeen] = {{0}}, c_side2[kNTrig][kNSeen] = {{0}};
    long c_opp_n[kNTrig][kNSeen] = {{0}};
    std::vector<std::vector<double> > trig_eta(kNTrig, std::vector<double>(kNEta, 0.0)), seen_eta(kNSeen, std::vector<double>(kNEta, 0.0));

    // The same opposite and sideways cone counts in bins of trigger pT, for the plot
    double pb_w[kNPB] = {0}, pb_opp[kNPB][kNSeen] = {{0}}, pb_opp2[kNPB][kNSeen] = {{0}}, pb_side[kNPB][kNSeen] = {{0}}, pb_side2[kNPB][kNSeen] = {{0}};
    // and for the band |dphi| < R about the opposite and sideways azimuths, at any eta in the range
    double bd_opp[kNPB][kNSeen] = {{0}}, bd_opp2[kNPB][kNSeen] = {{0}}, bd_side[kNPB][kNSeen] = {{0}}, bd_side2[kNPB][kNSeen] = {{0}};
    // the band counts again, in bins of trigger energy with thresholds on the energy of the jet seen
    double eb_w[kNPB] = {0}, eb_opp[kNPB][kNSeen] = {{0}}, eb_opp2[kNPB][kNSeen] = {{0}}, eb_side[kNPB][kNSeen] = {{0}}, eb_side2[kNPB][kNSeen] = {{0}};

    const fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    Long64_t n_events = 0;
    double sw_events = 0;

    for (size_t f = 0; f < simple_files.size(); f++)
    {
        if (file_weight[f] < 0) continue;
        TFile *fs = TFile::Open(simple_files[f].c_str());
        if (!fs || fs->IsZombie()) { std::cerr << "WARNING: cannot open " << simple_files[f] << std::endl; continue; }
        TTree *tree = (TTree *)fs->Get("data");
        if (!tree) { fs->Close(); delete fs; continue; }

        tree->SetBranchStatus("*", 0);
        const char *br[] = {"mcpart_num", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_geid",
                            "mcpart_idVtx", "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
        for (int i = 0; i < 10; i++) tree->SetBranchStatus(br[i], 1);
        tree->SetBranchAddress("mcpart_num", &mcpart_num);
        tree->SetBranchAddress("mcpart_px", mcpart_px);
        tree->SetBranchAddress("mcpart_py", mcpart_py);
        tree->SetBranchAddress("mcpart_pz", mcpart_pz);
        tree->SetBranchAddress("mcpart_E", mcpart_E);
        tree->SetBranchAddress("mcpart_geid", mcpart_geid);
        tree->SetBranchAddress("mcpart_idVtx", mcpart_idVtx);
        tree->SetBranchAddress("mcpart_Vtx_x", mcpart_Vtx_x);
        tree->SetBranchAddress("mcpart_Vtx_y", mcpart_Vtx_y);
        tree->SetBranchAddress("mcpart_Vtx_z", mcpart_Vtx_z);

        const double w = file_weight[f];
        const Long64_t n = tree->GetEntries();
        for (Long64_t ev = 0; ev < n; ev++)
        {
            if (tree->GetEntry(ev) <= 0) continue;
            n_events++; sw_events += w;

            // Primary vertex: creation point of the first idVtx == 1 particle
            bool found = false; double vx = 0, vy = 0, vz = 0;
            for (int i = 0; i < mcpart_num; i++)
                if (mcpart_idVtx[i] == 1) { vx = mcpart_Vtx_x[i]; vy = mcpart_Vtx_y[i]; vz = mcpart_Vtx_z[i]; found = true; break; }
            if (!found) continue;

            std::vector<fastjet::PseudoJet> parts;
            for (int i = 0; i < mcpart_num; i++)
            {
                if (mcpart_geid[i] != kGeantPhoton) continue;
                const double dx = mcpart_Vtx_x[i] - vx, dy = mcpart_Vtx_y[i] - vy, dz = mcpart_Vtx_z[i] - vz;
                if (std::sqrt(dx * dx + dy * dy + dz * dz) > kDecayDisplacementMax) continue;
                parts.push_back(fastjet::PseudoJet(mcpart_px[i], mcpart_py[i], mcpart_pz[i], mcpart_E[i]));
            }
            if (parts.empty()) continue;

            fastjet::ClusterSequence cs(parts, jet_def);
            const std::vector<fastjet::PseudoJet> all = cs.inclusive_jets(kClusterPtMin);
            std::vector<double> phi, pt, etas;
            std::vector<double> e_phi, e_E; // all jets down to kClusterPtMin, for the energy counts
            for (size_t j = 0; j < all.size(); j++)
            {
                if (all[j].pz() <= 0) continue;
                const double eta = all[j].eta();
                if (eta <= eta_min || eta >= eta_max) continue;
                e_phi.push_back(all[j].phi()); e_E.push_back(all[j].E());
                if (all[j].pt() < truth_jet_pt_min_em) continue;
                phi.push_back(all[j].phi()); pt.push_back(all[j].pt()); etas.push_back(eta);
            }

            // band probabilities in energy
            for (size_t a = 0; a < e_phi.size(); a++)
            {
                int eb = -1;
                for (int k = 0; k < kNPB; k++) if (e_E[a] >= kEB[k] && e_E[a] < kEB[k + 1]) eb = k;
                if (eb < 0) continue;
                eb_w[eb] += w;
                for (size_t b = 0; b < e_phi.size(); b++)
                {
                    if (b == a) continue;
                    for (int place = 0; place < 3; place++)
                    {
                        const double target = e_phi[a] + (place == 0 ? M_PI : (place == 1 ? 0.5 * M_PI : -0.5 * M_PI));
                        if (std::fabs(std::remainder(e_phi[b] - target, 2.0 * M_PI)) >= R) continue;
                        for (int q = 0; q < kNSeen; q++)
                        {
                            if (e_E[b] < kSeenE[q]) continue;
                            if (place == 0) { eb_opp[eb][q] += w; eb_opp2[eb][q] += w * w; }
                            else { eb_side[eb][q] += 0.5 * w; eb_side2[eb][q] += 0.25 * w * w; }
                        }
                    }
                }
            }

            // cone probabilities
            for (size_t a = 0; a < phi.size(); a++)
            {
                const int eb = (int)((etas[a] - kEtaLo) / (kEtaHi - kEtaLo) * kNEta);
                for (int q = 0; q < kNSeen; q++) if (pt[a] >= kSeen[q]) seen_eta[q][eb] += w;
                int pb = -1;
                for (int k = 0; k < kNPB; k++) if (pt[a] >= kPB[k] && pt[a] < kPB[k + 1]) pb = k;
                if (pb >= 0)
                {
                    pb_w[pb] += w;
                    for (size_t b = 0; b < phi.size(); b++)
                    {
                        if (b == a) continue;
                        const double deta = etas[b] - etas[a];
                        for (int place = 0; place < 3; place++)
                        {
                            const double target = phi[a] + (place == 0 ? M_PI : (place == 1 ? 0.5 * M_PI : -0.5 * M_PI));
                            const double dphi = std::remainder(phi[b] - target, 2.0 * M_PI);
                            if (std::fabs(dphi) >= R) continue;
                            const bool in_cone = deta * deta + dphi * dphi < R * R;
                            for (int q = 0; q < kNSeen; q++)
                            {
                                if (pt[b] < kSeen[q]) continue;
                                if (place == 0) { bd_opp[pb][q] += w; bd_opp2[pb][q] += w * w; }
                                else { bd_side[pb][q] += 0.5 * w; bd_side2[pb][q] += 0.25 * w * w; }
                                if (!in_cone) continue;
                                if (place == 0) { pb_opp[pb][q] += w; pb_opp2[pb][q] += w * w; }
                                else { pb_side[pb][q] += 0.5 * w; pb_side2[pb][q] += 0.25 * w * w; }
                            }
                        }
                    }
                }
                for (int g = 0; g < kNTrig; g++)
                {
                    if (pt[a] < kTrig[g]) continue;
                    c_w[g] += w; c_w2[g] += w * w; c_n[g]++; trig_eta[g][eb] += w;
                    for (size_t b = 0; b < phi.size(); b++)
                    {
                        if (b == a) continue;
                        const double deta = etas[b] - etas[a];
                        if (std::fabs(deta) >= R) continue;
                        for (int place = 0; place < 3; place++) // 0: opposite, 1 and 2: sideways
                        {
                            const double target = phi[a] + (place == 0 ? M_PI : (place == 1 ? 0.5 * M_PI : -0.5 * M_PI));
                            const double dphi = std::remainder(phi[b] - target, 2.0 * M_PI);
                            if (deta * deta + dphi * dphi >= R * R) continue;
                            for (int q = 0; q < kNSeen; q++)
                            {
                                if (pt[b] < kSeen[q]) continue;
                                if (place == 0) { c_opp[g][q] += w; c_opp2[g][q] += w * w; c_opp_n[g][q]++; }
                                else { c_side[g][q] += 0.5 * w; c_side2[g][q] += 0.25 * w * w; }
                            }
                        }
                    }
                }
            }

            for (size_t s = 0; s < sels.size(); s++)
            {
                Sel &S = sels[s];
                int n_sel = 0;
                for (size_t a = 0; a < phi.size(); a++)
                {
                    if (pt[a] < S.pt_min) continue;
                    n_sel++;
                    S.sw += w; S.sw2 += w * w; S.n_jets++;
                    for (size_t b = 0; b < phi.size(); b++)
                    {
                        if (b == a || pt[b] < S.pt_min) continue;
                        double d = std::fabs(phi[a] - phi[b]) * 180.0 / M_PI;
                        if (d > 180.0) d = 360.0 - d;
                        const int k = std::min(kNDphi - 1, (int)(d / 5.0));
                        S.pair[k] += w; S.pair_w2[k] += w * w; S.pair_n[k]++;
                    }
                }
                if (n_sel >= 2) S.sw_events_ge2 += w;
            }
        }
        fs->Close();
        delete fs;
    }
    std::cout << "Read " << n_events << " events; truth jets with " << eta_min << " < eta < " << eta_max << ", all azimuth" << std::endl;

    // ---------------- numbers ----------------
    for (size_t s = 0; s < sels.size(); s++)
    {
        const Sel &S = sels[s];
        double tot = 0, side = 0; int n_side = 0;
        for (int k = 0; k < kNDphi; k++)
        {
            tot += S.pair[k];
            if (k >= 12 && k < 24) { side += S.pair[k]; n_side++; } // 60-120 deg
        }
        const int kb = kNDphi - 1;
        const double p_back = S.pair[kb] / S.sw, e_back = std::sqrt(S.pair_w2[kb]) / S.sw, p_side = side / n_side / S.sw;
        std::cout << "\n" << S.name << ": " << S.n_jets << " truth jets (" << S.sw * S.sw / S.sw2 << " effective), "
                  << S.sw / sw_events << " per event; events with >= 2: " << 100.0 * S.sw_events_ge2 / sw_events << " %" << std::endl;
        std::cout << "  other truth jets per truth jet, all angles          = " << 100.0 * tot / S.sw << " %" << std::endl;
        std::cout << "  in 175-180 deg (back-to-back)                       = " << 100.0 * p_back << " +- " << 100.0 * e_back
                  << " %  (" << S.pair_n[kb] << " raw pairs)" << std::endl;
        std::cout << "  average per 5 deg bin over 60-120 deg (sideways)    = " << 100.0 * p_side << " %" << std::endl;
        std::cout << "  back-to-back / sideways                             = " << p_back / p_side << " +- " << e_back / p_side << std::endl;
    }

    // ---------------- cone probabilities ----------------
    std::cout << "\nProbability of a truth-jet axis inside an R = " << R << " cone, per trigger jet:" << std::endl;
    for (int g = 0; g < kNTrig; g++)
    {
        std::cout << "\ntrigger truth jet pT >= " << kTrig[g] << " GeV: " << c_n[g] << " jets (" << c_w[g] * c_w[g] / c_w2[g] << " effective)" << std::endl;
        for (int q = 0; q < kNSeen; q++)
        {
            // Unrelated event: a jet at eta_j is inside a cone at (eta_t, random phi) for a fraction sqrt(R^2 - deta^2) / pi of azimuth
            double unrel = 0.0;
            for (int bt = 0; bt < kNEta; bt++)
            {
                if (trig_eta[g][bt] <= 0) continue;
                const double eta_t = kEtaLo + (bt + 0.5) * (kEtaHi - kEtaLo) / kNEta;
                double p = 0.0;
                for (int bj = 0; bj < kNEta; bj++)
                {
                    if (seen_eta[q][bj] <= 0) continue;
                    const double deta = kEtaLo + (bj + 0.5) * (kEtaHi - kEtaLo) / kNEta - eta_t;
                    if (std::fabs(deta) < R) p += seen_eta[q][bj] * std::sqrt(R * R - deta * deta) / M_PI;
                }
                unrel += trig_eta[g][bt] * p / sw_events;
            }
            unrel /= c_w[g];
            const double opp = c_opp[g][q] / c_w[g], e_opp = std::sqrt(c_opp2[g][q]) / c_w[g];
            const double side = c_side[g][q] / c_w[g], e_side = std::sqrt(c_side2[g][q]) / c_w[g];
            std::cout << "  jet seen with pT >= " << kSeen[q] << " GeV:  opposite cone " << 100.0 * opp << " +- " << 100.0 * e_opp << " % (" << c_opp_n[g][q]
                      << " raw)   sideways cone " << 100.0 * side << " +- " << 100.0 * e_side << " %   unrelated event " << 100.0 * unrel << " %" << std::endl;
            std::cout << "      opposite / sideways = " << (side > 0 ? opp / side : 0) << "   opposite / unrelated = " << opp / unrel
                      << " +- " << e_opp / unrel << std::endl;
        }
    }

    // ---------------- save the 5 deg histograms and plot ----------------
    TFile fout((out_base + ".root").c_str(), "RECREATE");
    std::vector<TH1D *> h;
    for (size_t s = 0; s < sels.size(); s++)
    {
        const Sel &S = sels[s];
        TH1D *hh = new TH1D(Form("dphi_%d", (int)s), "", kNDphi, 0, 180);
        for (int k = 0; k < kNDphi; k++)
        {
            hh->SetBinContent(k + 1, S.pair[k] / S.sw);
            hh->SetBinError(k + 1, std::sqrt(S.pair_w2[k]) / S.sw);
        }
        hh->Write();
        h.push_back(hh);
    }
    Plot(h, out_base, eta_min, eta_max);

    std::vector<TH1D *> hc;
    for (int q = 0; q < kNSeen; q++)
        for (int place = 0; place < 2; place++)
        {
            TH1D *hh = new TH1D(Form("cone_%s_%d", place == 0 ? "opp" : "side", q), "", kNPB, kPB);
            for (int k = 0; k < kNPB; k++)
            {
                if (pb_w[k] <= 0) continue;
                hh->SetBinContent(k + 1, 100.0 * (place == 0 ? pb_opp[k][q] : pb_side[k][q]) / pb_w[k]);
                hh->SetBinError(k + 1, 100.0 * std::sqrt(place == 0 ? pb_opp2[k][q] : pb_side2[k][q]) / pb_w[k]);
            }
            hh->Write();
            hc.push_back(hh);
        }
    PlotCone(hc, out_base + "_cone_prob", eta_min, eta_max, false);

    std::vector<TH1D *> hb;
    for (int q = 0; q < kNSeen; q++)
        for (int place = 0; place < 2; place++)
        {
            TH1D *hh = new TH1D(Form("band_%s_%d", place == 0 ? "opp" : "side", q), "", kNPB, kPB);
            for (int k = 0; k < kNPB; k++)
            {
                if (pb_w[k] <= 0) continue;
                hh->SetBinContent(k + 1, 100.0 * (place == 0 ? bd_opp[k][q] : bd_side[k][q]) / pb_w[k]);
                hh->SetBinError(k + 1, 100.0 * std::sqrt(place == 0 ? bd_opp2[k][q] : bd_side2[k][q]) / pb_w[k]);
            }
            hh->Write();
            hb.push_back(hh);
        }
    PlotCone(hb, out_base + "_band_prob", eta_min, eta_max, true);
    PlotIncrease(hb, out_base + "_band_increase", eta_min, eta_max, false);

    std::vector<TH1D *> he;
    for (int q = 0; q < kNSeen; q++)
        for (int place = 0; place < 2; place++)
        {
            TH1D *hh = new TH1D(Form("eband_%s_%d", place == 0 ? "opp" : "side", q), "", kNPB, kEB);
            for (int k = 0; k < kNPB; k++)
            {
                if (eb_w[k] <= 0) continue;
                hh->SetBinContent(k + 1, 100.0 * (place == 0 ? eb_opp[k][q] : eb_side[k][q]) / eb_w[k]);
                hh->SetBinError(k + 1, 100.0 * std::sqrt(place == 0 ? eb_opp2[k][q] : eb_side2[k][q]) / eb_w[k]);
            }
            hh->Write();
            he.push_back(hh);
        }
    PlotIncrease(he, out_base + "_band_increase_energy", eta_min, eta_max, true);
    fout.Close();
    return 0;
}
