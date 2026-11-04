// Yield and sanity plots from AnHistograms outputs (one file per job, one run
// per job), summed over the files given:
//   cutflow.png            jets at each cut
//   qa_jets.png            jet eta, phi, energies, pT, x_F, points: jets of
//                          the events passing the event cuts and the jets
//                          entering the yields
//   qa_xy.png              jet positions (both stages, and the selected jets'
//                          ratio to all jets)
//   qa_events.png          jets per event, vertex z and source, listed EM
//                          triggers, Spin_config
//   pt_levels.png          uncorrected and corrected pT by highest listed EM
//                          level, with the 1.15 x threshold cut
//   yields_xF_<level>.png  yields vs x_F per side and spin state, Poisson
//                          errors (levels EM0, EM1, EM2 and all)
//   yields_xF_classes.png  yields vs x_F per photon class and per EM level
//   runs.png, runs.txt     per run: events read, events passing, selected
//                          jets per passing event, x>0 / x<0 jets, spin-state
//                          fractions, EM level mix, mean corrected E
// plus all pages in plots.pdf.
//
// Usage: ./AnPlots <yield file list> <out dir>
// Build: g++ -m64 -O2 AnPlots.cpp -o AnPlots `root-config --cflags --libs` -std=c++11

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "AnConfig.h"

static std::map<std::string, TH1 *> gSum;
static std::string gOut;
static TCanvas *gPdfCanvas = 0;

static void AddFile(TFile &f)
{
    // files without spin_convention used the old Spin_config mapping: the
    // yields of state p belong to physical state 3 - p (AnConfig.h)
    const bool old_labels = !f.Get("spin_convention");
    TIter next(f.GetListOfKeys());
    while (TObject *k = next())
    {
        TObject *o = f.Get(k->GetName());
        if (!o || !o->InheritsFrom(TH1::Class())) continue;
        std::string name = k->GetName();
        if (old_labels && name.compare(0, 5, "h_xF_") == 0)
            for (int p = 0; p < kNSpinStates; p++)
            {
                const std::string from = std::string("_") + kSpinStateName[p];
                if (name.size() > from.size() && name.compare(name.size() - from.size(), from.size(), from) == 0)
                {
                    name = name.substr(0, name.size() - from.size()) + "_" + kSpinStateName[kNSpinStates - 1 - p];
                    break;
                }
            }
        std::map<std::string, TH1 *>::iterator it = gSum.find(name);
        if (it != gSum.end()) it->second->Add((TH1 *)o);
        else
        {
            TH1 *c = (TH1 *)o->Clone(name.c_str());
            c->SetDirectory(0);
            gSum[name] = c;
        }
    }
}

static TH1 *Get(const std::string &n)
{
    std::map<std::string, TH1 *>::iterator it = gSum.find(n);
    if (it == gSum.end()) { std::cerr << "WARNING: no " << n << std::endl; return 0; }
    return it->second;
}

static void Save(TCanvas *c, const char *name)
{
    c->SaveAs((gOut + "/" + name + ".png").c_str());
    c->Print((gOut + "/plots.pdf").c_str(), "pdf");
}

static const int kCol[8] = {kBlue + 1, kOrange + 7, kGreen + 2, kMagenta + 1, kCyan + 2, kRed + 1, kGray + 2, kViolet - 6};

static void Style(TH1 *h, int i, int marker = 20)
{
    h->SetLineColor(kCol[i % 8]);
    h->SetMarkerColor(kCol[i % 8]);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.7);
    h->SetLineWidth(2);
}

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <yield file list> <out dir>" << std::endl; return 1; }
    gOut = argv[2];
    gSystem->mkdir(gOut.c_str(), true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    gStyle->SetTitleFontSize(0.05);

    // per-run quantities, filled file by file
    struct RunRow { double ev_read, ev_pass, jets, xpos, xneg, spin[kNSpinStates], lvl[kNEmLevels], esum, en; };
    std::map<int, RunRow> rows;
    std::ifstream list(argv[1]);
    std::string fn;
    int nfiles = 0;
    while (list >> fn)
    {
        TFile f(fn.c_str());
        if (f.IsZombie()) { std::cerr << "WARNING: cannot open " << fn << std::endl; continue; }
        AddFile(f);
        nfiles++;
        TNamed *rn = (TNamed *)f.Get("runs");
        const int run = rn ? std::atoi(rn->GetTitle()) : 0;
        RunRow &r = rows[run];
        TH1 *cut = (TH1 *)f.Get("h_cutflow");
        if (cut) r.ev_read += cut->GetBinContent(cut->GetNbinsX());
        for (int l = 0; l < kNEmLevels; l++)
        {
            TH1 *he = (TH1 *)f.Get(TString::Format("h_events_%s", kEmLevelName[l]));
            if (!he) continue;
            // files without spin_convention used the old mapping: state p is physical 3 - p
            const bool old_labels = !f.Get("spin_convention");
            for (int p = 0; p < kNSpinStates; p++)
            {
                r.spin[old_labels ? kNSpinStates - 1 - p : p] += he->GetBinContent(p + 1);
                r.ev_pass += he->GetBinContent(p + 1);
            }
            for (int k = 0; k < kNPhotonClasses; k++)
                for (int s = 0; s < kNSides; s++)
                    for (int p = 0; p < kNSpinStates; p++)
                    {
                        TH1 *hy = (TH1 *)f.Get(TString::Format("h_xF_%s_%s_%s_%s", kEmLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]));
                        if (!hy) continue;
                        const double n = hy->Integral();
                        r.jets += n; r.lvl[l] += n;
                        (s ? r.xpos : r.xneg) += n;
                    }
        }
        TH1 *ec = (TH1 *)f.Get("qa_Ecorr_sel");
        if (ec) { r.esum += ec->GetMean() * ec->GetEntries(); r.en += ec->GetEntries(); }
    }
    if (nfiles == 0) { std::cerr << "ERROR: no files" << std::endl; return 1; }
    printf("AnPlots: %d files, %zu runs\n", nfiles, rows.size());

    TCanvas *c0 = new TCanvas("c0", "", 10, 10);
    c0->Print((gOut + "/plots.pdf[").c_str());
    gPdfCanvas = c0;

    // cut flow
    {
        TCanvas *c = new TCanvas("c_cut", "", 900, 600);
        c->SetLogy(); c->SetBottomMargin(0.15);
        TH1 *h = Get("h_cutflow");
        if (h) { h->SetTitle(TString::Format("cut flow (%d files, %zu runs)", nfiles, rows.size())); h->SetFillColor(kBlue - 9); h->Draw("hist text0"); }
        Save(c, "cutflow");
    }
    // jet distributions, two stages
    {
        const char *v[10] = {"eta", "phi", "E", "ueE", "Ecorr", "pt", "ptcorr", "xF", "np", "njet"};
        TCanvas *c = new TCanvas("c_qa", "", 1600, 1000);
        c->Divide(4, 3);
        TLegend *lg = 0;
        for (int i = 0; i < 10; i++)
        {
            c->cd(i + 1); gPad->SetLogy();
            TH1 *a = Get(std::string("qa_") + v[i] + "_evt"), *b = Get(std::string("qa_") + v[i] + "_sel");
            if (!a || !b) continue;
            Style(a, 0); Style(b, 1);
            a->SetTitle(""); a->SetMinimum(0.5); a->Draw("hist"); b->Draw("hist same");
            if (!lg)
            {
                lg = new TLegend(0.40, 0.80, 0.89, 0.89);
                lg->SetTextSize(0.035);
                lg->AddEntry(a, "jets, events passing event cuts", "l");
                lg->AddEntry(b, "jets in the A_{N} yields", "l");
                lg->SetBorderSize(0); lg->SetFillStyle(0);
            }
            lg->Draw();
            // uncorrected pT: the A_N trigger cut of each EM level
            if (std::string(v[i]) == "pt")
            {
                gPad->Update();
                const double y0 = std::pow(10, gPad->GetUymin()), y1 = std::pow(10, gPad->GetUymax() - 0.25 * (gPad->GetUymax() - gPad->GetUymin()));
                const int lc[kNEmLevels] = {kGreen + 2, kMagenta + 1, kRed + 1};
                for (int l = 0; l < kNEmLevels; l++)
                {
                    const double xc = kTriggerPtFactor * kEmThresholdGeV[l];
                    TLine *ln = new TLine(xc, y0, xc, y1);
                    ln->SetLineStyle(2); ln->SetLineWidth(2); ln->SetLineColor(lc[l]); ln->Draw();
                    TLatex tl; tl.SetTextSize(0.045); tl.SetTextColor(lc[l]); tl.SetTextAlign(11);
                    tl.DrawLatex(xc + 0.15, 3 * std::pow(10.0, l) * y0, kEmLevelName[l]);
                }
            }
        }
        // the two stages' cuts, in the free pads (A_N values from AnConfig.h; the
        // jet-finding ones are JetParameters.h's R, point_energy_min,
        // reco_jet_pt_min_em and reco_fiducial_buffer_ecal)
        {
            TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.050);
            const double dy = 0.068;
            c->cd(11);
            double y = 0.92;
            l.SetTextColor(kBlue + 1); l.DrawLatex(0.02, y, "#bf{Event cuts} (blue: every jet of these events)"); y -= dy; l.SetTextColor(kBlack);
            l.DrawLatex(0.04, y, "#bullet no FCS LED (fcs_led = 0)"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet an EM trigger recorded the event: fcsEM0/1/2"); y -= dy;
            l.DrawLatex(0.08, y, "(or its _tpc id) listed, i.e. after prescale"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet a spin state: Spin_config 5, 6, 9 or 10"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet run not flagged bad (sample drawn from good runs)"); y -= dy;
            y -= 0.03;
            l.DrawLatex(0.02, y, "#bf{Jet finding} (both curves)"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet anti-k_{T}, R = 0.7, FCS ECal points with E > 1 GeV"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet uncorrected p_{T} #geq 2 GeV"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet axis inside the ECal, 20 cm from its edges"); y -= dy;
            c->cd(12);
            y = 0.92;
            l.SetTextColor(kOrange + 7); l.DrawLatex(0.02, y, "#bf{A_{N} yield cuts} (orange), per jet, in addition"); y -= dy; l.SetTextColor(kBlack);
            l.DrawLatex(0.04, y, "#bullet calibrated (UE subtracted, JES applied): E_{corr} > 0"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet E_{corr} #leq #sqrt{s}/2 and 0 < x_{F} #leq 1"); y -= dy;
            l.DrawLatex(0.04, y, TString::Format("#bullet corrected p_{T} > %.0f GeV", kJetPtCorrMin)); y -= dy;
            l.DrawLatex(0.04, y, TString::Format("#bullet uncorrected p_{T} > %.2f #times the threshold of an EM", kTriggerPtFactor)); y -= dy;
            l.DrawLatex(0.08, y, "trigger listed in the event (dashed lines):"); y -= dy;
            l.DrawLatex(0.08, y, TString::Format("EM0 > %.2f, EM1 > %.2f, EM2 > %.2f GeV",
                                                 kTriggerPtFactor * kEmThresholdGeV[0], kTriggerPtFactor * kEmThresholdGeV[1],
                                                 kTriggerPtFactor * kEmThresholdGeV[2])); y -= dy;
        }
        Save(c, "qa_jets");
    }
    // positions
    {
        TCanvas *c = new TCanvas("c_xy", "", 1650, 500);
        c->Divide(3, 1);
        TH2D *a = (TH2D *)Get("qa_xy_evt"), *b = (TH2D *)Get("qa_xy_sel");
        if (a && b)
        {
            c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.14); a->SetTitle("jets, events passing event cuts"); a->Draw("colz");
            c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.14); b->SetTitle("jets in the A_{N} yields"); b->Draw("colz");
            TH2D *r = (TH2D *)b->Clone("qa_xy_ratio");
            r->Divide(a);
            r->SetTitle("selected / all jets");
            r->SetMaximum(1);
            c->cd(3); gPad->SetRightMargin(0.14); r->Draw("colz");
        }
        Save(c, "qa_xy");
    }
    // event quantities
    {
        TCanvas *c = new TCanvas("c_ev", "", 1500, 900);
        c->Divide(3, 2);
        const char *n[6] = {"qa_njet_evt", "qa_njet_sel", "qa_vz", "qa_vz_source", "qa_em_listed", "qa_spin_config"};
        const char *t[6] = {"jets per event (event cuts passed)", "selected jets per event", "vertex z", "vertex source",
                            "listed EM triggers (no LED)", "Spin_config (EM-triggered)"};
        for (int i = 0; i < 6; i++)
        {
            c->cd(i + 1);
            TH1 *h = Get(n[i]);
            if (!h) continue;
            if (i != 2) gPad->SetLogy();
            Style(h, 0); h->SetTitle(t[i]); h->SetFillColor(kBlue - 9); h->Draw("hist");
            if (i == 2) { TLatex l; l.SetNDC(); l.SetTextSize(0.045); l.DrawLatex(0.14, 0.84, TString::Format("mean %.1f, RMS %.1f cm", h->GetMean(), h->GetRMS())); }
        }
        Save(c, "qa_events");
    }
    // pT by level, with the threshold cut
    {
        TCanvas *c = new TCanvas("c_pt", "", 1400, 550);
        c->Divide(2, 1);
        const char *pre[2] = {"h_ptraw_", "h_ptcorr_"};
        const char *tit[2] = {"uncorrected p_{T}, jets passing the corrected p_{T} > 2 GeV cut", "corrected p_{T}, same jets"};
        for (int i = 0; i < 2; i++)
        {
            c->cd(i + 1); gPad->SetLogy();
            TLegend *lg = new TLegend(0.55, 0.70, 0.88, 0.88);
            lg->SetBorderSize(0); lg->SetFillStyle(0);
            for (int l = 0; l < kNEmLevels; l++)
            {
                TH1 *h = Get(std::string(pre[i]) + kEmLevelName[l]);
                if (!h) continue;
                Style(h, l); h->SetTitle(tit[i]); h->SetMinimum(0.5);
                h->Draw(l ? "hist same" : "hist");
                lg->AddEntry(h, TString::Format("highest listed %s", kEmLevelName[l]), "l");
                if (i == 0)
                {
                    const double xc = kTriggerPtFactor * kEmThresholdGeV[l];
                    TLine *ln = new TLine(xc, 0.5, xc, h->GetMaximum());
                    ln->SetLineColor(kCol[l]); ln->SetLineStyle(2); ln->Draw();
                }
            }
            lg->Draw();
        }
        Save(c, "pt_levels");
    }
    // yields vs x_F per side and spin state
    for (int L = 0; L <= kNEmLevels; L++)
    {
        const int l0 = L < kNEmLevels ? L : 0, l1 = L < kNEmLevels ? L + 1 : kNEmLevels;
        const char *lname = L < kNEmLevels ? kEmLevelName[L] : "all";
        TCanvas *c = new TCanvas(TString::Format("c_y%d", L), "", 1000, 700);
        c->SetLogy();
        TLegend *lg = new TLegend(0.62, 0.55, 0.88, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        bool first = true;
        for (int s = 0; s < kNSides; s++)
            for (int p = 0; p < kNSpinStates; p++)
            {
                TH1D *y = 0;
                for (int l = l0; l < l1; l++)
                    for (int k = 0; k < kNPhotonClasses; k++)
                    {
                        TH1 *h = Get(TString::Format("h_xF_%s_%s_%s_%s", kEmLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]).Data());
                        if (!h) continue;
                        if (!y) { y = (TH1D *)h->Clone(TString::Format("y_%d_%d_%d", L, s, p)); y->Sumw2(); }
                        else y->Add(h);
                    }
                if (!y) continue;
                Style(y, p, s ? 20 : 24);
                y->SetTitle(TString::Format("yields, %s levels, all photon classes (Poisson errors);corrected x_{F};jets per bin", lname));
                y->SetMinimum(0.5);
                y->GetXaxis()->SetRangeUser(0, 0.8);
                y->Draw(first ? "e1" : "e1 same");
                first = false;
                lg->AddEntry(y, TString::Format("%s, %s", kSideName[s], kSpinStateName[p]), "pe");
            }
        lg->Draw();
        Save(c, TString::Format("yields_xF_%s", lname).Data());
    }
    // yields vs x_F per photon class and per level
    {
        TCanvas *c = new TCanvas("c_cls", "", 1400, 600);
        c->Divide(2, 1);
        for (int mode = 0; mode < 2; mode++)
        {
            c->cd(mode + 1); gPad->SetLogy();
            TLegend *lg = new TLegend(0.62, 0.6, 0.88, 0.88);
            lg->SetBorderSize(0); lg->SetFillStyle(0);
            const int nc = mode == 0 ? kNPhotonClasses : kNEmLevels;
            for (int i = 0; i < nc; i++)
            {
                TH1D *y = 0;
                for (int l = 0; l < kNEmLevels; l++)
                    for (int k = 0; k < kNPhotonClasses; k++)
                    {
                        if ((mode == 0 && k != i) || (mode == 1 && l != i)) continue;
                        for (int s = 0; s < kNSides; s++)
                            for (int p = 0; p < kNSpinStates; p++)
                            {
                                TH1 *h = Get(TString::Format("h_xF_%s_%s_%s_%s", kEmLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]).Data());
                                if (!h) continue;
                                if (!y) { y = (TH1D *)h->Clone(TString::Format("yc_%d_%d", mode, i)); y->Sumw2(); }
                                else y->Add(h);
                            }
                    }
                if (!y) continue;
                Style(y, i);
                y->SetTitle(mode == 0 ? "yields by points in the jet;corrected x_{F};jets per bin" : "yields by EM level;corrected x_{F};jets per bin");
                y->SetMinimum(0.5);
                y->GetXaxis()->SetRangeUser(0, 0.8);
                y->Draw(i ? "e1 same" : "e1");
                lg->AddEntry(y, mode == 0 ? PhotonClassName(i) : kEmLevelName[i], "pe");
            }
            lg->Draw();
        }
        Save(c, "yields_xF_classes");
    }
    // per-run trends
    {
        FILE *t = std::fopen((gOut + "/runs.txt").c_str(), "w");
        std::fprintf(t, "# run events_read events_pass sel_jets jets_per_event xpos/xneg BuYu BuYd BdYu BdYd EM0 EM1 EM2 <Ecorr>\n");
        const int nr = rows.size();
        TGraphErrors *g[6];
        const char *gt[6] = {"events read (EM-filtered SimpleTree events)", "fraction passing the event cuts",
                             "selected jets per passing event", "jets at x > 0 / x < 0",
                             "spin-state fraction of passing events", "mean corrected E of selected jets (GeV)"};
        for (int i = 0; i < 6; i++) { g[i] = new TGraphErrors(); g[i]->SetTitle(TString::Format("%s;run - 23000000", gt[i])); }
        TGraphErrors *gs[kNSpinStates], *gl[kNEmLevels];
        for (int p = 0; p < kNSpinStates; p++) gs[p] = new TGraphErrors();
        for (int l = 0; l < kNEmLevels; l++) gl[l] = new TGraphErrors();
        for (std::map<int, RunRow>::const_iterator it = rows.begin(); it != rows.end(); ++it)
        {
            const RunRow &r = it->second;
            const double x = it->first - 23000000;
            const int i = g[0]->GetN();
            g[0]->SetPoint(i, x, r.ev_read);
            g[1]->SetPoint(i, x, r.ev_read > 0 ? r.ev_pass / r.ev_read : 0);
            g[2]->SetPoint(i, x, r.ev_pass > 0 ? r.jets / r.ev_pass : 0);
            g[2]->SetPointError(i, 0, r.ev_pass > 0 ? std::sqrt(r.jets) / r.ev_pass : 0);
            const double lr = r.xneg > 0 ? r.xpos / r.xneg : 0;
            g[3]->SetPoint(i, x, lr);
            g[3]->SetPointError(i, 0, r.xpos > 0 && r.xneg > 0 ? lr * std::sqrt(1 / r.xpos + 1 / r.xneg) : 0);
            g[5]->SetPoint(i, x, r.en > 0 ? r.esum / r.en : 0);
            for (int p = 0; p < kNSpinStates; p++) gs[p]->SetPoint(i, x, r.ev_pass > 0 ? r.spin[p] / r.ev_pass : 0);
            for (int l = 0; l < kNEmLevels; l++) gl[l]->SetPoint(i, x, r.jets > 0 ? r.lvl[l] / r.jets : 0);
            std::fprintf(t, "%d %.0f %.0f %.0f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.2f\n", it->first, r.ev_read, r.ev_pass, r.jets,
                         r.ev_pass > 0 ? r.jets / r.ev_pass : 0, lr, r.spin[0] / std::max(1.0, r.ev_pass), r.spin[1] / std::max(1.0, r.ev_pass),
                         r.spin[2] / std::max(1.0, r.ev_pass), r.spin[3] / std::max(1.0, r.ev_pass), r.lvl[0] / std::max(1.0, r.jets),
                         r.lvl[1] / std::max(1.0, r.jets), r.lvl[2] / std::max(1.0, r.jets), r.en > 0 ? r.esum / r.en : 0);
        }
        std::fclose(t);
        TCanvas *c = new TCanvas("c_runs", "", 1600, 1000);
        c->Divide(3, 2);
        for (int i = 0; i < 6; i++)
        {
            c->cd(i + 1);
            gPad->SetGridy();
            if (i == 4)
            {
                TLegend *lg = new TLegend(0.15, 0.15, 0.85, 0.3);
                lg->SetNColumns(4); lg->SetBorderSize(0); lg->SetFillStyle(0);
                for (int p = 0; p < kNSpinStates; p++)
                {
                    gs[p]->SetMarkerStyle(20); gs[p]->SetMarkerSize(0.7); gs[p]->SetMarkerColor(kCol[p]); gs[p]->SetLineColor(kCol[p]);
                    gs[p]->SetTitle(TString::Format("%s;run - 23000000", gt[4]));
                    gs[p]->SetMinimum(0); gs[p]->SetMaximum(0.5);
                    gs[p]->Draw(p ? "P same" : "AP");
                    lg->AddEntry(gs[p], kSpinStateName[p], "p");
                }
                lg->Draw();
                continue;
            }
            g[i]->SetMarkerStyle(20); g[i]->SetMarkerSize(0.7); g[i]->SetMarkerColor(kCol[0]); g[i]->SetLineColor(kCol[0]);
            g[i]->Draw("AP");
        }
        Save(c, "runs");
        TCanvas *c2 = new TCanvas("c_runs_lvl", "", 900, 600);
        c2->SetGridy();
        TLegend *lg = new TLegend(0.15, 0.75, 0.6, 0.88);
        lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0);
        for (int l = 0; l < kNEmLevels; l++)
        {
            gl[l]->SetMarkerStyle(20); gl[l]->SetMarkerSize(0.7); gl[l]->SetMarkerColor(kCol[l]); gl[l]->SetLineColor(kCol[l]);
            gl[l]->SetTitle("EM level of the selected jets, fraction per run;run - 23000000;fraction");
            gl[l]->SetMinimum(0); gl[l]->SetMaximum(1);
            gl[l]->Draw(l ? "P same" : "AP");
            lg->AddEntry(gl[l], kEmLevelName[l], "p");
        }
        lg->Draw();
        Save(c2, "runs_em_levels");
    }
    c0->Print((gOut + "/plots.pdf]").c_str());
    printf("wrote %s/*.png, plots.pdf, runs.txt\n", gOut.c_str());
    return 0;
}
