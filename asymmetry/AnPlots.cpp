// Yield and sanity plots from AnHistograms outputs (one file per job, one run
// per job), summed over the files given:
//   cutflow.png            jets at each cut
//   qa_jets.png            jet eta, phi, energies, pT, x_F, points: jets of
//                          the events passing the event cuts and the jets
//                          entering the yields (grey lines: sqrt(s)/2 and
//                          x_F = 1, past which jets are kept apart)
//   qa_jets_trigger.png    the same per EM level (EM0, EM1, EM2)
//   yields_hiE.png         yields vs x_F per EM level, E <= sqrt(s)/2 and
//                          above it (the _hiE histograms)
//   companions_pairs.png   companion jets (the _comp histograms: below every
//                          listed threshold, in events with a jet above one)
//                          beside the triggering jets; two-photon jets' pair
//                          mass and orientation (grey: the pi0 mass)
//   qa_xy.png              jet positions (both stages, and the selected jets'
//                          ratio to all jets)
//   qa_events.png          jets per event, vertex z and source, listed EM
//                          triggers, Spin_config
//   pt_levels.png          uncorrected and corrected pT by the event's level
//   trigger_values.png     per level, the jets' trigger value / nominal
//                          threshold (slice edges dotted, the default
//                          factor CutFactor() solid)
//   yields_xF_<level>.png  yields vs x_F per side and spin state, Poisson
//                          errors (levels EM0, EM1, EM2 and all)
//   yields_xF_classes.png  yields vs x_F per photon class and per EM level
//   runs.png, runs.txt     per run: events read, events passing, selected
//                          jets per passing event, x>0 / x<0 jets, spin-state
//                          fractions, EM level mix, mean corrected E
// plus all pages in plots.pdf.
//
// The runs of AnConfig.h kSuppBadRuns (left out of A_N) are left out of the
// summed pages too, and kept in the per-run rows (runs.png, runs.txt).
//
// Two steps: reading the yields (the slow part, ~40 min over a full report's
// files) writes everything the pages use to <out dir>/qa_sums.root (the
// summed histograms, the per-run rows as TTree "runs", the file count); -s
// draws the pages again from that file alone, in seconds, for cosmetic
// changes.
//
// Usage: ./AnPlots <yield file list> <out dir>      read the yields, save qa_sums.root, draw
//        ./AnPlots -s <qa_sums.root> <out dir>      draw from a saved qa_sums.root
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
#include <TTree.h>

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

// The default-yield name of a triggering-slice histogram: the slice suffix
// (and a JP level's EM tag) removed, if the slice is at or above the level's default factor (CutFactor); "" otherwise (or if
// the name has no slice suffix). AddFile also sums such slices under that
// name, so the pages below read the default yields by their plain names.
static std::string DefaultName(const std::string &name)
{
    for (int g = 0; g < kNTrgSlices; g++)
    {
        const std::string sfx = kTrgSliceSuffix[g];
        const size_t at = name.find(sfx);
        if (at == std::string::npos) continue;
        std::string rest = name.substr(at + sfx.size());
        // JP levels: the EM tag after the slice (every tag summed)
        for (int t = 0; t < kNEmTags; t++)
            if (rest.compare(0, std::string(kEmTagSuffix[t]).size(), kEmTagSuffix[t]) == 0) { rest = rest.substr(std::string(kEmTagSuffix[t]).size()); break; }
        if (!rest.empty() && rest != kHiESuffix) continue;
        const bool jp = name.find("_JP") != std::string::npos;
        if (kTrgSliceEdge[g] < (jp ? kJpCutFactor : kEmCutFactor) - 1e-9) return "";
        return name.substr(0, at) + rest;
    }
    return "";
}

static void AddTo(const std::string &name, TH1 *o)
{
    std::map<std::string, TH1 *>::iterator it = gSum.find(name);
    if (it != gSum.end()) it->second->Add(o);
    else
    {
        TH1 *c = (TH1 *)o->Clone(name.c_str());
        c->SetDirectory(0);
        gSum[name] = c;
    }
}

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
        AddTo(name, (TH1 *)o);
        const std::string dn = DefaultName(name);
        if (!dn.empty()) AddTo(dn, (TH1 *)o);
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

// sqrt(s)/2 (GeV) from qa_sqrt_s, 0 if the files have none (made before 2026-11-04)
static double HalfSqrtS()
{
    std::map<std::string, TH1 *>::iterator it = gSum.find("qa_sqrt_s");
    return it != gSum.end() && it->second->GetEntries() > 0 ? 0.5 * it->second->GetMean() : 0;
}

// on the current pad: a grey line at sqrt(s)/2 on the corrected-energy panel
// and at x_F = 1 on the x_F panel (jets past them are the _hiE ones)
static void HalfSqrtSLine(const std::string &var)
{
    double xc = 0;
    if (var == "Ecorr") xc = HalfSqrtS();
    else if (var == "xF") xc = 1.0;
    if (!(xc > 0)) return;
    gPad->Update();
    const double y0 = std::pow(10, gPad->GetUymin()), y1 = std::pow(10, gPad->GetUymax());
    TLine *ln = new TLine(xc, y0, xc, y1);
    ln->SetLineColor(kGray + 1); ln->SetLineWidth(2); ln->Draw();
}

int main(int argc, char **argv)
{
    const bool from_sums = argc > 3 && std::string(argv[1]) == "-s";
    if (argc < 3 || (std::string(argv[1]) == "-s" && argc < 4))
    {
        std::cerr << "Usage: " << argv[0] << " <yield file list> <out dir>\n       " << argv[0] << " -s <qa_sums.root> <out dir>" << std::endl;
        return 1;
    }
    const char *in_arg = from_sums ? argv[2] : argv[1];
    gOut = from_sums ? argv[3] : argv[2];
    gSystem->mkdir(gOut.c_str(), true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    gStyle->SetTitleFontSize(0.05);

    // per-run quantities, filled file by file
    struct RunRow { double ev_read, ev_pass, jets, xpos, xneg, spin[kNSpinStates], lvl[kNTrgLevels], esum, en; };
    std::map<int, RunRow> rows;
    std::ifstream list(from_sums ? "" : in_arg);
    std::string fn;
    int nfiles = 0, n_supp = 0;
    if (from_sums)
    {
        // everything the pages need, as saved by an earlier run on the yields
        TFile in(in_arg);
        if (in.IsZombie()) { std::cerr << "ERROR: cannot open " << in_arg << std::endl; return 1; }
        TIter next(in.GetListOfKeys());
        while (TObject *k = next())
        {
            TObject *o = in.Get(k->GetName());
            if (o && o->InheritsFrom(TH1::Class())) { TH1 *h = (TH1 *)o; h->SetDirectory(0); gSum[k->GetName()] = h; }
        }
        TNamed *nf = (TNamed *)in.Get("n_files");
        nfiles = nf ? std::atoi(nf->GetTitle()) : 0;
        TTree *t = (TTree *)in.Get("runs");
        if (t)
        {
            int run; RunRow r;
            t->SetBranchAddress("run", &run);
            t->SetBranchAddress("ev_read", &r.ev_read); t->SetBranchAddress("ev_pass", &r.ev_pass);
            t->SetBranchAddress("jets", &r.jets); t->SetBranchAddress("xpos", &r.xpos); t->SetBranchAddress("xneg", &r.xneg);
            t->SetBranchAddress("spin", r.spin); t->SetBranchAddress("lvl", r.lvl);
            t->SetBranchAddress("esum", &r.esum); t->SetBranchAddress("en", &r.en);
            for (Long64_t i = 0; i < t->GetEntries(); i++) { t->GetEntry(i); rows[run] = r; }
        }
        printf("AnPlots: from %s (%d files, %zu runs, %zu histograms)\n", in_arg, nfiles, rows.size(), gSum.size());
    }
    while (!from_sums && list >> fn)
    {
        TFile f(fn.c_str());
        if (f.IsZombie()) { std::cerr << "WARNING: cannot open " << fn << std::endl; continue; }
        TNamed *rn = (TNamed *)f.Get("runs");
        const int run = rn ? std::atoi(rn->GetTitle()) : 0;
        // runs AnAsymmetry always leaves out (AnConfig.h kSuppBadRuns): kept
        // in the per-run rows (runs.png, runs.txt), not in the summed pages
        if (SuppBadRunReason(run)) n_supp++;
        else AddFile(f);
        nfiles++;
        RunRow &r = rows[run];
        TH1 *cut = (TH1 *)f.Get("h_cutflow");
        if (cut) r.ev_read += cut->GetBinContent(cut->GetNbinsX());
        for (int l = 0; l < kNTrgLevels; l++)
        {
            TH1 *he = (TH1 *)f.Get(TString::Format("h_events_%s", kTrgLevelName[l]));
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
                        double n = 0;
                        for (int g = 0; g < kNTrgSlices; g++)
                        {
                            if (kTrgSliceEdge[g] < CutFactor(l) - 1e-9) continue;
                            // JP levels: only jets not firing EM at the default factor (each jet once)
                            for (int tg = 0; tg < (IsJpLevel(l) ? kNEmTags : 1); tg++)
                            {
                                if (IsJpLevel(l) && tg > 0 && kTrgSliceEdge[tg - 1] >= kEmCutFactor - 1e-9) continue;
                                TH1 *hy = (TH1 *)f.Get(TString::Format("h_xF_%s_%s_%s_%s%s%s", kTrgLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p],
                                                                       kTrgSliceSuffix[g], IsJpLevel(l) ? kEmTagSuffix[tg] : ""));
                                if (hy) n += hy->Integral();
                            }
                        }
                        r.jets += n; r.lvl[l] += n;
                        (s ? r.xpos : r.xneg) += n;
                    }
        }
        TH1 *ec = (TH1 *)f.Get("qa_Ecorr_sel");
        if (ec) { r.esum += ec->GetMean() * ec->GetEntries(); r.en += ec->GetEntries(); }
    }
    if (nfiles == 0) { std::cerr << "ERROR: no files" << std::endl; return 1; }
    if (!from_sums)
    {
        printf("AnPlots: %d files, %zu runs; %d files of AnConfig.h kSuppBadRuns runs only in the per-run rows\n", nfiles, rows.size(), n_supp);
        // save what the pages use, so they can be redrawn with -s
        TFile out((gOut + "/qa_sums.root").c_str(), "RECREATE");
        for (std::map<std::string, TH1 *>::const_iterator it = gSum.begin(); it != gSum.end(); ++it) it->second->Write(it->first.c_str());
        TNamed("n_files", TString::Format("%d", nfiles).Data()).Write();
        TTree t("runs", "per-run rows of AnPlots");
        int run; RunRow r;
        t.Branch("run", &run, "run/I");
        t.Branch("ev_read", &r.ev_read, "ev_read/D"); t.Branch("ev_pass", &r.ev_pass, "ev_pass/D");
        t.Branch("jets", &r.jets, "jets/D"); t.Branch("xpos", &r.xpos, "xpos/D"); t.Branch("xneg", &r.xneg, "xneg/D");
        t.Branch("spin", r.spin, TString::Format("spin[%d]/D", kNSpinStates));
        t.Branch("lvl", r.lvl, TString::Format("lvl[%d]/D", kNTrgLevels));
        t.Branch("esum", &r.esum, "esum/D"); t.Branch("en", &r.en, "en/D");
        for (std::map<int, RunRow>::const_iterator it = rows.begin(); it != rows.end(); ++it) { run = it->first; r = it->second; t.Fill(); }
        t.Write();
        out.Close();
        printf("AnPlots: wrote %s/qa_sums.root\n", gOut.c_str());
    }

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
            HalfSqrtSLine(v[i]);
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
            l.DrawLatex(0.04, y, "#bullet an EM or JP trigger recorded the event (fcsEM0/1/2,"); y -= dy;
            l.DrawLatex(0.08, y, "_tpc ids, fcsJP*): listed, i.e. after prescale"); y -= dy;
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
            l.DrawLatex(0.04, y, "#bullet x_{F} > 0; E_{corr} > #sqrt{s}/2 (grey line) kept apart,"); y -= dy;
            l.DrawLatex(0.08, y, "shown here, in A_{N} only with AnAsymmetry -E"); y -= dy;
            l.DrawLatex(0.04, y, TString::Format("#bullet corrected p_{T} > %.0f GeV", kJetPtCorrMin)); y -= dy;
            l.DrawLatex(0.04, y, TString::Format("#bullet a trigger level fired: the jet's trigger value > %.2f (EM) / %.2f (JP) #times", kEmCutFactor, kJpCutFactor)); y -= dy;
            l.DrawLatex(0.08, y, "the nominal threshold of a listed level (EM: 4x4 window,"); y -= dy;
            l.DrawLatex(0.08, y, "HCal veto margin; JP: patch, listed region); trigger_values"); y -= dy;
        }
        Save(c, "qa_jets");
    }
    // the same jet distributions per EM level: dashed, every jet of the events
    // whose highest listed trigger is that level; solid, the jets in that
    // level's yields
    {
        const char *v[10] = {"eta", "phi", "E", "ueE", "Ecorr", "pt", "ptcorr", "xF", "np", "njet"};
        const int lc[kNTrgLevels] = {kGreen + 2, kMagenta + 1, kRed + 1, kCyan + 2, kBlue + 1, kOrange + 7};
        TCanvas *c = new TCanvas("c_qa_trg", "", 1600, 1000);
        c->Divide(4, 3);
        TLegend *lg = new TLegend(0.02, 0.30, 0.98, 0.95);
        lg->SetTextSize(0.055); lg->SetBorderSize(0); lg->SetFillStyle(0);
        for (int i = 0; i < 10; i++)
        {
            c->cd(i + 1); gPad->SetLogy();
            double ymax = 0;
            TH1 *hs[kNTrgLevels][2] = {};
            for (int l = 0; l < kNTrgLevels; l++)
                for (int g = 0; g < 2; g++)
                {
                    std::map<std::string, TH1 *>::iterator it = gSum.find(std::string("qa_") + v[i] + (g ? "_sel_" : "_evt_") + kTrgLevelName[l]);
                    hs[l][g] = it == gSum.end() ? 0 : it->second;
                    if (hs[l][g]) ymax = std::max(ymax, hs[l][g]->GetMaximum());
                }
            bool first = true;
            for (int l = 0; l < kNTrgLevels; l++)
                for (int g = 0; g < 2; g++)
                {
                    TH1 *h = hs[l][g];
                    if (!h) continue;
                    h->SetLineColor(lc[l]); h->SetLineWidth(2); h->SetLineStyle(g ? 1 : 2);
                    h->SetTitle(""); h->SetMinimum(0.5); h->SetMaximum(2 * ymax);
                    h->Draw(first ? "hist" : "hist same");
                    first = false;
                    if (i == 0)
                        lg->AddEntry(h, g ? TString::Format("%s: jets in the %s yields", kTrgLevelName[l], kTrgLevelName[l]).Data()
                                          : TString::Format("%s: every jet, highest listed trigger %s", kTrgLevelName[l], kTrgLevelName[l]).Data(), "l");
                }
            if (first)
            {
                if (i == 0) std::cerr << "WARNING: no per-level QA histograms (yield files made before 2026-11-04); qa_jets_trigger is empty" << std::endl;
                continue;
            }
            HalfSqrtSLine(v[i]);
        }
        c->cd(11);
        lg->Draw();
        c->cd(12);
        {
            TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.050);
            double y = 0.92;
            const double dy = 0.068;
            l.DrawLatex(0.02, y, "#bf{Per EM trigger level}"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet dashed: every jet of the events passing the"); y -= dy;
            l.DrawLatex(0.08, y, "event cuts, by the highest listed EM trigger"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet solid: jets in the A_{N} yields, by their level (the"); y -= dy;
            l.DrawLatex(0.08, y, "level the jet fired: EM first, lowest listed)"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet jets per event: by the highest listed trigger"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet grey: #sqrt{s}/2, x_{F} = 1"); y -= dy;
            l.DrawLatex(0.04, y, "#bullet all levels together: qa_jets"); y -= dy;
        }
        Save(c, "qa_jets_trigger");
    }
    // jets above sqrt(s)/2 (the _hiE yields, left out of A_N unless -E), per
    // EM level, beside the default yields
    {
        TCanvas *c = new TCanvas("c_hie", "", 1000, 700);
        c->SetLogy();
        TLegend *lg = new TLegend(0.50, 0.65, 0.88, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        bool first = true, any_hie = false;
        for (int e = 0; e < 2; e++)
            for (int l = 0; l < kNTrgLevels; l++)
            {
                TH1D *y = 0;
                for (int k = 0; k < kNPhotonClasses; k++)
                    for (int s = 0; s < kNSides; s++)
                        for (int p = 0; p < kNSpinStates; p++)
                        {
                            std::map<std::string, TH1 *>::iterator it = gSum.find(TString::Format("h_xF_%s_%s_%s_%s%s", kTrgLevelName[l], PhotonClassName(k),
                                                                                                  kSideName[s], kSpinStateName[p], e ? kHiESuffix : "").Data());
                            if (it == gSum.end()) continue;
                            if (!y) { y = (TH1D *)it->second->Clone(TString::Format("yh_%d_%d", e, l)); y->Sumw2(); }
                            else y->Add(it->second);
                        }
                if (!y) continue;
                if (e) any_hie = true;
                Style(y, l, e ? 24 : 20);
                y->SetLineStyle(e ? 2 : 1);
                y->SetTitle(TString::Format("yields with E_{corr} #leq #sqrt{s}/2 (solid) and above it (dashed, %s);corrected x_{F};jets per bin",
                                            "left out of A_{N} unless AnAsymmetry -E"));
                y->SetMinimum(0.5);
                y->Draw(first ? "e1" : "e1 same");
                first = false;
                lg->AddEntry(y, TString::Format("%s, %s", kTrgLevelName[l], e ? "E_{corr} > #sqrt{s}/2" : "E_{corr} #leq #sqrt{s}/2"), "pe");
            }
        if (!any_hie) std::cerr << "WARNING: no _hiE yields (yield files made before 2026-11-04); yields_hiE shows the default yields only" << std::endl;
        lg->Draw();
        HalfSqrtSLine("xF");
        Save(c, "yields_hiE");
    }
    // companion jets (below every listed threshold, in events with a jet above
    // one; the _comp histograms) beside the triggering jets, and the two-photon
    // jets' pair mass and orientation
    {
        const bool have_comp = gSum.count("qa_pt_comp") > 0, have_pair = gSum.count("qa_mgg_sel") > 0;
        if (!have_comp) std::cerr << "WARNING: no companion histograms (yield files made before 2026-11-05); companions_pairs shows what exists" << std::endl;
        TCanvas *c = new TCanvas("c_comp", "", 1800, 900);
        c->Divide(4, 2);
        const int lc[kNTrgLevels] = {kGreen + 2, kMagenta + 1, kRed + 1, kCyan + 2, kBlue + 1, kOrange + 7};
        // 1: yields vs x_F, triggering (solid) and companion (dashed) jets per level
        c->cd(1); gPad->SetLogy();
        {
            TLegend *lg = new TLegend(0.45, 0.62, 0.89, 0.89);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
            bool first = true;
            for (int cmp = 0; cmp < 2; cmp++)
                for (int l = 0; l < kNTrgLevels; l++)
                {
                    TH1D *y = 0;
                    for (int k = 0; k < kNPhotonClasses; k++)
                        for (int s = 0; s < kNSides; s++)
                            for (int p = 0; p < kNSpinStates; p++)
                            {
                                std::map<std::string, TH1 *>::iterator it = gSum.find(TString::Format("h_xF_%s_%s_%s_%s%s", kTrgLevelName[l], PhotonClassName(k),
                                                                                                      kSideName[s], kSpinStateName[p], cmp ? kCompSuffix : "").Data());
                                if (it == gSum.end()) continue;
                                if (!y) { y = (TH1D *)it->second->Clone(TString::Format("yc_%d_%d", cmp, l)); y->Sumw2(); }
                                else y->Add(it->second);
                            }
                    if (!y) continue;
                    y->SetLineColor(lc[l]); y->SetMarkerColor(lc[l]); y->SetLineWidth(2); y->SetLineStyle(cmp ? 2 : 1);
                    y->SetMarkerStyle(cmp ? 24 : 20); y->SetMarkerSize(0.6);
                    y->SetTitle("yields, E_{corr} #leq #sqrt{s}/2;corrected x_{F};jets per bin");
                    y->SetMinimum(0.5);
                    y->Draw(first ? "e1" : "e1 same");
                    first = false;
                    lg->AddEntry(y, TString::Format("%s %s", kTrgLevelName[l], cmp ? "companion" : "triggering"), "pe");
                }
            lg->Draw();
        }
        // 2-3: uncorrected pT and points, triggering vs companion jets (all levels)
        const char *v2[2] = {"pt", "np"};
        for (int i = 0; i < 2; i++)
        {
            c->cd(2 + i); gPad->SetLogy();
            TH1 *a = gSum.count(std::string("qa_") + v2[i] + "_sel") ? gSum[std::string("qa_") + v2[i] + "_sel"] : 0;
            TH1 *b = gSum.count(std::string("qa_") + v2[i] + "_comp") ? gSum[std::string("qa_") + v2[i] + "_comp"] : 0;
            if (!a) continue;
            Style(a, 1); a->SetTitle("triggering (orange) and companion (blue) jets"); a->SetMinimum(0.5);
            a->Draw("hist");
            if (b) { Style(b, 0); b->Draw("hist same"); }
        }
        // 4-5: pair mass and orientation, three stages
        const char *v3[2] = {"mgg", "pairang"};
        const char *st[3] = {"evt", "sel", "comp"};
        for (int i = 0; i < 2; i++)
        {
            c->cd(4 + i);
            TLegend *lg = new TLegend(0.50, 0.72, 0.89, 0.89);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
            double ymax = 0;
            for (int g = 0; g < 3; g++)
                if (gSum.count(std::string("qa_") + v3[i] + "_" + st[g])) ymax = std::max(ymax, gSum[std::string("qa_") + v3[i] + "_" + st[g]]->GetMaximum());
            bool first = true;
            for (int g = 0; g < 3; g++)
            {
                const std::string n = std::string("qa_") + v3[i] + "_" + st[g];
                if (!gSum.count(n)) continue;
                TH1 *h = gSum[n];
                Style(h, g == 0 ? 6 : g == 1 ? 1 : 0);
                h->SetTitle(i ? "two-photon jets: pair orientation (0 = along x, 90 = along y)" : "two-photon jets: M_{#gamma#gamma}");
                h->SetMinimum(0); h->SetMaximum(1.15 * ymax);
                h->Draw(first ? "hist" : "hist same");
                first = false;
                lg->AddEntry(h, g == 0 ? "every jet (event cuts)" : g == 1 ? "triggering jets" : "companion jets", "l");
            }
            if (!first) lg->Draw();
            if (i == 0 && !first)
            {
                gPad->Update();
                TLine *ln = new TLine(0.135, gPad->GetUymin(), 0.135, gPad->GetUymax());
                ln->SetLineColor(kGray + 1); ln->SetLineStyle(2); ln->Draw();
            }
        }
        if (!have_pair) std::cerr << "WARNING: no pair histograms (yield files made before 2026-11-05)" << std::endl;
        // 6: mass per level (triggering jets); 7-8: mass vs E and vs orientation
        c->cd(6);
        {
            bool first = true;
            TLegend *lg = new TLegend(0.60, 0.72, 0.89, 0.89);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
            for (int l = 0; l < kNTrgLevels; l++)
            {
                const std::string n = std::string("qa_mgg_sel_") + kTrgLevelName[l];
                if (!gSum.count(n) || gSum[n]->Integral() <= 0) continue;
                TH1 *h = (TH1 *)gSum[n]->Clone((n + "_norm").c_str());
                h->Scale(1.0 / h->Integral());
                h->SetLineColor(lc[l]); h->SetLineWidth(2);
                h->SetTitle("triggering two-photon jets per level (unit area);M_{#gamma#gamma} (GeV);fraction");
                h->Draw(first ? "hist" : "hist same");
                first = false;
                lg->AddEntry(h, kTrgLevelName[l], "l");
            }
            if (!first) lg->Draw();
        }
        c->cd(7); gPad->SetLogz();
        if (gSum.count("qa_mggE_sel")) { gSum["qa_mggE_sel"]->SetTitle("triggering two-photon jets;corrected jet E (GeV);M_{#gamma#gamma} (GeV)"); gSum["qa_mggE_sel"]->Draw("colz"); }
        c->cd(8); gPad->SetLogz();
        if (gSum.count("qa_mggang_sel")) { gSum["qa_mggang_sel"]->SetTitle("triggering two-photon jets;pair orientation (deg);M_{#gamma#gamma} (GeV)"); gSum["qa_mggang_sel"]->Draw("colz"); }
        Save(c, "companions_pairs");
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
            for (int l = 0; l < kNTrgLevels; l++)
            {
                TH1 *h = Get(std::string(pre[i]) + kTrgLevelName[l]);
                if (!h) continue;
                Style(h, l); h->SetTitle(tit[i]); h->SetMinimum(0.5);
                h->Draw(l ? "hist same" : "hist");
                lg->AddEntry(h, TString::Format("highest listed %s", kTrgLevelName[l]), "l");
            }
            lg->Draw();
        }
        Save(c, "pt_levels");
    }
    // trigger value / nominal threshold per level (h_trgratio_*: triggering
    // and _veto jets, before the slice cut), slice edges dotted, the default
    // factor solid
    {
        TCanvas *c = new TCanvas("c_trg", "", 1500, 900);
        c->Divide(3, 2);
        const int lc[kNTrgLevels] = {kGreen + 2, kMagenta + 1, kRed + 1, kCyan + 2, kBlue + 1, kOrange + 7};
        bool any = false;
        for (int l = 0; l < kNTrgLevels; l++)
        {
            c->cd(l + 1); gPad->SetLogy();
            std::map<std::string, TH1 *>::iterator it = gSum.find(std::string("h_trgratio_") + kTrgLevelName[l]);
            if (it == gSum.end()) continue;
            any = true;
            TH1 *h = it->second;
            h->SetLineColor(lc[l]); h->SetLineWidth(2); h->SetMinimum(0.5);
            h->SetTitle(TString::Format("%s jets: trigger value / %.1f GeV (%s)", kTrgLevelName[l], kTrgThresholdGeV[l],
                                        IsJpLevel(l) ? "largest JP patch with the jet" : "largest ECal 4x4 in the jet"));
            h->Draw("hist");
            gPad->Update();
            const double y0 = std::pow(10, gPad->GetUymin()), y1 = std::pow(10, gPad->GetUymax());
            for (int g = 0; g < kNTrgSlices; g++)
            {
                TLine *ln = new TLine(kTrgSliceEdge[g], y0, kTrgSliceEdge[g], y1);
                const bool d = std::fabs(kTrgSliceEdge[g] - CutFactor(l)) < 1e-9;
                ln->SetLineColor(d ? kBlack : kGray + 1); ln->SetLineStyle(d ? 1 : 3); ln->SetLineWidth(2); ln->Draw();
            }
        }
        if (!any) std::cerr << "WARNING: no h_trgratio_* (yield files made before 2026-11-05); trigger_values is empty" << std::endl;
        Save(c, "trigger_values");
    }
    // yields vs x_F per side and spin state
    for (int L = 0; L <= kNTrgLevels; L++)
    {
        const int l0 = L < kNTrgLevels ? L : 0, l1 = L < kNTrgLevels ? L + 1 : kNTrgLevels;
        const char *lname = L < kNTrgLevels ? kTrgLevelName[L] : "all";
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
                        TH1 *h = Get(TString::Format("h_xF_%s_%s_%s_%s", kTrgLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]).Data());
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
            const int nc = mode == 0 ? kNPhotonClasses : kNTrgLevels;
            for (int i = 0; i < nc; i++)
            {
                TH1D *y = 0;
                for (int l = 0; l < kNTrgLevels; l++)
                    for (int k = 0; k < kNPhotonClasses; k++)
                    {
                        if ((mode == 0 && k != i) || (mode == 1 && l != i)) continue;
                        for (int s = 0; s < kNSides; s++)
                            for (int p = 0; p < kNSpinStates; p++)
                            {
                                TH1 *h = Get(TString::Format("h_xF_%s_%s_%s_%s", kTrgLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]).Data());
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
                lg->AddEntry(y, mode == 0 ? PhotonClassName(i) : kTrgLevelName[i], "pe");
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
        TGraphErrors *gs[kNSpinStates], *gl[kNTrgLevels];
        for (int p = 0; p < kNSpinStates; p++) gs[p] = new TGraphErrors();
        for (int l = 0; l < kNTrgLevels; l++) gl[l] = new TGraphErrors();
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
            for (int l = 0; l < kNTrgLevels; l++) gl[l]->SetPoint(i, x, r.jets > 0 ? r.lvl[l] / r.jets : 0);
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
        for (int l = 0; l < kNTrgLevels; l++)
        {
            gl[l]->SetMarkerStyle(20); gl[l]->SetMarkerSize(0.7); gl[l]->SetMarkerColor(kCol[l]); gl[l]->SetLineColor(kCol[l]);
            gl[l]->SetTitle("EM level of the selected jets, fraction per run;run - 23000000;fraction");
            gl[l]->SetMinimum(0); gl[l]->SetMaximum(1);
            gl[l]->Draw(l ? "P same" : "AP");
            lg->AddEntry(gl[l], kTrgLevelName[l], "p");
        }
        lg->Draw();
        Save(c2, "runs_em_levels");
    }
    c0->Print((gOut + "/plots.pdf]").c_str());
    printf("wrote %s/*.png, plots.pdf, runs.txt\n", gOut.c_str());
    return 0;
}
