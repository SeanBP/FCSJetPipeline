// A_N yield histograms from calibrated jet files (after jet finding, the
// underlying-event step and ApplyCorrections_em_only), so the jet and
// SimpleTree files need not be kept.
//
// Per event: no FCS LED (fcs_led), the highest EM trigger level (AnConfig.h)
// and the spin state (Spin_config); events failing any are dropped. Per jet:
// a valid calibration (reco_E_corr > 0), the physical cut
// (pass_physical_jet_cut: E <= sqrt(s)/2, 0 < x_F <= 1), corrected pT =
// reco_E_corr / cosh(eta) above kJetPtCorrMin and above kTriggerPtFactor
// times the event's EM threshold. Each passing jet
// fills x_F (reco_x_F_corr) in
//   h_xF_<EM level>_<photon class>_<side>_<spin state>
// e.g. h_xF_EM1_np2_xpos_BuYd. Also written:
//   h_cutflow                    jets (and events) at each cut
//   h_events_<EM level>          events per spin state, all EM-level events
//   h_ptcorr_<EM level>          corrected pT before the threshold cut
//   runs                         TNamed, the run numbers of the input (from
//                                the file names, jet_output_<run>_<file>.root)
// Output files of one run, or of a fill, add with hadd.
//
// Usage: ./AnHistograms <jet file list> <output.root>
// Build: g++ -m64 -O2 AnHistograms.cpp VectorDict.cxx -o AnHistograms
//            `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "VectorDict.h"
#include "AnConfig.h"

// The run number in a jet_output_<run>_<file>.root name, 0 if none
static int RunFromFileName(const std::string &path)
{
    const std::string base = path.substr(path.find_last_of('/') + 1);
    int run = 0, file_number = 0;
    return std::sscanf(base.c_str(), "jet_output_%d_%d.root", &run, &file_number) == 2 ? run : 0;
}

enum Cut { kAllJets, kNoLed, kEmLevel, kSpin, kCalibrated, kPhysical, kPtMin, kTriggerPt, kNCuts };
static const char *const kCutName[kNCuts] = {"all jets", "no LED", "EM trigger", "spin state", "calibrated", "physical", "pT > 2", "pT > 1.15 thr"};

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <jet file list> <output.root>" << std::endl;
        return 1;
    }
    TChain t("jetTree");
    {
        std::ifstream in(argv[1]);
        std::string f;
        while (in >> f) t.Add(f.c_str());
    }
    if (t.GetNtrees() == 0) { std::cerr << "ERROR: no files in " << argv[1] << std::endl; return 1; }
    if (!t.GetBranch("reco_E_corr")) { std::cerr << "ERROR: no reco_E_corr; run ApplyCorrections_em_only first" << std::endl; return 1; }

    Int_t spin_config = -1, led = 0, em[kNEmLevels] = {0};
    Double_t sqrt_s = 0;
    std::vector<float> *eta = 0, *x = 0, *e_corr = 0, *xf_corr = 0;
    std::vector<int> *nphotons = 0;
    t.SetBranchAddress("Spin_config", &spin_config);
    t.SetBranchAddress("fcs_led", &led);
    t.SetBranchAddress("event_sqrt_s_gev", &sqrt_s);
    t.SetBranchAddress("fcsEM0", &em[0]);
    t.SetBranchAddress("fcsEM1", &em[1]);
    t.SetBranchAddress("fcsEM2", &em[2]);
    t.SetBranchAddress("reco_eta", &eta);
    t.SetBranchAddress("reco_x", &x);
    t.SetBranchAddress("reco_E_corr", &e_corr);
    t.SetBranchAddress("reco_x_F_corr", &xf_corr);
    t.SetBranchAddress("reco_nphotons", &nphotons);

    TFile out(argv[2], "RECREATE");
    TH1D *h[kNEmLevels][kNPhotonClasses][kNSides][kNSpinStates];
    TH1D *h_events[kNEmLevels], *h_pt[kNEmLevels];
    for (int l = 0; l < kNEmLevels; l++)
    {
        for (int k = 0; k < kNPhotonClasses; k++)
            for (int s = 0; s < kNSides; s++)
                for (int p = 0; p < kNSpinStates; p++)
                {
                    TString name = TString::Format("h_xF_%s_%s_%s_%s", kEmLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]);
                    h[l][k][s][p] = new TH1D(name, name + ";x_{F};jets", kNXfBins, kXfMin, kXfMax);
                }
        h_events[l] = new TH1D(TString::Format("h_events_%s", kEmLevelName[l]), ";spin state;events", kNSpinStates, 0, kNSpinStates);
        for (int p = 0; p < kNSpinStates; p++) h_events[l]->GetXaxis()->SetBinLabel(p + 1, kSpinStateName[p]);
        h_pt[l] = new TH1D(TString::Format("h_ptcorr_%s", kEmLevelName[l]), ";corrected jet p_{T} (GeV);jets", 100, 0, 20);
    }
    TH1D *h_cut = new TH1D("h_cutflow", ";;jets", kNCuts + 1, 0, kNCuts + 1);
    for (int c = 0; c < kNCuts; c++) h_cut->GetXaxis()->SetBinLabel(c + 1, kCutName[c]);
    h_cut->GetXaxis()->SetBinLabel(kNCuts + 1, "events read");

    std::set<int> runs;
    int tree_number = -1;
    const Long64_t n = t.GetEntries();
    for (Long64_t i = 0; i < n; i++)
    {
        t.GetEntry(i);
        if (t.GetTreeNumber() != tree_number)
        {
            tree_number = t.GetTreeNumber();
            runs.insert(RunFromFileName(t.GetFile()->GetName()));
        }
        h_cut->Fill(kNCuts);
        const int level = em[2] ? 2 : em[1] ? 1 : em[0] ? 0 : -1;
        const int spin = SpinState(spin_config);
        const size_t nj = eta->size();
        h_cut->Fill(kAllJets, nj);
        if (led) continue;
        h_cut->Fill(kNoLed, nj);
        if (level < 0) continue;
        h_cut->Fill(kEmLevel, nj);
        if (spin < 0) continue;
        h_cut->Fill(kSpin, nj);
        h_events[level]->Fill(spin);
        const double pt_min = kTriggerPtFactor * kEmThresholdGeV[level];
        for (size_t j = 0; j < nj; j++)
        {
            const double E = (*e_corr)[j];
            if (!(E > 0)) continue;
            h_cut->Fill(kCalibrated);
            if (!pass_physical_jet_cut(E, (*eta)[j], sqrt_s)) continue;
            h_cut->Fill(kPhysical);
            const double pt = E / std::cosh((*eta)[j]);
            if (pt <= kJetPtCorrMin) continue;
            h_cut->Fill(kPtMin);
            h_pt[level]->Fill(pt);
            if (pt <= pt_min) continue;
            h_cut->Fill(kTriggerPt);
            const int side = computeSide((*x)[j]) > 0 ? 1 : 0;
            h[level][PhotonClass((*nphotons)[j])][side][spin]->Fill((*xf_corr)[j]);
        }
    }

    std::ostringstream rs;
    for (std::set<int>::const_iterator r = runs.begin(); r != runs.end(); ++r) rs << (r == runs.begin() ? "" : " ") << *r;
    TNamed("runs", rs.str().c_str()).Write();
    out.Write();
    printf("AnHistograms: %lld events, runs %s\n", n, rs.str().c_str());
    for (int c = 0; c < kNCuts; c++) printf("  %-14s %.0f jets\n", kCutName[c], h_cut->GetBinContent(c + 1));
    out.Close();
    return 0;
}
