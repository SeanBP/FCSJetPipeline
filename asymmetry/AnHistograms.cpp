// A_N yield histograms from calibrated jet files (after jet finding, the
// underlying-event step and ApplyCorrections_em_only), so the jet and
// SimpleTree files need not be kept.
//
// Per event: no FCS LED (fcs_led), at least one EM trigger listed
// (AnConfig.h) and a spin state (Spin_config); events failing any are
// dropped. Per jet: a valid calibration (reco_E_corr > 0), the physical cut
// (pass_physical_jet_cut: E <= sqrt(s)/2, 0 < x_F <= 1), corrected pT =
// reco_E_corr / cosh(eta) above kJetPtCorrMin, and, as the FMS EM-jet
// analysis, uncorrected pT (reco_pt: no UE subtraction, no JES, what the
// hardware trigger sees) above kTriggerPtFactor times the threshold of at
// least one EM trigger listed in the event. In data the listed triggers are the ones that
// recorded the event (MuDst nominal trigger ids: condition fired and the
// prescale accepted it); a trigger whose condition fired but was prescaled
// away is not listed and admits no jet. The jet's EM level is the lowest
// listed level whose threshold it passes. Each passing jet
// fills x_F (reco_x_F_corr) in
//   h_xF_<EM level>_<photon class>_<side>_<spin state>
// e.g. h_xF_EM1_np2_xpos_BuYd, and with weight |cos phi| (phi = atan2(reco_y,
// reco_x), the jet's azimuth; the same for both beams, whose left is +x or -x)
//   h_abscos_<EM level>_<photon class>
// summed over sides and spin states: divided by the summed yields it is the
// bin's mean |cos phi|, which a left-right asymmetry is diluted by
// (AnAsymmetry.cpp); likewise, with weight x_F,
//   h_xfsum_<EM level>_<photon class>
// for the bin's mean corrected x_F, where its A_N point is placed. Also written:
//   h_cutflow                    jets (and events) at each cut
//   h_events_<EM level>          events per spin state, by highest listed level
//   h_ptcorr_<EM level>          corrected pT before the threshold cut, by
//                                highest listed level
//   h_ptraw_<EM level>           uncorrected pT before the threshold cut, the
//                                same jets
//   runs                         TNamed, the run numbers of the input (from
//                                the file names, jet_output_<run>_<file>.root)
//   spin_convention              TNamed, the Spin_config -> spin state
//                                mapping used (AnConfig.h)
// and sanity distributions, qa_<variable>_<stage>, at two stages: "evt",
// every jet of the events passing the event cuts (no LED, an EM trigger, a
// spin state), and "sel", the jets entering the yields: eta, phi, energy
// (uncorrected reco_E, UE reco_ue_E, corrected), pT (uncorrected, corrected),
// corrected x_F, points in the jet, jet position (reco_x, reco_y), and per
// event the jets and selected jets, the vertex z (event_vertex_z, -9999 =
// none, and its source) and the listed EM triggers. qa_xy_sel_<side> is
// split by side for the left/right acceptance.
// Output files of one run, or of a fill, add with hadd.
//
// Usage: ./AnHistograms <jet file list> <output.root>
// Build: g++ -m64 -O2 AnHistograms.cpp VectorDict.cxx -o AnHistograms
//            `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
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

    Int_t spin_config = -1, led = 0, em[kNEmLevels] = {0}, em_tpc[kNEmLevels] = {0};
    Double_t sqrt_s = 0;
    Float_t vz = -9999;
    Int_t vz_source = 0;
    std::vector<float> *eta = 0, *phi = 0, *e_raw = 0, *e_ue = 0, *x = 0, *y = 0, *pt_raw = 0, *e_corr = 0, *xf_corr = 0;
    std::vector<int> *nphotons = 0;
    t.SetBranchAddress("Spin_config", &spin_config);
    t.SetBranchAddress("fcs_led", &led);
    t.SetBranchAddress("event_sqrt_s_gev", &sqrt_s);
    t.SetBranchAddress("fcsEM0", &em[0]);
    t.SetBranchAddress("fcsEM1", &em[1]);
    t.SetBranchAddress("fcsEM2", &em[2]);
    // the _tpc copies: the same EM conditions under their own trigger ids
    for (int l = 0; l < kNEmLevels; l++)
    {
        const TString b = TString::Format("fcsEM%d_tpc", l);
        if (t.GetBranch(b)) t.SetBranchAddress(b, &em_tpc[l]);
    }
    if (t.GetBranch("event_vertex_z")) t.SetBranchAddress("event_vertex_z", &vz);
    if (t.GetBranch("event_vertex_source")) t.SetBranchAddress("event_vertex_source", &vz_source);
    t.SetBranchAddress("reco_eta", &eta);
    t.SetBranchAddress("reco_phi", &phi);
    t.SetBranchAddress("reco_E", &e_raw);
    t.SetBranchAddress("reco_ue_E", &e_ue);
    t.SetBranchAddress("reco_x", &x);
    t.SetBranchAddress("reco_y", &y);
    t.SetBranchAddress("reco_pt", &pt_raw);
    t.SetBranchAddress("reco_E_corr", &e_corr);
    t.SetBranchAddress("reco_x_F_corr", &xf_corr);
    t.SetBranchAddress("reco_nphotons", &nphotons);

    TFile out(argv[2], "RECREATE");
    TH1D *h[kNEmLevels][kNPhotonClasses][kNSides][kNSpinStates];
    TH1D *h_cos[kNEmLevels][kNPhotonClasses], *h_xfs[kNEmLevels][kNPhotonClasses];
    TH1D *h_events[kNEmLevels], *h_pt[kNEmLevels], *h_ptraw[kNEmLevels];
    for (int l = 0; l < kNEmLevels; l++)
    {
        for (int k = 0; k < kNPhotonClasses; k++)
            for (int s = 0; s < kNSides; s++)
                for (int p = 0; p < kNSpinStates; p++)
                {
                    TString name = TString::Format("h_xF_%s_%s_%s_%s", kEmLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]);
                    h[l][k][s][p] = new TH1D(name, name + ";x_{F};jets", kNXfBins, kXfMin, kXfMax);
                }
        for (int k = 0; k < kNPhotonClasses; k++)
        {
            TString name = TString::Format("h_abscos_%s_%s", kEmLevelName[l], PhotonClassName(k));
            h_cos[l][k] = new TH1D(name, name + ";x_{F};sum of |cos #phi|", kNXfBins, kXfMin, kXfMax);
            name = TString::Format("h_xfsum_%s_%s", kEmLevelName[l], PhotonClassName(k));
            h_xfs[l][k] = new TH1D(name, name + ";x_{F};sum of x_{F}", kNXfBins, kXfMin, kXfMax);
        }
        h_events[l] = new TH1D(TString::Format("h_events_%s", kEmLevelName[l]), ";spin state;events", kNSpinStates, 0, kNSpinStates);
        for (int p = 0; p < kNSpinStates; p++) h_events[l]->GetXaxis()->SetBinLabel(p + 1, kSpinStateName[p]);
        h_pt[l] = new TH1D(TString::Format("h_ptcorr_%s", kEmLevelName[l]), ";corrected jet p_{T} (GeV);jets", 100, 0, 20);
        h_ptraw[l] = new TH1D(TString::Format("h_ptraw_%s", kEmLevelName[l]), ";uncorrected jet p_{T} (GeV);jets", 100, 0, 20);
    }
    // sanity distributions: [0] jets of events passing the event cuts, [1] selected jets
    const char *const stage[2] = {"evt", "sel"};
    TH1D *q_eta[2], *q_phi[2], *q_e[2], *q_ue[2], *q_ecorr[2], *q_pt[2], *q_ptcorr[2], *q_xf[2], *q_np[2], *q_njet[2];
    TH2D *q_xy[2];
    for (int g = 0; g < 2; g++)
    {
        const char *st = stage[g];
        q_eta[g] = new TH1D(TString::Format("qa_eta_%s", st), ";jet #eta;jets", 60, 2.0, 5.0);
        q_phi[g] = new TH1D(TString::Format("qa_phi_%s", st), ";jet #phi (reco_phi, 0 to 2#pi);jets", 72, 0, 2 * M_PI);
        q_e[g] = new TH1D(TString::Format("qa_E_%s", st), ";uncorrected jet E (GeV);jets", 125, 0, 125);
        q_ue[g] = new TH1D(TString::Format("qa_ueE_%s", st), ";UE energy in the jet (GeV);jets", 100, 0, 10);
        q_ecorr[g] = new TH1D(TString::Format("qa_Ecorr_%s", st), ";corrected jet E (GeV);jets", 125, 0, 125);
        q_pt[g] = new TH1D(TString::Format("qa_pt_%s", st), ";uncorrected jet p_{T} (GeV);jets", 100, 0, 20);
        q_ptcorr[g] = new TH1D(TString::Format("qa_ptcorr_%s", st), ";corrected jet p_{T} (GeV);jets", 100, 0, 20);
        q_xf[g] = new TH1D(TString::Format("qa_xF_%s", st), ";corrected x_{F};jets", kNXfBins, kXfMin, kXfMax);
        q_np[g] = new TH1D(TString::Format("qa_np_%s", st), ";points in the jet;jets", 20, 0, 20);
        q_njet[g] = new TH1D(TString::Format("qa_njet_%s", st), ";jets per event;events", 10, 0, 10);
        q_xy[g] = new TH2D(TString::Format("qa_xy_%s", st), ";jet x (cm);jet y (cm)", 60, -150, 150, 44, -110, 110);
    }
    TH1D *q_vz = new TH1D("qa_vz", ";event vertex z (cm), events with a vertex;events", 80, -80, 80);
    TH1D *q_vzsrc = new TH1D("qa_vz_source", ";vertex source;events", 4, 0, 4);
    const char *const vz_name[4] = {"none", "TPC", "VPD", "BBC"};
    for (int i = 0; i < 4; i++) q_vzsrc->GetXaxis()->SetBinLabel(i + 1, vz_name[i]);
    TH1D *q_trig = new TH1D("qa_em_listed", ";listed EM triggers;events", 8, 0, 8);
    for (int i = 0; i < 8; i++)
        q_trig->GetXaxis()->SetBinLabel(i + 1, i == 0 ? "none" : TString::Format("%s%s%s", i & 1 ? "0" : "", i & 2 ? "1" : "", i & 4 ? "2" : ""));
    TH1D *q_spin = new TH1D("qa_spin_config", ";Spin_config;events", 17, -1, 16);
    TH2D *q_xy_side[kNSides];
    for (int s = 0; s < kNSides; s++)
        q_xy_side[s] = new TH2D(TString::Format("qa_xy_sel_%s", kSideName[s]), ";jet x (cm);jet y (cm)", 60, -150, 150, 44, -110, 110);
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
        bool listed[kNEmLevels];
        for (int l = 0; l < kNEmLevels; l++) listed[l] = em[l] || em_tpc[l];
        const int level = listed[2] ? 2 : listed[1] ? 1 : listed[0] ? 0 : -1;
        const int spin = SpinState(spin_config);
        const size_t nj = eta->size();
        h_cut->Fill(kAllJets, nj);
        if (led) continue;
        h_cut->Fill(kNoLed, nj);
        q_trig->Fill((listed[0] ? 1 : 0) + (listed[1] ? 2 : 0) + (listed[2] ? 4 : 0));
        if (level < 0) continue;
        h_cut->Fill(kEmLevel, nj);
        q_spin->Fill(spin_config);
        if (spin < 0) continue;
        h_cut->Fill(kSpin, nj);
        h_events[level]->Fill(spin);
        q_vzsrc->Fill(vz_source >= 0 && vz_source < 4 ? vz_source : 0);
        if (vz > -9000) q_vz->Fill(vz);
        q_njet[0]->Fill(nj);
        int nsel = 0;
        for (size_t j = 0; j < nj; j++)
        {
            const double E = (*e_corr)[j];
            // sanity distributions of every jet here; the selected jets' below
            q_eta[0]->Fill((*eta)[j]); q_phi[0]->Fill((*phi)[j]); q_e[0]->Fill((*e_raw)[j]); q_ue[0]->Fill((*e_ue)[j]);
            q_ecorr[0]->Fill(E); q_pt[0]->Fill((*pt_raw)[j]); q_ptcorr[0]->Fill(E > 0 ? E / std::cosh((*eta)[j]) : -1); q_xf[0]->Fill((*xf_corr)[j]);
            q_np[0]->Fill((*nphotons)[j]); q_xy[0]->Fill((*x)[j], (*y)[j]);
            if (!(E > 0)) continue;
            h_cut->Fill(kCalibrated);
            if (!pass_physical_jet_cut(E, (*eta)[j], sqrt_s)) continue;
            h_cut->Fill(kPhysical);
            const double pt = E / std::cosh((*eta)[j]);
            if (pt <= kJetPtCorrMin) continue;
            h_cut->Fill(kPtMin);
            h_pt[level]->Fill(pt);
            h_ptraw[level]->Fill((*pt_raw)[j]);
            // the lowest listed EM level whose threshold the jet's uncorrected pT passes (FMS: any listed trigger)
            int jl = -1;
            for (int l = 0; l < kNEmLevels && jl < 0; l++)
                if (listed[l] && (*pt_raw)[j] > kTriggerPtFactor * kEmThresholdGeV[l]) jl = l;
            if (jl < 0) continue;
            h_cut->Fill(kTriggerPt);
            const int side = computeSide((*x)[j]) > 0 ? 1 : 0;
            const int cls = PhotonClass((*nphotons)[j]);
            h[jl][cls][side][spin]->Fill((*xf_corr)[j]);
            h_xfs[jl][cls]->Fill((*xf_corr)[j], (*xf_corr)[j]);
            h_cos[jl][cls]->Fill((*xf_corr)[j], std::fabs((*x)[j]) / std::hypot((*x)[j], (*y)[j]));
            nsel++;
            q_eta[1]->Fill((*eta)[j]); q_phi[1]->Fill((*phi)[j]); q_e[1]->Fill((*e_raw)[j]); q_ue[1]->Fill((*e_ue)[j]);
            q_ecorr[1]->Fill(E); q_pt[1]->Fill((*pt_raw)[j]); q_ptcorr[1]->Fill(pt); q_xf[1]->Fill((*xf_corr)[j]);
            q_np[1]->Fill((*nphotons)[j]); q_xy[1]->Fill((*x)[j], (*y)[j]); q_xy_side[side]->Fill((*x)[j], (*y)[j]);
        }
        q_njet[1]->Fill(nsel);
    }

    std::ostringstream rs;
    for (std::set<int>::const_iterator r = runs.begin(); r != runs.end(); ++r) rs << (r == runs.begin() ? "" : " ") << *r;
    TNamed("runs", rs.str().c_str()).Write();
    TNamed("spin_convention", kSpinConventionTag).Write();
    out.Write();
    printf("AnHistograms: %lld events, runs %s\n", n, rs.str().c_str());
    for (int c = 0; c < kNCuts; c++) printf("  %-14s %.0f jets\n", kCutName[c], h_cut->GetBinContent(c + 1));
    out.Close();
    return 0;
}
