// JetMatcher_safe_v6_fixed.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <algorithm>
#include <sstream>
#include "JetParameters.h"
#include "FcsTriggerDefs.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TDatabasePDG.h"
#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;


// ------------------ Main ------------------

int main(int argc, char** argv) {

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " file_list.list [trigger_filter]" << endl;
        cerr << "  trigger_filter: comma-separated FCS trigger flag names (see FcsTriggerDefs.h)." << endl;
        cerr << "                  An event is kept if ANY of them fired. Omit/empty = keep all events." << endl;
        return 1;
    }

    ifstream infile(argv[1]);
    vector<string> file_names;
    string line;
    while (getline(infile, line))
        if (!line.empty()) file_names.push_back(line);

    // ------------------ Trigger filter ------------------
    // Comma-separated list of FcsTriggerDefs.h flag names; an event is
    // kept if any of them fired. Empty/omitted -> no filtering, keep
    // every event (unchanged behavior). Same convention as RecoJets.cpp
    // (data_to_jet's jet finder).
    vector<int> filterFlagIndices;
    if (argc >= 3 && string(argv[2]).size() > 0) {
        string filterArg = argv[2];
        stringstream ss(filterArg);
        string name;
        while (getline(ss, name, ',')) {
            if (name.empty()) continue;
            int idx = -1;
            for (int i = 0; i < kNTrigFlags; i++) {
                if (name == kTrigFlagName[i]) { idx = i; break; }
            }
            if (idx < 0) {
                cerr << "Error: unknown trigger filter name '" << name
                     << "' (not in FcsTriggerDefs.h)" << endl;
                return 1;
            }
            filterFlagIndices.push_back(idx);
        }
        cout << "Trigger filter active: " << filterFlagIndices.size() << " trigger(s) from '"
             << filterArg << "'" << endl;
    } else {
        cout << "No trigger filter: keeping all events" << endl;
    }

    TFile* fout = new TFile("jet_output.root", "RECREATE");
    TTree* outTree = new TTree("jetTree", "Jets");

    // ------------------ branches ------------------

    vector<float> truth_pt, truth_eta, truth_phi, truth_E;
    vector<float> reco_pt, reco_eta, reco_phi, reco_E;

    vector<float> truth_EMF, reco_EMF;

    vector<int> match_truth_idx, match_reco_idx;

    vector<float> truth_tau1, truth_tau2;
    vector<float> reco_tau1, reco_tau2;

    vector<float> reco_const_px, reco_const_py, reco_const_pz, reco_const_E;
    vector<int> reco_const_jetIndex, reco_const_detid, reco_nconst;

    vector<float> truth_const_px, truth_const_py, truth_const_pz, truth_const_E;
    vector<int> truth_const_jetIndex, truth_const_pid, truth_nconst;

    vector<float> truth_x, truth_y;
    vector<float> reco_x, reco_y;

    vector<float> truth_x_F, reco_x_F;
    vector<int> reco_side;

    outTree->Branch("truth_pt", &truth_pt);
    outTree->Branch("truth_eta", &truth_eta);
    outTree->Branch("truth_phi", &truth_phi);
    outTree->Branch("truth_E", &truth_E);
    outTree->Branch("truth_EMF", &truth_EMF);
    outTree->Branch("truth_x", &truth_x);
    outTree->Branch("truth_y", &truth_y);
    outTree->Branch("truth_x_F", &truth_x_F);

    outTree->Branch("reco_pt", &reco_pt);
    outTree->Branch("reco_eta", &reco_eta);
    outTree->Branch("reco_phi", &reco_phi);
    outTree->Branch("reco_E", &reco_E);
    outTree->Branch("reco_x", &reco_x);
    outTree->Branch("reco_y", &reco_y);
    outTree->Branch("reco_x_F", &reco_x_F);
    outTree->Branch("reco_side", &reco_side);

    outTree->Branch("reco_EMF", &reco_EMF);

    outTree->Branch("match_truth_idx", &match_truth_idx);
    outTree->Branch("match_reco_idx", &match_reco_idx);

    outTree->Branch("truth_tau1", &truth_tau1);
    outTree->Branch("truth_tau2", &truth_tau2);

    outTree->Branch("reco_tau1", &reco_tau1);
    outTree->Branch("reco_tau2", &reco_tau2);

    outTree->Branch("reco_const_px", &reco_const_px);
    outTree->Branch("reco_const_py", &reco_const_py);
    outTree->Branch("reco_const_pz", &reco_const_pz);
    outTree->Branch("reco_const_E", &reco_const_E);
    outTree->Branch("reco_const_jetIndex", &reco_const_jetIndex);
    outTree->Branch("reco_const_detid", &reco_const_detid);
    outTree->Branch("reco_nconst", &reco_nconst);

    outTree->Branch("truth_const_px", &truth_const_px);
    outTree->Branch("truth_const_py", &truth_const_py);
    outTree->Branch("truth_const_pz", &truth_const_pz);
    outTree->Branch("truth_const_E", &truth_const_E);
    outTree->Branch("truth_const_jetIndex", &truth_const_jetIndex);
    outTree->Branch("truth_const_pid", &truth_const_pid);
    outTree->Branch("truth_nconst", &truth_nconst);

    // Per-event FCS trigger flags, passed through from the SimpleTree
    // (see FcsTriggerDefs.h), same convention as RecoJets.cpp
    // (data_to_jet's jet finder) so both pipelines' JetTrees carry the
    // same branches.
    Int_t Trig_flag[kNTrigFlags];
    for (int i = 0; i < kNTrigFlags; i++){
        outTree->Branch(kTrigFlagName[i], &Trig_flag[i], TString::Format("%s/I", kTrigFlagName[i]));
    }

    // Per-event spin configuration, same branch as RecoJets.cpp/
    // data_to_jet for schema consistency -- always -1 ("no data") here,
    // since sim_to_jet's readMudst.C never calls StSimpleReaderMaker::
    // SetSpinDb() (simulated events have no real polarization pattern to
    // report; see StSimpleReaderMaker.h's Spin_config comment).
    Int_t Spin_config;
    outTree->Branch("Spin_config", &Spin_config, "Spin_config/I");

    TDatabasePDG* pdgDB = TDatabasePDG::Instance();

    // -------- editable EM PID list --------
    set<int> EM_PIDS = {11, -11, 22};

    const int MAX = 10000;

    Int_t mcpart_num, Cal_nhits;

    Float_t mcpart_px[MAX], mcpart_py[MAX], mcpart_pz[MAX], mcpart_E[MAX];
    Int_t mcpart_geid[MAX];
    Int_t mcpart_idVtx[MAX];

    Float_t Cal_hit_energy[MAX], Cal_hit_posx[MAX], Cal_hit_posy[MAX], Cal_hit_posz[MAX];
    Int_t Cal_detid[MAX];

    struct Candidate { int t, r; double dr; };

    for (auto &fname : file_names) {

        TFile* f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) continue;

        TTree* tree = (TTree*)f->Get("data");
        if (!tree) { f->Close(); continue; }

        tree->SetBranchAddress("mcpart_num", &mcpart_num);
        tree->SetBranchAddress("mcpart_px", mcpart_px);
        tree->SetBranchAddress("mcpart_py", mcpart_py);
        tree->SetBranchAddress("mcpart_pz", mcpart_pz);
        tree->SetBranchAddress("mcpart_E", mcpart_E);
        tree->SetBranchAddress("mcpart_geid", mcpart_geid);
        tree->SetBranchAddress("mcpart_idVtx", mcpart_idVtx);

        tree->SetBranchAddress("Cal_nhits", &Cal_nhits);
        tree->SetBranchAddress("Cal_hit_energy", Cal_hit_energy);
        tree->SetBranchAddress("Cal_hit_posx", Cal_hit_posx);
        tree->SetBranchAddress("Cal_hit_posy", Cal_hit_posy);
        tree->SetBranchAddress("Cal_hit_posz", Cal_hit_posz);
        tree->SetBranchAddress("Cal_detid", Cal_detid);

        for (int i = 0; i < kNTrigFlags; i++){
            tree->SetBranchAddress(kTrigFlagName[i], &Trig_flag[i]);
        }
        tree->SetBranchAddress("Spin_config", &Spin_config);

        Long64_t nentries = tree->GetEntries();

        for (Long64_t ev = 0; ev < nentries; ev++) {

            if (tree->GetEntry(ev) <= 0) continue;

            if (!filterFlagIndices.empty()) {
                bool passed = false;
                for (int idx : filterFlagIndices) {
                    if (Trig_flag[idx]) { passed = true; break; }
                }
                if (!passed) continue;
            }

            truth_pt.clear(); truth_eta.clear(); truth_phi.clear(); truth_E.clear();
            reco_pt.clear(); reco_eta.clear(); reco_phi.clear(); reco_E.clear();
            reco_x.clear(); reco_y.clear();
            truth_x_F.clear(); reco_x_F.clear(); reco_side.clear();

            truth_EMF.clear(); reco_EMF.clear();

            truth_tau1.clear(); truth_tau2.clear();
            reco_tau1.clear(); reco_tau2.clear();

            truth_nconst.clear(); reco_nconst.clear();

            reco_const_px.clear(); reco_const_py.clear(); reco_const_pz.clear(); reco_const_E.clear();
            reco_const_jetIndex.clear(); reco_const_detid.clear();

            truth_const_px.clear(); truth_const_py.clear(); truth_const_pz.clear(); truth_const_E.clear();
            truth_const_jetIndex.clear(); truth_const_pid.clear();
            truth_x.clear(); truth_y.clear();

            match_truth_idx.clear();
            match_reco_idx.clear();

            // ---------------- truth ----------------

            vector<PseudoJet> truth_particles;
            vector<int> truth_pid_particles;

            for (int i = 0; i < mcpart_num; i++) {
                if (mcpart_idVtx[i] != 1) continue;

                int pdg = pdgDB->ConvertGeant3ToPdg(mcpart_geid[i]);

                truth_particles.emplace_back(mcpart_px[i], mcpart_py[i], mcpart_pz[i], mcpart_E[i]);
                truth_pid_particles.push_back(pdg);
                truth_particles.back().set_user_index(truth_pid_particles.size() - 1);
            }

            ClusterSequence cs_truth(truth_particles, JetDefinition(antikt_algorithm, R));
            auto truth_all = cs_truth.inclusive_jets();

            vector<PseudoJet> truth_jets_selected;

            for (auto &jet : truth_all) {
                float jetXE = z_proj * jet.px() / jet.pz();
                float jetYE = z_proj * jet.py() / jet.pz();
                if (pass_jet_scale_cut(jetXE, jetYE, reco_fiducial_buffer, -R/2.0f)) truth_jets_selected.push_back(jet);
            }

            // ---------------- reco ----------------

            vector<PseudoJet> reco_particles;
            vector<int> reco_detid;

            for (int i = 0; i < Cal_nhits; i++) {

                float e = Cal_hit_energy[i];
                int det = Cal_detid[i];

                if (det == 4 || det == 5) continue;

                if (!(((det == 0 || det == 1) && e > mip_threshold * ecal_mip) ||
                      ((det == 2 || det == 3) && e > mip_threshold * hcal_mip)))
                    continue;

                float x = Cal_hit_posx[i];
                float y = Cal_hit_posy[i];
                float z = Cal_hit_posz[i];

                float norm = sqrt(x*x + y*y + z*z);

                reco_particles.emplace_back(e*x/norm, e*y/norm, e*z/norm, e);
                reco_detid.push_back(det);
                reco_particles.back().set_user_index(reco_detid.size() - 1);
            }

            ClusterSequence cs_reco(reco_particles, JetDefinition(antikt_algorithm, R));
            vector<PseudoJet> reco_jets_selected;

            for (auto &jet : cs_reco.inclusive_jets()) {
                float jetXE = z_proj * jet.px() / jet.pz();
                float jetYE = z_proj * jet.py() / jet.pz();
                if (pass_fiducial_cut(jetXE, jetYE, reco_fiducial_buffer)) reco_jets_selected.push_back(jet);
            }

            // ---------------- matching ----------------

            vector<Candidate> cands;

            for (int ir = 0; ir < (int)reco_jets_selected.size(); ir++) {
                for (int it = 0; it < (int)truth_jets_selected.size(); it++) {

                    double deta = reco_jets_selected[ir].eta() - truth_jets_selected[it].eta();
                    double dphi = reco_jets_selected[ir].phi() - truth_jets_selected[it].phi();

                    if (dphi > M_PI) dphi -= 2*M_PI;
                    if (dphi < -M_PI) dphi += 2*M_PI;

                    double dr = sqrt(deta*deta + dphi*dphi);

                    cands.push_back({it, ir, dr});
                }
            }

            sort(cands.begin(), cands.end(),
                 [](const Candidate &a, const Candidate &b){ return a.dr < b.dr; });

            set<int> used_t, used_r;

            for (auto &c : cands) {
                if (c.dr > R/2.) continue;
                if (used_t.count(c.t) || used_r.count(c.r)) continue;

                match_truth_idx.push_back(c.t);
                match_reco_idx.push_back(c.r);

                used_t.insert(c.t);
                used_r.insert(c.r);
            }

            // ---------------- fill truth jets ----------------

            for (int j = 0; j < (int)truth_jets_selected.size(); j++) {

                auto &jet = truth_jets_selected[j];
                auto consts = jet.constituents();

                double E_em = 0, E_tot = 0;

                truth_pt.push_back(jet.perp());
                truth_eta.push_back(jet.eta());
                truth_phi.push_back(jet.phi());
                truth_E.push_back(jet.E());
                truth_nconst.push_back(consts.size());

                truth_x.push_back(z_proj*jet.px()/jet.pz());
                truth_y.push_back(z_proj*jet.py()/jet.pz());
                truth_x_F.push_back(computeFeynmanX(jet.E(), jet.eta()));

                for (auto &c : consts) {
                    int idx = c.user_index();
                    int pid = truth_pid_particles[idx];

                    truth_const_px.push_back(c.px());
                    truth_const_py.push_back(c.py());
                    truth_const_pz.push_back(c.pz());
                    truth_const_E.push_back(c.E());
                    truth_const_jetIndex.push_back(j);
                    truth_const_pid.push_back(pid);

                    E_tot += c.E();
                    if (EM_PIDS.count(pid)) E_em += c.E();
                }

                double EMF = E_em / jet.E();
                truth_EMF.push_back(EMF);

                double t1 = computeTauN(consts, 1, R);
                double t2 = computeTauN(consts, 2, R);

                truth_tau1.push_back(t1);
                truth_tau2.push_back(t2);
            }

            // ---------------- fill reco jets ----------------

            for (int j = 0; j < (int)reco_jets_selected.size(); j++) {

                auto &jet = reco_jets_selected[j];
                auto consts = jet.constituents();

                double E_ecal = 0, E_hcal = 0;

                reco_pt.push_back(jet.perp());
                reco_eta.push_back(jet.eta());
                reco_phi.push_back(jet.phi());
                reco_E.push_back(jet.E());
                reco_nconst.push_back(consts.size());

                float rx = z_proj*jet.px()/jet.pz();
                reco_x.push_back(rx);
                reco_y.push_back(z_proj*jet.py()/jet.pz());
                reco_x_F.push_back(computeFeynmanX(jet.E(), jet.eta()));
                reco_side.push_back(computeSide(rx));

                for (auto &c : consts) {
                    int idx = c.user_index();
                    int det = reco_detid[idx];

                    reco_const_px.push_back(c.px());
                    reco_const_py.push_back(c.py());
                    reco_const_pz.push_back(c.pz());
                    reco_const_E.push_back(c.E());
                    reco_const_jetIndex.push_back(j);
                    reco_const_detid.push_back(det);

                    if (det == 0 || det == 1) E_ecal += c.E();
                    if (det == 2 || det == 3) E_hcal += c.E();
                }

                double EMF = E_ecal / jet.E();
                reco_EMF.push_back(EMF);

                double t1 = computeTauN(consts, 1, R);
                double t2 = computeTauN(consts, 2, R);

                reco_tau1.push_back(t1);
                reco_tau2.push_back(t2);
            }
            outTree->Fill();
        }

        f->Close();
    }

    fout->Write();
    fout->Close();
    return 0;
}