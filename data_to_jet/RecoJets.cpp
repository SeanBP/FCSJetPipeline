// RecoJets.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "JetParameters.h"
#include "FcsTriggerDefs.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
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

    string list_fname = argv[1];
    ifstream infile(list_fname);
    if (!infile) {
        cerr << "Error opening list file: " << list_fname << endl;
        return 1;
    }

    // Comma-separated FcsTriggerDefs.h flag names; kept if any fired, empty = keep all.
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

    vector<string> file_names;
    string line;
    while (getline(infile, line))
        if (!line.empty()) file_names.push_back(line);
    infile.close();

    // ------------------ Output ROOT file ------------------
    TFile* fout = new TFile("jet_output.root", "RECREATE");
    if (!fout || fout->IsZombie()) {
        cerr << "Error: could not create output file jet_output.root" << endl;
        return 1;
    }

    // Jet-level
    vector<float> reco_pt, reco_eta, reco_phi, reco_E;
    vector<float> reco_x, reco_y;
    vector<float> reco_x_F;
    vector<int> reco_side;
    vector<float> reco_EMF;
    vector<float> reco_tau1, reco_tau2;

    // Flat constituent storage
    vector<float> reco_const_px, reco_const_py, reco_const_pz, reco_const_E;
    vector<int>   reco_const_jetIndex, reco_const_detid, reco_nconst;

    TTree* outTree = new TTree("jetTree", "Reco jets");
    outTree->SetDirectory(fout);

    outTree->Branch("reco_pt", &reco_pt);
    outTree->Branch("reco_eta", &reco_eta);
    outTree->Branch("reco_phi", &reco_phi);
    outTree->Branch("reco_E", &reco_E);
    outTree->Branch("reco_x", &reco_x);
    outTree->Branch("reco_y", &reco_y);
    outTree->Branch("reco_x_F", &reco_x_F);
    outTree->Branch("reco_side", &reco_side);
    outTree->Branch("reco_hadfrac", &reco_EMF);
    outTree->Branch("reco_const_px", &reco_const_px);
    outTree->Branch("reco_const_py", &reco_const_py);
    outTree->Branch("reco_const_pz", &reco_const_pz);
    outTree->Branch("reco_const_E",  &reco_const_E);

    outTree->Branch("reco_const_jetIndex", &reco_const_jetIndex);
    outTree->Branch("reco_const_detid", &reco_const_detid); 
    outTree->Branch("reco_nconst", &reco_nconst);

    outTree->Branch("reco_tau1", &reco_tau1);
    outTree->Branch("reco_tau2", &reco_tau2);

    // Per-event FCS trigger flags, passed through from the SimpleTree (FcsTriggerDefs.h)
    Int_t Trig_flag[kNTrigFlags];
    for (int i = 0; i < kNTrigFlags; i++){
        outTree->Branch(kTrigFlagName[i], &Trig_flag[i], TString::Format("%s/I", kTrigFlagName[i]));
    }

    // Per-event spin config, passed through from the SimpleTree; -1 = no spin-DB coverage
    Int_t Spin_config;
    outTree->Branch("Spin_config", &Spin_config, "Spin_config/I");

    // MC cross-section weighting branches, passed through from the SimpleTree
    // for schema parity with sim_to_jet's jetTree (see JetMatcher.cpp /
    // StSimpleReaderMaker::SetMCXSec()). Always the -1 sentinel here --
    // data_to_jet's runMudst.C never calls SetMCXSec(), real data has no
    // generator cross section.
    Float_t mc_sigma_pb, mc_sigma_err_pb;
    Int_t mc_n_gen;
    outTree->Branch("mc_sigma_pb", &mc_sigma_pb, "mc_sigma_pb/F");
    outTree->Branch("mc_sigma_err_pb", &mc_sigma_err_pb, "mc_sigma_err_pb/F");
    outTree->Branch("mc_n_gen", &mc_n_gen, "mc_n_gen/I");

    const int MAX_HITS = 10000;

    Int_t Cal_nhits;
    Float_t Cal_hit_energy[MAX_HITS], Cal_hit_posx[MAX_HITS], Cal_hit_posy[MAX_HITS], Cal_hit_posz[MAX_HITS];
    Int_t Cal_detid[MAX_HITS];

    Long64_t cumulative_event = 0;

    for (auto &fname : file_names) {

        TFile* f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            cerr << "Error opening file " << fname << ", skipping..." << endl;
            continue;
        }

        TTree* tree = (TTree*)f->Get("data");
        if (!tree) {
            cerr << "Error: TTree 'data' not found in " << fname << endl;
            f->Close();
            continue;
        }

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

        mc_sigma_pb = -1; mc_sigma_err_pb = -1; mc_n_gen = -1;
        if ( tree->GetBranch("mc_sigma_pb") ) {
            tree->SetBranchAddress("mc_sigma_pb", &mc_sigma_pb);
            tree->SetBranchAddress("mc_sigma_err_pb", &mc_sigma_err_pb);
            tree->SetBranchAddress("mc_n_gen", &mc_n_gen);
        }

        Long64_t nentries = tree->GetEntries();
        cout << "Processing file: " << fname << " with " << nentries << " events." << endl;

        for (Long64_t ev = 0; ev < nentries; ev++) {
            tree->GetEntry(ev);
            if (Cal_nhits > MAX_HITS) continue;

            if (!filterFlagIndices.empty()) {
                bool passed = false;
                for (int idx : filterFlagIndices) {
                    if (Trig_flag[idx]) { passed = true; break; }
                }
                if (!passed) continue;
            }

            reco_const_px.clear(); 
            reco_const_py.clear(); 
            reco_const_pz.clear(); 
            reco_const_E.clear();
            reco_const_jetIndex.clear(); 
            reco_const_detid.clear();
            reco_pt.clear(); 
            reco_eta.clear(); 
            reco_phi.clear(); 
            reco_E.clear();
            reco_x.clear();
            reco_y.clear();
            reco_x_F.clear();
            reco_side.clear();
            reco_EMF.clear();
            reco_tau1.clear(); 
            reco_tau2.clear();
            reco_nconst.clear();

            vector<PseudoJet> reco_particles;
            vector<int> reco_detid;

            for (int i=0; i<Cal_nhits; i++) {

                float e = Cal_hit_energy[i];
                int detid = Cal_detid[i];
                if (detid == 4 || detid == 5) continue;

                if (!(((detid==0||detid==1)&&e>mip_threshold*ecal_mip) ||
                      ((detid==2||detid==3)&&e>mip_threshold*hcal_mip))) continue;

                float x=Cal_hit_posx[i];
                float y=Cal_hit_posy[i];
                float z=Cal_hit_posz[i];

                float norm = sqrt(x*x+y*y+z*z);

                reco_particles.emplace_back(
                    e*x/norm,
                    e*y/norm,
                    e*z/norm,
                    e
                );
                reco_detid.push_back(detid);
                reco_particles.back().set_user_index(reco_detid.size() - 1);
            }

            ClusterSequence cs_reco(reco_particles, JetDefinition(antikt_algorithm, R));
            vector<PseudoJet> reco_jets_selected;

            for (auto &jet : cs_reco.inclusive_jets()) {
                float jetXE = z_proj*jet.px()/jet.pz();
                float jetYE = z_proj*jet.py()/jet.pz();
                if (pass_fiducial_cut(jetXE, jetYE, reco_fiducial_buffer)) reco_jets_selected.push_back(jet);
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
            cumulative_event++;

            if (cumulative_event % 50000 == 0) {
                outTree->AutoSave("SaveSelf");
                cout << "Auto-saved at event " << cumulative_event << endl;
            }
        }

        f->Close();
    }

    fout->Write();
    fout->Close();

    cout << "Finished processing all files. Total events: " << cumulative_event << endl;
    return 0;
}