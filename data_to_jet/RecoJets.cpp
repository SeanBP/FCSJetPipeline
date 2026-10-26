// Clusters reco jets from real-data SimpleTree files.
// Usage: ./RecoJets file_list.list [trigger_filter] [em_only] [use_points]
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
#include "fastjet/ClusterSequenceArea.hh"

using namespace std;
using namespace fastjet;


// ------------------ Main ------------------
int main(int argc, char** argv) {

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " file_list.list [trigger_filter] [em_only] [use_points]" << endl;
        cerr << "  trigger_filter: comma-separated FCS trigger flag names (see FcsTriggerDefs.h)." << endl;
        cerr << "                  An event is kept if ANY of them fired. Omit/empty = keep all events." << endl;
        cerr << "  em_only: 1 = form reco jets from ECAL hits only (no HCAL), and use the" << endl;
        cerr << "           ECAL-only fiducial boundary. 0/omit = use both ECAL+HCAL." << endl;
        cerr << "  use_points: 1 = build reco jets from FCS points (StMuFcsPoint, the" << endl;
        cerr << "              single-photon-candidate level) instead of raw tower hits." << endl;
        cerr << "              0/omit = hits. HCAL produces no points, so this is" << endl;
        cerr << "              ECAL-only input regardless of em_only when enabled." << endl;
        return 1;
    }

    bool em_only = (argc >= 4) && (string(argv[3]) == "1");

    bool use_points = (argc >= 5) && (string(argv[4]) == "1");

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
    vector<int> reco_nphotons;

    // Flat constituent storage
    vector<float> reco_const_px, reco_const_py, reco_const_pz, reco_const_E;
    vector<int>   reco_const_jetIndex, reco_const_detid, reco_nconst;
    vector<float> reco_area; // jet areas (JetAreaDefinition)

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
    outTree->Branch("reco_area", &reco_area);

    outTree->Branch("reco_tau1", &reco_tau1);
    outTree->Branch("reco_tau2", &reco_tau2);

    // ECAL FCS points with energy > point_energy_min within deltaR < R of the
    // jet axis (the FMS photon definition)
    outTree->Branch("reco_nphotons", &reco_nphotons);

    // Per-event FCS trigger flags, passed through from the SimpleTree (FcsTriggerDefs.h)
    Int_t Trig_flag[kNTrigFlags];
    for (int i = 0; i < kNTrigFlags; i++){
        outTree->Branch(kTrigFlagName[i], &Trig_flag[i], TString::Format("%s/I", kTrigFlagName[i]));
    }

    // Per-event spin config, passed through from the SimpleTree; -1 = no spin-DB coverage
    Int_t Spin_config;
    outTree->Branch("Spin_config", &Spin_config, "Spin_config/I");

    // Event z-vertex from the SimpleTree: TPC, then VPD, then BBC, the first
    // with |z| < 80 cm (-9999 = none). event_vertex_source: 0 none, 1 TPC,
    // 2 VPD, 3 BBC.
    Float_t event_vertex_z;
    Int_t event_vertex_source;
    outTree->Branch("event_vertex_z", &event_vertex_z, "event_vertex_z/F");
    outTree->Branch("event_vertex_source", &event_vertex_source, "event_vertex_source/I");

    // Per-event sqrt(s), GeV, used for reco_x_F
    Double_t event_sqrt_s_gev;
    outTree->Branch("event_sqrt_s_gev", &event_sqrt_s_gev, "event_sqrt_s_gev/D");

    // Always -1 in data; kept for schema parity with sim_to_jet
    Double_t sigma_job;
    Int_t N_job;
    outTree->Branch("sigma_job", &sigma_job, "sigma_job/D");
    outTree->Branch("N_job", &N_job, "N_job/I");

    // Per-file settings: em_only, use_points
    Int_t reco_em_only = em_only ? 1 : 0;
    outTree->Branch("reco_em_only", &reco_em_only, "reco_em_only/I");

    Int_t reco_from_points = use_points ? 1 : 0;
    outTree->Branch("reco_from_points", &reco_from_points, "reco_from_points/I");

    // Entry of the source event in its SimpleTree "data" tree (used by the
    // underlying-event step to read the event back)
    Long64_t simpletree_entry = -1;
    outTree->Branch("simpletree_entry", &simpletree_entry, "simpletree_entry/L");

    const int MAX_HITS = 10000;
    const int MAX_PTS = 200; // matches StSimpleReaderMaker.h's Cal_pt_* array size

    Int_t Cal_nhits;
    Float_t Cal_hit_energy[MAX_HITS], Cal_hit_posx[MAX_HITS], Cal_hit_posy[MAX_HITS], Cal_hit_posz[MAX_HITS];
    Int_t Cal_detid[MAX_HITS];

    Int_t Cal_npts;
    Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
    Int_t Cal_pt_detid[MAX_PTS];

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

        sigma_job = -1; N_job = -1;
        if ( tree->GetBranch("sigma_job") ) {
            tree->SetBranchAddress("sigma_job", &sigma_job);
            tree->SetBranchAddress("N_job", &N_job);
        }

        // Defaults for files without these branches
        event_vertex_z = -9999; event_vertex_source = 0;
        if ( tree->GetBranch("event_vertex_z") ) {
            tree->SetBranchAddress("event_vertex_z", &event_vertex_z);
            tree->SetBranchAddress("event_vertex_source", &event_vertex_source);
        }
        event_sqrt_s_gev = kSqrtS_MainRun22_GeV;
        if ( tree->GetBranch("event_sqrt_s_gev") ) {
            tree->SetBranchAddress("event_sqrt_s_gev", &event_sqrt_s_gev);
        }
        Cal_npts = 0;
        if ( tree->GetBranch("Cal_npts") ) {
            tree->SetBranchAddress("Cal_npts", &Cal_npts);
            tree->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
            tree->SetBranchAddress("Cal_pt_x", Cal_pt_x);
            tree->SetBranchAddress("Cal_pt_y", Cal_pt_y);
            tree->SetBranchAddress("Cal_pt_z", Cal_pt_z);
            tree->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);
        } else if (use_points) {
            cerr << "Warning: use_points=1 but file has no Cal_npts branch -- "
                 << "0 reco jets will be built for this file." << endl;
        }

        Long64_t nentries = tree->GetEntries();
        cout << "Processing file: " << fname << " with " << nentries << " events." << endl;

        for (Long64_t ev = 0; ev < nentries; ev++) {
            tree->GetEntry(ev);
            simpletree_entry = ev;
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
            reco_area.clear();
            reco_nphotons.clear();

            vector<PseudoJet> reco_particles;
            vector<int> reco_detid;

            // Tower hits (default) or FCS points (use_points) feed the same
            // clustering. Hits use the MIP-fraction threshold, points
            // pass_point_energy_cut().
            int nSource = use_points ? Cal_npts : Cal_nhits;

            for (int i=0; i<nSource; i++) {

                float e = use_points ? Cal_pt_energy[i] : Cal_hit_energy[i];
                int detid = use_points ? Cal_pt_detid[i] : Cal_detid[i];
                if (detid == 4 || detid == 5) continue;
                if (em_only && (detid == 2 || detid == 3)) continue;

                if (use_points) {
                    if (!pass_point_energy_cut(e, event_sqrt_s_gev)) continue;
                } else if (!(((detid==0||detid==1)&&e>mip_threshold*ecal_mip) ||
                             ((detid==2||detid==3)&&e>mip_threshold*hcal_mip))) continue;

                float x = use_points ? Cal_pt_x[i] : Cal_hit_posx[i];
                float y = use_points ? Cal_pt_y[i] : Cal_hit_posy[i];
                float z = use_points ? Cal_pt_z[i] : Cal_hit_posz[i];

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

            ClusterSequenceArea cs_reco(reco_particles, JetDefinition(antikt_algorithm, R), JetAreaDefinition());
            vector<PseudoJet> reco_jets_selected;

            // reco_jet_pt_min_em for EM-only reco jets; ECAL+HCAL jets have none
            for (auto &jet : cs_reco.inclusive_jets(em_only ? reco_jet_pt_min_em : 0.0)) {
                float jetXE = z_proj*jet.px()/jet.pz();
                float jetYE = z_proj*jet.py()/jet.pz();
                bool passesFiducial = em_only
                    ? pass_fiducial_cut(jetXE, jetYE, reco_fiducial_buffer_ecal, kFiducialRectEcal)
                    : pass_fiducial_cut(jetXE, jetYE, reco_fiducial_buffer);
                if (passesFiducial) reco_jets_selected.push_back(jet);
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
                reco_area.push_back(jet.area());

                float rx = z_proj*jet.px()/jet.pz();
                reco_x.push_back(rx);
                reco_y.push_back(z_proj*jet.py()/jet.pz());
                reco_x_F.push_back(computeFeynmanX(jet.E(), jet.eta(), event_sqrt_s_gev));
                reco_side.push_back(computeSide(rx));

                int nphotons = 0;
                for (int p = 0; p < Cal_npts; p++) {
                    if (Cal_pt_detid[p] != 0 && Cal_pt_detid[p] != 1) continue;
                    if (!pass_point_energy_cut(Cal_pt_energy[p], event_sqrt_s_gev)) continue;
                    double pEta, pPhi;
                    computeEtaPhi(Cal_pt_x[p], Cal_pt_y[p], pEta, pPhi);
                    double dPhi = fabs(pPhi - jet.phi());
                    if (dPhi > M_PI) dPhi = 2 * M_PI - dPhi;
                    double dEta = pEta - jet.eta();
                    if (sqrt(dEta * dEta + dPhi * dPhi) < R) nphotons++;
                }
                reco_nphotons.push_back(nphotons);

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