// Clusters truth and reco jets from SimpleTree files and matches them.
// Usage: ./JetMatcher file_list.list [trigger_filter] [em_only] [use_points] [jet_pt_cuts]
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
        cerr << "Usage: " << argv[0] << " file_list.list [trigger_filter] [em_only] [use_points] [jet_pt_cuts]" << endl;
        cerr << "  jet_pt_cuts: 0 = no minimum pT on EM-only reco and truth jets (for a" << endl;
        cerr << "              jet-energy-scale calibration sample). 1/omit = apply them." << endl;
        cerr << "  trigger_filter: comma-separated FCS trigger flag names (see FcsTriggerDefs.h)." << endl;
        cerr << "                  An event is kept if ANY of them fired. Omit/empty = keep all events." << endl;
        cerr << "  em_only: 1 = form reco jets from ECAL hits only (no HCAL), and use the" << endl;
        cerr << "           ECAL-only fiducial boundary for both reco and truth jets." << endl;
        cerr << "           0/omit = use both ECAL+HCAL." << endl;
        cerr << "  use_points: 1 = build reco jets from FCS points (StMuFcsPoint, the" << endl;
        cerr << "              single-photon-candidate level) instead of raw tower hits." << endl;
        cerr << "              Truth jets are unaffected (they're generator-level already)." << endl;
        cerr << "              0/omit = hits. HCAL currently" << endl;
        cerr << "              produces no points at all, so this is ECAL-only input" << endl;
        cerr << "              regardless of em_only when enabled." << endl;
        return 1;
    }

    bool em_only = (argc >= 4) && (string(argv[3]) == "1");

    bool use_points = (argc >= 5) && (string(argv[4]) == "1");

    // The EM-only jet pT minimums are analysis cuts; a JES calibration sample
    // must not have them, since a reco threshold biases the response.
    bool jet_pt_cuts = !((argc >= 6) && (string(argv[5]) == "0"));
    const double reco_pt_min = (em_only && jet_pt_cuts) ? reco_jet_pt_min_em : 0.0;
    const double truth_pt_min = (em_only && jet_pt_cuts) ? truth_jet_pt_min_em : 0.0;
    cout << "Jet pT minimums: reco " << reco_pt_min << " GeV, truth " << truth_pt_min << " GeV" << endl;

    ifstream infile(argv[1]);
    vector<string> file_names;
    string line;
    while (getline(infile, line))
        if (!line.empty()) file_names.push_back(line);

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
    vector<int> reco_nphotons;
    vector<int> truth_nphotons;

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

    // FCS points with energy > point_energy_min within deltaR < R of the jet
    // axis (as data_to_jet/RecoJets.cpp)
    outTree->Branch("reco_nphotons", &reco_nphotons);

    // Photon constituents of the truth jet, with no energy cut, as the FMS
    // analysis (reco_nphotons needs 1 GeV; how to cut truth photons to match
    // the detector photons is an open question, as it was for FMS)
    outTree->Branch("truth_nphotons", &truth_nphotons);

    outTree->Branch("truth_const_px", &truth_const_px);
    outTree->Branch("truth_const_py", &truth_const_py);
    outTree->Branch("truth_const_pz", &truth_const_pz);
    outTree->Branch("truth_const_E", &truth_const_E);
    outTree->Branch("truth_const_jetIndex", &truth_const_jetIndex);
    outTree->Branch("truth_const_pid", &truth_const_pid);
    outTree->Branch("truth_nconst", &truth_nconst);

    // Per-event FCS trigger flags, passed through from the SimpleTree (FcsTriggerDefs.h)
    Int_t Trig_flag[kNTrigFlags];
    for (int i = 0; i < kNTrigFlags; i++){
        outTree->Branch(kTrigFlagName[i], &Trig_flag[i], TString::Format("%s/I", kTrigFlagName[i]));
    }

    // Always -1 in simulation; kept for schema consistency with data_to_jet
    Int_t Spin_config;
    outTree->Branch("Spin_config", &Spin_config, "Spin_config/I");

    // Per-job cross-section bookkeeping from the SimpleTree
    // (StSimpleReaderMaker::SetMCJobStats): PYTHIA's cross-section estimate
    // (mb) and number of trials. -1 if the input lacks the branches.
    Double_t sigma_job;
    Int_t N_job;
    outTree->Branch("sigma_job", &sigma_job, "sigma_job/D");
    outTree->Branch("N_job", &N_job, "N_job/I");

    // Generator ptHat range of the job, -1 if absent
    Float_t mc_pthatmin_gev;
    Float_t mc_pthatmax_gev;
    outTree->Branch("mc_pthatmin_gev", &mc_pthatmin_gev, "mc_pthatmin_gev/F");
    outTree->Branch("mc_pthatmax_gev", &mc_pthatmax_gev, "mc_pthatmax_gev/F");

    // Per-file settings: em_only, jet pT minimums used (0 = none), use_points
    Int_t reco_em_only = em_only ? 1 : 0;
    outTree->Branch("reco_em_only", &reco_em_only, "reco_em_only/I");
    Float_t reco_jet_pt_min_used = (Float_t)reco_pt_min, truth_jet_pt_min_used = (Float_t)truth_pt_min;
    outTree->Branch("reco_jet_pt_min", &reco_jet_pt_min_used, "reco_jet_pt_min/F");
    outTree->Branch("truth_jet_pt_min", &truth_jet_pt_min_used, "truth_jet_pt_min/F");

    Int_t reco_from_points = use_points ? 1 : 0;
    outTree->Branch("reco_from_points", &reco_from_points, "reco_from_points/I");

    // Entry of the source event in its SimpleTree "data" tree (used by the
    // underlying-event step to read the event back)
    Long64_t simpletree_entry = -1;
    outTree->Branch("simpletree_entry", &simpletree_entry, "simpletree_entry/L");

    TDatabasePDG* pdgDB = TDatabasePDG::Instance();

    // Species counted as electromagnetic in the EMF fractions
    set<int> EM_PIDS = {11, -11, 22};

    // em_only truth particles: photons only, as the FMS EM-jet particle level
    // (arXiv:2012.11428, Sec. II.D)
    set<int> EM_ONLY_TRUTH_PIDS = {22};

    const int MAX = 10000;

    Int_t mcpart_num, Cal_nhits;

    Float_t mcpart_px[MAX], mcpart_py[MAX], mcpart_pz[MAX], mcpart_E[MAX];
    Int_t mcpart_geid[MAX];
    Int_t mcpart_idVtx[MAX];
    Float_t mcpart_Vtx_x[MAX], mcpart_Vtx_y[MAX], mcpart_Vtx_z[MAX];

    // em_only truth: maximum displacement (cm) of a photon's creation point
    // from the primary vertex. 100 um keeps pi0/eta/Sigma0 decay photons
    // (decay lengths of microns) and excludes K0s/Lambda/K0L decays and GEANT
    // secondaries (cm and beyond). See JetParameters.h.
    const double DECAY_DISPLACEMENT_MAX = em_truth_displacement_max;

    Float_t Cal_hit_energy[MAX], Cal_hit_posx[MAX], Cal_hit_posy[MAX], Cal_hit_posz[MAX];
    Int_t Cal_detid[MAX];

    const int MAX_PTS = 200; // matches StSimpleReaderMaker.h's Cal_pt_* array size
    Int_t Cal_npts;
    Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
    Int_t Cal_pt_detid[MAX_PTS];

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
        tree->SetBranchAddress("mcpart_Vtx_x", mcpart_Vtx_x);
        tree->SetBranchAddress("mcpart_Vtx_y", mcpart_Vtx_y);
        tree->SetBranchAddress("mcpart_Vtx_z", mcpart_Vtx_z);

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

        mc_pthatmin_gev = -1;
        if ( tree->GetBranch("mc_pthatmin_gev") ) {
            tree->SetBranchAddress("mc_pthatmin_gev", &mc_pthatmin_gev);
        }

        mc_pthatmax_gev = -1;
        if ( tree->GetBranch("mc_pthatmax_gev") ) {
            tree->SetBranchAddress("mc_pthatmax_gev", &mc_pthatmax_gev);
        }

        // Files without points give Cal_npts = 0
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

        for (Long64_t ev = 0; ev < nentries; ev++) {

            if (tree->GetEntry(ev) <= 0) continue;
            simpletree_entry = ev;

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
            reco_nphotons.clear(); truth_nphotons.clear();

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

            // The collision vertex is smeared event by event, so displacement
            // is measured from this event's primary vertex, not the origin.
            float primVtx_x = 0, primVtx_y = 0, primVtx_z = 0;
            bool foundPrimVtx = false;
            if (em_only)
                foundPrimVtx = FindEmTruthPrimaryVertex(mcpart_num, mcpart_idVtx, mcpart_Vtx_x, mcpart_Vtx_y, mcpart_Vtx_z,
                                                        primVtx_x, primVtx_y, primVtx_z);

            for (int i = 0; i < mcpart_num; i++) {
                if (em_only) {
                    // Not idVtx == 1: every decay in flight (including
                    // pi0 -> gamma gamma) gets a new vertex ID, so that would
                    // drop the decay photons. A photon that later converts
                    // is counted with its energy at creation; at most one
                    // generation of a lineage lies within the displacement
                    // cut, so nothing is double counted.
                    if (!foundPrimVtx) continue;
                    if (!IsEmTruthParticle(mcpart_geid[i], mcpart_Vtx_x[i], mcpart_Vtx_y[i], mcpart_Vtx_z[i],
                                           primVtx_x, primVtx_y, primVtx_z)) continue;
                } else {
                    // ECAL+HCAL truth: particles from the first vertex only
                    // (this also drops pi0 decay photons)
                    if (mcpart_idVtx[i] != 1) continue;
                }

                int pdg = pdgDB->ConvertGeant3ToPdg(mcpart_geid[i]);

                // EM-only truth jets are clustered from photons only, within
                // the truth particle eta band (full azimuth)
                if (em_only && !EM_ONLY_TRUTH_PIDS.count(pdg)) continue;
                if (em_only)
                {
                    const double ptp = std::sqrt(mcpart_px[i] * mcpart_px[i] + mcpart_py[i] * mcpart_py[i]);
                    if (!(ptp > 0) || !InTruthParticleEtaRange(std::asinh(mcpart_pz[i] / ptp))) continue;
                }

                truth_particles.emplace_back(mcpart_px[i], mcpart_py[i], mcpart_pz[i], mcpart_E[i]);
                truth_pid_particles.push_back(pdg);
                truth_particles.back().set_user_index(truth_pid_particles.size() - 1);
            }

            // truth_pt_min is 0 for ECAL+HCAL jets or when jet_pt_cuts = 0
            ClusterSequence cs_truth(truth_particles, JetDefinition(antikt_algorithm, R));
            auto truth_all = cs_truth.inclusive_jets(truth_pt_min);

            vector<PseudoJet> truth_jets_selected;

            for (auto &jet : truth_all) {
                // Backward jets would project onto the forward plane spuriously
                if (jet.pz() <= 0) continue;

                // EM-only: the truth jet eta band (full azimuth); ECAL+HCAL:
                // the fiducial rectangle dilated by R/2
                float jetXE = z_proj * jet.px() / jet.pz();
                float jetYE = z_proj * jet.py() / jet.pz();
                bool passesFiducial = em_only
                    ? InTruthJetEtaRange(jet.eta())
                    : pass_jet_scale_cut(jetXE, jetYE, reco_fiducial_buffer, -R/2.0f);
                if (passesFiducial) truth_jets_selected.push_back(jet);
            }

            // ---------------- reco ----------------

            vector<PseudoJet> reco_particles;
            vector<int> reco_detid;

            // Tower hits (default) or FCS points (use_points) feed the same
            // clustering. Hits use the MIP-fraction threshold, points
            // pass_point_energy_cut().
            int nSource = use_points ? Cal_npts : Cal_nhits;

            for (int i = 0; i < nSource; i++) {

                float e = use_points ? Cal_pt_energy[i] : Cal_hit_energy[i];
                int det = use_points ? Cal_pt_detid[i] : Cal_detid[i];

                if (det == 4 || det == 5) continue;
                if (em_only && (det == 2 || det == 3)) continue;

                if (use_points) {
                    if (!pass_point_energy_cut(e, kSqrtS_MainRun22_GeV)) continue;
                } else if (!(((det == 0 || det == 1) && e > mip_threshold * ecal_mip) ||
                             ((det == 2 || det == 3) && e > mip_threshold * hcal_mip)))
                    continue;

                float x = use_points ? Cal_pt_x[i] : Cal_hit_posx[i];
                float y = use_points ? Cal_pt_y[i] : Cal_hit_posy[i];
                float z = use_points ? Cal_pt_z[i] : Cal_hit_posz[i];

                float norm = sqrt(x*x + y*y + z*z);

                reco_particles.emplace_back(e*x/norm, e*y/norm, e*z/norm, e);
                reco_detid.push_back(det);
                reco_particles.back().set_user_index(reco_detid.size() - 1);
            }

            ClusterSequence cs_reco(reco_particles, JetDefinition(antikt_algorithm, R));
            vector<PseudoJet> reco_jets_selected;

            // reco_pt_min is 0 for ECAL+HCAL jets or when jet_pt_cuts = 0
            for (auto &jet : cs_reco.inclusive_jets(reco_pt_min)) {
                float jetXE = z_proj * jet.px() / jet.pz();
                float jetYE = z_proj * jet.py() / jet.pz();
                bool passesFiducial = em_only
                    ? pass_fiducial_cut(jetXE, jetYE, reco_fiducial_buffer_ecal, kFiducialRectEcal)
                    : pass_fiducial_cut(jetXE, jetYE, reco_fiducial_buffer);
                if (passesFiducial) reco_jets_selected.push_back(jet);
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
                int nphotons_truth = 0;

                truth_pt.push_back(jet.perp());
                truth_eta.push_back(jet.eta());
                truth_phi.push_back(jet.phi());
                truth_E.push_back(jet.E());
                truth_nconst.push_back(consts.size());

                truth_x.push_back(z_proj*jet.px()/jet.pz());
                truth_y.push_back(z_proj*jet.py()/jet.pz());
                truth_x_F.push_back(computeFeynmanX(jet.E(), jet.eta(), kSqrtS_MainRun22_GeV));

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
                    if (pid == 22) nphotons_truth++;
                }
                truth_nphotons.push_back(nphotons_truth);

                // ~1 in em_only mode by construction
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
                reco_x_F.push_back(computeFeynmanX(jet.E(), jet.eta(), kSqrtS_MainRun22_GeV));
                reco_side.push_back(computeSide(rx));

                int nphotons = 0;
                for (int p = 0; p < Cal_npts; p++) {
                    if (Cal_pt_detid[p] != 0 && Cal_pt_detid[p] != 1) continue;
                    if (!pass_point_energy_cut(Cal_pt_energy[p], kSqrtS_MainRun22_GeV)) continue;
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
        }

        f->Close();
    }

    fout->Write();
    fout->Close();
    return 0;
}