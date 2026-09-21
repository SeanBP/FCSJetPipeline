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
        cerr << "Usage: " << argv[0] << " file_list.list [trigger_filter] [em_only]" << endl;
        cerr << "  trigger_filter: comma-separated FCS trigger flag names (see FcsTriggerDefs.h)." << endl;
        cerr << "                  An event is kept if ANY of them fired. Omit/empty = keep all events." << endl;
        cerr << "  em_only: 1 = form reco jets from ECAL hits only (no HCAL), and use the" << endl;
        cerr << "           ECAL-only fiducial boundary for both reco and truth jets." << endl;
        cerr << "           0/omit = use both ECAL+HCAL, the previous default behavior." << endl;
        return 1;
    }

    // Optional 3rd arg: 1 = EM-only reco jets (ECAL hits only, ECAL-only
    // fiducial boundary for reco and truth jets alike). Default 0 keeps
    // the previous ECAL+HCAL behavior unchanged.
    bool em_only = (argc >= 4) && (string(argv[3]) == "1");

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

    // Per-event FCS trigger flags, passed through from the SimpleTree (FcsTriggerDefs.h)
    Int_t Trig_flag[kNTrigFlags];
    for (int i = 0; i < kNTrigFlags; i++){
        outTree->Branch(kTrigFlagName[i], &Trig_flag[i], TString::Format("%s/I", kTrigFlagName[i]));
    }

    // Schema consistency with data_to_jet; always -1 here (sim has no real spin data)
    Int_t Spin_config;
    outTree->Branch("Spin_config", &Spin_config, "Spin_config/I");

    // Per-event MC job-level weighting, passed through from the SimpleTree
    // (see StSimpleReaderMaker::SetMCJobStats()): sigma_job is PYTHIA's raw
    // final cross-section estimate (mb, unconverted) and N_job its true
    // total generated trial count, for the job that produced this event.
    // Constant for the whole job/file, same convention as mc_pthatmin_gev
    // below -- bound directly to each input file's own branch in the loop,
    // so a jetTree built from several SimpleTree files (different ptHatMin
    // bins/generation jobs) still carries the right value per event. -1
    // sentinel if a source SimpleTree predates this bookkeeping.
    Double_t sigma_job;
    Int_t N_job;
    outTree->Branch("sigma_job", &sigma_job, "sigma_job/D");
    outTree->Branch("N_job", &N_job, "N_job/I");

    // Generator-level ptHatMin/ptHatMax cut for this job (see
    // StSimpleReaderMaker::SetMCPtHatMin()/SetMCPtHatMax()), same
    // passthrough/-1-sentinel convention as sigma_job above.
    Float_t mc_pthatmin_gev;
    Float_t mc_pthatmax_gev;
    outTree->Branch("mc_pthatmin_gev", &mc_pthatmin_gev, "mc_pthatmin_gev/F");
    outTree->Branch("mc_pthatmax_gev", &mc_pthatmax_gev, "mc_pthatmax_gev/F");

    // 1 if this file's reco jets were built from ECAL hits only (no HCAL,
    // ECAL-only fiducial boundary), 0 for the previous ECAL+HCAL default.
    // Constant for the whole job/file -- see the em_only CLI arg above.
    Int_t reco_em_only = em_only ? 1 : 0;
    outTree->Branch("reco_em_only", &reco_em_only, "reco_em_only/I");

    TDatabasePDG* pdgDB = TDatabasePDG::Instance();

    // -------- editable EM PID list --------
    // Used for the (ECAL+HCAL) truth/reco EMF fraction computation below,
    // the standard calorimetry definition of "electromagnetic" energy
    // (photons + e+/e-).
    set<int> EM_PIDS = {11, -11, 22};

    // em_only truth-level pre-clustering filter (see its use below):
    // photons only, matching STAR's own FMS EM-jet particle-level
    // definition verbatim (arXiv:2012.11428, Sec. II.D): "We define the
    // 'particle level' as the stable particles (photons here) produced in
    // a proton-proton event in PYTHIA prior to the GEANT simulation of
    // detector responses." Electrons/positrons are never mentioned in
    // that paper's truth-jet definition -- deliberately narrower than the
    // EM_PIDS used for the EMF fraction above.
    set<int> EM_ONLY_TRUTH_PIDS = {22};

    const int MAX = 10000;

    Int_t mcpart_num, Cal_nhits;

    Float_t mcpart_px[MAX], mcpart_py[MAX], mcpart_pz[MAX], mcpart_E[MAX];
    Int_t mcpart_geid[MAX];
    Int_t mcpart_idVtx[MAX];
    Float_t mcpart_Vtx_x[MAX], mcpart_Vtx_y[MAX], mcpart_Vtx_z[MAX];

    // em_only truth filter: max displacement (cm) of a particle's creation
    // vertex from the primary IP (fixed at (0,0,0), see
    // starsim_pythia8_filter.C's _primary->SetVertex(0,0,0)) for it to
    // still count as "prompt" -- see the em_only branch below for why this
    // exists instead of the plain idVtx==1 cut used for the default mode.
    //
    // Value chosen from the actual decay-length hierarchy, not a round
    // guess: the resonances we actually want to catch here (pi0, eta,
    // Sigma0 -- i.e. the ones producing a jet's genuine prompt EM content)
    // have c*tau of 25 nm, 0.15 nm, and 22 pm respectively, so even a
    // boosted decay length (gamma*beta*c*tau) stays at the few-micron
    // level for the highest constituent energies this analysis reaches
    // (~40 GeV/c => gamma ~ a few hundred for pi0 => ~microns). Long-lived
    // species that must stay excluded (K_short c*tau=2.68 cm, Lambda
    // c*tau=7.89 cm, K_long c*tau=15.5 m) sit 6+ orders of magnitude
    // higher and need momentum p < ~2 MeV/c (i.e. carrying negligible
    // energy, since decay length is proportional to p) to spuriously
    // decay within this cut -- so the residual contamination from those
    // is both rare and bounded to a vanishingly small energy contribution.
    // 100 microns (0.01 cm) sits comfortably above the prompt-decay scale
    // and comfortably below where K_short/Lambda/K_long become an issue.
    const double DECAY_DISPLACEMENT_MAX = 0.01;

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

        // -1 sentinel by default; overridden below only if this particular
        // input file actually has the branch (older SimpleTree files won't).
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

            // em_only's vertex-displacement cut below needs the position
            // of THIS event's own primary vertex, not the coordinate
            // origin -- StarPrimaryMaker smears the actual collision point
            // event-by-event (confirmed directly: one test event's idVtx==1
            // particles all shared vtx=(0.105,-0.148,7.053), ~7 cm from
            // the origin along z, well within RHIC's normal luminous-region
            // spread) despite starsim_pythia8_filter.C's
            // _primary->SetVertex(0,0,0) call, which only sets a nominal
            // reference/mean, not a fixed per-event position. Measuring
            // displacement from the origin instead of this position would
            // reject essentially every particle in the event, prompt or
            // not, once the smearing exceeds DECAY_DISPLACEMENT_MAX.
            float primVtx_x = 0, primVtx_y = 0, primVtx_z = 0;
            bool foundPrimVtx = false;
            if (em_only) {
                for (int i = 0; i < mcpart_num; i++) {
                    if (mcpart_idVtx[i] == 1) {
                        primVtx_x = mcpart_Vtx_x[i];
                        primVtx_y = mcpart_Vtx_y[i];
                        primVtx_z = mcpart_Vtx_z[i];
                        foundPrimVtx = true;
                        break;
                    }
                }
            }

            for (int i = 0; i < mcpart_num; i++) {
                if (em_only) {
                    // idVtx==1 (StMcTrack::IdVx(), the StMcVertex ID of a
                    // particle's creation point) only tags particles made
                    // at the very first vertex; STAR/GEANT hands out a new
                    // vertex ID for every decay-in-flight, including
                    // pi0/eta -> gamma gamma, whose decay length is
                    // microns even at high energy. A plain idVtx==1 cut
                    // therefore silently drops the dominant EM content of
                    // a real jet (pi0-decay photons) from the truth
                    // object, even though a real ECAL sees that energy
                    // fully -- this is what was producing reco_E > truth_E
                    // for low-truth_E matched jets (single stray
                    // primary-vertex photons loosely matched to much
                    // richer real reco jets), not genuine EM leakage.
                    //
                    // Use a vertex-displacement cut instead: accept a
                    // particle if it was created within
                    // DECAY_DISPLACEMENT_MAX of THIS event's own primary
                    // vertex (see above), which keeps prompt hadronic decay
                    // products (pi0/eta/Sigma0, all sub-mm decay lengths)
                    // while still excluding genuine GEANT detector-material
                    // shower secondaries created at real detector radii
                    // (cm to m scale). Deliberately NOT requiring
                    // idVtxEnd==0 (i.e. not requiring the particle be
                    // "final/stable"): a prompt decay photon that later
                    // pair-converts in real material (confirmed in test
                    // data: a 21.6 GeV prompt photon converting to e+e- at
                    // 346 cm from the primary vertex) is still counted here
                    // using its energy at creation, which is the correct
                    // truth-level value -- requiring stability would drop
                    // that energy from truth entirely (neither the
                    // pre-conversion photon nor its far-displaced daughters
                    // would pass), reintroducing a reco_E > truth_E bias
                    // for any jet with an early conversion, which is not
                    // rare (photon conversion probability before reaching
                    // a calorimeter is typically ~10-20% given realistic
                    // tracker/beampipe material budgets). This can't double
                    // count: a real decay/interaction chain only moves
                    // farther from the primary vertex at each successive
                    // step (nonzero flight distance is required to reach
                    // the next vertex), so at most one generation along any
                    // lineage can ever sit within a tight micron-scale
                    // radius of it -- confirmed in the same test data
                    // (this photon's own creation vertex sits ~10 microns
                    // from the primary vertex, matching pi0's real decay
                    // length, while its conversion vertex is 346 cm away).
                    if (!foundPrimVtx) continue;
                    double dx = mcpart_Vtx_x[i] - primVtx_x;
                    double dy = mcpart_Vtx_y[i] - primVtx_y;
                    double dz = mcpart_Vtx_z[i] - primVtx_z;
                    double vtx_r = sqrt(dx*dx + dy*dy + dz*dz);
                    if (vtx_r > DECAY_DISPLACEMENT_MAX) continue;
                } else {
                    // ECAL+HCAL (default) truth-jet definition: unchanged.
                    // Note this same idVtx==1 limitation (excludes prompt
                    // hadron-decay photons like pi0->gg) likely also
                    // affects this mode's truth jets, just proportionally
                    // less since photons are only one contribution among
                    // many hadronic constituents here -- not yet revisited.
                    if (mcpart_idVtx[i] != 1) continue;
                }

                int pdg = pdgDB->ConvertGeant3ToPdg(mcpart_geid[i]);

                // EM-only mode: restrict truth-level clustering input to
                // photons before jet-finding (see EM_ONLY_TRUTH_PIDS above
                // for the literature citation), rather than clustering all
                // truth particles and computing an EM fraction after the
                // fact. This matches STAR's own FMS "EM-jet" truth
                // definition: the truth-level object is built only from
                // the same particle species the EM-only detector can ever
                // see, so reco/truth are the same kind of object by
                // construction instead of a full hadronic jet the
                // detector has no hope of fully capturing.
                if (em_only && !EM_ONLY_TRUTH_PIDS.count(pdg)) continue;

                truth_particles.emplace_back(mcpart_px[i], mcpart_py[i], mcpart_pz[i], mcpart_E[i]);
                truth_pid_particles.push_back(pdg);
                truth_particles.back().set_user_index(truth_pid_particles.size() - 1);
            }

            ClusterSequence cs_truth(truth_particles, JetDefinition(antikt_algorithm, R));
            auto truth_all = cs_truth.inclusive_jets();

            vector<PseudoJet> truth_jets_selected;

            for (auto &jet : truth_all) {
                // Truth isn't eta-restricted before clustering, so a backward-going
                // (pz<=0) jet's forward-z projection below would be meaningless but
                // could spuriously pass the fiducial cut. Reco jets don't need this
                // guard (calorimeter hits are always at fixed forward z).
                if (jet.pz() <= 0) continue;

                float jetXE = z_proj * jet.px() / jet.pz();
                float jetYE = z_proj * jet.py() / jet.pz();
                bool passesFiducial = em_only
                    ? pass_jet_scale_cut(jetXE, jetYE, reco_fiducial_buffer_ecal, -R/2.0f, kFiducialRectEcal)
                    : pass_jet_scale_cut(jetXE, jetYE, reco_fiducial_buffer, -R/2.0f);
                if (passesFiducial) truth_jets_selected.push_back(jet);
            }

            // ---------------- reco ----------------

            vector<PseudoJet> reco_particles;
            vector<int> reco_detid;

            for (int i = 0; i < Cal_nhits; i++) {

                float e = Cal_hit_energy[i];
                int det = Cal_detid[i];

                if (det == 4 || det == 5) continue;
                if (em_only && (det == 2 || det == 3)) continue;

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

                // In em_only mode every truth constituent is already EM by
                // construction (see the truth-particle filter above), so
                // this is trivially ~1.0 there -- same situation as
                // reco_EMF below, kept for schema consistency across modes.
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