// Truth-labelled point jets for the per-jet pi0 probability (shared/Pi0Prob.h).
//
// Builds the EM-only reco jets from FCS points exactly as
// sim_to_jet/JetMatcher.cpp does (em_only = 1, use_points = 1, jet pT cuts on:
// pass_point_energy_cut, anti-kt R, ECal fiducial rectangle with its 20 cm
// margin, pT > reco_jet_pt_min_em) and labels every ECal point with the truth
// photon that made it:
//   - candidate photons: the EM-only truth photons (IsEmTruthParticle: prompt
//     and pi0/eta-decay photons, created within 100 um of the primary vertex)
//     projected from their creation point along their momentum to the point's z;
//   - points and photons are paired greedily by that transverse distance,
//     closest first, each used once, up to kMatchMaxCm;
//   - a photon's parent is the particle whose end vertex is the photon's
//     creation vertex (mcpart_idVtxEnd == mcpart_idVtx); pi0 (Geant 7) and
//     eta (17) decay in GEANT, so their photons are labelled this way.
// A pi0 is "resolved" in a jet when both its photons matched points of that
// jet (the signal pairs of the method). A pi0 with one photon matched to a jet
// point and the other not is "merged" when the other photon projects within
// kMergeCm of a point of the jet (its energy went into that point), else
// "lost" (below the 1 GeV point floor, outside the jet, or not reconstructed).
//
// Output tree "jets", one entry per reco jet:
//   file_idx, ev, pthat_min, pthat_max, vz_true,
//   jet_E, jet_pt, jet_eta, jet_phi, jet_xF (raw, no UE/JES),
//   n (ECal points in the jet), pt_E/x/y/z[n] (point energy and position),
//   pt_pi0[n]   index of the parent pi0 in this jet's pi0 list, -1 if none,
//   pt_par[n]   Geant id of the matched photon's parent (0 prompt / unknown,
//               -1 point not matched to a photon),
//   pt_dist[n]  match distance, cm (-1 unmatched),
//   pt_gE[n]    matched photon energy (-1 unmatched),
//   pt_ngeid[n], pt_ndist[n]  unmatched points only: Geant id and distance (cm)
//               of the nearest particle above 0.3 GeV created before the ECal
//               (any species, GEANT secondaries included; 0 / -1 otherwise),
//   k_res, n_merged, n_lost, n_eta_res (resolved eta pairs).
// Tree "files": file_idx, path, and the SimpleTree's sigma_job, N_job,
// pthat_min/max (first entry), so the fitter can cross-section weight without
// the SimpleTrees (Pi0DumpWeights.h); sigma_job = -1 if the file had none.
//
// Build:
//   g++ -O2 pi0_jet_dump.cpp -o pi0_jet_dump -I../../shared \
//       `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
// Run:  ./pi0_jet_dump <simpletree.list> <out.root> [first_file_idx] [matchcm] [recorded path]
//   (recorded path: stored as "path" for a one-file list, e.g. the name the
//   SimpleTree gets in the sim_to_jet job, where it is read as SimpleTree_mudst.root)
#include "JetParameters.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "fastjet/ClusterSequence.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;
using namespace fastjet;

static double kMatchMaxCm = 4.0;  // ~0.7 ECal tower (argv[4] overrides, for checks)
static const double kMergeCm = 6.0;     // ~1 ECal tower
static const int kMaxPt = 200, kMaxMc = 10000;

int main(int argc, char **argv)
{
    if (argc < 3) { cerr << "Usage: " << argv[0] << " <simpletree.list> <out.root> [first_file_idx]" << endl; return 1; }
    vector<string> files;
    { ifstream in(argv[1]); string l; while (getline(in, l)) if (!l.empty()) files.push_back(l); }
    const int idx0 = argc > 3 ? atoi(argv[3]) : 0;
    if (argc > 4 && atof(argv[4]) > 0) kMatchMaxCm = atof(argv[4]);
    const string recPath = argc > 5 ? argv[5] : "";

    TFile fout(argv[2], "RECREATE");
    TTree *tj = new TTree("jets", "truth-labelled point jets");
    TTree *tf = new TTree("files", "input files");
    Int_t file_idx, ev; Float_t pthat_min, pthat_max, vz_true;
    Float_t jet_E, jet_pt, jet_eta, jet_phi, jet_xF;
    Int_t n, k_res, n_merged, n_lost, n_eta_res;
    Float_t pE[kMaxPt], px_[kMaxPt], py_[kMaxPt], pz_[kMaxPt], pdist[kMaxPt], pgE[kMaxPt];
    Int_t ppi0[kMaxPt], ppar[kMaxPt], pngeid[kMaxPt]; Float_t pndist[kMaxPt];
    tj->Branch("file_idx", &file_idx, "file_idx/I"); tj->Branch("ev", &ev, "ev/I");
    tj->Branch("pthat_min", &pthat_min, "pthat_min/F"); tj->Branch("pthat_max", &pthat_max, "pthat_max/F");
    tj->Branch("vz_true", &vz_true, "vz_true/F");
    tj->Branch("jet_E", &jet_E, "jet_E/F"); tj->Branch("jet_pt", &jet_pt, "jet_pt/F");
    tj->Branch("jet_eta", &jet_eta, "jet_eta/F"); tj->Branch("jet_phi", &jet_phi, "jet_phi/F");
    tj->Branch("jet_xF", &jet_xF, "jet_xF/F");
    tj->Branch("n", &n, "n/I");
    tj->Branch("pt_E", pE, "pt_E[n]/F"); tj->Branch("pt_x", px_, "pt_x[n]/F");
    tj->Branch("pt_y", py_, "pt_y[n]/F"); tj->Branch("pt_z", pz_, "pt_z[n]/F");
    tj->Branch("pt_pi0", ppi0, "pt_pi0[n]/I"); tj->Branch("pt_par", ppar, "pt_par[n]/I");
    tj->Branch("pt_dist", pdist, "pt_dist[n]/F"); tj->Branch("pt_gE", pgE, "pt_gE[n]/F");
    tj->Branch("pt_ngeid", pngeid, "pt_ngeid[n]/I"); tj->Branch("pt_ndist", pndist, "pt_ndist[n]/F");
    tj->Branch("k_res", &k_res, "k_res/I"); tj->Branch("n_merged", &n_merged, "n_merged/I");
    tj->Branch("n_lost", &n_lost, "n_lost/I"); tj->Branch("n_eta_res", &n_eta_res, "n_eta_res/I");
    TString path;
    Double_t f_sigma; Int_t f_N; Float_t f_tmin, f_tmax;
    tf->Branch("file_idx", &file_idx, "file_idx/I"); tf->Branch("path", &path);
    tf->Branch("sigma_job", &f_sigma, "sigma_job/D"); tf->Branch("N_job", &f_N, "N_job/I");
    tf->Branch("pthat_min", &f_tmin, "pthat_min/F"); tf->Branch("pthat_max", &f_tmax, "pthat_max/F");

    static Int_t mc_n, mc_geid[kMaxMc], mc_vtx[kMaxMc], mc_vtxEnd[kMaxMc];
    static Float_t mc_px[kMaxMc], mc_py[kMaxMc], mc_pz[kMaxMc], mc_E[kMaxMc], mc_vx[kMaxMc], mc_vy[kMaxMc], mc_vz[kMaxMc];
    Int_t npts, pdet[kMaxPt]; Float_t pe[kMaxPt], pxx[kMaxPt], pyy[kMaxPt], pzz[kMaxPt];
    Float_t tmin, tmax;
    long nPts = 0, nPtsMatched = 0, nJets = 0;

    for (size_t fi = 0; fi < files.size(); fi++) {
        file_idx = idx0 + (int)fi;
        path = (files.size() == 1 && !recPath.empty()) ? recPath.c_str() : files[fi].c_str();
        f_sigma = -1; f_N = -1; f_tmin = f_tmax = -1;
        TFile *f = TFile::Open(files[fi].c_str());
        if (!f || f->IsZombie()) { cerr << "skip " << files[fi] << endl; tf->Fill(); continue; }
        TTree *t = (TTree *)f->Get("data");
        if (t && t->GetEntries() > 0 && t->GetBranch("sigma_job") && t->GetBranch("N_job"))
        {
            t->SetBranchStatus("*", 0);
            for (auto b : {"sigma_job", "N_job", "mc_pthatmin_gev", "mc_pthatmax_gev"}) if (t->GetBranch(b)) t->SetBranchStatus(b, 1);
            t->SetBranchAddress("sigma_job", &f_sigma); t->SetBranchAddress("N_job", &f_N);
            if (t->GetBranch("mc_pthatmin_gev")) t->SetBranchAddress("mc_pthatmin_gev", &f_tmin);
            if (t->GetBranch("mc_pthatmax_gev")) t->SetBranchAddress("mc_pthatmax_gev", &f_tmax);
            t->GetEntry(0);
            t->ResetBranchAddresses();
        }
        tf->Fill();
        if (!t || !t->GetBranch("Cal_npts") || !t->GetBranch("mcpart_idVtxEnd")) { cerr << "no tree/branches " << files[fi] << endl; f->Close(); continue; }
        t->SetBranchStatus("*", 0);
        const char *br[] = {"mcpart_num", "mcpart_geid", "mcpart_idVtx", "mcpart_idVtxEnd", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E",
                            "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z", "Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z",
                            "Cal_pt_detid", "mc_pthatmin_gev", "mc_pthatmax_gev"};
        for (auto b : br) t->SetBranchStatus(b, 1);
        t->SetBranchAddress("mcpart_num", &mc_n); t->SetBranchAddress("mcpart_geid", mc_geid);
        t->SetBranchAddress("mcpart_idVtx", mc_vtx); t->SetBranchAddress("mcpart_idVtxEnd", mc_vtxEnd);
        t->SetBranchAddress("mcpart_px", mc_px); t->SetBranchAddress("mcpart_py", mc_py); t->SetBranchAddress("mcpart_pz", mc_pz);
        t->SetBranchAddress("mcpart_E", mc_E);
        t->SetBranchAddress("mcpart_Vtx_x", mc_vx); t->SetBranchAddress("mcpart_Vtx_y", mc_vy); t->SetBranchAddress("mcpart_Vtx_z", mc_vz);
        t->SetBranchAddress("Cal_npts", &npts); t->SetBranchAddress("Cal_pt_energy", pe);
        t->SetBranchAddress("Cal_pt_x", pxx); t->SetBranchAddress("Cal_pt_y", pyy); t->SetBranchAddress("Cal_pt_z", pzz);
        t->SetBranchAddress("Cal_pt_detid", pdet);
        tmin = tmax = -1;
        if (t->GetBranch("mc_pthatmin_gev")) t->SetBranchAddress("mc_pthatmin_gev", &tmin);
        if (t->GetBranch("mc_pthatmax_gev")) t->SetBranchAddress("mc_pthatmax_gev", &tmax);

        const Long64_t N = t->GetEntries();
        for (Long64_t e = 0; e < N; e++) {
            if (t->GetEntry(e) <= 0) continue;
            if (mc_n > kMaxMc) mc_n = kMaxMc;
            ev = (int)e; pthat_min = tmin; pthat_max = tmax;
            float pvx = 0, pvy = 0, pvz = 0;
            if (!FindEmTruthPrimaryVertex(mc_n, mc_vtx, mc_vx, mc_vy, mc_vz, pvx, pvy, pvz)) continue;
            vz_true = pvz;

            // reco: ECal points as JetMatcher (em_only, use_points)
            vector<PseudoJet> parts; vector<int> src;
            for (int i = 0; i < npts && i < kMaxPt; i++) {
                const int d = pdet[i];
                if (d != 0 && d != 1) continue;
                if (!pass_point_energy_cut(pe[i])) continue;
                const double nrm = sqrt(pxx[i] * pxx[i] + pyy[i] * pyy[i] + pzz[i] * pzz[i]);
                parts.emplace_back(pe[i] * pxx[i] / nrm, pe[i] * pyy[i] / nrm, pe[i] * pzz[i] / nrm, pe[i]);
                src.push_back(i); parts.back().set_user_index((int)src.size() - 1);
            }
            if (parts.size() < 1) continue;
            ClusterSequence cs(parts, JetDefinition(antikt_algorithm, R));
            vector<PseudoJet> jets;
            for (auto &j : cs.inclusive_jets(reco_jet_pt_min_em)) {
                const float X = z_proj * j.px() / j.pz(), Y = z_proj * j.py() / j.pz();
                if (pass_fiducial_cut(X, Y, reco_fiducial_buffer_ecal, kFiducialRectEcal)) jets.push_back(j);
            }
            if (jets.empty()) continue;

            // truth photons, their parents, and the point-photon matching (all ECal points above the floor)
            map<int, int> endOf; // vertex id -> particle ending there
            for (int i = 0; i < mc_n; i++) if (mc_vtxEnd[i] > 0) endOf[mc_vtxEnd[i]] = i;
            vector<int> gam;
            for (int i = 0; i < mc_n; i++)
                if (IsEmTruthParticle(mc_geid[i], mc_vx[i], mc_vy[i], mc_vz[i], pvx, pvy, pvz) && mc_pz[i] > 0) gam.push_back(i);
            struct C { int p, g; double d; };
            vector<C> cand;
            auto proj = [&](int g, double z, double &x, double &y) {
                const double s = (z - mc_vz[g]) / mc_pz[g];
                x = mc_vx[g] + s * mc_px[g]; y = mc_vy[g] + s * mc_py[g];
            };
            for (size_t a = 0; a < src.size(); a++) {
                const int i = src[a];
                for (int g : gam) {
                    double x, y; proj(g, pzz[i], x, y);
                    const double d = hypot(x - pxx[i], y - pyy[i]);
                    if (d < kMatchMaxCm) cand.push_back({(int)a, g, d});
                }
            }
            sort(cand.begin(), cand.end(), [](const C &a, const C &b) { return a.d < b.d; });
            vector<int> ptG(src.size(), -1); vector<double> ptD(src.size(), -1);
            map<int, int> gPt;
            for (auto &c : cand) {
                if (ptG[c.p] >= 0 || gPt.count(c.g)) continue;
                ptG[c.p] = c.g; ptD[c.p] = c.d; gPt[c.g] = c.p;
            }
            auto parentOf = [&](int g) { auto it = endOf.find(mc_vtx[g]); return it == endOf.end() ? -1 : it->second; };

            for (auto &j : jets) {
                auto cons = j.constituents();
                n = 0;
                vector<int> jp; // point slots (index into parts) of this jet
                for (auto &c : cons) jp.push_back(c.user_index());
                sort(jp.begin(), jp.end());
                if ((int)jp.size() > kMaxPt) continue;
                map<int, int> pi0Local; // mcpart index of the pi0 -> local id
                map<int, vector<int>> pi0Pts, etaPts;
                for (int a : jp) {
                    const int i = src[a];
                    pE[n] = pe[i]; px_[n] = pxx[i]; py_[n] = pyy[i]; pz_[n] = pzz[i];
                    const int g = ptG[a];
                    ppi0[n] = -1; ppar[n] = g < 0 ? -1 : 0; pdist[n] = (float)ptD[a]; pgE[n] = g < 0 ? -1 : mc_E[g];
                    pngeid[n] = 0; pndist[n] = -1;
                    if (g < 0) { // nearest particle of any kind, for what made an unmatched point
                        double best = 1e9;
                        for (int q = 0; q < mc_n; q++) {
                            if (!(mc_E[q] > 0.3) || !(mc_pz[q] > 0) || mc_vz[q] > pzz[i] - 10) continue;
                            double x, y; proj(q, pzz[i], x, y);
                            const double dd = hypot(x - pxx[i], y - pyy[i]);
                            if (dd < best) { best = dd; pngeid[n] = mc_geid[q]; }
                        }
                        if (best < 1e9) pndist[n] = (float)best;
                    }
                    nPts++;
                    if (g >= 0) {
                        nPtsMatched++;
                        const int par = parentOf(g);
                        if (par >= 0) {
                            ppar[n] = mc_geid[par];
                            if (mc_geid[par] == 7) {
                                if (!pi0Local.count(par)) { int id = (int)pi0Local.size(); pi0Local[par] = id; }
                                ppi0[n] = pi0Local[par]; pi0Pts[par].push_back(n);
                            } else if (mc_geid[par] == 17) etaPts[par].push_back(n);
                        }
                    }
                    n++;
                }
                k_res = n_merged = n_lost = n_eta_res = 0;
                for (auto &kv : etaPts) if (kv.second.size() >= 2) n_eta_res++;
                for (auto &kv : pi0Pts) {
                    if (kv.second.size() >= 2) { k_res++; continue; }
                    // the other daughter photon: merged into a jet point, or lost
                    const int par = kv.first, mine = ptG[jp[kv.second[0]]];
                    int other = -1;
                    for (int g : gam) if (g != mine && parentOf(g) == par) { other = g; break; }
                    bool merged = false;
                    if (other >= 0)
                        for (int q = 0; q < n && !merged; q++) {
                            double x, y; proj(other, pz_[q], x, y);
                            if (hypot(x - px_[q], y - py_[q]) < kMergeCm) merged = true;
                        }
                    if (merged) n_merged++; else n_lost++;
                }
                jet_E = j.E(); jet_pt = j.perp(); jet_eta = j.eta(); jet_phi = j.phi_std();
                jet_xF = computeFeynmanX(j.E(), j.eta(), kSqrtS_MainRun22_GeV);
                tj->Fill(); nJets++;
            }
        }
        f->Close();
    }
    fout.cd(); tj->Write(); tf->Write(); fout.Close();
    cout << "jets " << nJets << ", jet points " << nPts << ", matched to a truth photon " << nPtsMatched
         << " (" << (nPts ? 100.0 * nPtsMatched / nPts : 0) << "%)" << endl;
    return 0;
}
