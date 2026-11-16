// Underlying-event (UE) estimate for every jet of a set of files, by the
// in-event eta-weighted side wedges (WeightedWedgeUE.h): the jet's own event, its eta band from
// 2R out to 150 deg in azimuth, over the whole ECal (not only the jets'
// fiducial region), weighted to the jet cone's eta distribution. Runs inside
// sim_to_jet and data_to_jet jobs; each file is independent.
//
// Usage: ./UnderlyingEvent pairs.list
//   pairs.list  one line per file: "<jet_output file> <SimpleTree file>"
//
// Each jet_output file is rewritten in place with, per reco jet (ECal points,
// the same selection as the reco jets' inputs):
//   reco_ue_rho    UE density, GeV per unit eta-phi area (live)
//   reco_ue_area   the jet cone's live eta-phi area on the ECal
//   reco_ue_E      reco_ue_rho x reco_ue_area: the UE to subtract
//   reco_ue_pt     the same in pT
//   reco_ue_cov    share of the jet's live area in eta bins the wedges reach
// and in simulation, per truth jet (EM-only truth photons, full azimuth):
//   truth_ue_rho, truth_ue_E (= truth_ue_rho x pi R^2), truth_ue_pt
// The energies and densities are -1 where the wedges reach none of the jet's
// eta bins (no estimate; the jet must be dropped where a UE is needed).
// Earlier reco_ue_* / truth_ue_* branches are replaced.
//
// simpletree_entry ties jetTree entries to SimpleTree entries; files without
// it need equal entry counts.
//
// Event vertex (2026-11-16): the reco wedges are seen from the vertex the reco
// jets were built from (jetTree reco_vz_used, FCS_JET_VERTEX in RecoJets /
// JetMatcher; absent or 0 = the origin): point directions and the tile
// geometry use z - vz, so the wedges and the jet cone's area match the jet
// axis. Truth needs nothing (momenta, exact areas).
//
// Full (ECAL+HCAL) jets (2026-11-16): a jet file with reco_em_only = 0 (sim_to_jet
// -J jet_full_*.root, or RecoJets in full mode) gets the same wedges over the
// full-jet rectangle (WeightedWedgeUE.h), with the reco jets' own inputs
// (reco_from_points = 1: ECal points + HCal clusters; 0: ECal and HCal tower
// hits above the jet finders' MIP thresholds) and, in simulation, the full-jet truth particles (first-vertex
// particles in the full-jet eta band, as JetMatcher). Not yet validated for full jets.
//
// Build (job folder):
//   g++ -m64 -O2 UnderlyingEvent.cpp VectorDict.cxx -o UnderlyingEvent `root-config --cflags --libs` -std=c++11
// from this folder:
//   g++ -m64 -O2 UnderlyingEvent.cpp ../shared/VectorDict.cxx -o UnderlyingEvent -I../shared `root-config --cflags --libs` -std=c++11

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "WeightedWedgeUE.h"

// One level's output branches, one value per jet
struct UeBranches
{
    std::vector<float> E, pt, rho, area, cov;
    void Attach(TTree *t, const std::string &level, bool with_area)
    {
        t->Branch((level + "_ue_E").c_str(), &E);
        t->Branch((level + "_ue_pt").c_str(), &pt);
        t->Branch((level + "_ue_rho").c_str(), &rho);
        if (with_area) { t->Branch((level + "_ue_area").c_str(), &area); t->Branch((level + "_ue_cov").c_str(), &cov); }
    }
    void Clear(size_t n) { E.assign(n, -1.f); pt.assign(n, -1.f); rho.assign(n, -1.f); area.assign(n, -1.f); cov.assign(n, 0.f); }
    void Set(size_t j, const UeResult &r)
    {
        E[j] = (float)r.E(); pt[j] = (float)r.pt(); rho[j] = (float)r.rho; area[j] = (float)r.area; cov[j] = (float)r.cov;
    }
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2)
    {
        std::cerr << "usage: ./UnderlyingEvent pairs.list" << std::endl;
        return 1;
    }
    std::vector<std::pair<std::string, std::string> > files;
    {
        std::ifstream in(argv[1]);
        std::string line;
        while (std::getline(in, line))
        {
            std::istringstream ss(line);
            std::string a, b;
            if (ss >> a >> b) files.push_back(std::make_pair(a, b));
        }
    }
    if (files.empty())
    {
        std::cerr << "ERROR: no file pairs in " << argv[1] << std::endl;
        return 1;
    }

    UeGeometry geo;
    geo.Build();
    const int MAX_HITS = 5000; // matches StSimpleReaderMaker.h's Cal_hit_* array size
    const TruthWeightedWedges truth_wedges;
    const int MAX_MC = 10000; // matches JetMatcher.cpp's mcpart_* array size
    const int MAX_PTS = 200;  // matches StSimpleReaderMaker.h's Cal_pt_* array size
    std::vector<Int_t> mc_geid(MAX_MC), mc_idVtx(MAX_MC);
    std::vector<Float_t> mc_px(MAX_MC), mc_py(MAX_MC), mc_pz(MAX_MC), mc_E(MAX_MC), mc_vx(MAX_MC), mc_vy(MAX_MC), mc_vz(MAX_MC);
    long n_bad = 0, n_reco = 0, n_reco_none = 0, n_truth = 0, n_truth_none = 0;
    double sum_reco = 0, sum_truth = 0;

    for (size_t f = 0; f < files.size(); f++)
    {
        const std::string &jet_file = files[f].first, &simple_file = files[f].second;
        TFile *fin = TFile::Open(jet_file.c_str());
        TFile *fs = TFile::Open(simple_file.c_str());
        TTree *jt = (fin && !fin->IsZombie()) ? (TTree *)fin->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        const bool has_entry = jt && jt->GetBranch("simpletree_entry");
        if (!jt || !st || (!st->GetBranch("Cal_npts") && !st->GetBranch("Cal_nhits")) || (!has_entry && jt->GetEntries() != st->GetEntries()))
        {
            std::cerr << "ERROR: cannot use " << jet_file << " with " << simple_file
                      << " (unreadable, no points, or no simpletree_entry branch and differing entry counts); left as it is" << std::endl;
            n_bad++;
            if (fin) { fin->Close(); delete fin; }
            if (fs) { fs->Close(); delete fs; }
            continue;
        }
        // earlier UE branches stay out of the copy
        TIter next_branch(jt->GetListOfBranches());
        while (TObject *b = next_branch())
        {
            const std::string name = b->GetName();
            if (name.compare(0, 8, "reco_ue_") == 0 || name.compare(0, 9, "truth_ue_") == 0) jt->SetBranchStatus(name.c_str(), 0);
        }

        std::vector<float> *reco_eta = 0, *reco_phi = 0, *truth_eta = 0, *truth_phi = 0;
        Long64_t simpletree_entry = -1;
        Double_t event_sqrt_s_gev = kSqrtS_MainRun22_GeV;
        jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi);
        if (has_entry) jt->SetBranchAddress("simpletree_entry", &simpletree_entry);
        // the point energy cut depends on sqrt(s), as in the jet finders
        if (jt->GetBranch("event_sqrt_s_gev")) jt->SetBranchAddress("event_sqrt_s_gev", &event_sqrt_s_gev);
        // the vertex the reco jets were built from (0 = origin, as before 2026-11-16)
        Float_t reco_vz_used = 0, reco_vx_used = 0, reco_vy_used = 0;
        if (jt->GetBranch("reco_vz_used")) jt->SetBranchAddress("reco_vz_used", &reco_vz_used);
        if (jt->GetBranch("reco_vx_used")) jt->SetBranchAddress("reco_vx_used", &reco_vx_used);
        if (jt->GetBranch("reco_vy_used")) jt->SetBranchAddress("reco_vy_used", &reco_vy_used);
        // jet type and reco inputs of this file (absent: EM-only point jets, as before)
        Int_t file_em_only = 1, file_from_points = 1;
        if (jt->GetEntries() > 0)
        {
            if (TBranch *b = jt->GetBranch("reco_em_only")) { b->GetEntry(0); file_em_only = (Int_t)jt->GetLeaf("reco_em_only")->GetValue(); }
            if (TBranch *b = jt->GetBranch("reco_from_points")) { b->GetEntry(0); file_from_points = (Int_t)jt->GetLeaf("reco_from_points")->GetValue(); }
        }
        const bool em_jets = file_em_only != 0, from_points = file_from_points != 0;
        if (geo.em_only != em_jets) { geo.em_only = em_jets; geo.Build(geo.vz, geo.vx, geo.vy); }
        const bool has_truth = jt->GetBranch("truth_eta") && st->GetBranch("mcpart_num");
        if (has_truth) { jt->SetBranchAddress("truth_eta", &truth_eta); jt->SetBranchAddress("truth_phi", &truth_phi); }

        // reco inputs: ECal points, or (from_points = 0) tower hits
        Int_t Cal_npts = 0, mcpart_num = 0, Cal_nhits = 0;
        static Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        static Int_t Cal_pt_detid[MAX_PTS];
        static Float_t Cal_hit_energy[MAX_HITS], Cal_hit_posx[MAX_HITS], Cal_hit_posy[MAX_HITS], Cal_hit_posz[MAX_HITS];
        static Int_t Cal_detid[MAX_HITS];
        const int MAX_CLUS = 100; // StSimpleReaderMaker.h Cal_clus_* size
        Int_t Cal_nclus = 0;
        static Float_t Cal_clus_energy[MAX_CLUS], Cal_clus_x[MAX_CLUS], Cal_clus_y[MAX_CLUS], Cal_clus_z[MAX_CLUS];
        static Int_t Cal_clus_detid[MAX_CLUS];
        const bool with_clusters = from_points && !em_jets && st->GetBranch("Cal_nclus");
        st->SetBranchStatus("*", 0);
        if (from_points)
        {
            const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
            for (int i = 0; i < 6; i++) st->SetBranchStatus(pb[i], 1);
            st->SetBranchAddress("Cal_npts", &Cal_npts);
            st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
            st->SetBranchAddress("Cal_pt_x", Cal_pt_x);
            st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
            st->SetBranchAddress("Cal_pt_z", Cal_pt_z);
            st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);
            if (with_clusters)
            {
                const char *cb[] = {"Cal_nclus", "Cal_clus_energy", "Cal_clus_x", "Cal_clus_y", "Cal_clus_z", "Cal_clus_detid"};
                for (int i = 0; i < 6; i++) st->SetBranchStatus(cb[i], 1);
                st->SetBranchAddress("Cal_nclus", &Cal_nclus);
                st->SetBranchAddress("Cal_clus_energy", Cal_clus_energy);
                st->SetBranchAddress("Cal_clus_x", Cal_clus_x);
                st->SetBranchAddress("Cal_clus_y", Cal_clus_y);
                st->SetBranchAddress("Cal_clus_z", Cal_clus_z);
                st->SetBranchAddress("Cal_clus_detid", Cal_clus_detid);
            }
        }
        else
        {
            const char *hb[] = {"Cal_nhits", "Cal_hit_energy", "Cal_hit_posx", "Cal_hit_posy", "Cal_hit_posz", "Cal_detid"};
            for (int i = 0; i < 6; i++) st->SetBranchStatus(hb[i], 1);
            st->SetBranchAddress("Cal_nhits", &Cal_nhits);
            st->SetBranchAddress("Cal_hit_energy", Cal_hit_energy);
            st->SetBranchAddress("Cal_hit_posx", Cal_hit_posx);
            st->SetBranchAddress("Cal_hit_posy", Cal_hit_posy);
            st->SetBranchAddress("Cal_hit_posz", Cal_hit_posz);
            st->SetBranchAddress("Cal_detid", Cal_detid);
        }
        if (has_truth)
        {
            const char *mb[] = {"mcpart_num", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_geid", "mcpart_idVtx",
                                "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
            for (int i = 0; i < 10; i++) st->SetBranchStatus(mb[i], 1);
            st->SetBranchAddress("mcpart_num", &mcpart_num);
            st->SetBranchAddress("mcpart_px", &mc_px[0]); st->SetBranchAddress("mcpart_py", &mc_py[0]); st->SetBranchAddress("mcpart_pz", &mc_pz[0]);
            st->SetBranchAddress("mcpart_E", &mc_E[0]); st->SetBranchAddress("mcpart_geid", &mc_geid[0]); st->SetBranchAddress("mcpart_idVtx", &mc_idVtx[0]);
            st->SetBranchAddress("mcpart_Vtx_x", &mc_vx[0]); st->SetBranchAddress("mcpart_Vtx_y", &mc_vy[0]); st->SetBranchAddress("mcpart_Vtx_z", &mc_vz[0]);
        }

        const std::string tmp_name = jet_file + ".ue_tmp";
        TFile *fout = TFile::Open(tmp_name.c_str(), "RECREATE");
        TTree *out = jt->CloneTree(0);
        UeBranches reco_br, truth_br;
        reco_br.Attach(out, "reco", true);
        if (has_truth) truth_br.Attach(out, "truth", false);

        std::vector<float> p_eta, p_phi, p_E;
        std::vector<int> p_tile;
        const Long64_t n = jt->GetEntries();
        for (Long64_t ev = 0; ev < n; ev++)
        {
            jt->GetEntry(ev);
            const size_t nr = reco_eta->size(), nt = has_truth ? truth_eta->size() : 0;
            reco_br.Clear(nr);
            if (has_truth) truth_br.Clear(nt);
            if ((nr > 0 || nt > 0) && st->GetEntry(has_entry ? simpletree_entry : ev) > 0)
            {
                // reco: ECal points, as the reco jets' own inputs (RecoJets.cpp / JetMatcher.cpp, use_points),
                // seen from the jets' vertex; the tile geometry rebuilt when the vertex changes
                if (nr > 0 && (std::fabs(reco_vz_used - geo.vz) > 1e-3 || std::fabs(reco_vx_used - geo.vx) > 1e-5 || std::fabs(reco_vy_used - geo.vy) > 1e-5))
                    geo.Build(reco_vz_used, reco_vx_used, reco_vy_used);
                p_eta.clear(); p_phi.clear(); p_E.clear(); p_tile.clear();
                const int nin = nr == 0 ? 0 : from_points ? std::min((int)Cal_npts, MAX_PTS) : std::min((int)Cal_nhits, MAX_HITS);
                for (int i = 0; i < nin; i++)
                {
                    float e, x, y, z;
                    if (from_points)
                    {
                        if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                        if (!pass_point_energy_cut(Cal_pt_energy[i])) continue;
                        e = Cal_pt_energy[i]; x = Cal_pt_x[i]; y = Cal_pt_y[i]; z = Cal_pt_z[i];
                    }
                    else
                    {
                        // the jet finders' hit selection (JetMatcher.cpp / RecoJets.cpp)
                        const int det = Cal_detid[i];
                        if (em_jets && (det == 2 || det == 3)) continue;
                        e = Cal_hit_energy[i];
                        if (!(((det == 0 || det == 1) && e > mip_threshold * ecal_mip) ||
                              ((det == 2 || det == 3) && e > mip_threshold * hcal_mip))) continue;
                        x = Cal_hit_posx[i]; y = Cal_hit_posy[i]; z = Cal_hit_posz[i];
                    }
                    const double px = x - geo.vx, py = y - geo.vy;
                    const double rho = std::sqrt(px * px + py * py);
                    p_eta.push_back(std::asinh((z - geo.vz) / rho)); p_phi.push_back(std::atan2(py, px));
                    p_E.push_back(e); p_tile.push_back(geo.TileOf(x, y));
                }
                // full point jets: plus the HCal clusters, as in the jet finders
                for (int i = 0; with_clusters && nr > 0 && i < std::min((int)Cal_nclus, MAX_CLUS); i++)
                {
                    if (Cal_clus_detid[i] != 2 && Cal_clus_detid[i] != 3) continue;
                    const float e = Cal_clus_energy[i];
                    if (!(e > 0)) continue;
                    const double px = Cal_clus_x[i] - geo.vx, py = Cal_clus_y[i] - geo.vy;
                    const double rho = std::sqrt(px * px + py * py);
                    p_eta.push_back(std::asinh((Cal_clus_z[i] - geo.vz) / rho)); p_phi.push_back(std::atan2(py, px));
                    // tiles live at the ECal front plane (z_proj): project the deeper HCal cluster there
                    const double sc = (z_proj - geo.vz) / (Cal_clus_z[i] - geo.vz);
                    p_E.push_back(e); p_tile.push_back(geo.TileOf(geo.vx + px * sc, geo.vy + py * sc));
                }
                for (size_t j = 0; j < nr; j++)
                {
                    const UeResult r = RecoWeightedWedges(geo, reco_eta->at(j), reco_phi->at(j), p_eta, p_phi, p_E, p_tile);
                    reco_br.Set(j, r);
                    n_reco++;
                    if (r.rho < 0) n_reco_none++; else sum_reco += r.E();
                }
                // truth: the particles the truth jets are built from, full azimuth
                // (EM-only: the EM truth photons; full: first-vertex particles, as JetMatcher)
                if (nt > 0)
                {
                    p_eta.clear(); p_phi.clear(); p_E.clear();
                    float vx = 0, vy = 0, vz = 0;
                    const int nmc = std::min((int)mcpart_num, MAX_MC);
                    const bool have_vtx = !em_jets || FindEmTruthPrimaryVertex(nmc, &mc_idVtx[0], &mc_vx[0], &mc_vy[0], &mc_vz[0], vx, vy, vz);
                    if (have_vtx)
                        for (int i = 0; i < nmc; i++)
                        {
                            if (em_jets ? !IsEmTruthParticle(mc_geid[i], mc_vx[i], mc_vy[i], mc_vz[i], vx, vy, vz) : mc_idVtx[i] != 1) continue;
                            const double pt = std::sqrt((double)mc_px[i] * mc_px[i] + (double)mc_py[i] * mc_py[i]);
                            if (!(pt > 0)) continue;
                            const double eta = std::asinh(mc_pz[i] / pt);
                            if (!InTruthParticleEtaRange(eta, em_jets)) continue;
                            p_eta.push_back(eta); p_phi.push_back(std::atan2(mc_py[i], mc_px[i])); p_E.push_back(mc_E[i]);
                        }
                    for (size_t j = 0; j < nt; j++)
                    {
                        const UeResult r = truth_wedges.Estimate(truth_eta->at(j), truth_phi->at(j), p_eta, p_phi, p_E);
                        truth_br.Set(j, r);
                        n_truth++;
                        if (r.rho < 0) n_truth_none++; else sum_truth += r.E();
                    }
                }
            }
            out->Fill();
        }
        fout->cd();
        out->Write("", TObject::kOverwrite);
        fout->Close(); delete fout;
        fin->Close(); delete fin;
        fs->Close(); delete fs;
        if (gSystem->Rename(tmp_name.c_str(), jet_file.c_str()) != 0)
            std::cerr << "ERROR: could not replace " << jet_file << " with " << tmp_name << std::endl;
    }

    std::cout << "UnderlyingEvent: " << files.size() << " file pairs (" << n_bad << " unusable); " << n_reco << " reco jets, "
              << n_reco_none << " without wedges, mean UE " << sum_reco / std::max(1L, n_reco - n_reco_none) << " GeV";
    if (n_truth > 0)
        std::cout << "; " << n_truth << " truth jets, " << n_truth_none << " without wedges, mean UE "
                  << sum_truth / std::max(1L, n_truth - n_truth_none) << " GeV";
    std::cout << std::endl;
    return 0;
}
