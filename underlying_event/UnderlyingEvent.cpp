// Underlying-event (UE) estimate for every jet of a set of files, by the
// in-event pixel ring (PixelRingUE.h): the jet's own event, its eta band from
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
//   reco_ue_cov    share of the jet's live area in eta bins the ring reaches
// and in simulation, per truth jet (EM-only truth photons, full azimuth):
//   truth_ue_rho, truth_ue_E (= truth_ue_rho x pi R^2), truth_ue_pt
// The energies and densities are -1 where the ring reaches none of the jet's
// eta bins (no estimate; the jet must be dropped where a UE is needed).
// Earlier reco_ue_* / truth_ue_* branches are replaced.
//
// simpletree_entry ties jetTree entries to SimpleTree entries; files without
// it need equal entry counts.
//
// Build (job folder):
//   g++ -m64 -O2 UnderlyingEvent.cpp VectorDict.cxx -o UnderlyingEvent `root-config --cflags --libs` -std=c++11
// from this folder:
//   g++ -m64 -O2 UnderlyingEvent.cpp ../shared/VectorDict.cxx -o UnderlyingEvent -I../shared `root-config --cflags --libs` -std=c++11

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "PixelRingUE.h"

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
    const TruthPixelRing truth_ring;
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
        if (!jt || !st || !st->GetBranch("Cal_npts") || (!has_entry && jt->GetEntries() != st->GetEntries()))
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
        const bool has_truth = jt->GetBranch("truth_eta") && st->GetBranch("mcpart_num");
        if (has_truth) { jt->SetBranchAddress("truth_eta", &truth_eta); jt->SetBranchAddress("truth_phi", &truth_phi); }

        Int_t Cal_npts = 0, mcpart_num = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        Int_t Cal_pt_detid[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pb[i], 1);
        st->SetBranchAddress("Cal_npts", &Cal_npts);
        st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x);
        st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z);
        st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);
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
        std::vector<int> p_pix;
        const Long64_t n = jt->GetEntries();
        for (Long64_t ev = 0; ev < n; ev++)
        {
            jt->GetEntry(ev);
            const size_t nr = reco_eta->size(), nt = has_truth ? truth_eta->size() : 0;
            reco_br.Clear(nr);
            if (has_truth) truth_br.Clear(nt);
            if ((nr > 0 || nt > 0) && st->GetEntry(has_entry ? simpletree_entry : ev) > 0)
            {
                // reco: ECal points, as the reco jets' own inputs (RecoJets.cpp / JetMatcher.cpp, use_points)
                p_eta.clear(); p_phi.clear(); p_E.clear(); p_pix.clear();
                for (int i = 0; i < std::min((int)Cal_npts, MAX_PTS) && nr > 0; i++)
                {
                    if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                    if (!pass_point_energy_cut(Cal_pt_energy[i], event_sqrt_s_gev)) continue;
                    const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                    p_eta.push_back(std::asinh(Cal_pt_z[i] / rho)); p_phi.push_back(std::atan2(Cal_pt_y[i], Cal_pt_x[i]));
                    p_E.push_back(Cal_pt_energy[i]); p_pix.push_back(geo.PixelOf(Cal_pt_x[i], Cal_pt_y[i]));
                }
                for (size_t j = 0; j < nr; j++)
                {
                    const UeResult r = RecoPixelRing(geo, reco_eta->at(j), reco_phi->at(j), p_eta, p_phi, p_E, p_pix);
                    reco_br.Set(j, r);
                    n_reco++;
                    if (r.rho < 0) n_reco_none++; else sum_reco += r.E();
                }
                // truth: the EM-only photons the truth jets are built from, full azimuth
                if (nt > 0)
                {
                    p_eta.clear(); p_phi.clear(); p_E.clear();
                    float vx = 0, vy = 0, vz = 0;
                    const int nmc = std::min((int)mcpart_num, MAX_MC);
                    if (FindEmTruthPrimaryVertex(nmc, &mc_idVtx[0], &mc_vx[0], &mc_vy[0], &mc_vz[0], vx, vy, vz))
                        for (int i = 0; i < nmc; i++)
                        {
                            if (!IsEmTruthParticle(mc_geid[i], mc_vx[i], mc_vy[i], mc_vz[i], vx, vy, vz)) continue;
                            const double pt = std::sqrt((double)mc_px[i] * mc_px[i] + (double)mc_py[i] * mc_py[i]);
                            if (!(pt > 0)) continue;
                            const double eta = std::asinh(mc_pz[i] / pt);
                            if (!InTruthParticleEtaRange(eta)) continue;
                            p_eta.push_back(eta); p_phi.push_back(std::atan2(mc_py[i], mc_px[i])); p_E.push_back(mc_E[i]);
                        }
                    for (size_t j = 0; j < nt; j++)
                    {
                        const UeResult r = truth_ring.Estimate(truth_eta->at(j), truth_phi->at(j), p_eta, p_phi, p_E);
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
              << n_reco_none << " without a ring, mean UE " << sum_reco / std::max(1L, n_reco - n_reco_none) << " GeV";
    if (n_truth > 0)
        std::cout << "; " << n_truth << " truth jets, " << n_truth_none << " without a ring, mean UE "
                  << sum_truth / std::max(1L, n_truth - n_truth_none) << " GeV";
    std::cout << std::endl;
    return 0;
}
