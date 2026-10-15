// Underlying-event (UE) estimate for every jet of a set of files, by the near
// off-axis proxy-event method (ProxySelection.h, ProxySearch.h). Runs inside
// sim_to_jet and data_to_jet jobs, and over whole simulation samples
// (submit_ue_by_pthat.sh).
//
// Usage: ./UnderlyingEvent pairs.list [n_max] [minbias]
//   pairs.list  one line per file: "<jet_output file> <SimpleTree file>"
//   n_max       most proxy events per reco jet, 0 (default) = all found
//   minbias     for samples without trigger or jet pT cuts (JES calibration):
//               every fiducial reco jet is estimated, untriggered events
//               match untriggered ones
//
// Each jet_output file is rewritten in place with, per reco jet:
//   reco_ue_E          mean ECal point energy (GeV) in an R cone at the jet
//                      position, over the proxy events
//   reco_ue_pt         the same for the points' summed pT
//   reco_ue_nproxy     proxy events used; 0 = none found, -1 = not estimated
//                      (no EM trigger, or jet outside the fiducial region)
//   reco_ue_pt_window  pT window half-width used, as a fraction of jet pT
// reco_ue_E and reco_ue_pt are -1 where reco_ue_nproxy <= 0.
//
// In simulation every truth jet also gets a same-event estimate, as the FMS
// analysis: the mean energy (and pT) of the EM-only truth photons in two R
// cones at the jet's eta and phi +- 90 deg, whatever the cones contain. The
// photons are those the truth jets are built from (IsEmTruthParticle, truth
// particle eta band), in momentum eta-phi, with full azimuth:
//   truth_ue_E, truth_ue_pt   per truth jet, GeV
//   truth_ue_ncone            2 (cones used); 0 = not estimated
//
// Reco groups: data files jet_output_<run>_<file>.root are grouped by run and
// ordered by the last five digits of <file> (time); simulation files
// ..._pt<ptHatMin>_<job>.root by ptHat bin (events are unweighted within a
// bin). Other names count as run 0, in list order.
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

#include "JetParameters.h"
#include "ProxySelection.h"
#include "ProxySearch.h"

struct FilePair
{
    std::string jet_file, simple_file;
    int run, sequence, order;                    // order: line number in the list
    std::vector<int> pool_event;                 // per jetTree entry: its pool event, or -1
    std::vector<std::vector<float> > truth_E, truth_pt;  // per jetTree entry and truth jet: same-event cone UE
    bool has_truth;
};

// Counts for the end-of-job summary, per level
struct UeSummary
{
    long n_jets, n_tried, n_min_reached, n_widened, n_none;
    double sum_n;
    UeSummary() : n_jets(0), n_tried(0), n_min_reached(0), n_widened(0), n_none(0), sum_n(0.0) {}
    void Print(const char *level, const char *which, int n_min) const
    {
        const long d = std::max(1L, n_tried);
        std::cout << "UnderlyingEvent: " << n_jets << " " << level << " jets, " << n_tried << " estimated (" << which << "); of those "
                  << "the mean number of proxy events is " << sum_n / d << ", " << 100.0 * n_min_reached / d << " % have at least "
                  << n_min << ", " << 100.0 * n_widened / d << " % used a widened pT window, " << 100.0 * n_none / d << " % found none"
                  << std::endl;
    }
};

// The four output branches of one level (reco or truth), one value per jet
struct UeBranches
{
    std::vector<float> E, pt, window;
    std::vector<int> nproxy;
    void Attach(TTree *t, const std::string &level)
    {
        t->Branch((level + "_ue_E").c_str(), &E);
        t->Branch((level + "_ue_pt").c_str(), &pt);
        t->Branch((level + "_ue_nproxy").c_str(), &nproxy);
        t->Branch((level + "_ue_pt_window").c_str(), &window);
    }
    // pool_event: the entry's event in the pool, or -1 (no estimate for any jet)
    void Fill(const ProxyPool &pool, int pool_event, size_t nj, int n_max, int n_min, UeSummary &sum)
    {
        E.assign(nj, -1.f); pt.assign(nj, -1.f);
        nproxy.assign(nj, -1); window.assign(nj, (float)proxy_pt_window_fraction);
        for (size_t j = 0; j < nj && pool_event >= 0; j++)
        {
            const ProxyResult res = pool.Estimate(pool_event, (int)j, n_max);
            nproxy[j] = res.n;
            window[j] = (float)res.window;
            if (res.n > 0) { E[j] = (float)res.ue_E; pt[j] = (float)res.ue_pt; }

            if (res.n < 0) continue;
            sum.n_tried++;
            sum.sum_n += res.n;
            if (res.n >= n_min) sum.n_min_reached++;
            if (res.n_nominal < n_min && res.n > res.n_nominal) sum.n_widened++;
            if (res.n == 0) sum.n_none++;
        }
        sum.n_jets += nj;
    }
};

// Truth-level UE of a jet at (eta, phi): mean energy and pT of the truth
// photons in two R cones at the jet's eta and phi +- 90 deg (same event)
static void TruthOffAxisCones(double eta, double phi, const std::vector<float> &p_eta, const std::vector<float> &p_phi,
                              const std::vector<float> &p_E, float &ue_E, float &ue_pt)
{
    double e = 0.0, pt = 0.0;
    for (int c = 0; c < 2; c++)
    {
        const double cone_phi = phi + (c == 0 ? 0.5 : -0.5) * M_PI;
        for (size_t i = 0; i < p_E.size(); i++)
            if (ProxyDeltaR(p_eta[i], p_phi[i], eta, cone_phi) < R)
            {
                e += p_E[i];
                pt += p_E[i] / std::cosh(p_eta[i]);
            }
    }
    ue_E = (float)(0.5 * e);
    ue_pt = (float)(0.5 * pt);
}

static bool ByTime(const FilePair &a, const FilePair &b)
{
    if (a.run != b.run) return a.run < b.run;
    if (a.sequence != b.sequence) return a.sequence < b.sequence;
    return a.order < b.order;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2)
    {
        std::cerr << "usage: ./UnderlyingEvent pairs.list [n_max] [minbias]" << std::endl;
        return 1;
    }
    const int n_max = argc > 2 ? std::atoi(argv[2]) : proxy_n_max; // 0 = every proxy event found
    const bool minimum_bias = argc > 3 && std::string(argv[3]) == "minbias";

    std::vector<FilePair> files;
    {
        std::ifstream in(argv[1]);
        std::string line;
        while (std::getline(in, line))
        {
            std::istringstream ss(line);
            FilePair fp;
            if (!(ss >> fp.jet_file >> fp.simple_file)) continue;
            fp.run = 0; fp.sequence = 0; fp.order = (int)files.size(); fp.has_truth = false;
            const size_t slash = fp.jet_file.find_last_of('/');
            const std::string base = slash == std::string::npos ? fp.jet_file : fp.jet_file.substr(slash + 1);
            int run = 0, file_number = 0;
            if (std::sscanf(base.c_str(), "jet_output_%d_%d.root", &run, &file_number) == 2)
            {
                fp.run = run; fp.sequence = file_number % 100000;
            }
            else
            {
                // simulation: ..._pt<ptHatMin>_<job>.root; groups are negative so they never equal a run number
                const size_t p = base.rfind("_pt");
                int pthat = 0, job = 0;
                if (p != std::string::npos && std::sscanf(base.c_str() + p, "_pt%d_%d.root", &pthat, &job) == 2)
                {
                    fp.run = -1 - pthat; fp.sequence = job;
                }
            }
            files.push_back(fp);
        }
    }
    if (files.empty())
    {
        std::cerr << "ERROR: no file pairs in " << argv[1] << std::endl;
        return 1;
    }
    std::sort(files.begin(), files.end(), ByTime);

    // ---------------- pass 1: read every event into the pool ----------------
    ProxyPool pool;
    if (minimum_bias) pool.SetMinimumBias();
    const int MAX_MC = 10000; // matches JetMatcher.cpp's mcpart_* array size
    std::vector<Int_t> mc_geid(MAX_MC), mc_idVtx(MAX_MC);
    std::vector<Float_t> mc_px(MAX_MC), mc_py(MAX_MC), mc_pz(MAX_MC), mc_E(MAX_MC), mc_vx(MAX_MC), mc_vy(MAX_MC), mc_vz(MAX_MC);
    long n_truth_files = 0;
    const int MAX_PTS = 200; // matches StSimpleReaderMaker.h's Cal_pt_* array size
    int n_bad = 0;

    for (size_t f = 0; f < files.size(); f++)
    {
        FilePair &fp = files[f];
        TFile *fj = TFile::Open(fp.jet_file.c_str());
        TFile *fs = TFile::Open(fp.simple_file.c_str());
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        const bool has_entry = jt && jt->GetBranch("simpletree_entry");
        if (!jt || !st || !st->GetBranch("Cal_npts") || (!has_entry && jt->GetEntries() != st->GetEntries()))
        {
            std::cerr << "ERROR: cannot use " << fp.jet_file << " with " << fp.simple_file
                      << " (unreadable, no points, or no simpletree_entry branch and differing entry counts);"
                      << " its jets get no UE estimate" << std::endl;
            n_bad++;
            if (fj) { fj->Close(); delete fj; }
            if (fs) { fs->Close(); delete fs; }
            continue;
        }

        std::vector<float> *reco_pt = 0, *reco_eta = 0, *reco_phi = 0, *reco_x = 0, *reco_y = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        Long64_t simpletree_entry = -1;
        Double_t event_sqrt_s_gev = kSqrtS_MainRun22_GeV;
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int i = 0; i < 8; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &reco_pt);
        jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi);
        jt->SetBranchAddress("reco_x", &reco_x);
        jt->SetBranchAddress("reco_y", &reco_y);
        jt->SetBranchAddress("fcsEM0", &em0);
        jt->SetBranchAddress("fcsEM1", &em1);
        jt->SetBranchAddress("fcsEM2", &em2);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &simpletree_entry); }
        // The energy cut on points depends on sqrt(s), as in the jet finders
        if (jt->GetBranch("event_sqrt_s_gev")) { jt->SetBranchStatus("event_sqrt_s_gev", 1); jt->SetBranchAddress("event_sqrt_s_gev", &event_sqrt_s_gev); }

        Int_t Cal_npts = 0;
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

        // Truth level (simulation only)
        fp.has_truth = jt->GetBranch("truth_pt") && st->GetBranch("mcpart_num");
        std::vector<float> *truth_pt = 0, *truth_eta = 0, *truth_phi = 0;
        Int_t mcpart_num = 0;
        if (fp.has_truth)
        {
            n_truth_files++;
            const char *tb[] = {"truth_pt", "truth_eta", "truth_phi"};
            for (int i = 0; i < 3; i++) jt->SetBranchStatus(tb[i], 1);
            jt->SetBranchAddress("truth_pt", &truth_pt);
            jt->SetBranchAddress("truth_eta", &truth_eta);
            jt->SetBranchAddress("truth_phi", &truth_phi);
            const char *mb[] = {"mcpart_num", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_geid", "mcpart_idVtx",
                                "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
            for (int i = 0; i < 10; i++) st->SetBranchStatus(mb[i], 1);
            st->SetBranchAddress("mcpart_num", &mcpart_num);
            st->SetBranchAddress("mcpart_px", &mc_px[0]);
            st->SetBranchAddress("mcpart_py", &mc_py[0]);
            st->SetBranchAddress("mcpart_pz", &mc_pz[0]);
            st->SetBranchAddress("mcpart_E", &mc_E[0]);
            st->SetBranchAddress("mcpart_geid", &mc_geid[0]);
            st->SetBranchAddress("mcpart_idVtx", &mc_idVtx[0]);
            st->SetBranchAddress("mcpart_Vtx_x", &mc_vx[0]);
            st->SetBranchAddress("mcpart_Vtx_y", &mc_vy[0]);
            st->SetBranchAddress("mcpart_Vtx_z", &mc_vz[0]);
        }

        const Long64_t n = jt->GetEntries();
        fp.pool_event.assign(n, -1);
        if (fp.has_truth) { fp.truth_E.assign(n, std::vector<float>()); fp.truth_pt.assign(n, std::vector<float>()); }
        std::vector<ProxyJet> jets;
        std::vector<float> p_eta, p_phi, p_E;
        for (Long64_t ev = 0; ev < n; ev++)
        {
            if (jt->GetEntry(ev) <= 0) continue;
            const bool any_truth = fp.has_truth && !truth_pt->empty();
            if (reco_pt->empty() && !any_truth) continue;      // no jet to estimate for
            if (st->GetEntry(has_entry ? simpletree_entry : ev) <= 0) continue;
            const int em_level = ProxyEmLevel(em0, em1, em2);

            // Truth level: same-event off-axis cones over the EM-only truth
            // photons in the truth particle eta band, as the truth jets
            if (any_truth)
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
                std::vector<float> &tE = fp.truth_E[ev], &tpt = fp.truth_pt[ev];
                tE.resize(truth_pt->size()); tpt.resize(truth_pt->size());
                for (size_t j = 0; j < truth_pt->size(); j++)
                    TruthOffAxisCones(truth_eta->at(j), truth_phi->at(j), p_eta, p_phi, p_E, tE[j], tpt[j]);
            }
            if (reco_pt->empty()) continue;

            jets.clear(); p_eta.clear(); p_phi.clear(); p_E.clear();
            for (size_t j = 0; j < reco_pt->size(); j++)
            {
                ProxyJet jet = {reco_pt->at(j), reco_eta->at(j), reco_phi->at(j), reco_x->at(j), reco_y->at(j)};
                jets.push_back(jet);
            }
            // Points: the same selection and direction convention as the reco
            // jets' own inputs (RecoJets.cpp / JetMatcher.cpp, use_points)
            for (int i = 0; i < Cal_npts; i++)
            {
                if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                if (!pass_point_energy_cut(Cal_pt_energy[i], event_sqrt_s_gev)) continue;
                const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                p_eta.push_back(std::asinh(Cal_pt_z[i] / rho)); p_phi.push_back(std::atan2(Cal_pt_y[i], Cal_pt_x[i]));
                p_E.push_back(Cal_pt_energy[i]);
            }
            fp.pool_event[ev] = pool.AddEvent(fp.run, em_level, jets, p_eta, p_phi, p_E);
        }
        fj->Close(); delete fj;
        fs->Close(); delete fs;
    }
    std::cout << "UnderlyingEvent: " << files.size() << " file pairs (" << n_bad << " unusable), " << pool.NEvents()
              << " events with a reco jet, " << pool.NJets() << " jets, " << pool.MemoryMB() << " MB" << std::endl;

    // ---------------- pass 2: estimate, and rewrite each jet file ----------------
    const int n_min = n_max > 0 ? std::min(proxy_n_min, n_max) : proxy_n_min;
    UeSummary reco_sum;
    long n_truth_jets = 0;
    double truth_ue_sum = 0.0;

    for (size_t f = 0; f < files.size(); f++)
    {
        FilePair &fp = files[f];
        TFile *fin = TFile::Open(fp.jet_file.c_str());
        TTree *jt = (fin && !fin->IsZombie()) ? (TTree *)fin->Get("jetTree") : 0;
        if (!jt)
        {
            std::cerr << "ERROR: cannot reopen " << fp.jet_file << "; left as it is" << std::endl;
            if (fin) { fin->Close(); delete fin; }
            continue;
        }
        // A file already processed keeps its old UE branches out of the copy
        TIter next_branch(jt->GetListOfBranches());
        while (TObject *b = next_branch())
        {
            const std::string name = b->GetName();
            if (name.compare(0, 8, "reco_ue_") == 0 || name.compare(0, 9, "truth_ue_") == 0) jt->SetBranchStatus(name.c_str(), 0);
        }

        std::vector<float> *reco_pt = 0, *truth_pt = 0;
        jt->SetBranchAddress("reco_pt", &reco_pt);
        if (fp.has_truth) jt->SetBranchAddress("truth_pt", &truth_pt);

        const std::string tmp_name = fp.jet_file + ".ue_tmp";
        TFile *fout = TFile::Open(tmp_name.c_str(), "RECREATE");
        TTree *out = jt->CloneTree(0);
        UeBranches reco_br;
        reco_br.Attach(out, "reco");
        std::vector<float> truth_ue_E, truth_ue_pt;
        std::vector<int> truth_ue_ncone;
        if (fp.has_truth)
        {
            out->Branch("truth_ue_E", &truth_ue_E);
            out->Branch("truth_ue_pt", &truth_ue_pt);
            out->Branch("truth_ue_ncone", &truth_ue_ncone);
        }

        const Long64_t n = jt->GetEntries();
        for (Long64_t ev = 0; ev < n; ev++)
        {
            jt->GetEntry(ev);
            const int pe = ev < (Long64_t)fp.pool_event.size() ? fp.pool_event[ev] : -1;
            reco_br.Fill(pool, pe, reco_pt->size(), n_max, n_min, reco_sum);
            if (fp.has_truth)
            {
                const size_t nj = truth_pt->size();
                const bool done = ev < (Long64_t)fp.truth_E.size() && fp.truth_E[ev].size() == nj;
                truth_ue_E.assign(nj, -1.f); truth_ue_pt.assign(nj, -1.f); truth_ue_ncone.assign(nj, 0);
                for (size_t j = 0; j < nj && done; j++)
                {
                    truth_ue_E[j] = fp.truth_E[ev][j]; truth_ue_pt[j] = fp.truth_pt[ev][j]; truth_ue_ncone[j] = 2;
                    n_truth_jets++; truth_ue_sum += truth_ue_E[j];
                }
            }
            out->Fill();
        }
        fout->cd();
        out->Write("", TObject::kOverwrite);
        fout->Close(); delete fout;
        fin->Close(); delete fin;

        if (gSystem->Rename(tmp_name.c_str(), fp.jet_file.c_str()) != 0)
            std::cerr << "ERROR: could not replace " << fp.jet_file << " with " << tmp_name << std::endl;
    }

    reco_sum.Print("reco", minimum_bias ? "fiducial; minimum-bias mode" : "EM-triggered and fiducial", n_min);
    if (n_truth_files > 0)
        std::cout << "UnderlyingEvent: " << n_truth_jets << " truth jets estimated from same-event off-axis cones, mean UE "
                  << truth_ue_sum / std::max(1L, n_truth_jets) << " GeV" << std::endl;
    return 0;
}
