// Truth decomposition of the proxy-event UE cone (simulation): what the
// proxy method's UE cones contain, for the separation test's variants.
//
// Reads one ptHat bin of a sim_to_jet sample (pairs "<jet_output> <SimpleTree>")
// and fills ProxyPools as underlying_event/UnderlyingEvent.cpp does (reco jets
// and ECal points; one group; every event at EM level 0, so no trigger
// matching; reco jets from reco_jet_pt_min_em). Built with the default
// separation (2.5R) it runs two pools, with and without the jet veto
// (sep_veto, and sep = the pipeline method); built with -DPROXY_MIN_DR_R=0 it runs
// one (nosep: no separation, so no veto).
//
// Each event's EM-only truth photons (IsEmTruthParticle, in the truth
// particle eta band) are clustered into anti-kt R truth jets. For a sample of
// reco jets (pT >= reco_jet_pt_min_em, fiducial), each variant's proxy events
// (the nearest n_max) are revisited and the truth photon energy in the UE
// cone (R, at the jet's position) is split by the truth jet each photon
// belongs to:
//   t_proxy   the proxy jet's own truth jet (nearest truth jet axis within
//             R/2 of the proxy reco jet)
//   t_other2  other truth jets with pT >= 2 GeV
//   t_other1  other truth jets with 1 <= pT < 2 GeV
//   t_soft    truth jets below 1 GeV
// averaged over the proxy events (t_all = their sum). reco_cone is the
// method's own mean reco cone energy (ECal points); half_width the band
// half-width the jet's position gives (the window is sized to the 5% area).
// For reference, the jet's own event's truth cones at its eta and phi +- 90
// deg (the truth UE that jet_scale subtracts, as UnderlyingEvent.cpp) are
// split the same way, with the jet's own truth jet in place of the proxy
// jet's: own_all, own_jet, own_other2, own_other1, own_soft (mean of the two).
//
// Options (after n_max):
//   trig   the EM trigger level from the jet tree's fcsEM0/1/2 flags (emulated
//          in simulation), as in data: only EM-triggered events are proxies
//          or estimated, and proxies match the jet's level
//   mixed  the pairs of several ptHat bins in one pool, as data's natural mix:
//          each event is kept with probability (its bin's sigma / trials) /
//          (the largest), and the kept events are shuffled; every kept jet
//          then has the same weight
//
// Output: tree "jets" (one row per sampled jet and variant: variant 0
// sep_veto, 1 sep, 2 nosep) and tree "info" (the bin's ptHat range, cross
// section, trials and events, and how many jets were sampled of those that
// could be estimated, for weights: a sampled jet of the bin stands for
// sigma / trials x n_estimable / n_sampled).
//
// Usage: ./ue_truth_decomposition <pairs.list> <out.root> [max_jets=20000] [n_max=500] [trig] [mixed]
// Build (from this folder):
//   g++ -m64 -O2 ue_truth_decomposition.cpp ../../shared/VectorDict.cxx -o ue_truth_decomposition_sep \
//       -I../../shared -I../../underlying_event -I../../jet_scale `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//   (add -DPROXY_MIN_DR_R=0 -o ue_truth_decomposition_nosep for the no-separation build)

#include <TFile.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "ProxySearch.h"
#include "JetXSecWeights.h"

// Truth content of one pool event: photons with the truth jet each belongs to
struct TruthEvent
{
    int ph_first, ph_n, tj_first, tj_n;
};

// One event as read, before it enters the pools
struct EvRec
{
    std::vector<ProxyJet> jets;
    std::vector<float> p_eta, p_phi, p_E;
    int em;
    std::vector<float> ph_eta, ph_phi, ph_E, tj_pt, tj_eta, tj_phi;
    std::vector<int> ph_tj;
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3)
    {
        std::cerr << "usage: ./ue_truth_decomposition <pairs.list> <out.root> [max_jets=20000] [n_max=500] [trig] [mixed]" << std::endl;
        return 1;
    }
    const long max_jets = argc > 3 ? std::atol(argv[3]) : 20000;
    const int n_max = argc > 4 ? std::atoi(argv[4]) : 500;
    bool trig = false, mixed = false;
    for (int i = 5; i < argc; i++)
    {
        const std::string o = argv[i];
        if (o == "trig") trig = true;
        else if (o == "mixed") mixed = true;
        else { std::cerr << "unknown option " << o << std::endl; return 1; }
    }
    const bool with_sep = PROXY_MIN_DR_R > 0;
    const int n_var = with_sep ? 2 : 1;
    const int var_id[2] = {with_sep ? 0 : 2, 1};
    const char *var_name[3] = {"sep_veto", "sep", "nosep"};

    std::vector<std::pair<std::string, std::string> > pairs;
    {
        std::ifstream in(argv[1]);
        std::string a, b;
        while (in >> a >> b) pairs.push_back(std::make_pair(a, b));
    }
    if (pairs.empty()) { std::cerr << "ERROR: no pairs in " << argv[1] << std::endl; return 1; }

    ProxyPool pool[2];
    pool[0].SetJetVeto(with_sep);
    pool[1].SetJetVeto(false);

    // Truth content, indexed by pool event
    std::vector<TruthEvent> tev;
    std::vector<float> ph_eta, ph_phi, ph_E;
    std::vector<int> ph_tj;
    std::vector<float> tj_pt, tj_eta, tj_phi;

    float pthat_min = -1, pthat_max = -1;
    double sum_N = 0, sum_N_sigma = 0;
    long n_events = 0;

    // Mixed: each file's keep probability, from its bin's sigma / trials
    std::vector<double> keep_p(pairs.size(), 1.0);
    TRandom3 rng(4357);
    if (mixed)
    {
        typedef std::pair<float, float> Key;
        std::map<Key, std::pair<double, double> > bins; // sum N, sum N sigma
        std::vector<Key> key(pairs.size(), Key(-1, -1));
        for (size_t f = 0; f < pairs.size(); f++)
        {
            TFile fj(pairs[f].first.c_str());
            TTree *jt = (TTree *)fj.Get("jetTree");
            if (!jt || jt->GetEntries() == 0) continue;
            Double_t sigma_job = -1; Int_t N_job = -1; Float_t ptmin = -1, ptmax = -1;
            jt->SetBranchStatus("*", 0);
            const char *b[] = {"sigma_job", "N_job", "mc_pthatmin_gev", "mc_pthatmax_gev"};
            for (int i = 0; i < 4; i++) jt->SetBranchStatus(b[i], 1);
            jt->SetBranchAddress("sigma_job", &sigma_job); jt->SetBranchAddress("N_job", &N_job);
            jt->SetBranchAddress("mc_pthatmin_gev", &ptmin); jt->SetBranchAddress("mc_pthatmax_gev", &ptmax);
            jt->GetEntry(0);
            const double sg = DecodePythia6Sigma(sigma_job);
            key[f] = Key(ptmin, ptmax);
            if (sg > 0 && N_job > 0) { bins[key[f]].first += N_job; bins[key[f]].second += N_job * sg; }
        }
        double wmax = 0;
        std::map<Key, double> w;
        for (std::map<Key, std::pair<double, double> >::iterator it = bins.begin(); it != bins.end(); ++it)
        {
            w[it->first] = it->second.first > 0 ? it->second.second / it->second.first / it->second.first : 0;
            wmax = std::max(wmax, w[it->first]);
        }
        for (size_t f = 0; f < pairs.size(); f++) keep_p[f] = w.count(key[f]) && wmax > 0 ? w[key[f]] / wmax : 0;
        std::cout << "mixed: keep probability per ptHat bin:";
        for (std::map<Key, double>::iterator it = w.begin(); it != w.end(); ++it)
            std::cout << " [" << it->first.first << "," << it->first.second << ") " << it->second / wmax;
        std::cout << std::endl;
    }
    std::vector<EvRec> recs;

    const int MAX_MC = 10000, MAX_PTS = 200;
    std::vector<Int_t> mc_geid(MAX_MC), mc_idVtx(MAX_MC);
    std::vector<Float_t> mc_px(MAX_MC), mc_py(MAX_MC), mc_pz(MAX_MC), mc_E(MAX_MC), mc_vx(MAX_MC), mc_vy(MAX_MC), mc_vz(MAX_MC);

    for (size_t f = 0; f < pairs.size(); f++)
    {
        TFile *fj = TFile::Open(pairs[f].first.c_str());
        TFile *fs = TFile::Open(pairs[f].second.c_str());
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        const bool has_entry = jt && jt->GetBranch("simpletree_entry");
        if (!jt || !st || !st->GetBranch("Cal_npts") || !st->GetBranch("mcpart_num") || (!has_entry && jt->GetEntries() != st->GetEntries()))
        {
            std::cerr << "skip " << pairs[f].first << std::endl;
            if (fj) { fj->Close(); delete fj; }
            if (fs) { fs->Close(); delete fs; }
            continue;
        }
        std::vector<float> *reco_pt = 0, *reco_eta = 0, *reco_phi = 0, *reco_x = 0, *reco_y = 0;
        Long64_t simpletree_entry = -1;
        Double_t sqrt_s = kSqrtS_MainRun22_GeV, sigma_job = -1;
        Int_t N_job = -1, em0 = 0, em1 = 0, em2 = 0;
        Float_t ptmin = -1, ptmax = -1;
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_x", "reco_y"};
        for (int i = 0; i < 5; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &reco_pt);
        jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi);
        jt->SetBranchAddress("reco_x", &reco_x);
        jt->SetBranchAddress("reco_y", &reco_y);
        const char *ob[] = {"simpletree_entry", "event_sqrt_s_gev", "sigma_job", "N_job", "mc_pthatmin_gev", "mc_pthatmax_gev"};
        void *oa[] = {&simpletree_entry, &sqrt_s, &sigma_job, &N_job, &ptmin, &ptmax};
        for (int i = 0; i < 6; i++)
            if (jt->GetBranch(ob[i])) { jt->SetBranchStatus(ob[i], 1); jt->SetBranchAddress(ob[i], oa[i]); }
        if (trig)
        {
            const char *eb[] = {"fcsEM0", "fcsEM1", "fcsEM2"};
            Int_t *ea[] = {&em0, &em1, &em2};
            for (int i = 0; i < 3; i++) { jt->SetBranchStatus(eb[i], 1); jt->SetBranchAddress(eb[i], ea[i]); }
        }

        Int_t Cal_npts = 0, Cal_pt_detid[MAX_PTS], mcpart_num = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid", "mcpart_num", "mcpart_px",
                            "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_geid", "mcpart_idVtx", "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
        for (int i = 0; i < 16; i++) st->SetBranchStatus(pb[i], 1);
        st->SetBranchAddress("Cal_npts", &Cal_npts);
        st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x);
        st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z);
        st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);
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

        const Long64_t n = jt->GetEntries();
        if (n > 0)
        {
            jt->GetEntry(0);
            pthat_min = ptmin; pthat_max = ptmax;
            const double sg = DecodePythia6Sigma(sigma_job);
            if (sg > 0 && N_job > 0) { sum_N += N_job; sum_N_sigma += N_job * sg; }
        }
        n_events += n;
        std::vector<fastjet::PseudoJet> parts;
        for (Long64_t ev = 0; ev < n; ev++)
        {
            if (mixed && rng.Rndm() >= keep_p[f]) continue;
            if (jt->GetEntry(ev) <= 0 || reco_pt->empty()) continue;
            if (st->GetEntry(has_entry ? simpletree_entry : ev) <= 0) continue;
            EvRec rec;
            rec.em = trig ? ProxyEmLevel(em0, em1, em2) : 0;
            std::vector<ProxyJet> &jets = rec.jets;
            std::vector<float> &p_eta = rec.p_eta, &p_phi = rec.p_phi, &p_E = rec.p_E;
            for (size_t j = 0; j < reco_pt->size(); j++)
            {
                ProxyJet q = {reco_pt->at(j), reco_eta->at(j), reco_phi->at(j), reco_x->at(j), reco_y->at(j)};
                jets.push_back(q);
            }
            for (int i = 0; i < Cal_npts && i < MAX_PTS; i++)
            {
                if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                if (!pass_point_energy_cut(Cal_pt_energy[i], sqrt_s)) continue;
                const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                p_eta.push_back(std::asinh(Cal_pt_z[i] / rho)); p_phi.push_back(std::atan2(Cal_pt_y[i], Cal_pt_x[i]));
                p_E.push_back(Cal_pt_energy[i]);
            }

            // Truth photons and their truth jets
            parts.clear();
            float vx = 0, vy = 0, vz = 0;
            const int nmc = std::min((int)mcpart_num, MAX_MC);
            if (FindEmTruthPrimaryVertex(nmc, &mc_idVtx[0], &mc_vx[0], &mc_vy[0], &mc_vz[0], vx, vy, vz))
                for (int i = 0; i < nmc; i++)
                {
                    if (!IsEmTruthParticle(mc_geid[i], mc_vx[i], mc_vy[i], mc_vz[i], vx, vy, vz)) continue;
                    const double pt = std::sqrt((double)mc_px[i] * mc_px[i] + (double)mc_py[i] * mc_py[i]);
                    if (!(pt > 0) || !InTruthParticleEtaRange(std::asinh(mc_pz[i] / pt))) continue;
                    fastjet::PseudoJet p(mc_px[i], mc_py[i], mc_pz[i], mc_E[i]);
                    p.set_user_index((int)parts.size());
                    parts.push_back(p);
                }
            std::vector<int> owner(parts.size(), -1);
            if (!parts.empty())
            {
                fastjet::ClusterSequence cs(parts, fastjet::JetDefinition(fastjet::antikt_algorithm, R));
                const std::vector<fastjet::PseudoJet> tjs = cs.inclusive_jets(0.0);
                for (size_t k = 0; k < tjs.size(); k++)
                {
                    const std::vector<fastjet::PseudoJet> c = tjs[k].constituents();
                    for (size_t m = 0; m < c.size(); m++) owner[c[m].user_index()] = (int)k;
                    rec.tj_pt.push_back(tjs[k].perp()); rec.tj_eta.push_back(tjs[k].eta()); rec.tj_phi.push_back(tjs[k].phi_std());
                }
            }
            for (size_t i = 0; i < parts.size(); i++)
            {
                rec.ph_eta.push_back(parts[i].eta()); rec.ph_phi.push_back(parts[i].phi_std()); rec.ph_E.push_back(parts[i].E());
                rec.ph_tj.push_back(owner[i]);
            }
            recs.push_back(rec);
        }
        fj->Close(); delete fj;
        fs->Close(); delete fs;
    }
    // Into the pools: in file order, or shuffled (mixed)
    std::vector<size_t> order(recs.size());
    for (size_t i = 0; i < order.size(); i++) order[i] = i;
    if (mixed)
        for (size_t i = order.size(); i > 1; i--) std::swap(order[i - 1], order[rng.Integer((UInt_t)i)]);
    for (size_t o = 0; o < order.size(); o++)
    {
        EvRec &rec = recs[order[o]];
        int idx = -1;
        for (int v = 0; v < n_var; v++) idx = pool[v].AddEvent(0, rec.em, rec.jets, rec.p_eta, rec.p_phi, rec.p_E);
        if (idx != (int)tev.size()) { std::cerr << "ERROR: pool and truth event indices diverge" << std::endl; return 1; }
        TruthEvent te;
        te.ph_first = (int)ph_E.size(); te.ph_n = (int)rec.ph_E.size();
        te.tj_first = (int)tj_pt.size(); te.tj_n = (int)rec.tj_pt.size();
        ph_eta.insert(ph_eta.end(), rec.ph_eta.begin(), rec.ph_eta.end());
        ph_phi.insert(ph_phi.end(), rec.ph_phi.begin(), rec.ph_phi.end());
        ph_E.insert(ph_E.end(), rec.ph_E.begin(), rec.ph_E.end());
        ph_tj.insert(ph_tj.end(), rec.ph_tj.begin(), rec.ph_tj.end());
        tj_pt.insert(tj_pt.end(), rec.tj_pt.begin(), rec.tj_pt.end());
        tj_eta.insert(tj_eta.end(), rec.tj_eta.begin(), rec.tj_eta.end());
        tj_phi.insert(tj_phi.end(), rec.tj_phi.begin(), rec.tj_phi.end());
        tev.push_back(te);
        rec = EvRec(); // in the pools and truth arrays now
    }
    std::vector<EvRec>().swap(recs);
    std::cout << "ue_truth_decomposition: " << (trig ? "EM-trigger matching, " : "") << (mixed ? "mixed ptHat bins, " : "")
              << pairs.size() << " pairs, ptHat [" << pthat_min << ", " << pthat_max << "), "
              << pool[0].NEvents() << " pool events, " << ph_E.size() << " truth photons, " << tj_pt.size() << " truth jets" << std::endl;

    // Truth jet of event e nearest to (eta, phi) within R/2, -1 if none
    auto match_truth = [&](int e, double eta, double phi) {
        int best = -1;
        double bd = 0.5 * R;
        for (int k = 0; k < tev[e].tj_n; k++)
        {
            const int t = tev[e].tj_first + k;
            const double d = ProxyDeltaR(tj_eta[t], tj_phi[t], eta, phi);
            if (d < bd) { bd = d; best = k; }
        }
        return best;
    };
    // Truth photon energy of event e in an R cone at (eta, phi), split by owner: c[0] the given truth jet,
    // c[1] other jets >= 2 GeV, c[2] other jets 1-2 GeV, c[3] below 1 GeV
    auto split_cone = [&](int e, double eta, double phi, int own, double c[4]) {
        const TruthEvent &te = tev[e];
        for (int i = te.ph_first; i < te.ph_first + te.ph_n; i++)
        {
            if (ProxyDeltaR(ph_eta[i], ph_phi[i], eta, phi) >= R) continue;
            const int k = ph_tj[i];
            const double pt = k >= 0 ? tj_pt[te.tj_first + k] : 0.0;
            const int cat = (k >= 0 && k == own) ? 0 : (pt >= 2.0 ? 1 : (pt >= 1.0 ? 2 : 3));
            c[cat] += ph_E[i];
        }
    };

    // Jets to estimate: every estimable reco jet, then an even sample of max_jets
    std::vector<std::pair<int, int> > cand;
    for (int e = 0; e < pool[0].NEvents(); e++)
    {
        const ProxyPool::Event &pe = pool[0].GetEvent(e);
        for (int j = 0; j < pe.jet_n; j++)
        {
            const ProxyJet &q = pool[0].GetJet(pe.jet_first + j);
            if (pe.em_level >= 0 && q.pt >= reco_jet_pt_min_em && ProxyInFiducial(q.x, q.y)) cand.push_back(std::make_pair(e, j));
        }
    }
    const double stride = cand.size() > (size_t)max_jets && max_jets > 0 ? (double)cand.size() / max_jets : 1.0;

    TFile fout(argv[2], "RECREATE");
    TTree tj("jets", "proxy UE cone truth decomposition per jet and variant");
    Int_t variant, n_proxy;
    Float_t jet_pt, jet_eta, jet_E, reco_cone, t_all, t_proxy, t_other2, t_other1, t_soft, proxy_dphi, proxy_matched, half_width;
    Float_t own_all, own_jet, own_other2, own_other1, own_soft;
    tj.Branch("variant", &variant, "variant/I");
    tj.Branch("n_proxy", &n_proxy, "n_proxy/I");
    tj.Branch("jet_pt", &jet_pt, "jet_pt/F");
    tj.Branch("jet_eta", &jet_eta, "jet_eta/F");
    tj.Branch("jet_E", &jet_E, "jet_E/F");
    tj.Branch("reco_cone", &reco_cone, "reco_cone/F");
    tj.Branch("t_all", &t_all, "t_all/F");
    tj.Branch("t_proxy", &t_proxy, "t_proxy/F");
    tj.Branch("t_other2", &t_other2, "t_other2/F");
    tj.Branch("t_other1", &t_other1, "t_other1/F");
    tj.Branch("t_soft", &t_soft, "t_soft/F");
    tj.Branch("proxy_dphi", &proxy_dphi, "proxy_dphi/F");
    tj.Branch("proxy_matched", &proxy_matched, "proxy_matched/F");
    tj.Branch("half_width", &half_width, "half_width/F"); // the band half-width the jet's position gives (deg)
    tj.Branch("own_all", &own_all, "own_all/F");
    tj.Branch("own_jet", &own_jet, "own_jet/F");
    tj.Branch("own_other2", &own_other2, "own_other2/F");
    tj.Branch("own_other1", &own_other1, "own_other1/F");
    tj.Branch("own_soft", &own_soft, "own_soft/F");

    long n_done[2] = {0, 0}, n_sampled = 0;
    for (double s = 0; s < cand.size(); s += stride)
    {
        n_sampled++;
        const int e = cand[(size_t)s].first, j = cand[(size_t)s].second;
        const ProxyJet &q = pool[0].GetJet(pool[0].GetEvent(e).jet_first + j);
        jet_pt = q.pt; jet_eta = q.eta; jet_E = q.pt * std::cosh(q.eta);
        // The jet's own event: truth cones at phi +- 90 deg, the jet's own truth jet as "own"
        double oc[4] = {0, 0, 0, 0};
        const int own_tj = match_truth(e, q.eta, q.phi);
        for (int c = 0; c < 2; c++) split_cone(e, q.eta, q.phi + (c == 0 ? 0.5 : -0.5) * M_PI, own_tj, oc);
        own_jet = 0.5 * oc[0]; own_other2 = 0.5 * oc[1]; own_other1 = 0.5 * oc[2]; own_soft = 0.5 * oc[3];
        own_all = own_jet + own_other2 + own_other1 + own_soft;
        for (int v = 0; v < n_var; v++)
        {
            const ProxyResult res = pool[v].Estimate(e, j, n_max);
            if (res.n <= 0) continue;
            double c[4] = {0, 0, 0, 0};
            int matched = 0;
            for (int i = 0; i < res.n; i++)
            {
                const int pe = res.proxy_event[i];
                const ProxyJet &pj = pool[v].GetJet(res.proxy_jet[i]);
                const int ptj = match_truth(pe, pj.eta, pj.phi);
                if (ptj >= 0) matched++;
                split_cone(pe, q.eta, q.phi, ptj, c);
            }
            variant = var_id[v]; n_proxy = res.n;
            reco_cone = res.ue_E;
            t_proxy = c[0] / res.n; t_other2 = c[1] / res.n; t_other1 = c[2] / res.n; t_soft = c[3] / res.n;
            t_all = t_proxy + t_other2 + t_other1 + t_soft;
            proxy_dphi = res.jet_dphi; proxy_matched = (float)matched / res.n; half_width = res.half_width;
            tj.Fill();
            n_done[v]++;
        }
    }
    for (int v = 0; v < n_var; v++)
        std::cout << "  " << var_name[var_id[v]] << ": " << n_done[v] << " jets estimated (of " << cand.size() << " estimable)" << std::endl;

    TTree ti("info", "ptHat bin");
    // Mixed: the kept events already follow the cross section, so every jet has the same weight
    Float_t lo = mixed ? -1 : pthat_min, hi = mixed ? -1 : pthat_max;
    Double_t sigma = mixed ? 1.0 : (sum_N > 0 ? sum_N_sigma / sum_N : -1), trials = mixed ? 1.0 : sum_N;
    Long64_t nev = n_events, n_estimable = (Long64_t)cand.size(), n_samp = n_sampled;
    ti.Branch("ptHatMin", &lo, "ptHatMin/F");
    ti.Branch("ptHatMax", &hi, "ptHatMax/F");
    ti.Branch("sigma_mb", &sigma, "sigma_mb/D");
    ti.Branch("trials", &trials, "trials/D");
    ti.Branch("n_events", &nev, "n_events/L");
    ti.Branch("n_estimable", &n_estimable, "n_estimable/L"); // reco jets that could be estimated
    ti.Branch("n_sampled", &n_samp, "n_sampled/L");          // of them sampled: each stands for n_estimable / n_sampled
    ti.Fill();
    fout.Write();
    fout.Close();
    return 0;
}
