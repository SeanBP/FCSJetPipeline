// Proxy-jet contamination against distance, data and simulation alike: in
// the no-separation proxy method, how the UE cone energy depends on the
// distance dR between the proxy jet and the cone.
//
// Fills ProxyPools as underlying_event/UnderlyingEvent.cpp does (pairs
// "<jet_output> <SimpleTree>", EM-triggered events, proxies matched in run or
// ptHat bin and EM level) and, for a sample of reco jets, takes their
// no-separation proxies (build with -DPROXY_MIN_DR_R=0; no veto). For each
// jet the proxies whose proxy jet is at dR >= 2.5R from the cone (where the
// pipeline's separation would put it) are the reference: same cone, same
// band, same pool. The excess of each proxy's cone energy over that
// reference mean, against its dR, is the contamination the separation
// removes, independent of the cone position. Jets with fewer than 10
// reference proxies are skipped. The excess is also split by whether the
// proxy jet is in the same FCS half (north x < 0, south x > 0) as the cone:
// FcsJetFilter (simulation only) requires > 50 GeV in one half, which puts
// extra energy in the proxy jet's half.
//
// Simulation (files ..._pt<ptHat>_<job>.root) is weighted per ptHat bin by
// sigma / events; each jet counts once (its proxies share its weight).
//
// Output: TProfiles per jet pT bin, p_excess_<pt> (cone E - reference, GeV)
// and p_cone_<pt> (cone E) against dR, plus printed tables.
//
// Build and run from this folder:
//   g++ -m64 -O2 -DPROXY_MIN_DR_R=0 -DPROXY_JET_VETO=0 proxy_dr_profile.cpp ../../shared/VectorDict.cxx -o proxy_dr_profile \
//       -I../../shared -I../../underlying_event -I../../jet_scale `root-config --cflags --libs` -std=c++11
//   ./proxy_dr_profile <pairs.list> <out.root> [max jets per group = 3000]

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "ProxySelection.h"
#include "ProxySearch.h"
#include "JetXSecWeights.h"

static const int kNPt = 4;
static const double kPtEdges[kNPt + 1] = {2, 3, 4, 6, 10};

struct Target { int event, jet, group; };

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::cerr << "usage: ./proxy_dr_profile <pairs.list> <out.root> [max jets per group]" << std::endl; return 1; }
    if (PROXY_MIN_DR_R != 0) { std::cerr << "build with -DPROXY_MIN_DR_R=0" << std::endl; return 1; }
    const int max_per_group = argc > 3 ? std::atoi(argv[3]) : 3000;
    const double ref_dr = 2.5 * R;

    std::vector<std::pair<std::string, std::string> > pairs;
    {
        std::ifstream in(argv[1]);
        std::string a, b;
        while (in >> a >> b) pairs.push_back(std::make_pair(a, b));
    }
    if (pairs.empty()) { std::cerr << "ERROR: no pairs" << std::endl; return 1; }

    // Group: run (data, jet_output_<run>_<file>) or -1-ptHatMin (simulation)
    ProxyPool pool;
    std::vector<Target> targets;
    std::map<int, double> group_N, group_sigmaN; // simulation: events and sum of N * sigma per group
    const int MAX_PTS = 200;
    for (size_t f = 0; f < pairs.size(); f++)
    {
        const std::string &jf = pairs[f].first;
        const std::string base = jf.substr(jf.find_last_of('/') + 1);
        int group = 0, a = 0, b = 0;
        if (std::sscanf(base.c_str(), "jet_output_%d_%d.root", &a, &b) == 2) group = a;
        else { size_t p = base.rfind("_pt"); if (p != std::string::npos && std::sscanf(base.c_str() + p, "_pt%d_%d.root", &a, &b) == 2) group = -1 - a; }

        TFile *fj = TFile::Open(jf.c_str());
        TFile *fs = TFile::Open(pairs[f].second.c_str());
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        const bool has_entry = jt && jt->GetBranch("simpletree_entry");
        if (!jt || !st || (!has_entry && jt->GetEntries() != st->GetEntries())) { std::cerr << "skip " << jf << std::endl; delete fj; delete fs; continue; }

        std::vector<float> *reco_pt = 0, *reco_eta = 0, *reco_phi = 0, *reco_x = 0, *reco_y = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        Long64_t simpletree_entry = -1;
        Double_t sqrt_s = kSqrtS_MainRun22_GeV, sigma_job = -1;
        Int_t N_job = -1;
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int i = 0; i < 8; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &reco_pt); jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi); jt->SetBranchAddress("reco_x", &reco_x); jt->SetBranchAddress("reco_y", &reco_y);
        jt->SetBranchAddress("fcsEM0", &em0); jt->SetBranchAddress("fcsEM1", &em1); jt->SetBranchAddress("fcsEM2", &em2);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &simpletree_entry); }
        if (jt->GetBranch("event_sqrt_s_gev")) { jt->SetBranchStatus("event_sqrt_s_gev", 1); jt->SetBranchAddress("event_sqrt_s_gev", &sqrt_s); }
        if (group < 0 && jt->GetBranch("sigma_job"))
        {
            jt->SetBranchStatus("sigma_job", 1); jt->SetBranchStatus("N_job", 1);
            jt->SetBranchAddress("sigma_job", &sigma_job); jt->SetBranchAddress("N_job", &N_job);
        }
        Int_t Cal_npts = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        Int_t Cal_pt_detid[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pb[i], 1);
        st->SetBranchAddress("Cal_npts", &Cal_npts); st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x); st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z); st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);

        bool counted = false;
        std::vector<ProxyJet> jets;
        std::vector<float> p_eta, p_phi, p_E;
        for (Long64_t ev = 0; ev < jt->GetEntries(); ev++)
        {
            if (jt->GetEntry(ev) <= 0 || reco_pt->empty()) continue;
            if (!counted && group < 0)
            {
                const double sg = DecodePythia6Sigma(sigma_job);
                if (sg > 0 && N_job > 0) { group_N[group] += N_job; group_sigmaN[group] += N_job * sg; }
                counted = true;
            }
            if (st->GetEntry(has_entry ? simpletree_entry : ev) <= 0) continue;
            jets.clear(); p_eta.clear(); p_phi.clear(); p_E.clear();
            for (size_t j = 0; j < reco_pt->size(); j++)
            {
                ProxyJet jet = {reco_pt->at(j), reco_eta->at(j), reco_phi->at(j), reco_x->at(j), reco_y->at(j)};
                jets.push_back(jet);
            }
            for (int i = 0; i < std::min((int)Cal_npts, MAX_PTS); i++)
            {
                if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                if (!pass_point_energy_cut(Cal_pt_energy[i], sqrt_s)) continue;
                const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                p_eta.push_back(std::asinh(Cal_pt_z[i] / rho)); p_phi.push_back(std::atan2(Cal_pt_y[i], Cal_pt_x[i]));
                p_E.push_back(Cal_pt_energy[i]);
            }
            const int e = pool.AddEvent(group, ProxyEmLevel(em0, em1, em2), jets, p_eta, p_phi, p_E);
            for (size_t j = 0; j < jets.size(); j++)
            {
                Target t = {e, (int)j, group};
                targets.push_back(t);
            }
        }
        fj->Close(); delete fj; fs->Close(); delete fs;
    }
    std::printf("proxy_dr_profile: %zu pairs, %d pool events, %zu candidate jets\n", pairs.size(), pool.NEvents(), targets.size());

    // Per-group weight (simulation: sigma / events, relative); sample up to max_per_group jets per group
    std::map<int, double> w_group;
    double wmax = 0;
    for (std::map<int, double>::iterator it = group_N.begin(); it != group_N.end(); ++it)
    {
        w_group[it->first] = group_sigmaN[it->first] / it->second / it->second;
        wmax = std::max(wmax, w_group[it->first]);
    }
    TRandom3 rng(4357);
    std::map<int, std::vector<int> > by_group;
    for (size_t i = 0; i < targets.size(); i++) by_group[targets[i].group].push_back((int)i);

    TFile out(argv[2], "RECREATE");
    const int nb = 12;
    const double lo = 0.6, hi = 3.0;
    TProfile *p_ex[kNPt], *p_cone[kNPt], *p_ex_half[2][kNPt];
    for (int p = 0; p < kNPt; p++)
    {
        p_ex[p] = new TProfile(Form("p_excess_%d", p), ";#DeltaR(proxy jet, UE cone);cone E - reference (GeV)", nb, lo, hi);
        for (int h = 0; h < 2; h++)
            p_ex_half[h][p] = new TProfile(Form("p_excess_%s_%d", h ? "opposite" : "same", p), ";#DeltaR(proxy jet, UE cone);cone E - reference (GeV)", nb, lo, hi);
        p_cone[p] = new TProfile(Form("p_cone_%d", p), ";#DeltaR(proxy jet, UE cone);cone E (GeV)", nb, lo, hi);
    }
    long n_used = 0;
    for (std::map<int, std::vector<int> >::iterator it = by_group.begin(); it != by_group.end(); ++it)
    {
        std::vector<int> &v = it->second;
        for (size_t i = v.size(); i > 1; i--) std::swap(v[i - 1], v[rng.Integer(i)]);
        const double wg = it->first < 0 ? (wmax > 0 && w_group.count(it->first) ? w_group[it->first] / wmax : 0) : 1.0;
        if (wg <= 0) continue;
        // sampling correction: each sampled jet stands for (all jets / sampled) of its group
        int n_est = 0, n_done = 0;
        std::vector<std::pair<int, ProxyResult> > done;
        for (size_t k = 0; k < v.size() && n_done < max_per_group; k++)
        {
            const Target &t = targets[v[k]];
            const ProxyJet &jet = pool.GetJet(pool.GetEvent(t.event).jet_first + t.jet);
            if (jet.pt < kPtEdges[0] || jet.pt >= kPtEdges[kNPt]) continue;
            ProxyResult res = pool.Estimate(t.event, t.jet, 500);
            if (res.n < 0) continue;
            n_est++;
            if (res.n <= 0) continue;
            done.push_back(std::make_pair(v[k], res));
            n_done++;
        }
        const double wsamp = wg; // group weights are per event; sampling equal fraction per group not needed for shapes within a pT bin
        for (size_t d = 0; d < done.size(); d++)
        {
            const Target &t = targets[done[d].first];
            const ProxyResult &res = done[d].second;
            const ProxyJet &jet = pool.GetJet(pool.GetEvent(t.event).jet_first + t.jet);
            int pb = -1;
            for (int p = 0; p < kNPt; p++) if (jet.pt >= kPtEdges[p] && jet.pt < kPtEdges[p + 1]) pb = p;
            std::vector<double> dr(res.n);
            double ref = 0; int nref = 0;
            for (int i = 0; i < res.n; i++)
            {
                const ProxyJet &pj = pool.GetJet(res.proxy_jet[i]);
                dr[i] = ProxyDeltaR(pj.eta, pj.phi, jet.eta, jet.phi);
                if (dr[i] >= ref_dr) { ref += res.E_list[i]; nref++; }
            }
            const double w = wsamp / res.n;
            for (int i = 0; i < res.n; i++) p_cone[pb]->Fill(dr[i], res.E_list[i], w);
            if (nref < 10) continue;
            ref /= nref;
            n_used++;
            for (int i = 0; i < res.n; i++)
            {
                p_ex[pb]->Fill(dr[i], res.E_list[i] - ref, w);
                // FCS half of the proxy jet against that of the cone (north x < 0, south x > 0)
                const bool same = (pool.GetJet(res.proxy_jet[i]).x > 0) == (jet.x > 0);
                p_ex_half[same ? 0 : 1][pb]->Fill(dr[i], res.E_list[i] - ref, w);
            }
        }
    }
    std::printf("proxy_dr_profile: %ld jets with >= 10 reference proxies (proxy jet at dR >= %.2f)\n", n_used, ref_dr);
    for (int p = 0; p < kNPt; p++)
    {
        std::printf("  jet pT %.0f-%.0f GeV: dR  excess (GeV) +- err  [cone E]  share of proxies | excess, proxy jet in same FCS half (share) | opposite half (share)\n", kPtEdges[p], kPtEdges[p + 1]);
        double tot = 0;
        for (int b = 1; b <= nb; b++) tot += p_cone[p]->GetBinEntries(b);
        for (int b = 1; b <= nb; b++)
        {
            if (p_cone[p]->GetBinEntries(b) <= 0) continue;
            const double ns = p_ex_half[0][p]->GetBinEntries(b), no = p_ex_half[1][p]->GetBinEntries(b);
            std::printf("    %.1f-%.1f  %+6.3f +- %.3f  [%.3f]  %.3f | %+6.3f (%.2f) | %+6.3f (%.2f)\n", p_ex[p]->GetBinLowEdge(b), p_ex[p]->GetBinLowEdge(b + 1),
                        p_ex[p]->GetBinContent(b), p_ex[p]->GetBinError(b), p_cone[p]->GetBinContent(b), p_cone[p]->GetBinEntries(b) / tot,
                        p_ex_half[0][p]->GetBinContent(b), ns + no > 0 ? ns / (ns + no) : 0., p_ex_half[1][p]->GetBinContent(b), ns + no > 0 ? no / (ns + no) : 0.);
        }
    }
    out.Write();
    out.Close();
    return 0;
}
