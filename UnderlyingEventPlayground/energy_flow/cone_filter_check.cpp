// Emulates generator filters offline on sim_to_jet SimpleTrees, to choose a
// FcsJetFilter that does not bias ECAL-only jets near the analysis threshold.
//
// Generator-level particles: SimpleTree mcpart_* from the primary vertex
// (idVtx == 1), projected to the FCS face with FcsJetFilter's box (ECAL
// columns, 20 cm fiducial margin). Per event:
//   half  the current filter's quantity, the energy summed per FCS half
//   cone  the largest energy in an R = 0.7 (eta, phi) cone around any of the
//         particles, per half
// Per EM-triggered fiducial reco jet (pT bins): the generator energy in an
// R = 0.7 cone around the jet axis, in its half. A cone filter at E_min keeps
// every event with such a jet whose cone energy exceeds E_min.
//
// Prints the per-jet cone-energy quantiles (unbiasedness) and the fraction of
// events passing cone >= E_min and half >= E (cost).
//
// Build and run from this folder:
//   g++ -m64 -O2 cone_filter_check.cpp ../../shared/VectorDict.cxx -o cone_filter_check -I../../shared -I../../underlying_event \
//       `root-config --cflags --libs` -std=c++11
//   ./cone_filter_check <pairs.list "<jet_output> <SimpleTree>"> [label]

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "VectorDict.h"
#include "JetParameters.h"
#include "ProxySelection.h"

static const float ZFCS = 710.16 + 13.9 + 15.0, XFCSMin = 16.69, XFCSMax = 16.69 + 22.0 * 5.542, YFCS = 34.0 / 2.0 * 5.542, FVCUT = 20.0;
static const int kNPt = 4;
static const double kPtEdges[kNPt + 1] = {2, 3, 4, 6, 10};

int main(int argc, char **argv)
{
    if (argc < 2) { std::cerr << "usage: ./cone_filter_check <pairs.list> [label]" << std::endl; return 1; }
    const std::string label = argc > 2 ? argv[2] : "";
    std::ifstream in(argv[1]);
    std::string jf, sf;
    const int MAX_MC = 10000;
    std::vector<Int_t> idVtx(MAX_MC);
    std::vector<Float_t> px(MAX_MC), py(MAX_MC), pz(MAX_MC), E(MAX_MC), vx(MAX_MC), vy(MAX_MC), vz(MAX_MC);
    std::vector<double> jet_cone[kNPt], jet_cone_pt[kNPt];
    std::vector<double> ev_cone, ev_half, ev_cone_pt;
    while (in >> jf >> sf)
    {
        TFile fj(jf.c_str()), fs(sf.c_str());
        TTree *jt = (TTree *)fj.Get("jetTree"), *st = (TTree *)fs.Get("data");
        if (!jt || !st) continue;
        std::vector<float> *pt = 0, *eta = 0, *phi = 0, *x = 0, *y = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0, nmc = 0;
        Long64_t entry = -1;
        const bool has_entry = jt->GetBranch("simpletree_entry");
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int i = 0; i < 8; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &pt); jt->SetBranchAddress("reco_eta", &eta); jt->SetBranchAddress("reco_phi", &phi);
        jt->SetBranchAddress("reco_x", &x); jt->SetBranchAddress("reco_y", &y);
        jt->SetBranchAddress("fcsEM0", &em0); jt->SetBranchAddress("fcsEM1", &em1); jt->SetBranchAddress("fcsEM2", &em2);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &entry); }
        st->SetBranchStatus("*", 0);
        const char *mb[] = {"mcpart_num", "mcpart_idVtx", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
        for (int i = 0; i < 9; i++) st->SetBranchStatus(mb[i], 1);
        st->SetBranchAddress("mcpart_num", &nmc); st->SetBranchAddress("mcpart_idVtx", &idVtx[0]);
        st->SetBranchAddress("mcpart_px", &px[0]); st->SetBranchAddress("mcpart_py", &py[0]); st->SetBranchAddress("mcpart_pz", &pz[0]);
        st->SetBranchAddress("mcpart_E", &E[0]); st->SetBranchAddress("mcpart_Vtx_x", &vx[0]); st->SetBranchAddress("mcpart_Vtx_y", &vy[0]);
        st->SetBranchAddress("mcpart_Vtx_z", &vz[0]);
        if (!has_entry && jt->GetEntries() != st->GetEntries()) continue;

        // every SimpleTree event counts for the event-level quantities; jets need the jetTree
        std::vector<Long64_t> jet_entry_of_st(st->GetEntries(), -1);
        for (Long64_t e = 0; e < jt->GetEntries(); e++) { jt->GetEntry(e); const Long64_t s = has_entry ? entry : e; if (s >= 0 && s < (Long64_t)jet_entry_of_st.size()) jet_entry_of_st[s] = e; }
        for (Long64_t s = 0; s < st->GetEntries(); s++)
        {
            st->GetEntry(s);
            std::vector<double> g_eta, g_phi, g_E, g_pt; std::vector<int> g_half;
            double half[2] = {0, 0};
            for (int i = 0; i < std::min(nmc, MAX_MC); i++)
            {
                if (idVtx[i] != 1 || pz[i] <= 0) continue;
                const double fx = vx[i] + px[i] / pz[i] * (ZFCS - vz[i]);
                const double fy = std::fabs(vy[i] + py[i] / pz[i] * (ZFCS - vz[i]));
                if (std::fabs(fx) < XFCSMin + FVCUT || std::fabs(fx) > XFCSMax - FVCUT || fy > YFCS - FVCUT) continue;
                const int h = fx > 0 ? 1 : 0;
                half[h] += E[i];
                const double ptp = std::sqrt((double)px[i] * px[i] + (double)py[i] * py[i]);
                g_eta.push_back(std::asinh(pz[i] / ptp)); g_phi.push_back(std::atan2(py[i], px[i])); g_E.push_back(E[i]); g_pt.push_back(ptp); g_half.push_back(h);
            }
            double cmax = 0, cmax_pt = 0;
            for (size_t i = 0; i < g_E.size(); i++)
            {
                double c = 0, cp = 0;
                for (size_t k = 0; k < g_E.size(); k++)
                    if (g_half[k] == g_half[i] && ProxyDeltaR(g_eta[k], g_phi[k], g_eta[i], g_phi[i]) < R) { c += g_E[k]; cp += g_pt[k]; }
                cmax = std::max(cmax, c); cmax_pt = std::max(cmax_pt, cp);
            }
            ev_cone.push_back(cmax); ev_cone_pt.push_back(cmax_pt); ev_half.push_back(std::max(half[0], half[1]));
            const Long64_t e = jet_entry_of_st[s];
            if (e < 0) continue;
            jt->GetEntry(e);
            if (!(em0 || em1 || em2)) continue;
            for (size_t j = 0; j < pt->size(); j++)
            {
                if (!ProxyInFiducial(x->at(j), y->at(j))) continue;
                int b = -1;
                for (int k = 0; k < kNPt; k++) if (pt->at(j) >= kPtEdges[k] && pt->at(j) < kPtEdges[k + 1]) b = k;
                if (b < 0) continue;
                const int h = x->at(j) > 0 ? 1 : 0;
                double c = 0, cp = 0;
                for (size_t k = 0; k < g_E.size(); k++)
                    if (g_half[k] == h && ProxyDeltaR(g_eta[k], g_phi[k], eta->at(j), phi->at(j)) < R) { c += g_E[k]; cp += g_pt[k]; }
                jet_cone[b].push_back(c); jet_cone_pt[b].push_back(cp);
            }
        }
    }
    std::printf("cone_filter_check %s: %zu events\n", label.c_str(), ev_cone.size());
    std::printf("  generator cone energy (R=0.7, in the filter box) around triggered fiducial reco jets, GeV: quantiles 0.1%% 0.5%% 1%% 5%% 50%%\n");
    for (int b = 0; b < kNPt; b++)
    {
        std::vector<double> &v = jet_cone[b];
        if (v.empty()) continue;
        std::sort(v.begin(), v.end());
        const double q[5] = {0.001, 0.005, 0.01, 0.05, 0.5};
        std::printf("    jet pT %2.0f-%-2.0f (%6zu jets):", kPtEdges[b], kPtEdges[b + 1], v.size());
        for (int k = 0; k < 5; k++) std::printf(" %6.2f", v[std::min(v.size() - 1, (size_t)(q[k] * v.size()))]);
        std::printf("   below 6/8/10 GeV:");
        for (double t : {6.0, 8.0, 10.0}) std::printf(" %.4f", (double)(std::lower_bound(v.begin(), v.end(), t) - v.begin()) / v.size());
        std::printf("\n");
    }
    std::printf("  generator cone pT sum around the same jets, GeV: quantiles 0.1%% 0.5%% 1%% 5%% 50%%\n");
    for (int b = 0; b < kNPt; b++)
    {
        std::vector<double> &v = jet_cone_pt[b];
        if (v.empty()) continue;
        std::sort(v.begin(), v.end());
        const double q[5] = {0.001, 0.005, 0.01, 0.05, 0.5};
        std::printf("    jet pT %2.0f-%-2.0f (%6zu jets):", kPtEdges[b], kPtEdges[b + 1], v.size());
        for (int k = 0; k < 5; k++) std::printf(" %6.2f", v[std::min(v.size() - 1, (size_t)(q[k] * v.size()))]);
        std::printf("   below 1/1.5/2 GeV:");
        for (double t : {1.0, 1.5, 2.0}) std::printf(" %.4f", (double)(std::lower_bound(v.begin(), v.end(), t) - v.begin()) / v.size());
        std::printf("\n");
    }
    std::printf("  events passing cone pT >= 1 / 1.5 / 2 / 3 GeV:");
    for (double t : {1.0, 1.5, 2.0, 3.0}) std::printf(" %.4f", (double)std::count_if(ev_cone_pt.begin(), ev_cone_pt.end(), [&](double c) { return c >= t; }) / ev_cone_pt.size());
    std::printf("\n");
    std::printf("  events passing: cone >= 6 / 8 / 10 / 15 GeV | half >= 10 / 20 / 30 / 50 GeV\n   ");
    const double ct[4] = {6, 8, 10, 15}, ht[4] = {10, 20, 30, 50};
    for (int k = 0; k < 4; k++) std::printf(" %.4f", (double)std::count_if(ev_cone.begin(), ev_cone.end(), [&](double c) { return c >= ct[k]; }) / ev_cone.size());
    std::printf(" |");
    for (int k = 0; k < 4; k++) std::printf(" %.4f", (double)std::count_if(ev_half.begin(), ev_half.end(), [&](double c) { return c >= ht[k]; }) / ev_half.size());
    std::printf("\n");
    return 0;
}
