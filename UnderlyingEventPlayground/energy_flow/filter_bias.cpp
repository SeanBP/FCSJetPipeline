// Bias of a generator cone-pT filter (FcsJetFilter -q) on the jet energy
// scale, measured on an unfiltered sample (submit_sim_to_jet.sh -e 0).
//
// For every event, the largest summed pT of generator-level particles in an
// R = 0.7 cone within one FCS half (FcsJetFilter's box; SimpleTree mcpart_*
// from the primary vertex, as cone_filter_check.cpp) decides whether a cone
// filter at threshold X would have kept it. For every fiducial reco jet
// matched to a truth jet (match_reco_idx / match_truth_idx, parallel lists of
// pairs), binned in reco E and in truth E:
//   the fraction whose event fails X = 0.5 / 1.0 / 1.5 / 2.0 GeV
//   the mean response reco E / truth E, all jets and those passing each X
//   the median response and the fraction of jets with response > 2 (all jets)
// Weights: sigma / events per ptHat bin, so the sample mixes as in nature.
//
// Build and run from this folder:
//   g++ -m64 -O2 filter_bias.cpp ../../shared/VectorDict.cxx -o filter_bias -I../../shared -I../../underlying_event \
//       -I../../jet_scale `root-config --cflags --libs` -std=c++11
//   ./filter_bias <pairs.list "<jet_output> <SimpleTree>">

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "VectorDict.h"
#include "JetParameters.h"
#include "ProxySelection.h"
#include "JetXSecWeights.h"

static const float ZFCS = 710.16 + 13.9 + 15.0, XFCSMin = 16.69, XFCSMax = 16.69 + 22.0 * 5.542, YFCS = 34.0 / 2.0 * 5.542, FVCUT = 20.0;
static const int kNX = 4;
static const double kX[kNX] = {0.5, 1.0, 1.5, 2.0};
static const int kNE = 7;
static const double kE[kNE + 1] = {3, 5, 10, 15, 20, 30, 50, 100};

struct Acc
{
    double w, w_fail[kNX], wr, wr_pass[kNX], w_pass[kNX], n, w_out;
    std::vector<std::pair<double, double> > resp_w;
    Acc() : w(0), wr(0), n(0), w_out(0) { for (int k = 0; k < kNX; k++) { w_fail[k] = 0; wr_pass[k] = 0; w_pass[k] = 0; } }
};

int main(int argc, char **argv)
{
    if (argc < 2) { std::cerr << "usage: ./filter_bias <pairs.list>" << std::endl; return 1; }
    std::vector<std::pair<std::string, std::string> > pairs;
    {
        std::ifstream in(argv[1]);
        std::string a, b;
        while (in >> a >> b) pairs.push_back(std::make_pair(a, b));
    }
    // sigma / events per ptHat bin (mc_pthatmin_gev), from each file's sigma_job and N_job
    std::map<float, std::pair<double, double> > bin_N_sN;
    std::vector<float> file_bin(pairs.size(), -1);
    for (size_t f = 0; f < pairs.size(); f++)
    {
        TFile fj(pairs[f].first.c_str());
        TTree *t = (TTree *)fj.Get("jetTree");
        if (!t || t->GetEntries() == 0) continue;
        Double_t sg = -1; Int_t N = -1; Float_t ptmin = -1;
        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("sigma_job", 1); t->SetBranchStatus("N_job", 1); t->SetBranchStatus("mc_pthatmin_gev", 1);
        t->SetBranchAddress("sigma_job", &sg); t->SetBranchAddress("N_job", &N); t->SetBranchAddress("mc_pthatmin_gev", &ptmin);
        t->GetEntry(0);
        const double s = DecodePythia6Sigma(sg);
        if (!(s > 0) || N <= 0) continue;
        file_bin[f] = ptmin;
        bin_N_sN[ptmin].first += N; bin_N_sN[ptmin].second += N * s;
    }
    std::map<float, double> bin_w;
    double wmax = 0;
    for (std::map<float, std::pair<double, double> >::iterator it = bin_N_sN.begin(); it != bin_N_sN.end(); ++it)
    {
        bin_w[it->first] = it->second.second / it->second.first / it->second.first;
        wmax = std::max(wmax, bin_w[it->first]);
    }

    Acc by_reco[kNE], by_truth[kNE];
    const int MAX_MC = 10000;
    std::vector<Int_t> idVtx(MAX_MC);
    std::vector<Float_t> px(MAX_MC), py(MAX_MC), pz(MAX_MC), vx(MAX_MC), vy(MAX_MC), vz(MAX_MC);
    long n_ev = 0;
    for (size_t f = 0; f < pairs.size(); f++)
    {
        if (file_bin[f] < 0) continue;
        const double w = bin_w[file_bin[f]] / wmax;
        TFile fj(pairs[f].first.c_str()), fs(pairs[f].second.c_str());
        TTree *jt = (TTree *)fj.Get("jetTree"), *st = (TTree *)fs.Get("data");
        if (!jt || !st) continue;
        std::vector<float> *rE = 0, *rx = 0, *ry = 0, *tE = 0;
        std::vector<int> *mt = 0, *mr = 0;
        Long64_t entry = -1;
        Int_t nmc = 0;
        const bool has_entry = jt->GetBranch("simpletree_entry");
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_E", "reco_x", "reco_y", "truth_E", "match_truth_idx", "match_reco_idx"};
        for (int i = 0; i < 6; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_E", &rE); jt->SetBranchAddress("reco_x", &rx); jt->SetBranchAddress("reco_y", &ry);
        jt->SetBranchAddress("truth_E", &tE); jt->SetBranchAddress("match_truth_idx", &mt); jt->SetBranchAddress("match_reco_idx", &mr);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &entry); }
        st->SetBranchStatus("*", 0);
        const char *mb[] = {"mcpart_num", "mcpart_idVtx", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
        for (int i = 0; i < 8; i++) st->SetBranchStatus(mb[i], 1);
        st->SetBranchAddress("mcpart_num", &nmc); st->SetBranchAddress("mcpart_idVtx", &idVtx[0]);
        st->SetBranchAddress("mcpart_px", &px[0]); st->SetBranchAddress("mcpart_py", &py[0]); st->SetBranchAddress("mcpart_pz", &pz[0]);
        st->SetBranchAddress("mcpart_Vtx_x", &vx[0]); st->SetBranchAddress("mcpart_Vtx_y", &vy[0]); st->SetBranchAddress("mcpart_Vtx_z", &vz[0]);
        if (!has_entry && jt->GetEntries() != st->GetEntries()) continue;
        for (Long64_t e = 0; e < jt->GetEntries(); e++)
        {
            jt->GetEntry(e);
            if (rE->empty()) continue;
            if (st->GetEntry(has_entry ? entry : e) <= 0) continue;
            n_ev++;
            // event's largest same-half cone pT sum of generator particles in the filter box
            std::vector<double> g_eta, g_phi, g_pt; std::vector<int> g_half;
            for (int i = 0; i < std::min(nmc, MAX_MC); i++)
            {
                if (idVtx[i] != 1 || pz[i] <= 0) continue;
                const double fx = vx[i] + px[i] / pz[i] * (ZFCS - vz[i]);
                const double fy = std::fabs(vy[i] + py[i] / pz[i] * (ZFCS - vz[i]));
                if (std::fabs(fx) < XFCSMin + FVCUT || std::fabs(fx) > XFCSMax - FVCUT || fy > YFCS - FVCUT) continue;
                const double ptp = std::sqrt((double)px[i] * px[i] + (double)py[i] * py[i]);
                if (!(ptp > 0)) continue;
                g_eta.push_back(std::asinh(pz[i] / ptp)); g_phi.push_back(std::atan2(py[i], px[i])); g_pt.push_back(ptp); g_half.push_back(fx > 0);
            }
            double cmax = 0;
            for (size_t i = 0; i < g_pt.size(); i++)
            {
                double c = 0;
                for (size_t k = 0; k < g_pt.size(); k++)
                    if (g_half[k] == g_half[i] && ProxyDeltaR(g_eta[k], g_phi[k], g_eta[i], g_phi[i]) < R) c += g_pt[k];
                cmax = std::max(cmax, c);
            }
            for (size_t m = 0; m < mr->size() && m < mt->size(); m++)
            {
                const int j = mr->at(m), ti = mt->at(m);
                if (j < 0 || j >= (int)rE->size() || !ProxyInFiducial(rx->at(j), ry->at(j))) continue;
                if (ti < 0 || ti >= (int)tE->size() || !(tE->at(ti) > 0)) continue;
                const double resp = rE->at(j) / tE->at(ti);
                const double ev[2] = {rE->at(j), tE->at(ti)};
                for (int a = 0; a < 2; a++)
                {
                    int b = -1;
                    for (int k = 0; k < kNE; k++) if (ev[a] >= kE[k] && ev[a] < kE[k + 1]) b = k;
                    if (b < 0) continue;
                    Acc &A = a == 0 ? by_reco[b] : by_truth[b];
                    A.w += w; A.wr += w * resp; A.n++;
                    A.resp_w.push_back(std::make_pair(resp, w));
                    if (resp > 2) A.w_out += w;
                    for (int x = 0; x < kNX; x++)
                    {
                        if (cmax < kX[x]) A.w_fail[x] += w;
                        else { A.w_pass[x] += w; A.wr_pass[x] += w * resp; }
                    }
                }
            }
        }
    }
    std::printf("filter_bias: %zu pairs, %ld events with a reco jet\n", pairs.size(), n_ev);
    for (int a = 0; a < 2; a++)
    {
        std::printf("  binned in %s E (GeV): jets | fraction failing cone pT >= 0.5 / 1 / 1.5 / 2 GeV | mean reco/truth E: all, passing 0.5 / 1 / 1.5 / 2 | median, fraction > 2\n", a ? "truth" : "reco");
        for (int b = 0; b < kNE; b++)
        {
            Acc &A = a == 0 ? by_reco[b] : by_truth[b];
            if (A.n < 1) continue;
            std::sort(A.resp_w.begin(), A.resp_w.end());
            double cum = 0, med = 0;
            for (size_t i = 0; i < A.resp_w.size(); i++) { cum += A.resp_w[i].second; if (cum >= 0.5 * A.w) { med = A.resp_w[i].first; break; } }
            std::printf("    %3.0f-%-3.0f %7.0f |", kE[b], kE[b + 1], A.n);
            for (int x = 0; x < kNX; x++) std::printf(" %.4f", A.w_fail[x] / A.w);
            std::printf(" | %.4f", A.wr / A.w);
            for (int x = 0; x < kNX; x++) std::printf(" %.4f", A.w_pass[x] > 0 ? A.wr_pass[x] / A.w_pass[x] : 0.);
            std::printf(" | %.4f %.4f\n", med, A.w_out / A.w);
        }
    }
    return 0;
}
