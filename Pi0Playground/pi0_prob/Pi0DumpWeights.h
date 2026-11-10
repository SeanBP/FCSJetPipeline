#ifndef PI0_DUMP_WEIGHTS_H
#define PI0_DUMP_WEIGHTS_H

// Cross-section weight per file_idx of pi0_jet_dump output, for
// pi0_prob_derive / _validate / _datafit.
//
// Dumps made in the sim_to_jet job (pi0dump_*.root) carry sigma_job, N_job and
// the ptHat range in their "files" tree; the weights are then computed from it
// with the formula of jet_scale/JetXSecWeights.h (w_bin = sigma_bin /
// sum_jobs N_job, sigma_bin the N_job-weighted mean; normalised to the
// largest). Files with no valid sigma_job (starsim died early) get no weight
// and their jets are skipped. Older dumps (no sigma_job) fall back to reading
// the SimpleTrees named in the files tree, mapped by file name through
// `list` (scratch copies), with ComputeXSecWeights.

#include "JetXSecWeights.h"

#include <TChain.h>
#include <TString.h>

#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

inline bool Pi0DumpXSecWeights(const char *glob, const char *list, std::map<int, double> &wOf)
{
    wOf.clear();
    TChain chf("files");
    chf.Add(glob);
    if (chf.GetEntries() == 0) { fprintf(stderr, "no files tree in %s\n", glob); return false; }
    Int_t fi = -1;
    TString *p = 0;
    chf.SetBranchAddress("file_idx", &fi);
    chf.SetBranchAddress("path", &p);
    const bool own = chf.GetBranch("sigma_job") != 0;
    if (own)
    {
        Double_t sig = -1; Int_t nj = -1; Float_t tmin = -1, tmax = -1;
        chf.SetBranchAddress("sigma_job", &sig); chf.SetBranchAddress("N_job", &nj);
        chf.SetBranchAddress("pthat_min", &tmin); chf.SetBranchAddress("pthat_max", &tmax);
        typedef std::pair<float, float> Key;
        std::map<Key, std::pair<double, double> > bins;  // sum N, sum N sigma
        std::map<int, Key> keyOf;
        int bad = 0;
        for (Long64_t i = 0; i < chf.GetEntries(); i++)
        {
            chf.GetEntry(i);
            const double s = DecodePythia6Sigma(sig);
            if (!(s > 0) || nj <= 0) { bad++; continue; }
            if (keyOf.count(fi)) { fprintf(stderr, "duplicate file_idx %d in the dumps\n", fi); return false; }
            Key k(tmin, tmax);
            bins[k].first += nj; bins[k].second += (double)nj * s;
            keyOf[fi] = k;
        }
        if (bins.empty()) { fprintf(stderr, "no dump file with a valid sigma_job\n"); return false; }
        std::map<Key, double> wb;
        double wmax = 0;
        printf("cross-section weights from the dumps (%zu ptHat bins; %d files without sigma_job skipped)\n", bins.size(), bad);
        for (auto &b : bins) { const double w = (b.second.second / b.second.first) / b.second.first; wb[b.first] = w; if (w > wmax) wmax = w; }
        for (auto &b : bins) printf("  ptHat %g-%g  sum N %.0f  sigma %.4g mb  rel w %.4g\n", b.first.first, b.first.second, b.second.first, b.second.second / b.second.first, wb[b.first] / wmax);
        for (auto &k : keyOf) wOf[k.first] = wb[k.second] / wmax;
        return true;
    }
    std::map<int, std::string> pathOf;
    for (Long64_t i = 0; i < chf.GetEntries(); i++) { chf.GetEntry(i); pathOf[fi] = p->Data(); }
    if (list && *list && std::string(list) != "-")
    {
        std::map<std::string, std::string> orig;
        std::ifstream in(list); std::string l;
        while (std::getline(in, l)) if (!l.empty()) orig[l.substr(l.rfind('/') + 1)] = l;
        for (auto &kv : pathOf)
        {
            auto it = orig.find(kv.second.substr(kv.second.rfind('/') + 1));
            if (it != orig.end()) kv.second = it->second;
        }
    }
    std::vector<std::string> paths; std::vector<int> idx;
    for (auto &kv : pathOf) { idx.push_back(kv.first); paths.push_back(kv.second); }
    std::vector<float> wv;
    if (!ComputeXSecWeights(paths, "data", wv)) { fprintf(stderr, "cross-section weights failed\n"); return false; }
    for (size_t i = 0; i < idx.size(); i++) if (wv[i] > 0) wOf[idx[i]] = wv[i];
    return true;
}

#endif
