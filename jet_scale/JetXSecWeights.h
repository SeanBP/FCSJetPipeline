#ifndef JET_XSEC_WEIGHTS_H
#define JET_XSEC_WEIGHTS_H

// Per-file MC cross-section weights, shared by the JES steps.
//
// Files are grouped by ptHat bin (mc_pthatmin_gev, mc_pthatmax_gev) and every
// event in a bin gets w_bin = sigma_bin / sum_jobs(N_job), with sigma_bin the
// N_job-weighted mean of the bin's sigma_job (bins can hold different numbers
// of jobs). Exact only for disjoint bins; ComputeXSecWeights() warns about
// open-ended bins (ptcutmax = -1) that overlap others. Weights are normalised
// to the largest.

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <utility>
#include <cstring>

// Decodes a raw Pythia 6 sigma_job: StarPythia6::Stats() reads the cross
// section 4 bytes off, giving a denormal (~5e-315) whose low word is the true
// value's high word. readMudst.C already decodes before storing, so normal
// values are returned unchanged.
inline double DecodePythia6Sigma(double sigma_job)
{
    if (!(sigma_job > 0.0) || sigma_job > 1e-300)
        return sigma_job;

    unsigned long long bits = 0;
    std::memcpy(&bits, &sigma_job, sizeof(bits));
    bits = (bits & 0xffffffffULL) << 32;

    double decoded = 0.0;
    std::memcpy(&decoded, &bits, sizeof(decoded));
    return decoded;
}

struct XSecBin
{
    double sum_N;
    double sum_N_sigma;
    int n_files;

    XSecBin() : sum_N(0), sum_N_sigma(0), n_files(0) {}
};

// Fills weights[i] for root_files[i]; -1 marks a file to skip (unreadable,
// no tree, or zero entries -- the same files the main loop gets nothing
// from anyway). Returns false if any readable file lacks valid
// sigma_job/N_job: mixing weighted and unweighted files in one fit would be
// meaningless, so the caller should stop.
inline bool ComputeXSecWeights(const std::vector<std::string> &root_files,
                               const std::string &tree_name,
                               std::vector<float> &weights)
{
    typedef std::pair<float, float> PtHatKey;

    std::map<PtHatKey, XSecBin> bins;
    std::vector<PtHatKey> file_key(root_files.size());
    std::vector<bool> file_ok(root_files.size(), false);

    int n_empty = 0;
    int n_decoded = 0;

    for (size_t fidx = 0; fidx < root_files.size(); fidx++)
    {
        TFile *f = TFile::Open(root_files[fidx].c_str());
        if (!f || f->IsZombie())
            continue;

        TTree *tree = (TTree *)f->Get(tree_name.c_str());
        if (!tree)
        {
            f->Close();
            delete f;
            continue;
        }

        if (tree->GetEntries() == 0)
        {
            n_empty++;
            f->Close();
            delete f;
            continue;
        }

        if (!tree->GetBranch("sigma_job") || !tree->GetBranch("N_job") ||
            !tree->GetBranch("mc_pthatmin_gev") || !tree->GetBranch("mc_pthatmax_gev"))
        {
            std::cerr << "ERROR: " << root_files[fidx]
                      << " has no sigma_job/N_job/mc_pthatmin_gev/mc_pthatmax_gev branches;"
                      << " cannot cross-section weight it" << std::endl;
            f->Close();
            delete f;
            return false;
        }

        Double_t sigma_job = -1;
        Int_t N_job = -1;
        Float_t mc_pthatmin_gev = -1;
        Float_t mc_pthatmax_gev = -1;

        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("sigma_job", 1);
        tree->SetBranchStatus("N_job", 1);
        tree->SetBranchStatus("mc_pthatmin_gev", 1);
        tree->SetBranchStatus("mc_pthatmax_gev", 1);

        tree->SetBranchAddress("sigma_job", &sigma_job);
        tree->SetBranchAddress("N_job", &N_job);
        tree->SetBranchAddress("mc_pthatmin_gev", &mc_pthatmin_gev);
        tree->SetBranchAddress("mc_pthatmax_gev", &mc_pthatmax_gev);

        tree->GetEntry(0);

        f->Close();
        delete f;

        if (sigma_job > 0.0 && sigma_job <= 1e-300)
        {
            sigma_job = DecodePythia6Sigma(sigma_job);
            n_decoded++;
        }

        if (!(sigma_job > 0) || N_job <= 0)
        {
            std::cerr << "ERROR: " << root_files[fidx] << " has sigma_job=" << sigma_job
                      << ", N_job=" << N_job << " (unset sentinel);"
                      << " cannot cross-section weight it" << std::endl;
            return false;
        }

        PtHatKey key(mc_pthatmin_gev, mc_pthatmax_gev);

        XSecBin &b = bins[key];
        b.sum_N += N_job;
        b.sum_N_sigma += (double)N_job * sigma_job;
        b.n_files++;

        file_key[fidx] = key;
        file_ok[fidx] = true;
    }

    if (bins.empty())
    {
        std::cerr << "ERROR: no readable, non-empty JetTree files to weight" << std::endl;
        return false;
    }

    std::map<PtHatKey, double> bin_weight;
    double w_max = 0.0;

    for (std::map<PtHatKey, XSecBin>::const_iterator it = bins.begin(); it != bins.end(); ++it)
    {
        double sigma_bin = it->second.sum_N_sigma / it->second.sum_N;
        double w = sigma_bin / it->second.sum_N;

        bin_weight[it->first] = w;
        w_max = std::max(w_max, w);
    }

    std::cout << "Cross-section weights per ptHat bin (" << bins.size() << " bins):" << std::endl;
    std::cout << "  ptHatMin  ptHatMax  files  sum(N_job)  sigma (mb)  relative weight" << std::endl;

    for (std::map<PtHatKey, XSecBin>::const_iterator it = bins.begin(); it != bins.end(); ++it)
    {
        std::cout << "  " << it->first.first << "  " << it->first.second
                  << "  " << it->second.n_files
                  << "  " << it->second.sum_N
                  << "  " << it->second.sum_N_sigma / it->second.sum_N
                  << "  " << bin_weight[it->first] / w_max << std::endl;
    }

    if (n_decoded > 0)
        std::cout << "NOTE: sigma_job was decoded from the raw Pythia 6 framework value in " << n_decoded
                  << " file(s) (see DecodePythia6Sigma)" << std::endl;

    if (n_empty > 0)
        std::cout << "WARNING: " << n_empty << " file(s) have zero jetTree entries; their"
                  << " generated trials are not counted in sum(N_job)" << std::endl;

    // Overlap check: bins are sorted by ptHatMin, so it's enough to compare
    // each bin's upper edge (< 0 means no upper cut) with the next lower edge.
    for (std::map<PtHatKey, XSecBin>::const_iterator it = bins.begin(); it != bins.end(); ++it)
    {
        std::map<PtHatKey, XSecBin>::const_iterator nxt = it;
        ++nxt;
        if (nxt == bins.end())
            break;

        if (it->first.second < 0 || it->first.second > nxt->first.first)
        {
            std::cout << "WARNING: ptHat bins overlap (e.g. [" << it->first.first << ", "
                      << it->first.second << ") and [" << nxt->first.first << ", "
                      << nxt->first.second << "), -1 = no upper cut); the overlap region is"
                      << " counted once per sample covering it, so the weighted spectrum is"
                      << " not the physical one. Use disjoint bins (ptcutmax on every bin)."
                      << std::endl;
            break;
        }
    }

    weights.assign(root_files.size(), -1.0f);
    for (size_t fidx = 0; fidx < root_files.size(); fidx++)
        if (file_ok[fidx])
            weights[fidx] = (float)(bin_weight[file_key[fidx]] / w_max);

    return true;
}

// Kish effective sample size, (sum w)^2 / sum w^2: how many equal-weight
// entries would give the same statistical precision.
inline double EffectiveEntries(const std::vector<double> &w)
{
    double sw = 0.0, sw2 = 0.0;
    for (size_t i = 0; i < w.size(); i++)
    {
        sw += w[i];
        sw2 += w[i] * w[i];
    }
    return sw2 > 0.0 ? sw * sw / sw2 : 0.0;
}

#endif
