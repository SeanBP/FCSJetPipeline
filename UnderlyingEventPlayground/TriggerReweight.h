// Event weights that give simulation the data's EM trigger mix. Events are
// classed by the highest EM trigger fired (fcsEM2 > fcsEM1 > fcsEM0), as the
// proxy selection does: data flags are recorded after prescales, so a data
// event above the EM2 threshold usually carries only EM2, while simulation
// (no prescales, nested thresholds) carries all three. Over all EM-triggered
// events, f(l) is the fraction at highest level l (simulation cross-section
// weighted); a simulation event at level l is weighted by f_data(l) / f_sim(l)
// on top of its cross-section weight.
#ifndef TRIGGER_REWEIGHT_H
#define TRIGGER_REWEIGHT_H

#include <TFile.h>
#include <TTree.h>

#include <cstdio>
#include <string>
#include <vector>

static const int kNTrigLevel = 4; // 0 = untriggered, 1-3 = highest of EM0/EM1/EM2

inline int TriggerLevel(int em0, int em1, int em2) { return em2 ? 3 : em1 ? 2 : em0 ? 1 : 0; }

// Weighted fraction of EM-triggered events at each highest level
inline std::vector<double> TriggerFractions(const std::vector<std::string> &files, const std::vector<float> &weight)
{
    std::vector<double> sum(kNTrigLevel, 0.0);
    for (size_t f = 0; f < files.size(); f++)
    {
        if (weight[f] < 0) continue;
        TFile *fin = TFile::Open(files[f].c_str());
        TTree *t = (fin && !fin->IsZombie()) ? (TTree *)fin->Get("jetTree") : 0;
        if (t)
        {
            Int_t em0 = 0, em1 = 0, em2 = 0;
            t->SetBranchStatus("*", 0);
            t->SetBranchStatus("fcsEM0", 1); t->SetBranchStatus("fcsEM1", 1); t->SetBranchStatus("fcsEM2", 1);
            t->SetBranchAddress("fcsEM0", &em0); t->SetBranchAddress("fcsEM1", &em1); t->SetBranchAddress("fcsEM2", &em2);
            for (Long64_t ev = 0; ev < t->GetEntries(); ev++)
            {
                t->GetEntry(ev);
                sum[TriggerLevel(em0, em1, em2)] += weight[f];
            }
        }
        if (fin) { fin->Close(); delete fin; }
    }
    double tot = 0;
    for (int l = 1; l < kNTrigLevel; l++) tot += sum[l];
    std::vector<double> frac(kNTrigLevel, 0.0);
    for (int l = 1; l < kNTrigLevel; l++) frac[l] = tot > 0 ? sum[l] / tot : 0.0;
    return frac;
}

// Per-level simulation weight factors f_data / f_sim
inline std::vector<double> TriggerReweightFactors(const std::vector<double> &f_sim, const std::vector<double> &f_data)
{
    std::vector<double> w(kNTrigLevel, 0.0);
    const char *name[kNTrigLevel] = {"", "EM0", "EM1", "EM2"};
    printf("Highest EM trigger  sim fraction  data fraction  sim weight factor\n");
    for (int l = 1; l < kNTrigLevel; l++)
    {
        w[l] = f_sim[l] > 0 ? f_data[l] / f_sim[l] : 0.0;
        printf("  %-16s %12.4f  %13.4f  %17.3f\n", name[l], f_sim[l], f_data[l], w[l]);
    }
    return w;
}

#endif
