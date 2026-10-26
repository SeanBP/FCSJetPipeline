// Generator filter threshold check (FcsJetFilter -e): yield of EM-triggered,
// fiducial reco jets per sim_to_jet job, and per pb of cross section, for one
// ptHat bin made with two filter thresholds.
//
// Reads jet_output files (one list per sample). Per file: accepted events,
// trials N_job and sigma_job; per jet pT bin: EM-triggered (fcsEM0/1/2),
// fiducial reco jets. Prints per sample the jets per job and the
// cross-section-weighted rate sigma * jets / trials (mb), whose ratio between
// samples is the fraction of the true jet rate the higher threshold keeps.
//
// Build and run from this folder:
//   g++ -m64 -O2 filter_yield.cpp ../../shared/VectorDict.cxx -o filter_yield -I../../shared -I../../underlying_event \
//       -I../../jet_scale `root-config --cflags --libs` -std=c++11
//   ./filter_yield <label> <jet_output list> [<label> <list> ...]

#include <TFile.h>
#include <TTree.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "VectorDict.h"
#include "JetParameters.h"
#include "ProxySelection.h"
#include "JetXSecWeights.h"

static const int kNPt = 4;
static const double kPtEdges[kNPt + 1] = {2, 3, 4, 6, 10};

int main(int argc, char **argv)
{
    if (argc < 3 || argc % 2 == 0) { std::cerr << "usage: ./filter_yield <label> <list> [<label> <list> ...]" << std::endl; return 1; }
    std::printf("%-12s %5s %9s %12s | jets per job: %s | sigma x jets / trials (nb): %s\n", "sample", "files", "events", "trials/event",
                "2-3 3-4 4-6 6-10 GeV", "2-3 3-4 4-6 6-10 GeV");
    for (int a = 1; a + 1 < argc; a += 2)
    {
        std::ifstream in(argv[a + 1]);
        std::string fn;
        long n_files = 0, n_events = 0;
        double trials = 0, sigma_trials = 0, jets[kNPt] = {0, 0, 0, 0};
        while (in >> fn)
        {
            TFile f(fn.c_str());
            TTree *t = (TTree *)f.Get("jetTree");
            if (!t) continue;
            std::vector<float> *pt = 0, *x = 0, *y = 0;
            Int_t em0 = 0, em1 = 0, em2 = 0, N_job = -1;
            Double_t sigma_job = -1;
            t->SetBranchStatus("*", 0);
            const char *b[] = {"reco_pt", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2", "sigma_job", "N_job"};
            for (int i = 0; i < 8; i++) t->SetBranchStatus(b[i], 1);
            t->SetBranchAddress("reco_pt", &pt); t->SetBranchAddress("reco_x", &x); t->SetBranchAddress("reco_y", &y);
            t->SetBranchAddress("fcsEM0", &em0); t->SetBranchAddress("fcsEM1", &em1); t->SetBranchAddress("fcsEM2", &em2);
            t->SetBranchAddress("sigma_job", &sigma_job); t->SetBranchAddress("N_job", &N_job);
            const Long64_t n = t->GetEntries();
            if (n == 0) continue;
            t->GetEntry(0);
            const double sg = DecodePythia6Sigma(sigma_job);
            if (!(sg > 0) || N_job <= 0) { std::cerr << "no sigma/N_job in " << fn << std::endl; continue; }
            n_files++; trials += N_job; sigma_trials += sg * N_job;
            for (Long64_t e = 0; e < n; e++)
            {
                t->GetEntry(e);
                n_events++;
                if (!(em0 || em1 || em2)) continue;
                for (size_t j = 0; j < pt->size(); j++)
                {
                    if (!ProxyInFiducial(x->at(j), y->at(j))) continue;
                    for (int k = 0; k < kNPt; k++) if (pt->at(j) >= kPtEdges[k] && pt->at(j) < kPtEdges[k + 1]) jets[k]++;
                }
            }
        }
        if (n_files == 0) { std::printf("%-12s no usable files\n", argv[a]); continue; }
        // sigma x jets / trials, with sigma the trial-weighted mean of the files (mb -> nb)
        const double sigma = sigma_trials / trials;
        std::printf("%-12s %5ld %9ld %12.1f |", argv[a], n_files, n_events, trials / n_events);
        for (int k = 0; k < kNPt; k++) std::printf(" %7.1f", jets[k] / n_files);
        std::printf(" |");
        for (int k = 0; k < kNPt; k++) std::printf(" %9.3g", sigma * jets[k] / trials * 1e6);
        std::printf("\n");
    }
    return 0;
}
