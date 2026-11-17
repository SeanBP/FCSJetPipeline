// Per-run A_N jet counts from an AnMerge output (TTree "runs", AnRunTools.h),
// the weights for pol_uncertainty.py (2026-11-17). Writes "run jets" lines.
// Usage: root -b -q -l 'dump_run_jets.C+("all_masked.root", "run_jets.txt")'
#include "TFile.h"
#include "TTree.h"
#include <cstdio>

void dump_run_jets(const char *in, const char *out)
{
    TFile f(in);
    TTree *t = (TTree *)f.Get("runs");
    if (!t) { printf("dump_run_jets: no runs tree in %s\n", in); return; }
    Int_t run = 0;
    Double_t jets = 0;
    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("jets", &jets);
    FILE *o = fopen(out, "w");
    fprintf(o, "# run A_N_jets (from %s, TTree runs)\n", in);
    double tot = 0;
    for (Long64_t i = 0; i < t->GetEntries(); i++) { t->GetEntry(i); fprintf(o, "%d %.0f\n", run, jets); tot += jets; }
    fclose(o);
    printf("dump_run_jets: %lld runs, %.0f jets -> %s\n", t->GetEntries(), tot, out);
}
