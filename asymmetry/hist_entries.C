// Prints "ENTRIES <n>": the summed entries of every histogram in the files
// matching a pattern (TChain-style wildcard), as a check that a hadd merge
// kept everything. Used by merge_hists_job.xml.
// Usage: root -l -b -q 'hist_entries.C("part_*.root")'
//    or: env HE_PATTERN='part_*.root' root -l -b -q hist_entries.C
// (the second avoids nested quotes in tcsh)

void hist_entries(const char *pattern = 0)
{
    if (!pattern) pattern = gSystem->Getenv("HE_PATTERN");
    if (!pattern) { printf("no pattern\n"); return; }
    TString dir = gSystem->DirName(pattern), base = gSystem->BaseName(pattern);
    if (dir == "") dir = ".";
    TRegexp re(base, kTRUE);
    void *d = gSystem->OpenDirectory(dir);
    if (!d) return;
    double tot = 0;
    int nf = 0;
    const char *e;
    while ((e = gSystem->GetDirEntry(d)))
    {
        TString n = e;
        Ssiz_t len = 0;
        if (re.Index(n, &len) != 0 || len != n.Length()) continue;
        TFile *f = TFile::Open(dir + "/" + n);
        if (!f || f->IsZombie()) { printf("BAD %s\n", n.Data()); return; }
        nf++;
        TIter it(f->GetListOfKeys());
        TKey *k;
        while ((k = (TKey *)it()))
        {
            if (!TString(k->GetClassName()).BeginsWith("TH")) continue;
            TH1 *h = (TH1 *)k->ReadObj();
            tot += h->GetEntries();
            delete h;
        }
        f->Close();
    }
    gSystem->FreeDirectory(d);
    if (nf > 0) printf("ENTRIES %.0f\n", tot);
}
