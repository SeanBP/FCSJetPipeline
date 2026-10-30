// reemulate.C -- re-runs the FCS trigger emulation on simulation SimpleTrees
// (StFcsTrgReEmu) for a chosen run number, DB date and trigger version, and
// writes one friend file per input, <outdir>/<input basename>, with tree
// "trgemu" (fcsEM0-2 re-emulated, orig_fcsEM0-2 as stored).
//
// The DB date sets the gains used to rebuild each hit's ADC and the
// emulator's online gains; the run number sets StFcsDb's E_T gain mode
// (mode 1 for 20000000 < run < 23027048, else 2) and is passed to
// StFcsTriggerSimMaker::InitRun.
//
//   simulation's own settings (validation, must reproduce the stored flags):
//     root4star -b -q 'reemulate.C("sim.list", "out_sim", 1, 20230410, 0, 202209)'
//   run 23007048's settings:
//     root4star -b -q 'reemulate.C("sim.list", "out_23007048", 23007048, 20220107, 133502, 202209)'
// The last argument scales every ADC fed to the emulator (default 1).
//
// Needs stage_params.txt (the EM/JP thresholds) in the working directory, and
// "cons" run in this folder first.
void reemulate(const char *list, const char *outdir, int run = 1, int date = 20230410, int time = 0, int trgSelect = 202209, float adcScale = 1.0)
{
    gROOT->Macro("Load.C");
    gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    gSystem->Load("StFcsDbMaker");
    gSystem->Load("RTS");
    gSystem->Load("StFcsTriggerSimMaker");
    gSystem->Load("StFcsTrgReEmu");

    StChain *chain = new StChain("StChain");
    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
    dbMk->SetDateTime(date, time);
    StFcsDbMaker *fcsDbMk = new StFcsDbMaker();
    StFcsTriggerSimMaker *trg = new StFcsTriggerSimMaker();
    trg->setTrigger(trgSelect);
    trg->setSimMode(1);
    trg->setEtGain(1.0);
    trg->setThresholdFile((char *)"stage_params.txt");
    StFcsTrgReEmu *re = new StFcsTrgReEmu(trg);
    re->setAdcScale(adcScale);

    chain->Init();
    // no events drive InitRun here: load the run's FCS DB settings and the
    // emulator's per-channel gains by hand, in that order
    fcsDbMk->InitRun(run);
    trg->InitRun(run);
    printf("reemulate: run %d, DB date %d %06d, trgSelect %d, ADC scale %.3f\n", run, date, time, trgSelect, adcScale);

    gSystem->mkdir(outdir, kTRUE);
    TString all = gSystem->GetFromPipe(Form("cat %s", list));
    TObjArray *files = all.Tokenize("\n");
    long nfiles = 0, nev = 0;
    for (int i = 0; i < files->GetEntries(); i++)
    {
        TString f = ((TObjString *)files->At(i))->GetString();
        f = f.Strip(TString::kBoth);
        if (f.Length() == 0) continue;
        TString out = TString(outdir) + "/" + gSystem->BaseName(f.Data());
        int n = re->Process(f.Data(), out.Data());
        if (n >= 0) { nfiles++; nev += n; }
    }
    printf("reemulate: %ld files, %ld events -> %s\n", nfiles, nev, outdir);
}
