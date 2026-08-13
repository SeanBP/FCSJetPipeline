// Minimal, standalone example of retrieving the actual per-event spin
// state -- not just whether the DB covers a run (see
// run_number_to_spin_coverage.sh for that), but the real spin4bits value
// StSpinDbMaker resolves for each individual event, exactly the way
// production's StSimpleReaderMaker::Make() does it (see its Spin_config
// branch, FCSJetPipeline/StRoot/StSpinPool/StSimpleReaderMaker/).
//
// Per event:
//   1. Get the run number from StMuEvent, and (re-)InitRun() the
//      StSpinDbMaker only when the run number changes (InitRun() does a
//      real DB query -- no need to repeat it every event).
//   2. Get this event's 7-bit bunch crossing id from StEvent's raw
//      StL0Trigger (StMuEvent::l0Trigger().bunchCrossingId7bit(runNumber)).
//   3. If isValid() (all 3 spin DB tables resolved for this run --
//      doesn't by itself guarantee REAL non-fallback data, see
//      run_number_to_spin_coverage.sh's header for the known stale-entry
//      traps), look up spin4usingBX7(bx7): a combined 0-10 code --
//      0 = both beams unpolarized/empty, then bit0=yellow(1=up,0=down
//      when set), bit1=... (see StSpinDbMaker::optimizeTables() for the
//      exact legend) -- otherwise -1, matching production's "no data"
//      sentinel.
//
// Usage:
//   root4star -b -q read_spin_state.C\(\"<path-or-xrootd-URL-to-a-muDst>\",maxEvents\)
//
// No local muDst is shipped here -- the default below streams a known-good
// real Run 22 file (run 22359013, which run_number_to_spin_coverage.sh
// confirms has genuine spin coverage) directly via xrootd.

void read_spin_state(
    const char* file =
        "root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/"
        "production_pp500_2022/ReversedFullField/P24ia/2021/359/22359013/"
        "st_fwd_22359013_raw_1500006.MuDst.root",
    int maxEvents = 20){

    gROOT->Macro("Load.C");
    gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    gSystem->Load("StSpinDbMaker");
    gMessMgr->SetLimit("I", 0);
    gMessMgr->SetLimit("Q", 0);
    gMessMgr->SetLimit("W", 0);

    StChain* chain = new StChain("StChain");
    StMuDstMaker* muDstMaker = new StMuDstMaker(0, 0, "", file, ".", 1000, "MuDst");

    // StSpinDbMaker needs a St_db_Maker in the chain to actually query the
    // DB (it asserts if none is found) -- blacklist is just to keep Init()
    // fast, none of these subsystems' calibrations are needed here.
    St_db_Maker* dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
    dbMk->SetAttr("blacklist", "tpc");
    dbMk->SetAttr("blacklist", "svt");
    dbMk->SetAttr("blacklist", "ssd");
    dbMk->SetAttr("blacklist", "ist");
    dbMk->SetAttr("blacklist", "pxl");
    dbMk->SetAttr("blacklist", "pp2pp");
    dbMk->SetAttr("blacklist", "ftpc");
    dbMk->SetAttr("blacklist", "emc");
    dbMk->SetAttr("blacklist", "eemc");
    dbMk->SetAttr("blacklist", "mtd");
    dbMk->SetAttr("blacklist", "pmd");
    dbMk->SetAttr("blacklist", "tof");
    dbMk->SetAttr("blacklist", "etof");
    dbMk->SetAttr("blacklist", "rhicf");
    dbMk->SetAttr("blacklist", "fcs");

    StSpinDbMaker* spinDb = new StSpinDbMaker();

    chain->Init();

    int n = muDstMaker->tree()->GetEntries();
    if (maxEvents >= 0 && maxEvents < n) n = maxEvents;
    printf("Found %d entries in MuDst, reading %d\n", (int)muDstMaker->tree()->GetEntries(), n);

    int lastRunNumber = -1;

    printf("\n%8s %10s %6s %8s %6s %s\n", "event", "run", "bx7", "isValid", "spin4", "meaning");

    for (int i = 0; i < n; i++) {
        chain->Clear();
        int stat = chain->Make(i);
        if (stat) { printf("chain->Make(%d) status=%d, stopping\n", i, stat); break; }

        StMuEvent* ev = muDstMaker->muDst()->event();
        if (!ev) { printf("%8d: no StMuEvent\n", i); continue; }

        int runNumber = ev->runNumber();
        if (runNumber != lastRunNumber) {
            spinDb->InitRun(runNumber);
            lastRunNumber = runNumber;
        }

        unsigned int bx7 = ev->l0Trigger().bunchCrossingId7bit(runNumber);

        int spin4 = -1;
        if (spinDb->isValid()) spin4 = spinDb->spin4usingBX7((int)bx7);

        // Legend (StSpinDbMaker::optimizeTables()): yellow up = {1,5,9},
        // yellow down = {2,6,10}, blue up = {4,5,6}, blue down = {8,9,10};
        // a beam is "empty/unpolarized" in this code for any other value.
        const char* yellow = "empty/unpol";
        const char* blue   = "empty/unpol";
        if      (spin4==1||spin4==5||spin4==9)  yellow = "up";
        else if (spin4==2||spin4==6||spin4==10) yellow = "down";
        if      (spin4==4||spin4==5||spin4==6)  blue = "up";
        else if (spin4==8||spin4==9||spin4==10) blue = "down";

        char meaning[64];
        if (spin4 < 0) snprintf(meaning, sizeof(meaning), "no data");
        else snprintf(meaning, sizeof(meaning), "yellow=%s blue=%s", yellow, blue);

        printf("%8d %10d %6u %8s %6d %s\n", i, runNumber, bx7,
               spinDb->isValid() ? "true" : "false", spin4, meaning);
    }

    chain->Finish();
}
