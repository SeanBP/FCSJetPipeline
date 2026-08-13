// TrgSimTest.C
//
// Playground validation driver: sets up the chain and hands the whole
// per-event loop to chain->EventLoop(), which calls straight into
// TrgSimHelperMaker::Make() (compiled, cons-built -- see its header for
// why the per-event loop needed to live there and not in this
// interpreted macro).
//
// IMPORTANT interpretation note printed with the summary: real L1
// trigger bits are heavily PRESCALED (only a fraction of events where
// the physics condition was true actually get flagged "fired"), while
// the simulator reproduces the raw physics condition with no
// prescaling. So sim=1,real=0 is EXPECTED/benign for high-rate
// triggers. The real validation question is P(sim=1|real=1): whenever
// real=1 (actually recorded as fired), does sim also =1 (physics
// condition met)? That should be ~100% if trgSelect is correct; a low
// value there is a real red flag.
//
// Usage:
//   root4star -b -q TrgSimTest.C\(\"<path-or-xrootd-URL-to-a-muDst>\",\"TriggerIDs.txt\",202209,-1\)
//
// No local muDst is shipped in this playground -- StMuDstMaker/TFile::Open
// read root://... xrootd URLs directly, so the default below streams a
// known-good real Run 22 file (run 22359013) straight from STAR's xrootd
// pool with no download needed. See README.md for how that URL was found.

void TrgSimTest(const char* file =
                    "root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/"
                    "production_pp500_2022/ReversedFullField/P24ia/2021/359/22359013/"
                    "st_fwd_22359013_raw_1500006.MuDst.root",
                const char* triggerIdFile = "TriggerIDs.txt",
                int trgSelect = 202209,
                int maxEvents = -1){

    gROOT->Macro("Load.C");
    gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    gSystem->Load("StFcsDbMaker");
    gSystem->Load("RTS");
    gSystem->Load("StFcsTriggerSimMaker");
    gSystem->Load("StSpinPoolTrgSimHelperMaker");
    gMessMgr->SetLimit("I", 0);
    gMessMgr->SetLimit("Q", 0);
    gMessMgr->SetLimit("W", 0);

    StChain* chain = new StChain("StChain");
    StMuDstMaker* muDstMaker = new StMuDstMaker(0, 0, "", file, ".", 1000, "MuDst");

    St_db_Maker* dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb");
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

    StFcsDbMaker* fcsDbMkr = new StFcsDbMaker();

    StFcsTriggerSimMaker* trgSim = new StFcsTriggerSimMaker();
    trgSim->setTrigger(trgSelect);
    trgSim->setDebug(0);
    // setThresholdDb() is a no-op in the official maker (readThresholdDb()
    // is an empty stub, "to be implemented before run22 online DB moves
    // to offline" -- never was). setThresholdFile() is the real,
    // implemented mechanism; stage_params.txt is the calibrated file
    // shipped alongside StFcsTriggerSimMaker itself (copied from
    // StRoot/StSpinPool/StFcsTriggerSimMaker/files/), used in place of
    // fcs_trg_base's hardcoded default thresholds.
    trgSim->setThresholdFile((char*)"stage_params.txt");

    TrgSimHelperMaker* trgHelper = new TrgSimHelperMaker();
    trgHelper->setMuDstMaker(muDstMaker);
    trgHelper->setTrgSim(trgSim);
    trgHelper->setTriggerIdFile(triggerIdFile);
    trgHelper->setMaxPrintEvents(30);

    chain->Init();

    // trgSim's own Make() segfaults in its MuDst-fallback path when run
    // outside a full StEvent chain (see TrgSimHelperMaker.h); its
    // Init()/InitRun() (gains, pedestals, stage_version) already ran via
    // chain->Init() above. TrgSimHelperMaker::Make() drives the same
    // fcs_trg_base engine directly instead.
    trgSim->SetActive(kFALSE);

    int n = muDstMaker->tree()->GetEntries();
    if ( maxEvents >= 0 && maxEvents < n ) n = maxEvents;
    printf("Found %d entries in MuDst, reading %d, trgSelect=%d\n",
           (int)muDstMaker->tree()->GetEntries(), n, trgSelect);

    chain->EventLoop(0, n);
    chain->Finish();
}
