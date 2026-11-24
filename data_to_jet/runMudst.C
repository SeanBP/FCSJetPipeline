void runMudst(char* file="st_fwd_23080044_raw_1000019.MuDst.root", 
	      int ifile=-1, Int_t nevt=10, char* outdir=".", int readMuDst=1, int debug=0,
	      const char* trigFilter=""){  // trigFilter: comma-separated FCS trigger flag names; only events that fired one are written (empty = all)
    gROOT->Macro("Load.C");
    gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    gSystem->Load("StEventMaker");
    gSystem->Load("StFcsDbMaker");
    gSystem->Load("StFcsRawHitMaker");
    gSystem->Load("StFcsClusterMaker");
    gSystem->Load("libMinuit");
    gSystem->Load("StFcsPointMaker");
    gSystem -> Load("StSimpleReaderMaker");
    gSystem -> Load("StEpdUtil");
    gSystem -> Load("StEpdDbMaker");
    gSystem -> Load("StEpdHitMaker");
    gSystem -> Load("StSpinDbMaker");
    gMessMgr->SetLimit("I", 0);
    gMessMgr->SetLimit("Q", 0);
    gMessMgr->SetLimit("W", 0);

    StChain* chain = new StChain("StChain"); chain->SetDEBUG(0);
    StMuDstMaker* muDstMaker = new StMuDstMaker(0, 0, "", file,".", 1000, "MuDst");
    // Branches StSimpleReaderMaker never reads (2026-11-16, job efficiency):
    // sTGC raw hits alone are ~42% of an st_fwd MuDst; with FST and TOF hits
    // off, reading takes ~30% less time. (TOF header stays: VPD vertex.)
    muDstMaker->SetStatus("Ftt*", 0);
    muDstMaker->SetStatus("Fst*", 0);
    muDstMaker->SetStatus("BTofHit", 0);
    muDstMaker->SetStatus("BTofRawHit", 0);
    int n=muDstMaker->tree()->GetEntries();
    printf("Found %d entries in Mudst\n",n);
    int start=0, stop=n;
    if(ifile>=0){
	int start=ifile*nevt;
	int stop=(ifile+1)*nevt-1;
	if(n<start) {printf(" No event left. Exiting\n"); return;}
	if(n<stop)  {printf(" Overwriting end event# stop=%d\n",n); stop=n;}
    }else if(nevt>=0 && nevt<n){
	stop=nevt;
    }else if(nevt==-2){
	stop=2000000000; 
    }
    printf("Doing Event=%d to %d\n",start,stop);
    
    St_db_Maker* dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb"); 
    if(dbMk){
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
    }
    
    StFcsDbMaker *fcsDbMkr= new StFcsDbMaker();
    StFcsDb* fcsDb = (StFcsDb*) chain->GetDataSet("fcsDb");

    // Recluster the FCS hits with the pipeline's tower thresholds, the same
    // as readMudst.C does for simulation, so data and simulation share them
    // whatever library made the MuDst: Ecal towers below 0.4 GeV (6 sigma
    // above the worst-case Run 22 pedestal noise, with margin) are not
    // clustered; seeds at 1 GeV (P24ia).
    // Early trigger filter (2026-11-16): events the reader's trigger filter
    // would drop skip the makers below (StSimpleTrigFilterMaker.h)
    StSimpleTrigFilterMaker* trgFilter = new StSimpleTrigFilterMaker();
    StFcsMuToEventMaker* fcsMuToEvent = new StFcsMuToEventMaker();

    // EPD hits for the vertex (StSimpleReaderMaker): where the MuDst has none
    // (P24ia st_fwd before ~day 50, 2026-11-16), rebuild them from the
    // trigger data's QT data with the EPD DB into the StEvent made above;
    // where it has hits, StEpdHitMaker leaves them alone
    StEpdDbMaker* epdDbMkr = new StEpdDbMaker("epdDb");
    StEpdHitMaker* epdHitMkr = new StEpdHitMaker();
    epdHitMkr->setReadMuDst();
    StFcsClusterMaker* fcsClu = new StFcsClusterMaker();
    fcsClu->setTowerEThreSeed(1.0, 1.0);
    fcsClu->setTowerEThreshold(0.4, 0.01);
    StFcsPointMaker* fcsPoi = new StFcsPointMaker();

    // Per-event spin config; -1 outside the spin-DB coverage window
    StSpinDbMaker* spinDbMkr = new StSpinDbMaker();

    //Analysis Maker
    StSimpleReaderMaker* AnalysisCode  =  new StSimpleReaderMaker(muDstMaker) ;
    AnalysisCode -> SetSpinDb(spinDbMkr) ;
    AnalysisCode -> SetTriggerFilter(trigFilter) ;
    trgFilter -> SetReader(AnalysisCode) ;
    AnalysisCode -> SetFcsFromStEvent() ; // the reclustered clusters and points
    // Event vertex: EPD, else BBC TAC, as measured, else the nominal z (see
    // the StSimpleReaderMaker header)

    // HCAL gain correction; a placeholder of 1.0 until per-tower corrections exist
    AnalysisCode -> SetHcalRetroactiveGainCorr("hcal_gain_corrections") ;

    // Hot / dead channel mask (asymmetry/channel_qa, Sean 2026-11-24): with a
    // channel_mask.txt in the working directory (submit_an_job.sh -M), masked
    // hits are zeroed before reclustering and left out of the SimpleTree's hit
    // list and trigger sums, for the runs of their period
    if ( !gSystem->AccessPathName("channel_mask.txt") ) {
      fcsMuToEvent -> SetChannelMask("channel_mask.txt") ;
      AnalysisCode -> SetChannelMask("channel_mask.txt") ;
    }

    TString foriternum(file);
    Ssiz_t last_ = foriternum.Last('_');
    TString filename = "SimpleTree_mudst_";
    TString iternum = foriternum(last_+1,7);
    // The file name ends <run>_<iternum>.MuDst.root, or, for the st_fwd
    // stream, <run>_raw_<iternum>.MuDst.root; either way the output is named
    // with the full 8-digit run number.
    Ssiz_t run_end = last_;
    if ( run_end >= 4 && foriternum(run_end-4,4) == "_raw" ) run_end -= 4;
    TString runnum = foriternum(run_end-8,8);
    filename += runnum + "_" + iternum + ".root";
    AnalysisCode -> SetOutputFileName(filename.Data()) ; 
    if ( nevt == 0 )  nevt = 10000000 ; 

    chain->Init();
    chain->EventLoop(start,stop);
    chain->Finish();
    delete chain;
}
