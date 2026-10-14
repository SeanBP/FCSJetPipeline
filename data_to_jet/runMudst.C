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
    gSystem -> Load("StSpinDbMaker");
    gMessMgr->SetLimit("I", 0);
    gMessMgr->SetLimit("Q", 0);
    gMessMgr->SetLimit("W", 0);

    StChain* chain = new StChain("StChain"); chain->SetDEBUG(0);
    StMuDstMaker* muDstMaker = new StMuDstMaker(0, 0, "", file,".", 1000, "MuDst");
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

    // Per-event spin config; -1 outside the spin-DB coverage window
    StSpinDbMaker* spinDbMkr = new StSpinDbMaker();

    //Analysis Maker
    StSimpleReaderMaker* AnalysisCode  =  new StSimpleReaderMaker(muDstMaker) ;
    AnalysisCode -> SetSpinDb(spinDbMkr) ;
    AnalysisCode -> SetTriggerFilter(trigFilter) ;

    // HCAL gain correction; a placeholder of 1.0 until per-tower corrections exist
    AnalysisCode -> SetHcalRetroactiveGainCorr("hcal_gain_corrections") ;

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
