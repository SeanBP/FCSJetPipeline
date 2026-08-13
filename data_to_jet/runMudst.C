void runMudst(char* file="st_fwd_23080044_raw_1000019.MuDst.root", 
	      int ifile=-1, Int_t nevt=10, char* outdir=".", int readMuDst=1, int debug=0){  
    gROOT->Macro("Load.C");
    gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    gSystem->Load("StEventMaker");
    gSystem->Load("StFcsDbMaker");
    gSystem->Load("StFcsRawHitMaker");
    //gSystem->Load("StFcsWaveformFitMaker");
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

    // Real per-event spin configuration (see StSimpleReaderMaker.h's
    // Spin_config comment): resolved from Calibrations/rhic via bunch
    // crossing id, -1 ("no data") wherever that DB isn't covered -- see
    // SpinPlayground/run_number_to_spin_coverage.sh for the actual
    // coverage window. Not used by sim_to_jet's readMudst.C -- simulated
    // events have no real polarization pattern to report.
    StSpinDbMaker* spinDbMkr = new StSpinDbMaker();

    //Analysis Maker
    StSimpleReaderMaker* AnalysisCode  =  new StSimpleReaderMaker(muDstMaker) ;
    AnalysisCode -> SetSpinDb(spinDbMkr) ;

    // Retroactive HCAL gain correction (see StSimpleReaderMaker.h's
    // mUseHcalRetroCorr comment and hcal_gain_corrections/README.md):
    // currently a no-op (the placeholder file is all 1.0), until real
    // per-tower HCAL corrections are measured. Not used by sim_to_jet's
    // readMudst.C -- simulated events already get real, live gain
    // straight from StFcsDb.
    AnalysisCode -> SetHcalRetroactiveGainCorr("hcal_gain_corrections") ;

    // In order to speed up the analysis and eliminate IO, turn off unneeded branches

    // Miscellaneous things we need before starting the chain
    //TString Name = JobIdName ; //Name.Append(".histograms.root") ;
    TString foriternum(file);
    Ssiz_t last_ = foriternum.Last('_');
    TString filename = "SimpleTree_mudst_";
    TString iternum = foriternum(last_+1,7);
    TString runnum = foriternum(last_-8,8);  //Since I know file name will have the form xrd_runnum_iternum.MuDst.root I can use this this hack to get the runnumber
    filename += runnum + "_" + iternum + ".root";
    //std::cout << iternum.Atoi() << std::endl;
    AnalysisCode -> SetOutputFileName(filename.Data()) ; // Name the output file for histograms
    if ( nevt == 0 )  nevt = 10000000 ; 

    chain->Init();
    chain->EventLoop(start,stop);
    chain->Finish();
    delete chain;
}
