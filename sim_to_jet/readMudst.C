void readMudst(Int_t nEvents, Int_t nFiles, TString InputFileList, TString OutputDir =".", TString JobIdName = "", Int_t trgSelect = 202209 )
{

  // Load libraries
  gROOT   -> Macro("loadMuDst.C");
  gROOT   -> Macro("Load.C");
  gSystem -> Load("StFcsDbMaker");
  gSystem -> Load("RTS");
  gSystem -> Load("StFcsTriggerSimMaker");
  gSystem -> Load("StSimpleReaderMaker.so");
  gSystem -> Load("StEpdUtil");

  // List of member links in the chain
  StChain*                    chain  =  new StChain ;
  StMuDstMaker*          muDstMaker  =  new StMuDstMaker(0,0,"",InputFileList,"MuDst",nFiles) ;
 
  // Load the DB Maker
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
    
  StFcsDbMaker *fcsDbMkr = new StFcsDbMaker();

  // FCS trigger simulation: this muDst came from BFC's "cmudst" output
  // (simulated events, no real online trigger decision to read back), so
  // StSimpleReaderMaker runs fcs_trg_base itself to compute the same
  // Trig_flag[] branches data_to_jet's runMudst.C fills from real fired
  // trigger ids -- see StSimpleReaderMaker::SetTriggerSim() for how, and
  // why trgSim's own Make() must not run (SetActive(kFALSE) below).
  // setEtGain(1.0)/setThresholdFile(...) match the (previously
  // commented-out) template already in this pipeline's own runSimBfc.C.
  StFcsTriggerSimMaker* trgSim = new StFcsTriggerSimMaker();
  trgSim->setTrigger(trgSelect);
  trgSim->setSimMode(1);
  trgSim->setEtGain(1.0);
  trgSim->setThresholdFile((char*)"stage_params.txt");
  trgSim->setDebug(0);

  // Analysis Maker
  StSimpleReaderMaker* AnalysisCode  =  new StSimpleReaderMaker(muDstMaker) ;
  AnalysisCode -> SetTriggerSim(trgSim) ;

  // In order to speed up the analysis and eliminate IO, turn off unneeded branches
  //FIX ME: StMuMCTrack won't set status to 1 after being turned off
  //muDstMaker -> SetStatus("*",0) ;                // Turn off all branches
  //muDstMaker -> SetStatus("MuEvent",1) ;          // Turn on the Event data (esp. Event number)
  //muDstMaker -> SetStatus("StMuMcTrack",1) ;      // Turn on the MC Track data
  //muDstMaker -> SetStatus("FcsHit",1) ;           // Turn on the FCS Hit data

  // Miscellaneous things we need before starting the chain
  //TString Name = JobIdName ; 
  //Name.Append(".histograms.root") ;
  AnalysisCode -> SetOutputFileName("SimpleTree_mudst.root") ; // Name the output file
  if ( nEvents == 0 )  nEvents = 10000000 ;       // Take all events in nFiles if nEvents = 0

  // Loop over the links in the chain
  chain -> Init() ;

  // trgSim's Init()/InitRun() (gains, pedestals, stage_version) already
  // ran as part of chain->Init() above; its own Make() would segfault
  // reading a muDst with no StEvent present (see
  // StSimpleReaderMaker::SetTriggerSim() in the header), so keep the
  // chain from calling it -- StSimpleReaderMaker drives the same
  // fcs_trg_base engine directly instead.
  trgSim -> SetActive(kFALSE) ;

  chain -> EventLoop(1,nEvents) ;
  chain -> Finish() ;

  // Cleanup
  delete chain ;
}
