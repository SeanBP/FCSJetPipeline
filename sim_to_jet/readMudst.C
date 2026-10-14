// Reads the StarGenStats object written by the starsim stage
// (pythia8.genstats.root or pythia6.genstats.root, renamed by sim_to_jet.xml).
// Returns false, leaving the arguments untouched, if none is found.
Bool_t ReadMCJobStats(Double_t &sigmaGen, Int_t &nTried){
  TString candidates[2] = {"pythia8.genstats.root", "pythia6.genstats.root"};
  for (int i = 0; i < 2; i++){
    if ( gSystem->AccessPathName(candidates[i]) ) continue; // AccessPathName returns 0 on success
    TFile f(candidates[i]);
    if ( f.IsZombie() ) continue;
    TIter next(f.GetListOfKeys());
    TKey* key;
    while ( (key = (TKey*)next()) ){
      if ( TString(key->GetClassName()) != "StarGenStats" ) continue;
      StarGenStats* stats = (StarGenStats*)key->ReadObj();
      if ( !stats ) continue;
      // mb, as PYTHIA reports it
      sigmaGen = stats->sigmaGen;
      nTried = stats->nTried;

      // Pythia 6: StarPythia6::Stats() reads the cross section 4 bytes off
      // (its PyInt5_t struct omits NGENPD), giving a denormal whose low word
      // is the true value's high word. Shift it back. Same decoding as
      // jet_scale/JetXSecWeights.h's DecodePythia6Sigma().
      if ( sigmaGen > 0.0 && sigmaGen < 1e-300 ) {
        Double_t raw = sigmaGen;
        Double_t decoded = 0.0;
        UInt_t* rawWords = (UInt_t*)&raw;
        UInt_t* decodedWords = (UInt_t*)&decoded;
        decodedWords[0] = 0;
        decodedWords[1] = rawWords[0];
        cout << "ReadMCJobStats: decoded Pythia 6 sigmaGen " << raw << " -> " << decoded << " mb" << endl;
        sigmaGen = decoded;
      }

      cout << "ReadMCJobStats: " << candidates[i] << " -> sigmaGen=" << sigmaGen
           << " mb, nTried=" << stats->nTried << endl;
      return kTRUE;
    }
    cout << "ReadMCJobStats: " << candidates[i] << " found but no StarGenStats object in it" << endl;
  }
  cout << "ReadMCJobStats: no genstats file found -- sigma_job/N_job will stay at the -1 sentinel" << endl;
  return kFALSE;
}

void readMudst(Int_t nEvents, Int_t nFiles, TString InputFileList, TString OutputDir =".", TString JobIdName = "", Int_t trgSelect = 202209, Double_t ptHatMin = -1, Double_t ptHatMax = -1, TString trigFilter = "" )
{

  // Load libraries
  gROOT   -> Macro("loadMuDst.C");
  gROOT   -> Macro("Load.C");
  gSystem -> Load("StFcsDbMaker");
  gSystem -> Load("RTS");
  gSystem -> Load("StFcsTriggerSimMaker");
  gSystem -> Load("StSimpleReaderMaker.so");
  gSystem -> Load("StEpdUtil");
  gSystem -> Load("StarGeneratorEvent.so"); // StarGenStats, for ReadMCJobStats() above

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

  // Simulated events have no real trigger decision; StSimpleReaderMaker
  // computes Trig_flag[] itself via fcs_trg_base (see SetTriggerSim()).
  StFcsTriggerSimMaker* trgSim = new StFcsTriggerSimMaker();
  trgSim->setTrigger(trgSelect);
  trgSim->setSimMode(1);
  trgSim->setEtGain(1.0);
  trgSim->setThresholdFile((char*)"stage_params.txt");
  trgSim->setDebug(0);

  // Analysis Maker
  StSimpleReaderMaker* AnalysisCode  =  new StSimpleReaderMaker(muDstMaker) ;
  AnalysisCode -> SetTriggerSim(trgSim) ;
  // Comma-separated FCS trigger flag names: only events that fired one are written (empty = all)
  AnalysisCode -> SetTriggerFilter(trigFilter.Data()) ;

  Double_t sigmaGen; Int_t nTried;
  if ( ReadMCJobStats(sigmaGen, nTried) )
    AnalysisCode -> SetMCJobStats(sigmaGen, nTried) ;

  AnalysisCode -> SetMCPtHatMin(ptHatMin) ;
  AnalysisCode -> SetMCPtHatMax(ptHatMax) ;

  AnalysisCode -> SetOutputFileName("SimpleTree_mudst.root") ; // Name the output file
  if ( nEvents == 0 )  nEvents = 10000000 ;       // Take all events in nFiles if nEvents = 0

  // Loop over the links in the chain
  chain -> Init() ;

  // trgSim's own Make() segfaults without a real StEvent; StSimpleReaderMaker
  // drives fcs_trg_base directly instead (Init()/InitRun() already ran above).
  trgSim -> SetActive(kFALSE) ;

  chain -> EventLoop(1,nEvents) ;
  chain -> Finish() ;

  // Cleanup
  delete chain ;
}
