// Reads back the StarGenStats object StarPrimaryMaker::Finish() writes at
// the end of the starsim stage (pythia8.genstats.root / pythia6.genstats.root
// -- sim_to_jet.xml renames it there right after the starsim.C step, since
// runSimBfc.C's own output file would otherwise silently overwrite the
// original pythia8.starsim.root/pythia6.starsim.root name and destroy it
// before this macro ever runs). Scans by class rather than by hardcoded key
// name (StarPythia8's object is named "Pythia8", StarPythia6's is
// "Pythia6Stats") so this doesn't break if that naming ever changes.
// Returns false (leaves sigma_pb/sigmaErr_pb/nGen untouched) if no genstats
// file is found -- e.g. a non-sim_to_jet caller, or an older job predating
// this bookkeeping.
Bool_t ReadMCXSec(Double_t &sigma_pb, Double_t &sigmaErr_pb, Int_t &nGen){
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
      // Pythia convention: sigmaGen/sigmaErr in mb; 1 mb = 1e9 pb.
      sigma_pb = stats->sigmaGen * 1e9;
      sigmaErr_pb = stats->sigmaErr * 1e9;
      nGen = stats->nAccepted;
      cout << "ReadMCXSec: " << candidates[i] << " -> sigmaGen=" << stats->sigmaGen
           << " mb (+-" << stats->sigmaErr << "), nAccepted=" << stats->nAccepted
           << " -> sigma_pb=" << sigma_pb << ", nGen=" << nGen << endl;
      return kTRUE;
    }
    cout << "ReadMCXSec: " << candidates[i] << " found but no StarGenStats object in it" << endl;
  }
  cout << "ReadMCXSec: no genstats file found -- mc_sigma_pb/mc_n_gen will stay at the -1 sentinel" << endl;
  return kFALSE;
}

void readMudst(Int_t nEvents, Int_t nFiles, TString InputFileList, TString OutputDir =".", TString JobIdName = "", Int_t trgSelect = 202209, Double_t ptHatMin = -1 )
{

  // Load libraries
  gROOT   -> Macro("loadMuDst.C");
  gROOT   -> Macro("Load.C");
  gSystem -> Load("StFcsDbMaker");
  gSystem -> Load("RTS");
  gSystem -> Load("StFcsTriggerSimMaker");
  gSystem -> Load("StSimpleReaderMaker.so");
  gSystem -> Load("StEpdUtil");
  gSystem -> Load("StarGeneratorEvent.so"); // StarGenStats, for ReadMCXSec() above

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

  Double_t mc_sigma_pb, mc_sigma_err_pb; Int_t mc_n_gen;
  if ( ReadMCXSec(mc_sigma_pb, mc_sigma_err_pb, mc_n_gen) )
    AnalysisCode -> SetMCXSec(mc_sigma_pb, mc_sigma_err_pb, mc_n_gen) ;

  AnalysisCode -> SetMCPtHatMin(ptHatMin) ;

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

  // trgSim's own Make() segfaults without a real StEvent; StSimpleReaderMaker
  // drives fcs_trg_base directly instead (Init()/InitRun() already ran above).
  trgSim -> SetActive(kFALSE) ;

  chain -> EventLoop(1,nEvents) ;
  chain -> Finish() ;

  // Cleanup
  delete chain ;
}
