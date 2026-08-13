TString input_dir   = "./";
TString output_dir  = "./";
//TString input_chain = "sdt20161210.120000,fzin,geant,evout,y2015,FieldOn,logger,MakeEvent,McEvout,IdTruth,ReverseField,db,fcsSim,fcsCluster,fcsPoint,-tpcDB";
//TString input_chain = "sdt20211025.120000,fzin,geant,FieldOn,logger,MakeEvent,fcsSim,fcsWFF,fcsCluster,fcsPoint";
TString input_chain = "y2023,AgML,USExgeom,fzin,geant,FieldOn,logger,MakeEvent,fcsSim,fcsWFF,fcsCluster,fcsPoint,cmudst";
class StFmsSimulatorMaker;

// Fallback only -- see RunNumberToDateTime() below, which tries the real
// RunLog DB first and only falls back to this arithmetic if unreachable.
//
// Converts a Run 22 pp500 run number (YYDDDnnnnn: 2-digit year code, 3-digit
// day-of-year, then a segment number) to a YYYYMMDD date.
//
// IMPORTANT: the run number's 2-digit year code is NOT the real calendar
// year -- it's one less (STAR's "Run 22" label starts partway through
// real calendar year 2021). Verified empirically (Aug 2026) against a
// real muDst: run 22359013's actual StMuEvent timestamp is 2021-12-25,
// not 2022-12-25 -- i.e. year code "22" -> real year 2021, "23" -> 2022.
// Cross-checked against several more (run, date) pairs from Akio Ogawa's
// run22 gain page and, more authoritatively, directly against the real
// RunLog.runDescriptor DB (see RunNumberToDateTime()) for runs 22359013,
// 23101043, and 23276001: all consistent with this same -1 offset.
int RunNumberToDateArithmetic(int runNumber){
    int yy  = runNumber / 1000000;
    int ddd = (runNumber / 1000) % 1000;
    int year = 2000 + yy - 1;
    int daysInMonth[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    bool leap = (year%4==0 && year%100!=0) || (year%400==0);
    if(leap) daysInMonth[1]=29;
    int month=1, day=ddd;
    for(int m=0; m<12; m++){
        if(day <= daysInMonth[m]){ month=m+1; break; }
        day -= daysInMonth[m];
    }
    return year*10000 + month*100 + day;
}

// Looks up a run's real start time from the authoritative RunLog.runDescriptor
// DB (same mechanism StFcsDbMaker/macro/fcsTimedepGainCorr_db.C uses), for
// St_db_Maker::SetDateTime(), so a specific real run's FCS calibration era
// can be pinned for the simulation instead of BFC's own nominal "y2023"
// chain-flag date. Falls back to RunNumberToDateArithmetic() (with itim=0)
// if the DB is unreachable -- e.g. a SUMS batch worker node with more
// restricted outbound network access than the interactive session this was
// developed/tested in (verified working there, Aug 2026; batch-node
// behavior untested) -- so a network hiccup degrades to the
// previously-verified-correct heuristic rather than failing the job.
//
// Port selection: the formula in fcsTimedepGainCorr_db.C (3400+year-1,
// year=runNumber/1000000-1) undershoots by 1 at the earliest Run 22 runs
// (gives 3420 for run 22359013, which returns nothing; 3421 is the port
// that actually has that run's data, and 3421 also serves every later
// run tested). So port 3421 is tried first as the most likely hit, then
// the formula's port, then a couple of neighbors, before giving up.
void RunNumberToDateTime(int runNumber, int &idat, int &itim){
    int yearCode = runNumber/1000000;
    int formulaPort = 3400 + (yearCode-1) - 1;
    int candidatePorts[4] = {3421, formulaPort, formulaPort-1, formulaPort+1};

    for(int i=0; i<4; i++){
        int port = candidatePorts[i];
        bool dup = false;
        for(int j=0;j<i;j++) if(candidatePorts[j]==port) dup=true;
        if(dup) continue;

        TString cmd = Form("mysql -h db04.star.bnl.gov --port=%d --connect-timeout=10 -N -s -e "
                            "\"SELECT startRunTime FROM RunLog.runDescriptor WHERE runNumber=%d LIMIT 1\" 2>/dev/null",
                            port, runNumber);
        TString result = gSystem->GetFromPipe(cmd.Data());
        result = result.Strip(TString::kBoth);

        if(result.IsDigit() && result.Length() > 0){
            int starttime = result.Atoi();
            TString dateStr = gSystem->GetFromPipe(Form("date -u -d @%d +%%Y%%m%%d", starttime));
            TString timeStr = gSystem->GetFromPipe(Form("date -u -d @%d +%%H%%M%%S", starttime));
            idat = dateStr.Atoi();
            itim = timeStr.Atoi();
            cout << "RunNumberToDateTime: run=" << runNumber << " -> RunLog DB (port=" << port
                 << ") startRunTime=" << starttime << " -> date=" << idat << " time=" << itim << endl;
            return;
        }
    }

    idat = RunNumberToDateArithmetic(runNumber);
    itim = 0;
    cout << "*** RunNumberToDateTime: RunLog DB unreachable/empty for run=" << runNumber
         << " on all candidate ports -- falling back to day-of-year arithmetic: date=" << idat
         << " (verified consistent with RunLog for 3 other test runs, but this run itself was"
         << " not cross-checked against the DB just now) ***" << endl;
}

void runSimBfc( Int_t nEvents=1000, Int_t run=1, const char* pid="jet", int TrgVersion=202207,
		int debug=0, int e=0, float pt=1.5, float vz=0.0,
		char* epdmask="0.0100",
		int leakyHcal=0,
		int eventDisplay=0,
		TString myDir=input_dir, TString myOutDir=output_dir,
		TString myChain=input_chain, Int_t mnEvents=0,
		Int_t calibRun=0){
    
  gROOT->LoadMacro("bfc.C");
  //  gROOT->Macro("loadMuDst.C");

  TString myDat;
  TString proc(pid);

if(proc.Contains("dy") || proc.Contains("mb") || proc.Contains("jet") || proc.Contains("dybg")){
      myDat=Form("pythia.%s.vz%d.run%i.fzd",pid,(int)vz,run);
  }else if(proc.Contains("pythia8")){
      myDat = "pythia8.starsim.fzd";
  }else if(proc.Contains("pythia6")){
      myDat = "pythia6.starsim.fzd";
  }else if(proc.Contains("herwig")){
      myDat = "herwig6.starsim.fzd";
  }else if(e>0.0){
      myDat=Form("%s.e%d.vz%d.run%i.fzd",pid,e,(int)vz,run);
  }else{
      myDat=Form("%s.pt%3.1f.vz%d.run%i.fzd",pid,pt,(int)vz,run);
  }
  printf("Opening %s\n",(myDir+myDat).Data());

  bfc( -1, myChain, myDir+myDat );
 
  TString outfile = myOutDir + myDat.ReplaceAll(".fzd",".root");      
  cout << "output file=" <<outfile<<endl;
  chain->SetOutputFile(outfile);
  
  St_db_Maker *dbMk= (St_db_Maker*) chain->GetMaker("db");
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
  dbMk->SetAttr("blacklist", "Calibrations_rich");

  // Pin the conditions-DB lookup (FCS gain/gainCorrection, mainly -- see
  // the blacklist above) to a specific real run's calibration era,
  // instead of BFC's own "y2023" nominal chain-flag date. calibRun=0
  // (default) leaves St_db_Maker on its normal BFC-chosen date.
  if(calibRun != 0){
    int idat, itim;
    RunNumberToDateTime(calibRun, idat, itim);
    dbMk->SetDateTime(idat, itim);
    cout << "Pinned conditions DB to calibRun=" << calibRun
         << " -> date=" << idat << " time=" << itim << endl;
  }

  StFcsDbMaker* fcsdbmkr = (StFcsDbMaker*) chain->GetMaker("fcsDbMkr");
  cout << "fcsdbmkr="<<fcsdbmkr<<endl;
  //fcsdbmkr->setDbAccess(0);

  StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");
  cout << "fcsdb="<<fcsdb<<endl;
  // Left on default GAINMODE::DB (no forceUniformGain/forceUniformGainCorrection
  // override): StFcsDb pulls real calibrated gain/gain-correction from the
  // STAR conditions DB by default. Verified directly (TriggerPlayground,
  // Aug 2026) against both a real run and this same simulated BFC chain --
  // ECAL gainCorr comes back genuinely per-tower (not flat), consistent
  // with Xilin's pi0-mass-fit calibration (avg ~1.2-1.4, matching Akio's
  // run22 gain page); HCAL gainCorr comes back flat 1.0, consistent with
  // HCAL not having per-tower calibration yet -- and HCAL's base gain
  // still correctly reflects the known 1.3 electronics-gain factor (and,
  // for later run periods, the average 1.21 ECAL-derived correction
  // baked into the HCAL base gain itself). Forcing a flat 0.0053/1.0
  // override here, as this code previously did, discarded all of that
  // and used a uniform placeholder for both detectors instead.

  StFcsFastSimulatorMaker *fcssim = (StFcsFastSimulatorMaker*) chain->GetMaker("fcsSim");
  fcssim->setDebug(0);
  fcssim->setLeakyHcal(leakyHcal);
 
  StFcsWaveformFitMaker *wff=(StFcsWaveformFitMaker *)chain->GetMaker("StFcsWaveformFitMaker");
  wff->setDebug(0);
  wff->setEnergySelect(0,0,0);

  gSystem->Load("StEpdUtil");

/*
  StFcsClusterMaker *clu=(StFcsClusterMaker *)chain->GetMaker("StFcsClusterMaker");
  clu->setDebug(1);

  StFcsPointMaker *poi=(StFcsPointMaker *)chain->GetMaker("StFcsPointMaker");
  poi->setDebug(1);
  poi->setShowerShape(3);
  

  gSystem->Load("RTS");
  gSystem->Load("StFcsTriggerSimMaker");
  StFcsTriggerSimMaker* fcsTrgSim = new StFcsTriggerSimMaker(); 
  fcsTrgSim->setSimMode(1);
  fcsTrgSim->setTrigger(TrgVersion);
  fcsTrgSim->setDebug(debug);
  fcsTrgSim->setEtGain(1.0); //ET match
  //fcsTrgSim->setEtGain(0.5); //halfway
  //fcsTrgSim->setEtGain(0.0); //E match
  //fcsTrgSim->setReadPresMask(Form("mask/fcs_ecal_epd_mask.ele.pt0.6.vz0.thr%s.txt",epdmask));
  //TString txfile(outfile); txfile.ReplaceAll(".root",".event.txt");  fcsTrgSim->setWriteEventText(txfile.Data());
  TString qafile(outfile); qafile.ReplaceAll(".root",".qahist.root"); fcsTrgSim->setWriteQaHist(qafile.Data());
  fcsTrgSim->setThresholdFile("stage_params.txt");

  gSystem->Load("StFcsTrgQaMaker");
  StFcsTrgQaMaker* fcsTrgQa = new StFcsTrgQaMaker(); 
  TString tqafile(outfile); tqafile.ReplaceAll(".root",Form(".thr%s.trgqa.root",epdmask)); 
  fcsTrgQa->setFilename(tqafile.Data());
  fcsTrgQa->setEcalPtThr(pt*0.75);
*/
  if(eventDisplay>0){
      gSystem->Load("StEpdUtil");
      gSystem->Load("StFcsEventDisplay");
      StFcsEventDisplay* fcsed = new StFcsEventDisplay();
      fcsed->setMaxEvents(eventDisplay);
      outfile.ReplaceAll(".root",".eventDisplay.png");
      fcsed->setFileName(outfile.Data());
  }

  chain->Init();
  StMaker::lsMakers(chain);
  chain->EventLoop(mnEvents,nEvents);

  // Sanity-check log line: confirms which gain/gainCorrection actually
  // got used (a representative ECAL/HCAL channel each), so a job's log
  // shows at a glance whether calibRun took effect as expected.
  cout << "Gain check (calibRun=" << calibRun << ", det,id,gain,gainCorr): "
       << "ecal(0,50)=" << fcsdb->getGain(0,50) << "/" << fcsdb->getGainCorrection(0,50)
       << " hcal(2,50)=" << fcsdb->getGain(2,50) << "/" << fcsdb->getGainCorrection(2,50) << endl;

  chain->Finish();
}
