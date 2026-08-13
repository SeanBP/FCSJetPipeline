// Minimal BFC chain (DB + FcsDbMaker only, no GEANT/simulation/input file
// needed) to query StFcsDb's gain/gainCorrection for a given calibration
// era -- same minimal chain as GeometryPlayground/get_fcs_geometry.C
// (same StFcsDbMaker::InitRun() call pulls both Geometry/fcs and
// Calibrations/fcs tables together), just printing gain instead of
// detector position.
//
// UNLIKE geometry, gain/gainCorrection is NOT constant across Run 22 --
// see README.md. This macro's job is to demonstrate the retrieval
// mechanism for one run at a time, not to find "the" Run 22 value.
//
// calibRun=0 (default) leaves St_db_Maker on BFC's own nominal "y2023"
// chain-flag date, same meaning/convention as sim_to_jet's runSimBfc.C
// calibRun parameter (this macro reuses that exact same run->date
// resolution mechanism, see RunNumberToDateTime() below).

class StChain;
StChain *chain;

// Fallback only -- see RunNumberToDateTime() below, which tries the real
// RunLog DB first and only falls back to this arithmetic if unreachable.
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
// DB (same mechanism StFcsDbMaker/macro/fcsTimedepGainCorr_db.C uses, and
// the same one GeometryPlayground/get_fcs_geometry.C and sim_to_jet's
// runSimBfc.C use), rather than inferring it from the run number's
// day-of-year encoding. Falls back to RunNumberToDateArithmetic() (with
// itim=0) if the DB is unreachable.
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
         << " (verified consistent with RunLog for 3 test runs, but this run itself was not"
         << " cross-checked against the DB just now) ***" << endl;
}

void get_fcs_gain(Int_t calibRun=22359013){
  gROOT->LoadMacro("bfc.C");

  // No fzin in this chain, so no input file is ever opened -- gain comes
  // entirely from the conditions DB, no simulated/real event needed.
  bfc(-1, "y2023,db,fcsDb", "");

  St_db_Maker *dbMk = (St_db_Maker*) chain->GetMaker("db");
  if(calibRun != 0){
    int idat, itim;
    RunNumberToDateTime(calibRun, idat, itim);
    dbMk->SetDateTime(idat, itim);
    cout << "Pinned conditions DB to calibRun=" << calibRun
         << " -> date=" << idat << " time=" << itim << endl;
  }

  chain->Init();

  StFcsDbMaker* fcsdbmkr = (StFcsDbMaker*) chain->GetMaker("fcsDbMkr");
  // InitRun() isn't triggered by chain->Init() alone (no event source to
  // signal a run-number change), so drive it by hand -- same pattern used
  // throughout this project (GeometryPlayground, TriggerPlayground).
  fcsdbmkr->InitRun(calibRun != 0 ? calibRun : 1);

  StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");

  // Representative channel per detector (id=50, arbitrary mid-range
  // channel, same one runSimBfc.C's own sanity-check log line uses --
  // see sim_to_jet/README.md).
  int dets[2] = {0,2};
  const char* names[2] = {"ECAL","HCAL"};
  int id = 50;

  cout.precision(9);
  ofstream fout("fcs_gain.txt");
  fout << "# det name id gain gainCorrection  (calibRun=" << calibRun << ")" << endl;
  for(int i=0;i<2;i++){
    int det = dets[i];
    float gain = fcsdb->getGain(det, id);
    float corr = fcsdb->getGainCorrection(det, id);
    cout << names[i] << " (det" << det << ", id=" << id << "): gain=" << gain
         << " gainCorrection=" << corr << endl;
    fout << det << " " << names[i] << " " << id << " " << gain << " " << corr << endl;
  }
  fout.close();
  cout << "\nWrote fcs_gain.txt" << endl;
}
