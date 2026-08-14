// Minimal BFC chain (DB + FcsDbMaker only, no GEANT/simulation/input file
// needed) to query StFcsDb's detector-position geometry for a given
// calibration era, and print the 4 corner StarXYZ positions (same
// convention as StFcsHit/StFcsPoint reconstructed positions -- i.e. cell
// centers at shower-max depth, computed via getStarXYZ(det,col,row)) for
// each of ECAL north/south and HCAL north/south.
//
// calibRun=0 (default) leaves St_db_Maker on BFC's own nominal "y2023"
// chain-flag date, same meaning/convention as sim_to_jet's runSimBfc.C
// calibRun parameter.

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
// DB (same mechanism StFcsDbMaker/macro/fcsTimedepGainCorr_db.C uses),
// rather than inferring it from the run number's day-of-year encoding.
// Falls back to RunNumberToDateArithmetic() (with itim=0) if the DB is
// unreachable -- e.g. a SUMS batch worker node with more restricted
// outbound network access than this interactive session; verified
// empirically (Aug 2026) working here, but batch-node behavior is untested.
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
         << " (verified consistent with RunLog for 3 test runs, but this run itself was not"
         << " cross-checked against the DB just now) ***" << endl;
}

// Geometry/fcs/fcsDetectorPosition's real, final survey-based position
// only took effect starting this date. Empirically bisected (Aug 2026,
// GeometryPlayground): xoff was a placeholder value (round numbers,
// -67.399/67.399 for ECAL) through 2021-12-20, then jumped to the final
// calibrated value (-17.399/17.399) and stayed fixed there for the rest
// of Run 22 (checked through the y2023 nominal date). Confirmed directly
// against the raw DB table: the real transition is precise to the
// minute, 2021-12-20 16:30:00 (dataID 17/22 in Geometry_fcs.
// fcsDetectorPosition); this constant rounds that up to the next whole
// day, so it flags a few conservative extra hours on the 20th itself.
// (yoff/y-center settles even earlier, by 2021-12-02, so xoff is the
// binding constraint.) Requesting a run/date before this returns
// geometry that was never the real, final detector position.
const int kFcsSurveyValidFromDate = 20211221;

void get_fcs_geometry(Int_t calibRun=0){
  gROOT->LoadMacro("bfc.C");

  // No fzin in this chain, so no input file is ever opened -- geometry
  // comes entirely from the conditions DB, no simulated/real event needed.
  bfc(-1, "y2023,db,fcsDb", "");

  St_db_Maker *dbMk = (St_db_Maker*) chain->GetMaker("db");
  if(calibRun != 0){
    int idat, itim;
    RunNumberToDateTime(calibRun, idat, itim);
    dbMk->SetDateTime(idat, itim);
    cout << "Pinned conditions DB to calibRun=" << calibRun
         << " -> date=" << idat << " time=" << itim << endl;
    if(idat < kFcsSurveyValidFromDate){
      cout << "*** WARNING: date=" << idat << " is before " << kFcsSurveyValidFromDate
           << ", the date the FCS detector-position survey was finalized in the DB."
           << " Geometry returned for this run is a pre-survey placeholder, NOT the"
           << " real detector position -- do not trust this tool's fiducial output"
           << " for this run. ***" << endl;
    }
  }

  chain->Init();

  StFcsDbMaker* fcsdbmkr = (StFcsDbMaker*) chain->GetMaker("fcsDbMkr");
  // InitRun() isn't triggered by chain->Init() alone (no event source to
  // signal a run-number change), so drive it by hand -- same pattern used
  // for StFcsTriggerSimMaker in the trigger-sim work.
  fcsdbmkr->InitRun(calibRun != 0 ? calibRun : 1);

  StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");

  int dets[4] = {0,1,2,3};
  const char* names[4] = {"ECAL north","ECAL south","HCAL north","HCAL south"};

  // Plain-text dump for fiducial_from_corners.cpp (pure C++, no ROOT/STAR
  // dependency) to consume: one line per corner, "det x y z" in cm, 4
  // corners per detector in a fixed col1row1 -> colMax,row1 ->
  // colMax,rowMax -> col1,rowMax order (a consistent winding around the
  // detector's rectangular footprint).
  ofstream fout("fcs_corners.txt");
  fout << "# det x_cm y_cm z_cm  (calibRun=" << calibRun << ")" << endl;

  for(int i=0;i<4;i++){
    int det = dets[i];
    int nCol = fcsdb->nColumn(det);
    int nRow = fcsdb->nRow(det);
    cout << "\nDetector " << det << " (" << names[i] << "): nCol=" << nCol << " nRow=" << nRow << endl;
    int cols[4] = {1, nCol, nCol, 1};
    int rows[4] = {1, 1, nRow, nRow};
    const char* corner[4] = {"corner(col1,row1)","corner(colMax,row1)","corner(colMax,rowMax)","corner(col1,rowMax)"};
    for(int c=0;c<4;c++){
      StThreeVectorD xyz = fcsdb->getStarXYZ(det, cols[c], rows[c]);
      cout << "  " << corner[c] << ": x=" << xyz.x() << " y=" << xyz.y() << " z=" << xyz.z() << " [cm]" << endl;
      fout << det << " " << xyz.x() << " " << xyz.y() << " " << xyz.z() << endl;
    }
  }
  fout.close();
  cout << "\nWrote fcs_corners.txt" << endl;

  // Per-detector offset (StFcsDb::getDetectorOffset(det,zdepth), full
  // double precision -- default cout precision (6 sig figs) is not enough
  // to tell a real difference from round-off, see the north/south
  // symmetry investigation in README.md).
  cout.precision(17);
  ofstream foff("fcs_detector_offsets.txt");
  foff.precision(17);
  foff << "# det name xoff_cm yoff_cm  (calibRun=" << calibRun << ", full double precision)" << endl;
  for(int i=0;i<4;i++){
    int det = dets[i];
    StThreeVectorD off = fcsdb->getDetectorOffset(det);
    cout << "Detector " << det << " (" << names[i] << ") offset: xoff=" << off.x()
         << " yoff=" << off.y() << endl;
    foff << det << " " << names[i] << " " << off.x() << " " << off.y() << endl;
  }
  foff.close();
  cout << "\nWrote fcs_detector_offsets.txt" << endl;
}
