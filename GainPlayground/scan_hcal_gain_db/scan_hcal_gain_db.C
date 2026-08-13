// Scans the STAR conditions DB (Calibrations/fcs) for real HCAL
// per-tower gain-correction values across a date range, detects when the
// values actually change (a new real calibration period), and writes
// one correction file per detected period directly into
// data_to_jet/hcal_gain_corrections/ -- regenerating manifest.txt to
// match. This is the automated version of the manual, single-snapshot
// approach GainPlayground/get_fcs_gain.C and
// data_to_jet/hcal_gain_corrections/generate_placeholder.C use; see
// FCSJetPipeline/data_to_jet/hcal_gain_corrections/README.md for how the
// output files this writes get consumed.
//
// WHY THIS EXISTS: as of Aug 2026, HCAL's gainCorrection is a flat 1.0
// placeholder for all of Run 22 (see ../README.md) -- running this
// script today correctly reproduces that (one period, all 1.0). It's
// meant to be RERUN once STAR's FCS calibration group determines real
// per-tower HCAL corrections: at that point it will automatically detect
// the real period boundaries and produce real correction files, with no
// manual DB inspection needed.
//
// METHOD: samples getGainCorrection(det,id) for EVERY real HCAL channel
// at regular date intervals (stepDays) across [startDate, endDate), via
// the same live StFcsDb call GainPlayground/get_fcs_gain.C uses -- NOT
// raw SQL blob parsing. (StDbLib's exact multi-entry conflict-resolution
// rule for overlapping DB validity windows was investigated and left
// unresolved in SpinPlayground/README.md; always asking the real
// framework API directly, rather than re-deriving its resolution logic
// in SQL, sidesteps that whole question -- this is guaranteed to match
// what the production pipeline itself sees, by construction.) Repeated
// SetDateTime()+InitRun() calls within one process correctly refresh
// StFcsDb's cached values -- verified directly (Aug 2026) by scanning
// two dates known to give different ECAL gainCorrection and confirming
// both the change and the reversion back to the first date's value.
// Consecutive samples with an identical full per-channel map are
// coalesced into one period; a period boundary is drawn wherever ANY
// channel's value changes between consecutive samples.
//
// A note on step size: the reference STAR macro this project's DB
// mechanisms are cross-checked against (StFcsDbMaker/macro/
// fcsTimedepGainCorr_db.C) shows real Run 22 ECAL calibration periods
// spaced weeks to ~2 months apart -- stepDays=7 (the default) comfortably
// resolves that granularity without an excessive number of DB queries
// (~600 tables/channels queried per sampled date). Shrink it if a
// change is suspected to have happened faster than that.
//
// Usage:
//   root4star -b -q scan_hcal_gain_db.C(20211201,20220701,7)
//
// CAUTION: overwrites manifest.txt and adds/updates correction files
// directly in data_to_jet/hcal_gain_corrections/. Backs up any existing
// manifest.txt to manifest.txt.bak first (overwriting any previous
// backup) -- this repo has no version control, so that single backup is
// the only undo available. The correction-factor files themselves are
// NOT backed up (they're regenerated deterministically from the DB by
// this same script, so there's nothing to lose by overwriting them).

class StChain;
StChain *chain;

const int kMaxChannels = 600; // > 520 real HCAL channels (260/side), safety margin
const char* kOutDir = "../../data_to_jet/hcal_gain_corrections";

int gNCh = 0;
int gEhp[kMaxChannels], gNs[kMaxChannels], gDep[kMaxChannels], gCh[kMaxChannels];
int gDet[kMaxChannels], gId[kMaxChannels];
float gPrevValues[kMaxChannels];
float gCurValues[kMaxChannels];

// idat is YYYYMMDD; TDatime handles month/year rollovers correctly
// (plain day+=stepDays arithmetic would not, near month/year boundaries).
int NextDate(int idat, int stepDays){
  int year = idat/10000, month = (idat/100)%100, day = idat%100;
  TDatime dt(year, month, day, 0, 0, 0);
  UInt_t epoch = dt.Convert();
  epoch += stepDays*24*3600;
  TDatime dt2(epoch);
  return dt2.GetDate();
}

void WritePeriod(int pStart, int pEnd, float* values, ofstream &fman){
  char fname[256];
  sprintf(fname, "hcal_gain_corr_%d.txt", pStart);
  TString path = TString(kOutDir) + "/" + fname;
  ofstream fout(path.Data());
  fout << "# Real HCAL retroactive gain-correction factors for [" << pStart << "," << pEnd << ")" << endl;
  fout << "# Format matches StFcsDb's own gain-correction text format (StFcsDb::readGainCorrFromText()): ehp ns dep ch factor" << endl;
  fout << "# Generated via GainPlayground/scan_hcal_gain_db/scan_hcal_gain_db.C" << endl;
  for(int i=0;i<gNCh;i++){
    fout << gEhp[i] << " " << gNs[i] << " " << gDep[i] << " " << gCh[i] << " " << values[i] << endl;
  }
  fout.close();
  cout << "Wrote " << path << " (" << gNCh << " channels, period [" << pStart << "," << pEnd << "))" << endl;
  fman << pStart << " " << pEnd << " " << fname << endl;
}

void scan_hcal_gain_db(int startDate=20211201, int endDate=20220701, int stepDays=7){
  gROOT->LoadMacro("bfc.C");

  // No fzin in this chain, so no input file is ever opened -- gain comes
  // entirely from the conditions DB, no simulated/real event needed
  // (same minimal chain as get_fcs_gain.C/get_fcs_geometry.C).
  bfc(-1, "y2023,db,fcsDb", "");

  St_db_Maker *dbMk = (St_db_Maker*) chain->GetMaker("db");
  chain->Init();

  StFcsDbMaker* fcsdbmkr = (StFcsDbMaker*) chain->GetMaker("fcsDbMkr");
  StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");

  // InitRun() isn't triggered by chain->Init() alone (no event source to
  // signal a run-number change) -- drive it by hand once here so
  // StFcsDb's channel map (needed by getDepfromId() below) is actually
  // populated before we enumerate channels. The date this happens to use
  // doesn't matter for the channel MAP itself (stable across Run 22),
  // only for the gainCorrection VALUES sampled in the date loop further
  // down, which each get their own fresh InitRun() call anyway.
  fcsdbmkr->InitRun(1);

  // Enumerate every real HCAL channel once (the channel MAP itself, i.e.
  // which (ehp,ns,dep,ch) exist, comes from Geometry/fcs, stable across
  // all of Run 22 -- see GeometryPlayground/README.md -- so this doesn't
  // need to be redone per sampled date, only the gainCorrection VALUES
  // at each channel do).
  int dets[2] = {2, 3}; // kFcsHcalNorthDetId, kFcsHcalSouthDetId (StEnumerations.h)
  gNCh = 0;
  for(int i=0;i<2;i++){
    int det = dets[i];
    int n = fcsdb->maxId(det);
    for(int id=0; id<n; id++){
      int ehp, ns, crt, slt, dep, ch;
      fcsdb->getDepfromId(det, id, ehp, ns, crt, slt, dep, ch);
      gEhp[gNCh]=ehp; gNs[gNCh]=ns; gDep[gNCh]=dep; gCh[gNCh]=ch;
      gDet[gNCh]=det; gId[gNCh]=id;
      gNCh++;
    }
  }
  cout << "Scanning " << gNCh << " HCAL channels from " << startDate << " to " << endDate
       << " in " << stepDays << "-day steps" << endl;

  TString manifestPath = TString(kOutDir) + "/manifest.txt";
  TString backupPath = manifestPath + ".bak";
  if ( gSystem->AccessPathName(manifestPath.Data()) == 0 ){ // 0 == exists
    gSystem->CopyFile(manifestPath.Data(), backupPath.Data(), kTRUE);
    cout << "Backed up existing manifest.txt -> " << backupPath << endl;
  }
  ofstream fman(manifestPath.Data());
  fman << "# Auto-generated by GainPlayground/scan_hcal_gain_db/scan_hcal_gain_db.C" << endl;
  fman << "# startDate endDate filename  (YYYYMMDD, [startDate,endDate))" << endl;
  fman << "# Regenerate by rerunning that script -- do not hand-edit, it will be overwritten." << endl;

  bool havePrev = false;
  int periodStart = startDate;
  int idat = startDate;
  int nSamples = 0;
  while ( idat < endDate ){
    dbMk->SetDateTime(idat, 0);
    fcsdbmkr->InitRun(1); // run number arg only affects StFcsDb's own logging; the date set above is what matters
    for(int i=0;i<gNCh;i++) gCurValues[i] = fcsdb->getGainCorrection(gDet[i], gId[i]);
    nSamples++;

    if ( havePrev ){
      bool changed = false;
      for(int i=0;i<gNCh;i++){
        if ( gCurValues[i] != gPrevValues[i] ){ changed = true; break; }
      }
      if ( changed ){
        cout << "Change detected at date=" << idat << " -- closing period starting " << periodStart << endl;
        WritePeriod(periodStart, idat, gPrevValues, fman);
        periodStart = idat;
      }
    }
    for(int i=0;i<gNCh;i++) gPrevValues[i] = gCurValues[i];
    havePrev = true;
    idat = NextDate(idat, stepDays);
  }
  // Final period, extended to endDate regardless of the last step's exact landing date.
  WritePeriod(periodStart, endDate, gPrevValues, fman);
  fman.close();

  cout << "Sampled " << nSamples << " date(s). Wrote " << manifestPath << endl;
}
