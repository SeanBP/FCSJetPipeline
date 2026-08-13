// Generates a placeholder HCAL retroactive-gain-correction file, entirely
// filled with 1.0 (no-op), covering every real HCAL channel -- north and
// south, every (dep,ch) StFcsDb actually maps to a channel.
//
// Enumerates channels the same way StFcsDb itself does when filling its
// own gain-correction arrays from the DB (StFcsDb.cxx: `for(int id=0;
// id<maxId(det); id++) getDepfromId(det,id,...)`), rather than guessing
// dep/ch loop bounds by hand -- this guarantees the output exactly
// matches the real channel map for whatever geometry/DB era this is run
// against, with no invalid or missing rows.
//
// Minimal BFC chain (DB + StFcsDbMaker only, no simulation/input file
// needed), same pattern as GeometryPlayground/get_fcs_geometry.C and
// GainPlayground/get_fcs_gain.C.
//
// Usage: root4star -b -q generate_placeholder.C(<run_number>,"<outfile>")
// calibRun only matters in that it selects which StFcsDb channel map
// (nColumn/nRow/maxId) is active -- irrelevant in practice, since this
// map has been stable for all of Run 22 (see GeometryPlayground's
// symmetry/stability finding), but a real run number is passed for
// consistency with the rest of this project's tools.

class StChain;
StChain *chain;

void generate_placeholder(Int_t calibRun=22359013, const char* outfile="hcal_gain_corr_run22_placeholder.txt"){
  gROOT->LoadMacro("bfc.C");
  bfc(-1, "y2023,db,fcsDb", "");

  St_db_Maker *dbMk = (St_db_Maker*) chain->GetMaker("db");
  chain->Init();

  StFcsDbMaker* fcsdbmkr = (StFcsDbMaker*) chain->GetMaker("fcsDbMkr");
  fcsdbmkr->InitRun(calibRun);

  StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");

  ofstream fout(outfile);
  fout << "# Placeholder HCAL retroactive gain-correction factors -- ALL 1.0"
       << " (no-op) until real per-tower HCAL corrections are measured."
       << " See ../../GainPlayground/README.md and this folder's README.md." << endl;
  fout << "# Format matches StFcsDb's own gain-correction text format"
       << " (StFcsDb::readGainCorrFromText()): ehp ns dep ch factor" << endl;
  fout << "# Generated via generate_placeholder.C(" << calibRun << ")" << endl;

  int dets[2] = {2, 3}; // kFcsHcalNorthDetId, kFcsHcalSouthDetId (StEnumerations.h)
  int nWritten = 0;
  for(int i=0;i<2;i++){
    int det = dets[i];
    int n = fcsdb->maxId(det);
    for(int id=0; id<n; id++){
      int ehp, ns, crt, slt, dep, ch;
      fcsdb->getDepfromId(det, id, ehp, ns, crt, slt, dep, ch);
      fout << ehp << " " << ns << " " << dep << " " << ch << " 1.0" << endl;
      nWritten++;
    }
  }
  fout.close();
  cout << "Wrote " << nWritten << " channel rows to " << outfile << endl;
}
