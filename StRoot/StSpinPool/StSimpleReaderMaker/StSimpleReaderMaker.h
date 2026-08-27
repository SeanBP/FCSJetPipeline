//Class StSimpleReaderMaker

#ifndef StSimpleReaderMaker_def
#define StSimpleReaderMaker_def

#include "StMaker.h"
#include "TString.h"
#include <vector>

class StMuDstMaker ;
class TFile        ;
class TTree        ;
class StFcsDb      ;
class StEpdGeom    ;
class StFcsTriggerSimMaker ;
class fcs_trg_base ;
class StSpinDbMaker ;

class StSimpleReaderMaker : public StMaker
{
  
 private:

  StMuDstMaker* mMuDstMaker ;                      //  Make MuDst pointer available to member functions
  StFcsDb* mFcsDb = 0 ;
  StEpdGeom* mEpdgeo = 0 ;

  //Output file and tree
  TFile* out_file;
  TTree* out_tree;
  TString mOutputFileName;

  UInt_t        mEventsProcessed ;                 //  Number of Events read and processed

  //TTree Branch variables
  int Cal_nhits;
  int Cal_detid[5000];
  int Cal_hitid[5000];
  int Cal_adcsum[5000];
  float Cal_hit_energy[5000];
  float Cal_hit_posx[5000];
  float Cal_hit_posy[5000];
  float Cal_hit_posz[5000];

  int Cal_nclus;
  int Cal_clus_detid[100];
  int Cal_clus_ntowers[100];
  float Cal_clus_energy[100];
  float Cal_clus_loc_x[100];
  float Cal_clus_loc_y[100];
  float Cal_clus_x[100];
  float Cal_clus_y[100];
  float Cal_clus_z[100];

  int Trk_ntrks;
  float Trk_px[500];
  float Trk_py[500];
  float Trk_pz[500];
  int Trk_charge[500];
  float Trk_chi2[500];
  float Trk_ndf[500]; 
  int Trk_nseedpoints[500];
  int Trk_nfitpoints[500];
  float Trk_dca_x[500];
  float Trk_dca_y[500];
  float Trk_dca_z[500];
  int Trk_vtxindex[500];
  float Trk_proj_ecal_x[500];
  float Trk_proj_ecal_y[500];
  float Trk_proj_ecal_z[500];
  float Trk_proj_hcal_x[500];
  float Trk_proj_hcal_y[500];
  float Trk_proj_hcal_z[500];

  int mcpart_num;
  //int mcpart_index[1000]; //Counts particles sequentially
  int mcpart_geid[1000];
  int mcpart_idVtx[1000];
  int mcpart_idVtxEnd[1000];
  float mcpart_px[1000];
  float mcpart_py[1000];
  float mcpart_pz[1000];
  float mcpart_E[1000];
  int mcpart_charge[1000];
  float mcpart_Vtx_x[1000];
  float mcpart_Vtx_y[1000];
  float mcpart_Vtx_z[1000];
  float mcpart_VtxEnd_x[1000];
  float mcpart_VtxEnd_y[1000];
  float mcpart_VtxEnd_z[1000];

  // Boolean flags: one per known FCS trigger name from TriggerIDs.txt,
  // resolved event-by-event either from the fired trigger id(s) + run
  // number (real data, see FcsTriggerDefs.h) or, if SetTriggerSim() was
  // called, from the fcs_trg_base trigger-simulation engine run on this
  // event's FCS hits (for simulated events, which have no real fired
  // trigger id -- see SetTriggerSim()). kNTrigFlags must match
  // kNTrigFlags in FcsTriggerDefs.h.
  enum { kNTrigFlags = 64 };
  int Trig_flag[kNTrigFlags];

  // Simulated-trigger mode: mTrgSim's fcs_trg_base engine only computes
  // 20 of the 64 named flags (see StSimpleReaderMaker.cxx); this maps
  // each of those 20 to its index in Trig_flag[]/kTrigFlagName, built
  // once in Init(). The rest of Trig_flag[] stays 0 in this mode --
  // there is no way to simulate the other (mostly historical-naming-era)
  // triggers from a physics-condition simulation.
  StFcsTriggerSimMaker* mTrgSim = 0 ;
  fcs_trg_base* mTrgBase = 0 ;
  Int_t mTrgLastRunNumber = -1 ;
  enum { kNSimTrg = 21 };
  int mSimFlagIndex[kNSimTrg] ;

  // Per-event spin configuration for real data (StSpinDbMaker's spin4bits
  // combined code, see spin4usingBX7()): -1 means "no data" -- either no
  // StSpinDbMaker was attached (SetSpinDb() not called, e.g. sim_to_jet,
  // where simulated events have no real polarization to report), the
  // offline spin DB (Calibrations/rhic) has no entry covering this run
  // (real DB coverage gaps exist -- see SpinPlayground/
  // run_number_to_spin_coverage.sh), or this specific bunch crossing
  // wasn't resolved even though the run is otherwise covered (e.g. an
  // abort gap). Always -1 for simulated events -- there is no way to
  // simulate a real polarization pattern the way SetTriggerSim() computes
  // a simulated trigger decision, since spin state isn't a physics
  // condition of the generated event, it's a property of which real
  // bunch happened to collide.
  StSpinDbMaker* mSpinDb = 0 ;
  Int_t mSpinLastRunNumber = -1 ;
  int Spin_config ;

  // MC cross-section bookkeeping (sim_to_jet only -- see SetMCXSec()).
  // Constant for the whole job/file (one generator run = one ptHatMin bin =
  // one sigmaGen), so these are set once via SetMCXSec() and just repeated
  // on every Fill() rather than recomputed per event. -1 is the "not MC /
  // not set" sentinel, same convention as Spin_config for real data (there
  // SetMCXSec() is never called, so these three stay -1 for every event).
  float mc_sigma_pb = -1;
  float mc_sigma_err_pb = -1;
  int mc_n_gen = -1;

  // Generator-level ptHatMin cut for this job (sim_to_jet only -- see
  // SetMCPtHatMin()). GeV; matches the ckin3 argument passed to
  // Pythia8()/Pythia6() in starsim_pythia8_filter.C/starsim_pythia6_filter.C
  // (the ${ptcut} shell variable in sim_to_jet.xml). Constant for the whole
  // job/file, same convention/sentinel as mc_sigma_pb above.
  float mc_pthatmin_gev = -1;

  // Retroactive HCAL gain-correction machinery (real data only -- see
  // FCSJetPipeline/data_to_jet/hcal_gain_corrections/README.md).
  // StFcsDb's HCAL gainCorrection is currently a flat 1.0 placeholder (no
  // real per-tower HCAL calibration exists yet -- see
  // GainPlayground/README.md); this multiplies each HCAL hit's energy
  // (already read from the muDst -- real data's calibration was baked in
  // once by official production and is never touched otherwise, see
  // Cal_hit_energy's comment in Make()) by a per-channel factor looked up
  // from a manifest of date-ranged correction files. This means
  // already-produced SimpleTree files can be corrected for the real HCAL
  // calibration later just by rerunning with an updated manifest/file --
  // no muDst-level reprocessing needed. Not wired into sim_to_jet:
  // simulated events already get real, live gainCorrection straight from
  // StFcsDb (StFcsFastSimulatorMaker calls getGain()/getGainCorrection()
  // directly), so there is nothing to retroactively fix there.
  struct HcalRetroCorrPeriod { int startDate; int endDate; TString filename; };
  std::vector<HcalRetroCorrPeriod> mHcalCorrPeriods;
  bool mUseHcalRetroCorr = false;
  TString mHcalCorrDir;
  TString mHcalCorrLoadedFile;
  float mHcalCorrFactor[2][24][32]; // [ns][dep][ch] -- kFcsNorthSouth/kFcsMaxDepBd/kFcsMaxDepCh (StEnumerations.h)
  int mHcalCorrCurrentDate = -1;

  void loadHcalCorrManifest();
  void loadHcalCorrFile(const TString& filename);
  float getHcalRetroCorr(int det, int id);

 protected:

 public:

  StSimpleReaderMaker(StMuDstMaker* maker) ;       //  Constructor
  virtual          ~StSimpleReaderMaker( ) ;       //  Destructor

  Int_t Init    ( ) ;                       //  Initiliaze the analysis tools ... done once
  Int_t Make    ( ) ;                       //  The main analysis that is done on each event
  Int_t Finish  ( ) ;                       //  Finish the analysis, close files, and clean up.

  void SetOutputFileName(TString name) {mOutputFileName = name;} // Make name available to member functions

  // Simulated-events mode (sim_to_jet): use trgSim's fcs_trg_base engine
  // to compute Trig_flag[] from this event's FCS hits, instead of the
  // real-data trigger-id lookup. trgSim must already be constructed and
  // configured (setTrigger/setSimMode(1)/setThresholdFile/etc.) and
  // still be in the chain at Init() time (its own Init()/InitRun() do
  // real setup we rely on -- gains, pedestals, stage_version -- even
  // though we bypass its Make(), which segfaults reading a muDst with no
  // StEvent present; the caller must SetActive(kFALSE) on it after
  // chain->Init() for the same reason data_to_jet's runMudst.C does).
  void SetTriggerSim(StFcsTriggerSimMaker* trgSim) { mTrgSim = trgSim; }

  // Real-data mode (data_to_jet): attach an StSpinDbMaker (already
  // constructed and in the chain, same pattern as SetTriggerSim's
  // trgSim) to resolve Spin_config per event from the event's real
  // bunch-crossing id. Not called by sim_to_jet -- see Spin_config's
  // comment above for why simulated events always report -1 regardless.
  void SetSpinDb(StSpinDbMaker* spinDb) { mSpinDb = spinDb; }

  // Sim-only mode (sim_to_jet): record this job's generator cross-section
  // normalization on the SimpleTree, so it survives all the way to the
  // JetTree (see JetMatcher.cpp) for per-event MC weighting. sigma_pb/
  // sigmaErr_pb are the generator's total cross section (and its error) for
  // this job's ptHatMin/ptHatMax phase-space cut, in pb; nGen is the number
  // of hard-scattering trials that cross section corresponds to (PYTHIA's
  // nAccepted/sumWeightGen -- the full generated sample, NOT the smaller
  // count that survives FcsJetFilter: the filter only decides what gets
  // written to disk, it doesn't change how much luminosity each surviving
  // event represents). Not called by data_to_jet -- see mc_sigma_pb's
  // comment in the header for why real data keeps the -1 sentinel.
  void SetMCXSec(double sigma_pb, double sigmaErr_pb, int nGen) {
    mc_sigma_pb = sigma_pb; mc_sigma_err_pb = sigmaErr_pb; mc_n_gen = nGen;
  }

  // Sim-only mode (sim_to_jet): record this job's generator-level ptHatMin
  // cut (ckin3 in starsim_pythia8_filter.C/starsim_pythia6_filter.C), so it
  // survives to the JetTree (see JetMatcher.cpp) alongside mc_sigma_pb. Not
  // called by data_to_jet -- see mc_pthatmin_gev's comment in the header.
  void SetMCPtHatMin(double pthatmin_gev) { mc_pthatmin_gev = pthatmin_gev; }

  // Real-data mode (data_to_jet): enable retroactive HCAL gain
  // correction. manifestDir is the folder containing manifest.txt and
  // the correction-factor files it references (see
  // data_to_jet/hcal_gain_corrections/README.md for the format). Not
  // called by sim_to_jet -- see mUseHcalRetroCorr's comment above.
  void SetHcalRetroactiveGainCorr(const char* manifestDir) { mUseHcalRetroCorr = true; mHcalCorrDir = manifestDir; }

  ClassDef(StSimpleReaderMaker,1)                  //  Macro for CINT compatability
    
};

#endif

