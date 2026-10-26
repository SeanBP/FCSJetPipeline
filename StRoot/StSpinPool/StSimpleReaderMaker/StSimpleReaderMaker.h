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

  // FCS points (StMuFcsPoint): single-photon candidates, used for reco jets
  // (-F 1) and the per-jet photon count. Up to ~50 per event observed.
  int Cal_npts;
  int Cal_pt_detid[200];
  float Cal_pt_energy[200];
  float Cal_pt_loc_x[200];
  float Cal_pt_loc_y[200];
  float Cal_pt_x[200];
  float Cal_pt_y[200];
  float Cal_pt_z[200];
  int Cal_pt_parentClusterId[200];

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

  // One flag per FCS trigger name (FcsTriggerDefs.h): from the fired trigger
  // ids and run number in data, or from fcs_trg_base in simulation
  // (SetTriggerSim()).
  enum { kNTrigFlags = 64 };
  int Trig_flag[kNTrigFlags];

  // Simulation: fcs_trg_base computes 20 of the named flags; this maps them
  // to Trig_flag[] indices (built in Init()). The others stay 0.
  StFcsTriggerSimMaker* mTrgSim = 0 ;

  // See SetTriggerFilter(). mTrigFilterIndex holds the Trig_flag[] index
  // of each requested name, resolved once in Init().
  TString mTrigFilterNames ;
  std::vector<int> mTrigFilterIndex ;
  fcs_trg_base* mTrgBase = 0 ;
  Int_t mTrgLastRunNumber = -1 ;
  enum { kNSimTrg = 21 };
  int mSimFlagIndex[kNSimTrg] ;

  // Spin configuration (StSpinDbMaker spin4bits, see spin4usingBX7()); -1 if
  // no spin DB is attached (simulation), the run is not covered, or the
  // bunch crossing is unresolved.
  StSpinDbMaker* mSpinDb = 0 ;

  bool mFcsFromStEvent = false ; // see SetFcsFromStEvent()
  Int_t mSpinLastRunNumber = -1 ;
  int Spin_config ;

  // Event z-vertex, data only; -9999 = none. The first of TPC primary vertex
  // (ranking > 0), VPD, BBC with |z| < 80 cm. In P24ia (st_fwd) only BBC is
  // populated (18-34% of events); the ZDC vertex is not used.
  Float_t event_vertex_z = -9999;
  Int_t event_vertex_source = 0; // 0=none, 1=TPC, 2=VPD, 3=BBC

  // Per-run sqrt(s), GeV (computeSqrtSForRun() in JetParameters.h)
  Double_t event_sqrt_s_gev = -1;

  // Simulation only, constant per job (SetMCJobStats()): PYTHIA's cross
  // section (mb) and number of trials. -1 if unset (data).
  Double_t sigma_job = -1;
  Int_t N_job = -1;

  // Simulation only: the job's generator ptHat range, GeV (-1 = unset / no
  // upper cut)
  float mc_pthatmin_gev = -1;
  float mc_pthatmax_gev = -1;

  // Data only: per-channel HCAL gain correction applied to hit energies,
  // looked up from a manifest of date-ranged files
  // (data_to_jet/hcal_gain_corrections/README.md), so SimpleTrees can be
  // recalibrated without reprocessing the MuDst. Simulation takes its gains
  // from StFcsDb directly.
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

  // Simulation: compute Trig_flag[] with trgSim's fcs_trg_base engine.
  // trgSim must be configured and in the chain at Init() (its Init/InitRun
  // set up gains and pedestals); its Make() needs an StEvent, so the caller
  // deactivates it after chain->Init().
  void SetTriggerSim(StFcsTriggerSimMaker* trgSim) { mTrgSim = trgSim; }

  // Comma-separated trigger flag names: only events that fired one are
  // written. Empty (default) writes all; unknown names fail Init().
  void SetTriggerFilter(const char* names) { mTrigFilterNames = names ? names : ""; }

  // Data: StSpinDbMaker (in the chain) used to resolve Spin_config
  void SetSpinDb(StSpinDbMaker* spinDb) { mSpinDb = spinDb; }

  // Simulation: the job's cross section (mb) and total trial count (before
  // FcsJetFilter, which does not change the luminosity per event)
  void SetMCJobStats(double sigmaGen, int nTried) {
    sigma_job = sigmaGen; N_job = nTried;
  }

  // Simulation: the job's generator ptHat range
  void SetMCPtHatMin(double pthatmin_gev) { mc_pthatmin_gev = pthatmin_gev; }
  void SetMCPtHatMax(double pthatmax_gev) { mc_pthatmax_gev = pthatmax_gev; }

  // Write FCS clusters and points from the StEvent FCS collection (reclustered
  // in the chain: StFcsMuToEventMaker, StFcsClusterMaker, StFcsPointMaker)
  // instead of the MuDst's. Hits always come from the MuDst. Points'
  // Cal_pt_parentClusterId is then the cluster's index within its detector.
  void SetFcsFromStEvent(bool on = true) { mFcsFromStEvent = on; }

  // Data: enable the HCAL gain correction; manifestDir holds manifest.txt and
  // the correction files
  void SetHcalRetroactiveGainCorr(const char* manifestDir) { mUseHcalRetroCorr = true; mHcalCorrDir = manifestDir; }

  ClassDef(StSimpleReaderMaker,1)                  //  Macro for CINT compatability
    
};

#endif

