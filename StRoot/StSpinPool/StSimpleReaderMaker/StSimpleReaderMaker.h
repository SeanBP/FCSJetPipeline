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
class TH1          ;
class TRandom3     ;

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

  // Event z-vertex (Sean 2026-11-16; ChooseVertex()): the EPD vertex, else
  // the BBC TAC vertex, as measured (no rescaling to the TPC scale: the
  // methods read ~0.8 z_TPC, but rescaling widens the error per event more
  // than it removes bias; the simulation reproduces the same response), else
  // the nominal z (mNominalVz, the beam-spot mean; with the beamline x, y in
  // shared/JetParameters.h, the nominal IP); no |z| cut here. Data: from the MuDst EPD hits (or, where the MuDst has none,
  // the hits StEpdHitMaker rebuilt from the trigger data into StEvent) and the
  // BBC TAC difference. Simulation: epd_vz and bbc_vz_tac emulated from the
  // truth vertex by SetVertexSmear(), then the same choice. The old
  // TPC / VPD / StBbcTriggerDetector::zVertex() choice is gone: st_fwd has no
  // TPC or VPD, and zVertex() uses pre-2003 BBC constants.
  Float_t event_vertex_z = -9999;
  Int_t event_vertex_source = 0; // 0=not set (sim without SetVertexSmear()), 1=EPD, 2=BBC TAC, 3=nominal
  Int_t epd_source = 0;          // EPD hits used: 0=none, 1=MuDst, 2=rebuilt (StEvent)
  Float_t mc_vz = -9999;         // simulation: truth primary-vertex z (first MC vertex)
  Float_t mc_vx = -9999, mc_vy = -9999; // simulation: its x, y (the generated beamline)
  // Simulation response of each method, z_raw = a + s z_true: the peak
  // (Gaussian core) of z_method per 20 cm z_TPC slice, linear fit over
  // |z_TPC| < 100 cm (P25ib st_physics, 40 runs, 53k events;
  // vertex_tpc_v1/code/make_smear.C). [0] EPD a = -1.31 cm, s = 0.8136;
  // [1] BBC TAC a = +3.48 cm, s = 0.7884. Constructor defaults; the
  // SetVertexSmear() file's "response" vector overrides.
  double mRespA[2], mRespS[2];
  double mNominalVz;             // cm: beam-spot mean (TPC, Gaussian fit: -7.68 cm); SetNominalVz()
  // Simulation vertex emulation (SetVertexSmear()): joint fire probabilities
  // of (EPD, BBC TAC) and the measured method - TPC residuals, from data
  bool mVertexSmear = false;
  double mFireProb[4];                // both, EPD only, BBC only, neither (set in the constructor)
  TH1* mVtxRes[2];                    // EPD - TPC, BBC TAC - TPC (cm)
  TRandom3* mVtxRng = 0;
  double SampleResidual(int m);
  void ChooseVertex();

  // Vertex test (2026-11-09, after D. Kapukchyan's StFwdAnaVertex /
  // StFwdAnaEpdQaAndVert): the raw per-detector z, no |z| cut, -9999 = none.
  //  bbc_tac_diff: StTriggerData::bbcTimeDifference() - 4096 (-4096 = no TAC)
  //  bbc_vz_tac:   bbc_tac_diff * -0.2475 cm (0.0165 ns * 30 cm/ns / 2), when
  //                bbc_tac_diff != 0
  //  bbc_vz_det:   StMuEvent bbcTriggerDetector().zVertex() (the source of
  //                event_vertex_source 3)
  //  epd_vz:       (avg TAC west - avg TAC east) * 0.2475 cm over hits with
  //                1 < nMIP < 15 and TAC > 50, when both averages are > 50
  //  vpd_vz:       StMuEvent::vpdVz() when the VPD has hits
  Int_t bbc_tac_diff = -4096;
  Float_t bbc_vz_tac = -9999, bbc_vz_det = -9999, epd_vz = -9999, vpd_vz = -9999;
  Float_t epd_avgtac_e = 0, epd_avgtac_w = 0;
  Int_t epd_ntac_e = 0, epd_ntac_w = 0;

  // Per-run sqrt(s), GeV (computeSqrtSForRun() in JetParameters.h)
  Double_t event_sqrt_s_gev = -1;

  // Run and event number of the MuDst event (to match other per-event data)
  Int_t event_run = 0;
  Int_t event_id = -1;

  // FCS stage-3 output bits as the trigger control unit saw them: before
  // prescale, so a condition that fired but did not record the event is
  // still set (the Trig_flag[] list holds only triggers that recorded it).
  // Data: lastDSM(2) + (lastDSM(5) << 16) of the event's trigger data, as
  // StFcsTrgQaMaker; simulation (SetTriggerSim()): fcs_trg_base's emulated
  // output, the same layout. Bits 0-2 HAD0-2, 3-5 EM0-2, 6 JP2, 7-9
  // JPA1/JPBC1/JPDE1, 10-12 JPA0/JPBC0/JPDE0, 13 DiJP, 14 DiJPAsy, 15 DiELEA,
  // 16-23 north ELE0-2/EM3/EHT/HHT/ETOT/HTOT, 24-31 the same south.
  // fcs_tcu_source: 0 none, 1 hardware (data), 2 emulated (simulation).
  UInt_t fcs_tcu_bits = 0;
  Int_t fcs_tcu_source = 0;

  // FCS trigger 2x2 sums from the offline hit energies (data and simulation
  // alike), in the stage-2 layout of the trigger algorithm
  // (StRoot/RTS/src/TRG_FCS/stage_2_202207.cxx): trg_e2x2[ns][r][c] is the
  // ECal 2x2 tower sum at stage-2 row r (0-15), column c (0-9);
  // trg_h2x2[ns][r][c] the HCal one (rows 0-9, columns 0-5). E_T in GeV
  // (hit energy x sin theta from the origin, negative energies clipped at 0
  // as the stage-0 pedestal subtraction), not truncated: the stage-1 count
  // truncation (0.0316 GeV, saturating at 255) is applied downstream
  // (RecoJets). Towers outside the trigger (ECal DEP >= 20) are not in it.
  // The tree "trg_geom" (one entry, written at Finish) has each 2x2's centre:
  // the mean STAR x, y, z of its towers.
  Float_t trg_e2x2[2][16][10];
  Float_t trg_h2x2[2][10][6];
  Float_t trg_e2x2_x[2][16][10], trg_e2x2_y[2][16][10], trg_e2x2_z[2][16][10];
  Float_t trg_h2x2_x[2][10][6], trg_h2x2_y[2][10][6], trg_h2x2_z[2][10][6];
  bool mTrgGeomDone = false;
  void buildTrgGeom();
  // stage-2 position of a tower from its DEP board and channel; false if
  // the tower is not in the stage-2 sums
  static bool trg2x2Of(int ehp, int dep, int ch, int &r, int &c);

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

  // Event vertex choice (see event_vertex_z): the nominal z, cm, for events
  // with neither method (default: the TPC beam-spot mean)
  void SetNominalVz(double cm) { mNominalVz = cm; }
  // Simulation: emulate the EPD and BBC TAC vertices from the truth vertex.
  // The file holds TH1 "fire" (4 bins: both, EPD only, BBC only, neither;
  // probabilities, in the st_fwd analysis sample) and TH1 "res_epd", "res_bbc"
  // (the scatter of the raw method z about a + s z_TPC, cm) and TVectorD
  // "response" (a_epd, s_epd, a_bbc, s_bbc); each event draws which fire and
  // sets z_raw = a + s mc_vz + a sampled residual, then the data's choice. Returns
  // false (and leaves the emulation off) if the file or a histogram is missing.
  bool SetVertexSmear(const char* file, unsigned int seed);

  // Data: enable the HCAL gain correction; manifestDir holds manifest.txt and
  // the correction files
  void SetHcalRetroactiveGainCorr(const char* manifestDir) { mUseHcalRetroCorr = true; mHcalCorrDir = manifestDir; }

  ClassDef(StSimpleReaderMaker,1)                  //  Macro for CINT compatability
    
};

#endif

