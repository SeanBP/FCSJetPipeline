// Class StSimpleReaderMaker

#include "StSimpleReaderMaker.h"

#include <iostream>

#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuFcsHit.h"
#include "StMuDSTMaker/COMMON/StMuFcsCollection.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"

#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

#include "StRoot/StEpdUtil/StEpdGeom.h"

#include "TTree.h"
#include "TFile.h"
#include "TDatime.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TString.h"

#include "FcsTriggerDefs.h"

#include <cstring>
#include <stdint.h>
#include <fstream>
#include <sstream>
#include <string>
#include "StSpinPool/StFcsTriggerSimMaker/StFcsTriggerSimMaker.h"
#include "RTS/src/TRG_FCS/fcs_trg_base.h"
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"

ClassImp(StSimpleReaderMaker)                   // Macro for CINT compatibility

// Same 21-name/bit table StFcsTriggerSimMaker::Make() decodes from
// dsm_out (copied here, not modifying the official maker, since we
// bypass its Make() -- see SetTriggerSim() in the header for why).
// Names match "fcs"+ctrg[i] in FcsTriggerDefs.h's 64-flag set for the
// ~20 that have one; DYNoEpd has no counterpart there and is skipped.
static const int kNSimTrgLocal = 21;
static const char* kSimTrgName[kNSimTrgLocal] = {
    "fcsJP2","fcsJPA1","fcsJPA0","fcsJPBC1","fcsJPBC0","fcsJPDE1","fcsJPDE0",
    "fcsDiJP","fcsDiJPAsy","fcsDY","fcsJPsi","DYNoEpd","fcsDYAsy",
    "fcsHad2","fcsHad1","fcsHad0","fcsEM2","fcsEM1","fcsEM0","fcsELE2","fcsEM3"
};

static void DecodeSimTrg(unsigned int dsm_out, bool simTrg[kNSimTrgLocal]){
    for (int i = 0; i < kNSimTrgLocal; i++) simTrg[i] = false;
    simTrg[0]  = (dsm_out>>6)&1;                              // JP2
    simTrg[1]  = (dsm_out>>7)&1;                              // JPA1
    simTrg[2]  = (dsm_out>>10)&1;                             // JPA0
    simTrg[3]  = (dsm_out>>8)&1;                              // JPBC1
    simTrg[4]  = (dsm_out>>11)&1;                             // JPBC0
    simTrg[5]  = (dsm_out>>9)&1;                              // JPDE1
    simTrg[6]  = (dsm_out>>12)&1;                             // JPDE0
    simTrg[7]  = (dsm_out>>13)&1;                             // DiJP
    simTrg[8]  = (dsm_out>>14)&1;                             // DiJPAsy
    simTrg[9]  = ((dsm_out>>18)&1) && ((dsm_out>>26)&1);      // DY
    simTrg[10] = ((dsm_out>>17)&1) && ((dsm_out>>25)&1);      // JPsi
    simTrg[11] = ((dsm_out>>19)&1) && ((dsm_out>>27)&1);      // DYNoEpd
    simTrg[12] = (dsm_out>>15)&1;                             // DYAsy
    simTrg[13] = (dsm_out>>2)&1;                              // Had2
    simTrg[14] = (dsm_out>>1)&1;                              // Had1
    simTrg[15] = (dsm_out>>0)&1;                              // Had0
    simTrg[16] = (dsm_out>>5)&1;                              // EM2
    simTrg[17] = (dsm_out>>4)&1;                              // EM1
    simTrg[18] = (dsm_out>>3)&1;                              // EM0
    simTrg[19] = ((dsm_out>>18)&1) || ((dsm_out>>26)&1);      // ELE2
    simTrg[20] = ((dsm_out>>19)&1) || ((dsm_out>>27)&1);      // EM3
}

StSimpleReaderMaker::StSimpleReaderMaker( StMuDstMaker* maker ) : StMaker("StSimpleReaderMaker")
{ // Initialize and/or zero all public/private data members here.
  out_file = NULL;
  out_tree = NULL;
  mOutputFileName = "";
  mEventsProcessed = 0;

  mEpdgeo = new StEpdGeom;  // EPD geom.

  mMuDstMaker = maker ;     // Pass MuDst pointer to DstAnlysisMaker Class member functions
}

StSimpleReaderMaker::~StSimpleReaderMaker() 
{ // Destroy and/or zero out all public/private data members here.
}

Int_t StSimpleReaderMaker::Init( )
{ // Do once at the start of the analysis

  // FCS DB
  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
  // mFcsDb->setDbAccess(0);
  if (!mFcsDb) {
  	LOG_ERROR << "StSimpleReaderMaker::InitRun Failed to get StFcsDb" << endm;
        return kStFatal;
  }

  if ( mUseHcalRetroCorr ) loadHcalCorrManifest();

  out_file = new TFile(mOutputFileName,"RECREATE");
  out_tree = new TTree("data","Simple Data Tree");

  out_tree->Branch("Cal_nhits",&Cal_nhits,"Cal_nhits/I");
  out_tree->Branch("Cal_detid",Cal_detid,"Cal_detid[Cal_nhits]/I");
  out_tree->Branch("Cal_hitid",Cal_hitid,"Cal_hitid[Cal_nhits]/I");
  out_tree->Branch("Cal_adcsum",Cal_adcsum,"Cal_adcsum[Cal_nhits]/I");
  out_tree->Branch("Cal_hit_energy",Cal_hit_energy,"Cal_hit_energy[Cal_nhits]/F");
  out_tree->Branch("Cal_hit_posx",Cal_hit_posx,"Cal_hit_posx[Cal_nhits]/F");
  out_tree->Branch("Cal_hit_posy",Cal_hit_posy,"Cal_hit_posy[Cal_nhits]/F");
  out_tree->Branch("Cal_hit_posz",Cal_hit_posz,"Cal_hit_posz[Cal_nhits]/F");

  out_tree->Branch("Cal_nclus",&Cal_nclus,"Cal_nclus/I");
  out_tree->Branch("Cal_clus_detid",Cal_clus_detid,"Cal_clus_detid[Cal_nclus]/I");
  out_tree->Branch("Cal_clus_ntowers",Cal_clus_ntowers,"Cal_clus_ntowers[Cal_nclus]/I");
  out_tree->Branch("Cal_clus_energy",Cal_clus_energy,"Cal_clus_energy[Cal_nclus]/F");
  out_tree->Branch("Cal_clus_loc_x",Cal_clus_loc_x,"Cal_clus_loc_x[Cal_nclus]/F");
  out_tree->Branch("Cal_clus_loc_y",Cal_clus_loc_y,"Cal_clus_loc_y[Cal_nclus]/F");
  out_tree->Branch("Cal_clus_x",Cal_clus_x,"Cal_clus_x[Cal_nclus]/F");
  out_tree->Branch("Cal_clus_y",Cal_clus_y,"Cal_clus_y[Cal_nclus]/F");
  out_tree->Branch("Cal_clus_z",Cal_clus_z,"Cal_clus_z[Cal_nclus]/F");
	
  out_tree->Branch("Trk_ntrks",&Trk_ntrks,"Trk_ntrks/I");
  out_tree->Branch("Trk_px",Trk_px,"Trk_px[Trk_ntrks]/F");
  out_tree->Branch("Trk_py",Trk_py,"Trk_py[Trk_ntrks]/F");
  out_tree->Branch("Trk_pz",Trk_pz,"Trk_pz[Trk_ntrks]/F");
  out_tree->Branch("Trk_charge",Trk_charge,"Trk_charge[Trk_ntrks]/I");
  out_tree->Branch("Trk_chi2",Trk_chi2,"Trk_chi2[Trk_ntrks]/F");
  out_tree->Branch("Trk_ndf",Trk_ndf,"Trk_ndf[Trk_ntrks]/F");
  out_tree->Branch("Trk_nseedpoints",Trk_nseedpoints,"Trk_nseedpoints[Trk_ntrks]/I");
  out_tree->Branch("Trk_nfitpoints",Trk_nfitpoints,"Trk_nfitpoints[Trk_ntrks]/I");
  out_tree->Branch("Trk_dca_x",Trk_dca_x,"Trk_dca_x[Trk_ntrks]/F");
  out_tree->Branch("Trk_dca_y",Trk_dca_y,"Trk_dca_y[Trk_ntrks]/F");
  out_tree->Branch("Trk_dca_z",Trk_dca_z,"Trk_dca_z[Trk_ntrks]/F");
  out_tree->Branch("Trk_vtxindex",Trk_vtxindex,"Trk_vtxindex[Trk_ntrks]/I");
  out_tree->Branch("Trk_proj_ecal_x",Trk_proj_ecal_x,"Trk_proj_ecal_x[Trk_ntrks]/F");
  out_tree->Branch("Trk_proj_ecal_y",Trk_proj_ecal_y,"Trk_proj_ecal_y[Trk_ntrks]/F");
  out_tree->Branch("Trk_proj_ecal_z",Trk_proj_ecal_z,"Trk_proj_ecal_z[Trk_ntrks]/F");
  out_tree->Branch("Trk_proj_hcal_x",Trk_proj_hcal_x,"Trk_proj_hcal_x[Trk_ntrks]/F");
  out_tree->Branch("Trk_proj_hcal_y",Trk_proj_hcal_y,"Trk_proj_hcal_y[Trk_ntrks]/F");
  out_tree->Branch("Trk_proj_hcal_z",Trk_proj_hcal_z,"Trk_proj_hcal_z[Trk_ntrks]/F");

  out_tree->Branch("mcpart_num",&mcpart_num,"mcpart_num/I");
  out_tree->Branch("mcpart_geid",mcpart_geid,"mcpart_geid[mcpart_num]/I");
  out_tree->Branch("mcpart_idVtx",mcpart_idVtx,"mcpart_idVtx[mcpart_num]/I");
  out_tree->Branch("mcpart_idVtxEnd",mcpart_idVtxEnd,"mcpart_idVtxEnd[mcpart_num]/I");
  out_tree->Branch("mcpart_px",mcpart_px,"mcpart_px[mcpart_num]/F");
  out_tree->Branch("mcpart_py",mcpart_py,"mcpart_py[mcpart_num]/F");
  out_tree->Branch("mcpart_pz",mcpart_pz,"mcpart_pz[mcpart_num]/F");
  out_tree->Branch("mcpart_E",mcpart_E,"mcpart_E[mcpart_num]/F");
  out_tree->Branch("mcpart_charge",mcpart_charge,"mcpart_charge[mcpart_num]/I");
  out_tree->Branch("mcpart_Vtx_x",mcpart_Vtx_x,"mcpart_Vtx_x[mcpart_num]/F");
  out_tree->Branch("mcpart_Vtx_y",mcpart_Vtx_y,"mcpart_Vtx_y[mcpart_num]/F");
  out_tree->Branch("mcpart_Vtx_z",mcpart_Vtx_z,"mcpart_Vtx_z[mcpart_num]/F");
  out_tree->Branch("mcpart_VtxEnd_x",mcpart_VtxEnd_x,"mcpart_VtxEnd_x[mcpart_num]/F");
  out_tree->Branch("mcpart_VtxEnd_y",mcpart_VtxEnd_y,"mcpart_VtxEnd_y[mcpart_num]/F");
  out_tree->Branch("mcpart_VtxEnd_z",mcpart_VtxEnd_z,"mcpart_VtxEnd_z[mcpart_num]/F");

  if ( kNTrigFlags != ::kNTrigFlags ){
    LOG_ERROR << "StSimpleReaderMaker::Init kNTrigFlags (" << (int)kNTrigFlags
              << ") does not match FcsTriggerDefs.h (" << ::kNTrigFlags
              << "); regenerate and rebuild" << endm;
    return kStFatal;
  }
  for (int i = 0; i < kNTrigFlags; i++){
    out_tree->Branch(kTrigFlagName[i], &Trig_flag[i], TString::Format("%s/I", kTrigFlagName[i]));
  }

  out_tree->Branch("Spin_config", &Spin_config, "Spin_config/I");

  out_tree->Branch("sigma_job", &sigma_job, "sigma_job/D");
  out_tree->Branch("N_job", &N_job, "N_job/I");
  out_tree->Branch("mc_pthatmin_gev", &mc_pthatmin_gev, "mc_pthatmin_gev/F");
  out_tree->Branch("mc_pthatmax_gev", &mc_pthatmax_gev, "mc_pthatmax_gev/F");

  if ( mTrgSim ){
    mTrgBase = mTrgSim->getTriggerEmu();
    for (int i = 0; i < kNSimTrg; i++){
      mSimFlagIndex[i] = -1;
      for (int j = 0; j < kNTrigFlags; j++){
        if ( TString(kSimTrgName[i]) == kTrigFlagName[j] ){ mSimFlagIndex[i] = j; break; }
      }
      if ( mSimFlagIndex[i] < 0 && TString(kSimTrgName[i]) != "DYNoEpd" ){
        LOG_ERROR << "StSimpleReaderMaker::Init sim trigger name '" << kSimTrgName[i]
                  << "' not found in FcsTriggerDefs.h" << endm;
        return kStFatal;
      }
    }
  }

  return kStOK ;

}

Int_t StSimpleReaderMaker::Make( ) 
{ // Do each event

  // Reset variables
  Cal_nhits = 0; Cal_nclus = 0; mcpart_num = 0;
  for (int i = 0; i < kNTrigFlags; i++) Trig_flag[i] = 0;
  Spin_config = -1; // "no data" until proven otherwise below

  // -----------
  // Get 'event' data
  StMuEvent* muEvent = mMuDstMaker->muDst()->event() ;

  // -----------
  // Resolve fired trigger(s) to named FCS trigger flags: either the real
  // fired trigger id(s) (real data -- the same numeric id gets
  // reassigned to a different named trigger across the run period, so
  // matching needs (id, run number) together, see FcsTriggerDefs.h /
  // TriggerIDs.txt), or, if SetTriggerSim() was called, the simulated
  // physics-condition decision from fcs_trg_base run on this event's FCS
  // hits (simulated events have no real fired trigger id to look up).
  if ( mTrgSim ){

    int runNumber = muEvent ? muEvent->runNumber() : 0;
    if ( runNumber != mTrgLastRunNumber ){
      mTrgSim->InitRun(runNumber);
      mTrgLastRunNumber = runNumber;
    }

    mTrgBase->start_event();

    const int kFcsNDetLocal = 6; // StEnumerations.h kFcsNDet; loop det=0..kFcsNDet inclusive, matches StFcsTriggerSimMaker
    StMuFcsCollection* simFcsColl = mMuDstMaker->muDst()->muFcsCollection();
    if ( simFcsColl ){
      for (int det = 0; det <= kFcsNDetLocal; det++){
        int ns  = mFcsDb->northSouth(det);
        int ehp = mFcsDb->ecalHcalPres(det);

        int nh = simFcsColl->numberOfHits(det);
        int det_hit_index = simFcsColl->indexOfFirstHit(det);

        for (int i = 0; i < nh; i++){
          StMuFcsHit* hit = simFcsColl->getHit(i + det_hit_index);
          unsigned short ch = hit->channel();
          if ( ehp < 0 || ch >= 32 ) continue;

          // Simulated-event ADC packing (StFcsTriggerSimMaker::feedADC's
          // mSimMode==1 branch): fast-sim hits don't carry a real
          // per-timebin waveform, just a single energy-equivalent ADC.
          uint16_t data[8];
          memset(data, 0, sizeof(data));
          data[1] = hit->adc(0) - 1;
          data[6] = 1;
          mTrgBase->fill_event(ehp, ns, hit->dep(), ch, data, 8);
        }
      }
    }

    unsigned int dsm_out = (unsigned int)mTrgBase->end_event();
    bool simTrg[kNSimTrg];
    DecodeSimTrg(dsm_out, simTrg);
    for (int i = 0; i < kNSimTrg; i++){
      if ( simTrg[i] && mSimFlagIndex[i] >= 0 ) Trig_flag[mSimFlagIndex[i]] = 1;
    }

  } else if ( muEvent ){
    int runNumber = muEvent->runNumber();
    const StTriggerId& nominalTrig = muEvent->triggerIdCollection().nominal();
    std::vector<unsigned int> firedIds = nominalTrig.triggerIds();

    for (size_t k = 0; k < firedIds.size(); k++){
      unsigned int fid = firedIds[k];
      for (int i = 0; i < kNTrigDefs; i++){
        if ( (unsigned int)kTrigDef[i][1] == fid &&
             runNumber >= kTrigDef[i][2] && runNumber <= kTrigDef[i][3] ){
          Trig_flag[kTrigDef[i][0]] = 1;
        }
      }
    }
  }

  // -----------
  // Resolve this event's real spin configuration (real data only --
  // Spin_config stays -1, set above, if SetSpinDb() was never called).
  // See SpinPlayground/run_number_to_spin_coverage.sh for how narrow the
  // real DB coverage window is (2021-12-14 through 2022-04-18 only, as of
  // this investigation) -- most events outside that window will
  // legitimately end up with Spin_config=-1 even with a valid muEvent and
  // a correctly-attached mSpinDb, because isValid() will be false or the
  // specific bunch crossing won't resolve.
  if ( mSpinDb && muEvent ){
    int runNumber = muEvent->runNumber();
    if ( runNumber != mSpinLastRunNumber ){
      mSpinDb->InitRun(runNumber);
      mSpinLastRunNumber = runNumber;
    }
    if ( mSpinDb->isValid() ){
      unsigned int bx7 = muEvent->l0Trigger().bunchCrossingId7bit(runNumber);
      Spin_config = mSpinDb->spin4usingBX7((int)bx7); // itself -1 if this bx wasn't resolved
    }
  }

  // Resolved once per event (not per hit) for the retroactive HCAL
  // correction lookup below -- real event time straight from the muDst,
  // no DB query needed (unlike Spin_config's run-level DB resolution).
  // StEventInfo::time() returns a raw Unix-epoch int, not a TDatime --
  // TDatime's epoch constructor does the conversion explicitly here.
  if ( mUseHcalRetroCorr ){
    mHcalCorrCurrentDate = muEvent ? TDatime(muEvent->eventInfo().time()).GetDate() : -1;
  }

  // -----------
  // Get FCS data
  StMuFcsCollection* fcs_coll = mMuDstMaker->muDst()->muFcsCollection();  // Array containing the FCS Hits and Clusters

  // Loop over FCS hits
  for(UInt_t ihit = 0 ; ihit < fcs_coll->numberOfHits() ; ihit++){

  	StMuFcsHit* hit = fcs_coll->getHit(ihit); // Pointer to a hit

	int det = hit->detectorId();

        if( det < kFcsNDet ){ // Some entries in MuDST file have detid > 5

    		Cal_detid[Cal_nhits] = det;
   		Cal_hitid[Cal_nhits] = hit->id();
    		Cal_adcsum[Cal_nhits] = hit->adcSum();
    		Cal_hit_energy[Cal_nhits] = hit->energy();

    		// Retroactive HCAL gain correction -- no-op (factor 1.0) unless
    		// SetHcalRetroactiveGainCorr() was called and a manifest period
    		// covers this event's date. See mUseHcalRetroCorr's comment in
    		// StSimpleReaderMaker.h.
    		if ( mUseHcalRetroCorr && (det == kFcsHcalNorthDetId || det == kFcsHcalSouthDetId) ){
    			Cal_hit_energy[Cal_nhits] *= getHcalRetroCorr(det, hit->id());
    		}

		if( det <= kFcsHcalSouthDetId ){
        		StThreeVectorD xyz = mFcsDb->getStarXYZ(det,hit->id());

    			Cal_hit_posx[Cal_nhits] = xyz.x();
    			Cal_hit_posy[Cal_nhits] = xyz.y();
    			Cal_hit_posz[Cal_nhits] = xyz.z();
		}
	      	else if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){ // EPD as Pres.
			// Adapted from code StFcsEventDisplay.cxx
			double zepd=375.0;
			double zfcs=710.0+13.90+15.0;
			double zr=zfcs/zepd;
			int pp,tt,n;
			double x[5],y[5];
			double xsum(0), ysum(0);
	       	 	mFcsDb->getEPDfromId(det,hit->id(),pp,tt);
			mEpdgeo->GetCorners(100*pp+tt,&n,x,y);
 		
			// Get average of corner positions
			// N.B. Number of corners is usually 4, sometimes 5	
			for(int i=0; i<n; i++){
			    xsum += zr*x[i];
			    ysum += zr*y[i];
			}
			Cal_hit_posx[Cal_nhits] = xsum/n;
                	Cal_hit_posy[Cal_nhits] = ysum/n;
                	Cal_hit_posz[Cal_nhits] = zepd;
	    	}
		// Increment number of hits
		Cal_nhits++;
	}
  } // Loop over FCS hits

  // Loop over FCS clusters
  for(UInt_t iclus = 0 ; iclus < fcs_coll->numberOfClusters() ; iclus++){

        StMuFcsCluster* clus = fcs_coll->getCluster(iclus); // Pointer to a cluster

	int det = clus->detectorId();

        if( det < kFcsNDet ){ // Some entries in MuDST file have detid > 5

		Cal_clus_detid[Cal_nclus] = det;
		Cal_clus_ntowers[Cal_nclus] = clus->nTowers();
		Cal_clus_energy[Cal_nclus] = clus->energy();

		// Get cluster position in local coordinates
		Cal_clus_loc_x[Cal_nclus] = clus->x();
		Cal_clus_loc_y[Cal_nclus] = clus->y();

		// Get cluster global position
		StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(det,clus->x(),clus->y());
		Cal_clus_x[Cal_nclus] = xyz.x();
		Cal_clus_y[Cal_nclus] = xyz.y();
		Cal_clus_z[Cal_nclus] = xyz.z();

		// Increment number of clusters
		Cal_nclus++;
	}
  } // Loop over FCS clusters

  // -----------
  // Get Fwd track data
  StMuFwdTrackCollection * ftc = mMuDstMaker->muDst()->muFwdTrackCollection();
  Trk_ntrks = ftc->numberOfFwdTracks();
  
  for ( size_t iTrack = 0; iTrack < ftc->numberOfFwdTracks(); iTrack++ ){

	StMuFwdTrack * muFwdTrack = ftc->getFwdTrack( iTrack );

	Trk_px[iTrack] = muFwdTrack->momentum().Px();
	Trk_py[iTrack] = muFwdTrack->momentum().Py();
	Trk_pz[iTrack] = muFwdTrack->momentum().Pz();
 	Trk_charge[iTrack] = muFwdTrack->charge();
  	Trk_chi2[iTrack] = muFwdTrack->chi2();
  	Trk_ndf[iTrack] = muFwdTrack->ndf();
	Trk_nseedpoints[iTrack] = muFwdTrack->numberOfSeedPoints();
	Trk_nfitpoints[iTrack] = muFwdTrack->numberOfFitPoints();
	Trk_dca_x[iTrack] = muFwdTrack->dca().x();
	Trk_dca_y[iTrack] = muFwdTrack->dca().y();
	Trk_dca_z[iTrack] = muFwdTrack->dca().z();
	Trk_vtxindex[iTrack] = muFwdTrack->vertexIndex();

	// Set track projections to large negative values initially
	// in case track projection fails
	Trk_proj_ecal_x[iTrack] = -9999.;
	Trk_proj_ecal_y[iTrack] = -9999.;
	Trk_proj_ecal_z[iTrack] = -9999.;
	Trk_proj_hcal_x[iTrack] = -9999.;
        Trk_proj_hcal_y[iTrack] = -9999.;
        Trk_proj_hcal_z[iTrack] = -9999.;

	// Get track projections
	for ( auto proj : muFwdTrack->mProjections ) {
		// FCS ECal
		if (proj.mDetId == 41) { // See StEvent/StDetectorDefinitions.h
			Trk_proj_ecal_x[iTrack] = proj.mXYZ.x();
			Trk_proj_ecal_y[iTrack] = proj.mXYZ.y();
			Trk_proj_ecal_z[iTrack] = proj.mXYZ.z();
		}
		// FCS HCal
		if (proj.mDetId == 42) { // See StEvent/StDetectorDefinitions.h
                        Trk_proj_hcal_x[iTrack] = proj.mXYZ.x();
                        Trk_proj_hcal_y[iTrack] = proj.mXYZ.y();
                        Trk_proj_hcal_z[iTrack] = proj.mXYZ.z();
                }
	} // Loop over track projections

  } // Loop over Fwd tracks

  // -----------
  // Retrieve pointer to MC tracks
  TClonesArray *mcTracks = mMuDstMaker->muDst()->mcArray(1);
  
  // Loop over MC tracks
  for (Int_t iTrk=0; iTrk<mcTracks->GetEntriesFast(); iTrk++) {

	// Retrieve i-th MC tracks from MuDst
  	StMuMcTrack *mcTrack = (StMuMcTrack*)mcTracks->UncheckedAt(iTrk);
	if ( !mcTrack ) continue;

	//mcpart_index[mcpart_num] = mcTrack->Id();       // Counts particles sequentially
	mcpart_geid[mcpart_num] = mcTrack->GePid();
  	mcpart_idVtx[mcpart_num] = mcTrack->IdVx();       // ID of creation vertex
        mcpart_idVtxEnd[mcpart_num] = mcTrack->IdVxEnd(); // ID of stop (end) vertex
  	mcpart_px[mcpart_num] = mcTrack->Pxyz().x();
  	mcpart_py[mcpart_num] = mcTrack->Pxyz().y();
  	mcpart_pz[mcpart_num] = mcTrack->Pxyz().z();
  	mcpart_E[mcpart_num] = mcTrack->E();
  	mcpart_charge[mcpart_num] = mcTrack->Charge();

	// Find associated creation and stop vertex locations
	// Retrieve pointer to MC vertices
	TClonesArray *mcVertices = mMuDstMaker->muDst()->mcArray(0);
	
	// Loop over MC vertices
	for (Int_t iVtx=0; iVtx<mcVertices->GetEntriesFast(); iVtx++) {

    		// Retrieve i-th MC vertex from MuDst
    		StMuMcVertex *mcVertex = (StMuMcVertex*)mcVertices->UncheckedAt(iVtx);
	
		// Creation vertex
		if ( mcVertex->Id()==mcTrack->IdVx() ) {
			mcpart_Vtx_x[mcpart_num] = mcVertex->XyzV().x();	
			mcpart_Vtx_y[mcpart_num] = mcVertex->XyzV().y();
			mcpart_Vtx_z[mcpart_num] = mcVertex->XyzV().z();
		}
	
		// Stop vertex
		if ( mcTrack->IdVxEnd()==0 ){ // MC track does not have stop vertex
			mcpart_VtxEnd_x[mcpart_num] = -9999.;
                        mcpart_VtxEnd_y[mcpart_num] = -9999.;
                        mcpart_VtxEnd_z[mcpart_num] = -9999.;		
		}
		else if ( mcVertex->Id()==mcTrack->IdVxEnd() ) {
			mcpart_VtxEnd_x[mcpart_num] = mcVertex->XyzV().x();
                        mcpart_VtxEnd_y[mcpart_num] = mcVertex->XyzV().y();
                        mcpart_VtxEnd_z[mcpart_num] = mcVertex->XyzV().z();
		}

	} // Loop over MC vertices

	// Increment number of MC tracks
	mcpart_num++;

  } // Loop over MC tracks

  mEventsProcessed++ ;
  
  // Fill TTree
  out_tree->Fill();
  return kStOK ;
  
}

Int_t StSimpleReaderMaker::Finish( )
{ // Do once at the end the analysis

  out_file->Write();
  out_file->Close();

  cout << "Total Events Processed in DstMaker " << mEventsProcessed << endl ;

  return kStOk ;

}

// Reads manifest.txt (startDate endDate filename, one line per period,
// '#' comments/blank lines skipped) from mHcalCorrDir. See
// data_to_jet/hcal_gain_corrections/README.md for the format.
void StSimpleReaderMaker::loadHcalCorrManifest(){
  TString path = mHcalCorrDir + "/manifest.txt";
  std::ifstream fin(path.Data());
  if ( !fin ){
    LOG_ERROR << "StSimpleReaderMaker::loadHcalCorrManifest: could not open " << path << endm;
    return;
  }
  std::string line;
  while ( std::getline(fin, line) ){
    if ( line.empty() || line[0] == '#' ) continue;
    std::istringstream iss(line);
    HcalRetroCorrPeriod p;
    std::string fn;
    if ( !(iss >> p.startDate >> p.endDate >> fn) ) continue;
    p.filename = fn;
    mHcalCorrPeriods.push_back(p);
  }
  LOG_INFO << "StSimpleReaderMaker::loadHcalCorrManifest: loaded " << mHcalCorrPeriods.size()
           << " period(s) from " << path << endm;
}

// Reads one correction-factor file (ehp ns dep ch factor per line, same
// format StFcsDb::readGainCorrFromText() uses) into mHcalCorrFactor.
// Rows with ehp!=1 (not HCAL) are skipped -- this machinery is HCAL-only,
// see mUseHcalRetroCorr's comment in the header.
void StSimpleReaderMaker::loadHcalCorrFile(const TString& filename){
  for (int ns = 0; ns < 2; ns++)
    for (int dep = 0; dep < 24; dep++)
      for (int ch = 0; ch < 32; ch++)
        mHcalCorrFactor[ns][dep][ch] = 1.0;

  TString path = mHcalCorrDir + "/" + filename;
  std::ifstream fin(path.Data());
  if ( !fin ){
    LOG_ERROR << "StSimpleReaderMaker::loadHcalCorrFile: could not open " << path << endm;
    return;
  }
  std::string line;
  int ehp, ns, dep, ch;
  float factor;
  while ( std::getline(fin, line) ){
    if ( line.empty() || line[0] == '#' ) continue;
    std::istringstream iss(line);
    if ( !(iss >> ehp >> ns >> dep >> ch >> factor) ) continue;
    if ( ehp != 1 ) continue;
    if ( ns < 0 || ns >= 2 || dep < 0 || dep >= 24 || ch < 0 || ch >= 32 ) continue;
    mHcalCorrFactor[ns][dep][ch] = factor;
  }
  mHcalCorrLoadedFile = filename;
  LOG_INFO << "StSimpleReaderMaker::loadHcalCorrFile: loaded " << path << endm;
}

// Retroactive HCAL correction factor for one hit, resolved from
// mHcalCorrCurrentDate (set once per event in Make()). Returns 1.0
// (no-op) if no manifest period covers this event's date -- not an
// error, matches the "no data -> no-op / sentinel" convention used
// elsewhere in this project (e.g. Spin_config's -1).
float StSimpleReaderMaker::getHcalRetroCorr(int det, int id){
  TString file;
  for (size_t i = 0; i < mHcalCorrPeriods.size(); i++){
    if ( mHcalCorrCurrentDate >= mHcalCorrPeriods[i].startDate &&
         mHcalCorrCurrentDate <  mHcalCorrPeriods[i].endDate ){
      file = mHcalCorrPeriods[i].filename;
      break;
    }
  }
  if ( file.IsNull() ) return 1.0;
  if ( file != mHcalCorrLoadedFile ) loadHcalCorrFile(file);

  int ehp, ns, crt, slt, dep, ch;
  mFcsDb->getDepfromId(det, id, ehp, ns, crt, slt, dep, ch);
  if ( ns < 0 || ns >= 2 || dep < 0 || dep >= 24 || ch < 0 || ch >= 32 ) return 1.0;
  return mHcalCorrFactor[ns][dep][ch];
}
