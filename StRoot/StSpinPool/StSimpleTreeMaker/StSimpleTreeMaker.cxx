//class StSimpleTreeMaker
//

#include "StSimpleTreeMaker.h"

#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFcsHit.h"
#include "StEventTypes.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StSpinPool/StFcsQaMaker/StFcsQaMaker.h"
#include "StSpinPool/StFcsRawDaqReader/StFcsRawDaqReader.h"
#include "StThreeVectorF.hh"
#include "Stypes.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"

#ifndef SKIPDefImp
ClassImp(StSimpleTreeMaker)
#endif

//------------------------
StSimpleTreeMaker::StSimpleTreeMaker(const Char_t* name) : StMaker(name) {

   out_file = NULL;
   out_tree = NULL;
   mOutputFileName = "";
   mEventsProcessed = 0;

}

StSimpleTreeMaker::~StSimpleTreeMaker() {}

//-----------------------
Int_t StSimpleTreeMaker::Init() {

   //FCS DB
   mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
   //mFcsDb->setDbAccess(0);
   if (!mFcsDb) {
      LOG_ERROR << "StSimpleTreeMaker::InitRun Failed to get StFcsDbMaker" << endm;
      return kStFatal;
   }

   //EPD geom.
   mEpdgeo = new StEpdGeom;

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

   return kStOK;
}

//-----------------------
Int_t StSimpleTreeMaker::Finish() {
   
   out_file->Write();
   out_file->Close();
   
   return kStOK;
}

//----------------------
Int_t StSimpleTreeMaker::Make() {

   //Get StEvent and FCS Collection
   StEvent* event = (StEvent*)GetInputDS("StEvent");
   if (!event) {
      LOG_ERROR << "StSimpleTreeMaker::Make did not find StEvent" << endm;
      return kStErr;
   }
   mFcsColl = event->fcsCollection();
   if (!mFcsColl) {
      LOG_ERROR << "StSimpleTreeMaker::Make did not find StEvent->StFcsCollection" << endm;
      return kStErr;
   }

   //Reset variables
   Cal_nhits = 0;

   //---
   //FCS Hits
   //Note kFcsNDet is equal to 6
   //Note kFcsHcalSouthDetId is equal to 3
   //----
   for (int det = 0; det < kFcsNDet; det++) {
      //check_cal seems to equal 0 for ecal but non-zero for hcal
      //int check_cal = mFcsDb->ecalHcalPres(det);

      StSPtrVecFcsHit& hits = mFcsColl->hits(det);
      int nh = mFcsColl->numberOfHits(det);

      for (int ihit = 0; ihit < nh; ihit++) {
            StFcsHit* hit = hits[ihit];
	    //hit->print();

            Cal_detid[Cal_nhits] = hit->detectorId();
            Cal_hitid[Cal_nhits] = hit->id();
            Cal_adcsum[Cal_nhits] = hit->adcSum();
            Cal_hit_energy[Cal_nhits] = hit->energy();

	    if(det<=kFcsHcalSouthDetId){ //ecal & hcal
            	StThreeVectorF xyz = mFcsDb->getStarXYZ(hit);
           	Cal_hit_posx[Cal_nhits] = xyz.x();
            	Cal_hit_posy[Cal_nhits] = xyz.y();
            	Cal_hit_posz[Cal_nhits] = xyz.z();
	    }
	    else if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){//EPD as Pres.
		//Adapted from code StFcsEventDisplay.cxx
		double zepd=375.0;
		double zfcs=710.0+13.90+15.0;
		double zr=zfcs/zepd;
		int pp,tt,n;
		double x[5],y[5];
		double xsum(0), ysum(0);
	        mFcsDb->getEPDfromId(det,hit->id(),pp,tt);
		mEpdgeo->GetCorners(100*pp+tt,&n,x,y);
 		
		//Get average of corner positions
		//N.B. Number of corners is usually 4, sometimes 5
		for(int i=0; i<n; i++){
		    xsum += zr*x[i];
		    ysum += zr*y[i];
		}
		Cal_hit_posx[Cal_nhits] = xsum/n;
                Cal_hit_posy[Cal_nhits] = ysum/n;
                Cal_hit_posz[Cal_nhits] = zepd;
	    }

            //Increment total hits
            Cal_nhits++;

      } //Loop over hits
   } //Loop over detectors
   
   //Fill Tree
   mEventsProcessed++;
   out_tree->Fill();
   return kStOK;
}
