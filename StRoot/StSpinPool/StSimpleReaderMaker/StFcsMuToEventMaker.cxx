#include "StFcsMuToEventMaker.h"

#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFcsCollection.h"
#include "StMuDSTMaker/COMMON/StMuFcsUtil.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFcsHit.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StEvent/StEnumerations.h"

#include "FcsChannelMask.h"

ClassImp(StFcsMuToEventMaker)

StFcsMuToEventMaker::StFcsMuToEventMaker(const char* name) : StMaker(name), mMask(new FcsChannelMask()) {}
StFcsMuToEventMaker::~StFcsMuToEventMaker() { delete mMask; }

int StFcsMuToEventMaker::SetChannelMask(const char* file)
{
  const int n = mMask->Load(file);
  LOG_INFO << "StFcsMuToEventMaker: channel mask " << file << ": " << n << " rows" << endm;
  return n;
}

Int_t StFcsMuToEventMaker::Make()
{
  StMuDst* muDst = (StMuDst*)GetInputDS("MuDst");
  if ( !muDst || !muDst->muFcsCollection() ) {
    LOG_WARN << "StFcsMuToEventMaker: no MuDst FCS collection; no StEvent made" << endm;
    return kStOK;
  }
  // Hits only (StMuFcsUtil::getFcs copies no clusters or points); the chain
  // deletes the StEvent at the end of the event
  StMuFcsUtil util;
  StEvent* event = new StEvent();
  StFcsCollection* fcs = util.getFcs(muDst->muFcsCollection());
  // hot / dead channel mask: zero the masked hits, so no cluster or point uses them
  if ( fcs && !mMask->Empty() && muDst->event() ) {
    const int run = muDst->event()->runNumber();
    for ( int det = 0; det <= kFcsHcalSouthDetId; det++ ) {
      StSPtrVecFcsHit& hits = fcs->hits(det);
      for ( unsigned int i = 0; i < hits.size(); i++ )
        if ( mMask->Masked(run, det, hits[i]->id()) ) { hits[i]->setEnergy(0); hits[i]->setAdcSum(0); }
    }
  }
  event->setFcsCollection(fcs);
  AddData(event);
  return kStOK;
}
