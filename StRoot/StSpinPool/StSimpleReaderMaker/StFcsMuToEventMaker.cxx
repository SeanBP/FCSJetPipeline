#include "StFcsMuToEventMaker.h"

#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFcsCollection.h"
#include "StMuDSTMaker/COMMON/StMuFcsUtil.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCollection.h"

ClassImp(StFcsMuToEventMaker)

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
  event->setFcsCollection(util.getFcs(muDst->muFcsCollection()));
  AddData(event);
  return kStOK;
}
