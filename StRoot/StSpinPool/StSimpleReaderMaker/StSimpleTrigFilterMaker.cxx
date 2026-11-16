#include "StSimpleTrigFilterMaker.h"
#include "StSimpleReaderMaker.h"
#include <cstdio>

ClassImp(StSimpleTrigFilterMaker)

StSimpleTrigFilterMaker::StSimpleTrigFilterMaker(const char* name)
  : StMaker(name), mReader(0), mNPass(0), mNSkip(0)
{
  // allowed to end the event early (StMaker::Make, kStSKIP)
  SetAttr(".Privilege", 1);
}

Int_t StSimpleTrigFilterMaker::Make()
{
  if ( !mReader || mReader->PassesTriggerFilterData() ){ mNPass++; return kStOK; }
  mNSkip++;
  return kStSKIP;
}

Int_t StSimpleTrigFilterMaker::Finish()
{
  printf("StSimpleTrigFilterMaker: %lld events passed, %lld skipped before the FCS makers%s\n",
         mNPass, mNSkip, mReader ? "" : " (no reader set: filter off)");
  return kStOK;
}
