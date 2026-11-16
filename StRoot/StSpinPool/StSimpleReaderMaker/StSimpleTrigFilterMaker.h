// Early trigger filter (2026-11-16, A_N job efficiency): placed in the chain
// right after the DB makers, before StFcsMuToEventMaker / StEpdHitMaker /
// StFcsClusterMaker / StFcsPointMaker, so events that StSimpleReaderMaker's
// trigger filter would drop anyway skip the FCS reclustering and point fits
// (about 60% of the per-event time; ~70% of st_fwd events fail the A_N
// trigger list). The decision is the reader's own
// (StSimpleReaderMaker::PassesTriggerFilterData()), so the SimpleTree is
// unchanged. Data only: with no filter, or in trigger-simulation mode,
// every event passes. Returns kStSKIP with the ".Privilege" attribute, so
// StMaker::Make stops the event there and StChain moves to the next one.
// Usage (runMudst.C): create it after the DB makers, then SetReader() once
// the reader exists.

#ifndef StSimpleTrigFilterMaker_def
#define StSimpleTrigFilterMaker_def

#include "StMaker.h"

class StSimpleReaderMaker;

class StSimpleTrigFilterMaker : public StMaker
{
 public:
  StSimpleTrigFilterMaker(const char* name = "SimpleTrigFilter");
  virtual ~StSimpleTrigFilterMaker() {}
  void SetReader(StSimpleReaderMaker* reader) { mReader = reader; }
  Int_t Make();
  Int_t Finish();

 private:
  StSimpleReaderMaker* mReader;
  Long64_t mNPass;
  Long64_t mNSkip;

  ClassDef(StSimpleTrigFilterMaker, 1)
};

#endif
