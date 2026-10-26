// Puts an StEvent holding only the FCS hits of the MuDst into the chain, so
// StFcsClusterMaker and StFcsPointMaker can recluster them (with the
// pipeline's tower thresholds) the same way in data and simulation.
// StSimpleReaderMaker::SetFcsFromStEvent() then writes the new clusters and
// points. Lighter than StMuDst2StEventMaker, which converts the whole event.

#ifndef StFcsMuToEventMaker_def
#define StFcsMuToEventMaker_def

#include "StMaker.h"

class StFcsMuToEventMaker : public StMaker
{
 public:
  StFcsMuToEventMaker(const char* name = "FcsMuToEvent") : StMaker(name) {}
  virtual ~StFcsMuToEventMaker() {}
  Int_t Make();

  ClassDef(StFcsMuToEventMaker, 1)
};

#endif
