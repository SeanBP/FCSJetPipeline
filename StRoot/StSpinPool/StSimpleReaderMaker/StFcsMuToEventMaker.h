// Puts an StEvent holding only the FCS hits of the MuDst into the chain, so
// StFcsClusterMaker and StFcsPointMaker can recluster them (with the
// pipeline's tower thresholds) the same way in data and simulation.
// StSimpleReaderMaker::SetFcsFromStEvent() then writes the new clusters and
// points. Lighter than StMuDst2StEventMaker, which converts the whole event.

#ifndef StFcsMuToEventMaker_def
#define StFcsMuToEventMaker_def

#include "StMaker.h"

class FcsChannelMask;

class StFcsMuToEventMaker : public StMaker
{
 public:
  StFcsMuToEventMaker(const char* name = "FcsMuToEvent");
  virtual ~StFcsMuToEventMaker();
  Int_t Make();
  // Hot / dead channel mask (FcsChannelMask.h): masked hits get zero energy
  // and ADC sum before reclustering. Returns the rows read (-1: no file).
  int SetChannelMask(const char* file);

 private:
  FcsChannelMask* mMask; //! (FcsChannelMask.h, kept out of the dictionary)

 public:

  ClassDef(StFcsMuToEventMaker, 1)
};

#endif
