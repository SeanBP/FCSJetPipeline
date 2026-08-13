// Class TrgSimHelperMaker
//
// Playground-only maker: runs the fcs_trg_base trigger-simulation engine
// on real MuDst FCS hits and cross-checks its decision against the real
// fired trigger id(s) already in the muDst (StMuEvent's
// triggerIdCollection + TriggerIDs.txt).
//
// This is cons-built (not ACLiC), and does its ENTIRE per-event loop
// inside Make() (called by chain->EventLoop(), no CINT round-trip per
// event) -- an earlier version drove the per-event loop from an
// interpreted macro calling out to compiled helper functions, which
// reproducibly corrupted CINT's interpreter state (crash on the *next*
// interpreted statement after a successful compiled call, not inside
// it) regardless of whether the helper was ACLiC- or cons-built. Since
// StSimpleReaderMaker (chain->EventLoop() calling straight into a
// compiled Make(), no CINT round-trip mid-event) has been rock solid in
// this same environment, that's the pattern this follows too.
//
// Bypasses StFcsTriggerSimMaker::Make(), which segfaults in its
// MuDst-fallback path when run outside a full StEvent chain;
// StFcsTriggerSimMaker is still used (via setTrgSim()) for its working
// Init()/InitRun() setup of gains/pedestals/stage_version -- this maker
// just drives its fcs_trg_base engine directly instead of calling its
// Make().

#ifndef TrgSimHelperMaker_def
#define TrgSimHelperMaker_def

#include "StMaker.h"
#include <map>
#include <vector>
#include <string>

class StMuDstMaker;
class StFcsDb;
class StFcsTriggerSimMaker;

class TrgSimHelperMaker : public StMaker
{
 public:

  TrgSimHelperMaker(const char* name="TrgSimHelperMaker") : StMaker(name) {}
  virtual ~TrgSimHelperMaker() {}

  Int_t Init();
  Int_t Make();
  Int_t Finish();

  void setMuDstMaker(StMuDstMaker* m) { mMuDstMaker = m; }
  void setTrgSim(StFcsTriggerSimMaker* t) { mTrgSim = t; }
  void setTriggerIdFile(const char* f) { mTriggerIdFile = f; }
  void setMaxPrintEvents(int n) { mMaxPrintEvents = n; }

 private:

  StMuDstMaker* mMuDstMaker = 0;
  StFcsDb* mFcsDb = 0;
  StFcsTriggerSimMaker* mTrgSim = 0;
  void* mTrgBase = 0; // fcs_trg_base* -- kept opaque (void*) in the header; only the .cxx needs the real type

  const char* mTriggerIdFile = 0;
  int mMaxPrintEvents = 20;
  int mNPrinted = 0;
  int mLastRunNumber = -1;

  std::vector<std::string> mDefName;
  std::vector<int> mDefId;
  std::vector<int> mDefRunStart;
  std::vector<int> mDefRunEnd;

  void loadTriggerIds();
  void lookupTriggerNames(int id, int runNumber, std::vector<std::string>& names);

  std::map<std::string, long long> mReal1Sim1, mReal1Sim0, mReal0Sim1, mReal0Sim0;

  ClassDef(TrgSimHelperMaker,1)

};

#endif
