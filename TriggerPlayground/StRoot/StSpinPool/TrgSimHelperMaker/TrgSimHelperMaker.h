// Class TrgSimHelperMaker
//
// Playground maker: runs the fcs_trg_base trigger simulation on real MuDst
// FCS hits and compares its decision with the fired trigger ids
// (TriggerIDs.txt). The whole per-event loop is in the compiled Make():
// driving it from an interpreted macro corrupts CINT's state. Drives
// fcs_trg_base directly; StFcsTriggerSimMaker (setTrgSim()) is used only for
// its Init()/InitRun() setup, since its Make() needs a full StEvent chain.

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
