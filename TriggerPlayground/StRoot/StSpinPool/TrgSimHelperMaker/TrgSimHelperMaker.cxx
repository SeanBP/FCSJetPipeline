// Class TrgSimHelperMaker

#include "TrgSimHelperMaker.h"

#include <cstring>
#include <stdint.h>
#include <fstream>
#include <cstdio>

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuFcsCollection.h"
#include "StMuDSTMaker/COMMON/StMuFcsHit.h"
#include "StEvent/StTriggerId.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "StSpinPool/StFcsTriggerSimMaker/StFcsTriggerSimMaker.h"
#include "RTS/src/TRG_FCS/fcs_trg_base.h"

ClassImp(TrgSimHelperMaker)

// Same 21-name/bit table StFcsTriggerSimMaker::Make() decodes from
// dsm_out (copied here, not modifying the official maker, since we
// bypass its Make()). Names match "fcs"+ctrg[i] in the 64-flag set used
// by StSimpleReaderMaker/RecoJets for direct cross-check, except
// DYNoEpd which has no counterpart there.
static const int kNSimTrg = 21;
static const char* kSimTrgName[kNSimTrg] = {
    "fcsJP2","fcsJPA1","fcsJPA0","fcsJPBC1","fcsJPBC0","fcsJPDE1","fcsJPDE0",
    "fcsDiJP","fcsDiJPAsy","fcsDY","fcsJPsi","DYNoEpd","fcsDYAsy",
    "fcsHad2","fcsHad1","fcsHad0","fcsEM2","fcsEM1","fcsEM0","fcsELE2","fcsEM3"
};

static void DecodeSimTrg(unsigned int dsm_out, bool simTrg[kNSimTrg]){
    for (int i = 0; i < kNSimTrg; i++) simTrg[i] = false;
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

void TrgSimHelperMaker::loadTriggerIds(){
    std::ifstream in(mTriggerIdFile);
    std::string name;
    int id, runStart, runEnd;
    while ( in >> name >> id >> runStart >> runEnd ){
        mDefName.push_back(name);
        mDefId.push_back(id);
        mDefRunStart.push_back(runStart);
        mDefRunEnd.push_back(runEnd);
    }
    printf("TrgSimHelperMaker: loaded %d trigger definitions from %s\n",
           (int)mDefName.size(), mTriggerIdFile);
}

void TrgSimHelperMaker::lookupTriggerNames(int id, int runNumber, std::vector<std::string>& names){
    names.clear();
    for (size_t i = 0; i < mDefId.size(); i++){
        if ( mDefId[i] == id && runNumber >= mDefRunStart[i] && runNumber <= mDefRunEnd[i] ){
            names.push_back(mDefName[i]);
        }
    }
}

Int_t TrgSimHelperMaker::Init(){
    mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
    if (!mFcsDb){
        LOG_ERROR << "TrgSimHelperMaker::Init Failed to get StFcsDb" << endm;
        return kStFatal;
    }
    if (!mMuDstMaker){
        LOG_ERROR << "TrgSimHelperMaker::Init mMuDstMaker not set (call setMuDstMaker())" << endm;
        return kStFatal;
    }
    if (!mTrgSim){
        LOG_ERROR << "TrgSimHelperMaker::Init mTrgSim not set (call setTrgSim())" << endm;
        return kStFatal;
    }
    if (!mTriggerIdFile){
        LOG_ERROR << "TrgSimHelperMaker::Init mTriggerIdFile not set (call setTriggerIdFile())" << endm;
        return kStFatal;
    }
    loadTriggerIds();
    mTrgBase = (void*)mTrgSim->getTriggerEmu();
    return kStOK;
}

Int_t TrgSimHelperMaker::Make(){
    fcs_trg_base* tb = (fcs_trg_base*)mTrgBase;

    StMuEvent* ev = mMuDstMaker->muDst()->event();
    if (!ev) return kStOK;
    int runNumber = ev->runNumber();

    // trgSim is SetActive(kFALSE) (see header for why -- its own Make()
    // segfaults), which also means chain scheduling never calls its
    // InitRun() with the real run number, so the per-run gain/pedestal
    // corrections it loads there (from the DB, via StFcsDb) never get
    // set. Call it ourselves whenever the run number changes.
    if (runNumber != mLastRunNumber){
        mTrgSim->InitRun(runNumber);
        mLastRunNumber = runNumber;
    }

    // Real fired trigger names for this event.
    std::map<std::string,bool> realFired;
    const StTriggerId& nominal = ev->triggerIdCollection().nominal();
    std::vector<unsigned int> firedIds = nominal.triggerIds();
    std::vector<std::string> names;
    for (size_t k = 0; k < firedIds.size(); k++){
        lookupTriggerNames((int)firedIds[k], runNumber, names);
        for (size_t m = 0; m < names.size(); m++) realFired[names[m]] = true;
    }

    // Simulated decision: drive fcs_trg_base ourselves (mirrors
    // StFcsTriggerSimMaker::Make()'s MuDst-fallback loop and
    // feedADC()'s mSimMode==0 logic, without calling its crashing Make()).
    tb->start_event();

    StMuFcsCollection* muFcsColl = mMuDstMaker->muDst()->muFcsCollection();
    const int kFcsNDetLocal = 6; // StEnumerations.h kFcsNDet; loop det=0..kFcsNDet inclusive, matches StFcsTriggerSimMaker
    const int kTrgTimebin = 50;  // matches StFcsTriggerSimMaker's own mTrgTimebin default

    if (muFcsColl){
        for (int det = 0; det <= kFcsNDetLocal; det++){
            int ns  = mFcsDb->northSouth(det);
            int ehp = mFcsDb->ecalHcalPres(det);

            int nh = muFcsColl->numberOfHits(det);
            int det_hit_index = muFcsColl->indexOfFirstHit(det);

            for (int i = 0; i < nh; i++){
                StMuFcsHit* hit = muFcsColl->getHit(i + det_hit_index);
                unsigned short ch = hit->channel();
                if (ehp < 0 || ch >= 32) continue;

                uint16_t data[8];
                memset(data, 0, sizeof(data));
                unsigned int ntb = hit->nTimeBin();
                for (unsigned int t = 0; t < ntb; t++){
                    int tbin = hit->timebin(t);
                    if (tbin >= kTrgTimebin-3 && tbin <= kTrgTimebin+4){
                        data[tbin-kTrgTimebin+3] = hit->adc(t);
                    }
                }
                tb->fill_event(ehp, ns, hit->dep(), ch, data, 8);
            }
        }
    }

    unsigned int dsm_out = (unsigned int)tb->end_event();
    bool simTrg[kNSimTrg];
    DecodeSimTrg(dsm_out, simTrg);

    if (mNPrinted < mMaxPrintEvents){
        mNPrinted++;
        printf("Event: run=%d dsm_out=0x%03X  real:", runNumber, dsm_out);
        for (std::map<std::string,bool>::iterator it = realFired.begin(); it != realFired.end(); ++it)
            printf(" %s", it->first.c_str());
        printf("  sim:");
        for (int i = 0; i < kNSimTrg; i++) if (simTrg[i]) printf(" %s", kSimTrgName[i]);
        printf("\n");
    }

    for (int i = 0; i < kNSimTrg; i++){
        std::string nm = kSimTrgName[i];
        if (nm == "DYNoEpd") continue; // no real counterpart
        bool r = realFired.count(nm) > 0;
        bool s = simTrg[i];
        if (r && s) mReal1Sim1[nm]++;
        else if (r && !s) mReal1Sim0[nm]++;
        else if (!r && s) mReal0Sim1[nm]++;
        else mReal0Sim0[nm]++;
    }

    return kStOK;
}

Int_t TrgSimHelperMaker::Finish(){
    printf("\n===== Cross-check summary =====\n");
    printf("%-14s %10s %10s %10s %10s   %s\n", "name", "real1sim1", "real1sim0", "real0sim1", "real0sim0", "P(sim=1|real=1)");
    for (int i = 0; i < kNSimTrg; i++){
        std::string nm = kSimTrgName[i];
        if (nm == "DYNoEpd") continue;
        long long a = mReal1Sim1[nm], b = mReal1Sim0[nm], c = mReal0Sim1[nm], d = mReal0Sim0[nm];
        long long realTot = a + b;
        printf("%-14s %10lld %10lld %10lld %10lld   ", nm.c_str(), a, b, c, d);
        if (realTot > 0) printf("%.1f%% (n=%lld)\n", 100.0 * a / realTot, realTot);
        else printf("n/a (real never fired)\n");
    }
    return kStOK;
}
