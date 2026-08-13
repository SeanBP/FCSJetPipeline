// class StSimpleTreeMaker
//

#ifndef STSIMPLETREEMAKER_HH
#define STSIMPLETREEMAKER_HH

#include "StEnumerations.h"
#include "StMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

#include "stdio.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"

//class TApplication;
class StFcsDbMaker;
class StFcsCollection;
class StFcsDb;
class StEpdGeom;

class StSimpleTreeMaker : public StMaker {
  public:
   StSimpleTreeMaker(const Char_t* name = "StSimpleTreeMaker");
   ~StSimpleTreeMaker();
   Int_t Init();
   Int_t Make();
   Int_t Finish();

   void SetOutputFileName(TString name) {mOutputFileName = name;}
  
  private:
   StFcsDb* mFcsDb = 0;
   //StFcsDbMaker* mFcsDbMaker=0;
   StFcsCollection* mFcsColl = 0;
   StEpdGeom* mEpdgeo=0;

   //Output file and tree
   TFile* out_file;
   TTree* out_tree;
   TString mOutputFileName;

   UInt_t mEventsProcessed;

   //TTree Branch variables
   int Cal_nhits;
   
   int Cal_detid[5000];
   int Cal_hitid[5000];
   int Cal_adcsum[5000];
   float Cal_hit_energy[5000];
   float Cal_hit_posx[5000];
   float Cal_hit_posy[5000];
   float Cal_hit_posz[5000];

#ifndef SKIPDefImp
   ClassDef(StSimpleTreeMaker, 0)
#endif
};

#endif
