#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"

#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoFcsHit.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

static const int mDebug = 0;
int mMaxEvent = 400000;

// Function to extract run number from filename
int getrun(const char* fname){
    TString f(fname);
    int idxst = f.Index("st_");
    int idxbr = f.Index("_", idxst+3);
    int idxar = f.Index("_", idxbr+1);
    TString r(f(idxbr+1, idxar-idxbr));
    int run = r.Atoi();
    if(run == 0) run = 23055058;
    return run;
}

void pico_to_root(const Char_t* inFile = "infiles.lis", const Char_t* outFileName = "./output/pico_all.root") {
    int run = getrun(inFile);

    // Create output directory if it doesn't exist
    TString outDir(outFileName);
    int lastSlash = outDir.Last('/');
    if(lastSlash >= 0) {
        TString dir = outDir(0, lastSlash);
        gSystem->mkdir(dir, kTRUE);
    }

    // Initialize FCS database
    StFcsDbMaker *fcsDbMk = new StFcsDbMaker();
    fcsDbMk->Init();
    StFcsDb* fcsDb = dynamic_cast<StFcsDb*>(fcsDbMk->GetDataSet("fcsDb"));
    fcsDb->setDbAccess(0);
    fcsDb->InitRun(run);

    // Initialize PicoDst reader
    StPicoDstReader* picoReader = new StPicoDstReader(inFile);
    picoReader->Init();

    picoReader->SetStatus("*", 0);
    picoReader->SetStatus("Event", 1);  
    picoReader->SetStatus("FcsHits", 1);

    if(!picoReader->chain()) {
        std::cout << "No chain found." << std::endl;
        return;
    }

    Long64_t events2read = picoReader->chain()->GetEntries();
    if(events2read > mMaxEvent) events2read = mMaxEvent;
    std::cout << "Reading " << events2read << " events." << std::endl;

    // Open output ROOT file and TTree
    TFile* outFile = new TFile(outFileName, "RECREATE");
    TTree* outTree = new TTree("events", "Event Tree");

    // Event-level variables
    int evtID;
    float vtx_x, vtx_y, vtx_z;
    int nHits;

    // Hit-level variables (vectors)
    std::vector<float> hit_energy;
    std::vector<float> hit_x;
    std::vector<float> hit_y;
    std::vector<float> hit_z;
    std::vector<int> hit_detID;

    // Set branches
    outTree->Branch("evtID", &evtID);
    outTree->Branch("vtx_x", &vtx_x);
    outTree->Branch("vtx_y", &vtx_y);
    outTree->Branch("vtx_z", &vtx_z);
    outTree->Branch("nHits", &nHits);
    outTree->Branch("hit.energy", &hit_energy);
    outTree->Branch("hit.position.x", &hit_x);
    outTree->Branch("hit.position.y", &hit_y);
    outTree->Branch("hit.position.z", &hit_z);
    outTree->Branch("hit.detID", &hit_detID);

    // Event loop
    for(Long64_t iEvent = 0; iEvent < events2read; iEvent++) {
        if(iEvent % 1000 == 0)
            std::cout << "Processing event " << iEvent << "/" << events2read << std::endl;

        if(!picoReader->readPicoEvent(iEvent)) continue;
        StPicoDst* dst = picoReader->picoDst();
        if(!dst) continue;
        StPicoEvent* event = dst->event();
        if(!event) continue;

        evtID = event->eventId();
        auto vtx = event->primaryVertex();
        vtx_x = vtx.x();
        vtx_y = vtx.y();
        vtx_z = vtx.z();

        nHits = dst->numberOfFcsHits();

        // Clear hit vectors
        hit_energy.clear();
        hit_x.clear();
        hit_y.clear();
        hit_z.clear();
        hit_detID.clear();

        for(int ihit = 0; ihit < nHits; ihit++) {
            StPicoFcsHit* h = dst->fcsHit(ihit);
            hit_energy.push_back(h->energy());
            hit_detID.push_back(h->detectorId());

            // Get hit position from FCS DB
            StThreeVectorD pos = fcsDb->getStarXYZ(h->detectorId(), h->id());
            hit_x.push_back(pos.x());
            hit_y.push_back(pos.y());
            hit_z.push_back(pos.z());
        }

        outTree->Fill();
    }

    // Write and close file
    outTree->Write();
    outFile->Close();

    picoReader->Finish();

    std::cout << "All events and FCS hits saved to " << outFileName << std::endl;
}
