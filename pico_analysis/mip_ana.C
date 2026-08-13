#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TText.h"
#include "TStyle.h"
#include "TLegend.h"

#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoFcsHit.h"
#include "StPicoEvent/StPicoFcsCluster.h"
#include "StPicoEvent/StPicoFwdTrack.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

static const int mDebug=0;
int mMaxEvent=400000;

// Track cuts
int mChi2CutMin = 0;
int mChi2CutMax = 65000;
double mNHitCut = 4;
int mVtxIndex = 0;

double mEcalTowCut = 3; //Cut on track projection to Ecal (cm)
double mEcalFrac = 0.1; //Energy fraction of adjacent towers in ECal
int mEcalRows = 22;     //Number of rows in ECal

double mHcalTowCut = 5; //Cut on track projection to Hcal (cm)
double mHcalFrac = 0.1; //Energy fraction of adjacent towers in HCal
int mHcalRows = 13;     //Number of rows in HCal

// Get run number
int getrun(const char* fname){
  TString f(fname);
  int idxst = f.Index("st_");
  int idxbr = f.Index("_",idxst+3);
  int idxar = f.Index("_",idxbr+1);
  TString r(f(idxbr+1,idxar-idxbr));
  int run=r.Atoi();
  std::cout << Form("RunNumber=%d found from %s (st=%d b=%d e=%d l=%d)",
		    run,fname,idxst,idxbr,idxar,idxar-idxbr) << std::endl;
  if(run==0){
    run=23055058;
    std::cout << Form("Overwrite RunNumber = %d",run) << std::endl; 
  }
  return run;
}

// MIP analysis function
void mip_ana(const Char_t *inFile = "infiles.lis") {
  int run=getrun(inFile);

  StFcsDbMaker *fcsDbMk = new StFcsDbMaker();
  fcsDbMk->Init();
  StFcsDb* fcsDb = dynamic_cast<StFcsDb*>(fcsDbMk->GetDataSet("fcsDb"));
  fcsDb->setDbAccess(0);
  fcsDb->InitRun(run);
  std::cout << "FcsDb initialized for run " << run << std::endl;
  for(int det=0; det<4; det++){
    StThreeVectorD off=fcsDb->getDetectorOffset(det);
    printf("det=%d offset=%6.3f %6.3f %6.3f\n",det,off.x(),off.y(),off.z());
  }

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);  
  picoReader->SetStatus("FcsHits",1);
  picoReader->SetStatus("FwdTracks",1);
  std::cout << "Status has been set" << std::endl;

  if( !picoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "Events in Tree:  "  << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Events in chain: " << events2read << std::endl;
  if(events2read>mMaxEvent){
    events2read=mMaxEvent;
    std::cout << "Limit # of event to read to " << events2read << std::endl;
  }
  
  // Define Histograms

  // Tower energy spectrum -- FCS ECal
  TH1 *h1a = new TH1D("h1a","FCS ECal",200,0,3);
  h1a->GetXaxis()->SetTitle("Tower Energy [GeV]");h1a->GetXaxis()->CenterTitle();
  h1a->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h1a->GetYaxis()->CenterTitle();
  h1a->SetLineColor(kBlue);h1a->SetLineWidth(2);

  // Tower energy spectrum Zoomed in -- FCS ECal
  TH1 *h1b = new TH1D("h1b","FCS ECal",500,0.1,0.2);
  h1b->GetXaxis()->SetTitle("Tower Energy [GeV]");h1b->GetXaxis()->CenterTitle();
  h1b->SetLineColor(kBlue);h1b->SetLineWidth(2);
 
  // Tower energy spectrum with track matching -- FCS ECal
  TH1 *h1c = new TH1D("h1c","FCS ECal",200,0,3);
  h1c->GetXaxis()->SetTitle("Tower Energy [GeV]");h1c->GetXaxis()->CenterTitle();
  h1c->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h1c->GetYaxis()->CenterTitle();
  h1c->SetLineColor(kGreen);h1c->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation -- FCS ECal
  TH1 *h1d = new TH1D("h1d","FCS ECal",200,0,3);
  h1d->GetXaxis()->SetTitle("Tower Energy [GeV]");h1d->GetXaxis()->CenterTitle();
  h1d->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h1d->GetYaxis()->CenterTitle();
  h1d->SetLineColor(kRed);h1d->SetLineWidth(2);

  // Tower energy spectrum with isolation but no track matching -- FCS ECal
  TH1 *h1e = new TH1D("h1e","FCS ECal",200,0,3);
  h1e->GetXaxis()->SetTitle("Tower Energy [GeV]");h1e->GetXaxis()->CenterTitle();
  h1e->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h1e->GetYaxis()->CenterTitle();
  h1e->SetLineColor(kMagenta);h1e->SetLineWidth(2);

  // Tower energy spectrum -- FCS HCal
  TH1 *h2a = new TH1D("h2a","FCS HCal",200,0,3);
  h2a->GetXaxis()->SetTitle("Tower Energy [GeV]");h2a->GetXaxis()->CenterTitle();
  h2a->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h2a->GetYaxis()->CenterTitle();
  h2a->SetLineColor(kBlue);h2a->SetLineWidth(2);

  // Tower energy spectrum Zoomed in -- FCS HCal
  TH1 *h2b = new TH1D("h2b","FCS HCal",500,0.1,0.2);
  h2b->GetXaxis()->SetTitle("Tower Energy [GeV]");h2b->GetXaxis()->CenterTitle();
  h2b->SetLineColor(kBlue);h2b->SetLineWidth(2);

  // Tower energy spectrum with track matching -- FCS HCal
  TH1 *h2c = new TH1D("h2c","FCS HCal",200,0,3);
  h2c->GetXaxis()->SetTitle("Tower Energy [GeV]");h2c->GetXaxis()->CenterTitle();
  h2c->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h2c->GetYaxis()->CenterTitle();
  h2c->SetLineColor(kGreen);h2c->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation -- FCS HCal
  TH1 *h2d = new TH1D("h2d","FCS HCal",200,0,3);
  h2d->GetXaxis()->SetTitle("Tower Energy [GeV]");h2d->GetXaxis()->CenterTitle();
  h2d->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h2d->GetYaxis()->CenterTitle();
  h2d->SetLineColor(kRed);h2d->SetLineWidth(2); 
 
  // Tower energy spectrum with isolation but no track matching -- FCS HCal
  TH1 *h2e = new TH1D("h2e","FCS HCal",200,0,3);
  h2e->GetXaxis()->SetTitle("Tower Energy [GeV]");h2e->GetXaxis()->CenterTitle();
  h2e->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h2e->GetYaxis()->CenterTitle();
  h2e->SetLineColor(kMagenta);h2e->SetLineWidth(2);

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {
    if(iEvent%1000==0) std::cout << "Working on event #[" << iEvent<< "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..."<< std::endl;
      break;
    }

    // Retrieve picoDst
    StPicoDst *dst = picoReader->picoDst();

    // Retrieve event information
    StPicoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..."<< std::endl;
      break;
    }

    // Vector to hold indices of good tracks
    std::vector<int> good_trk_index;

    // Loop over reconstructed tracks and find set of 'good' tracks
    int nTrks = dst->numberOfFwdTracks();
    for(int iTrk = 0; iTrk < nTrks; iTrk++){
	StPicoFwdTrack* t = dst->fwdTrack(iTrk);

      	if( abs(t->numberOfFitPoints()) >= mNHitCut &&
	    t->chi2() > mChi2CutMin && t->chi2() < mChi2CutMax &&
	    t->vertexIndex() == 0){
		good_trk_index.push_back(iTrk);
	}
    } //track loop

    //loop over FCS hits
    int nHits = dst->numberOfFcsHits();
    for (int ihit = 0; ihit < nHits; ++ihit){

	StPicoFcsHit* h = dst->fcsHit(ihit);

	//FCS Ecal
	if( h->detectorId() == 0 || h->detectorId() == 1 ){
		
           // Get tower information
           float tower_energy = h->energy();
           int tower_det = h->detectorId();
           int tower_row = int(h->id()/mEcalRows);
           int tower_column = h->id()%mEcalRows;
	   StThreeVectorD xyz = fcsDb->getStarXYZ(tower_det,h->id());
           auto tower_x = xyz.x();
	   auto tower_y = xyz.y();

	   // Some FCS tower hit have energy = 0.
	   if( tower_energy < 0.00001 ) continue;

	   h1a->Fill(tower_energy);
           h1b->Fill(tower_energy);

	   // Require hit to be close to a track projection
	   // Full width of ECal towers ~5.57cm
	   bool matched_to_trk = false;
	   for(Size_t iTrk = 0; iTrk < good_trk_index.size(); iTrk++){

	      StPicoFwdTrack* t = dst->fwdTrack( good_trk_index[iTrk] ); 
	      TVector3 pe = t->ecalProjection();
	      double pex=pe.X();
      	      double pey=pe.Y();

	      if( fabs(pex-tower_x) < mEcalTowCut &&
	          fabs(pey-tower_y) < mEcalTowCut ){
			matched_to_trk = true;
			break;
		}
	   }

	   if(matched_to_trk) h1c->Fill(tower_energy);

	   // Require small relative energy deposit in adjacent towers
	   bool keep_tower = true;
	   for (int iAdj = 0; iAdj < nHits; ++iAdj){

	     // Don't use the same tower
	     if(iAdj == ihit) continue;

	     StPicoFcsHit* ha = dst->fcsHit(iAdj);

	     float adj_energy = ha->energy();
             int adj_det = ha->detectorId();
             int adj_row = int(ha->id()/mEcalRows);
             int adj_column = ha->id()%mEcalRows;

	     if( adj_det == tower_det &&
	         abs(adj_row - tower_row) <=1 &&
	         abs(adj_column - tower_column) <=1 &&
                 (adj_energy / tower_energy) > mEcalFrac ){
			keep_tower = false;
			break;
	     }
	   } //Adjacent tower loop

           if(matched_to_trk && keep_tower) h1d->Fill(tower_energy);
	   if(keep_tower) h1e->Fill(tower_energy);

	} //ECal if statement

	//FCS Hcal
	if( h->detectorId() == 2 || h->detectorId() == 3 ){
	   // Get tower information
           float tower_energy = h->energy();
           int tower_det = h->detectorId();
           int tower_row = int(h->id()/mHcalRows);
           int tower_column = h->id()%mHcalRows;
           StThreeVectorD xyz = fcsDb->getStarXYZ(tower_det,h->id());
           auto tower_x = xyz.x();
           auto tower_y = xyz.y();

           // Some FCS tower hit have energy = 0.
           if( tower_energy < 0.00001 ) continue;

	   h2a->Fill(tower_energy);
           h2b->Fill(tower_energy);

	   // Require hit to be close to a track projection
	   // Full width of HCal towers ~9.99cm
	   bool matched_to_trk = false;
           for(Size_t iTrk = 0; iTrk < good_trk_index.size(); iTrk++){

              StPicoFwdTrack* t = dst->fwdTrack( good_trk_index[iTrk] );
              TVector3 ph = t->hcalProjection();
              double phx=ph.X();
              double phy=ph.Y();

              if( fabs(phx-tower_x) < mHcalTowCut &&
                  fabs(phy-tower_y) < mHcalTowCut ){
                        matched_to_trk = true;
                        break;
                }
           }	

	   if(matched_to_trk) h2c->Fill(tower_energy);

	   // Require small relative energy deposit in adjacent towers
           bool keep_tower = true;
           for (int iAdj = 0; iAdj < nHits; ++iAdj){

             // Don't use the same tower
             if(iAdj == ihit) continue;

             StPicoFcsHit* ha = dst->fcsHit(iAdj);

             float adj_energy = ha->energy();
             int adj_det = ha->detectorId();
             int adj_row = int(ha->id()/mHcalRows);
             int adj_column = ha->id()%mHcalRows;

             if( adj_det == tower_det &&
                 abs(adj_row - tower_row) <=1 &&
                 abs(adj_column - tower_column) <=1 &&
                 (adj_energy / tower_energy) > mHcalFrac ){
                        keep_tower = false;
                        break;
             }
           } //Adjacent tower loop

           if(matched_to_trk && keep_tower) h2d->Fill(tower_energy);
           if(keep_tower) h2e->Fill(tower_energy);

	} //HCal if statement
	
    } //FCS hit loop

  } //event loop
  picoReader->Finish();

  // Set Style
  gStyle->SetOptStat(0);

  // Make plots
  TCanvas *c1a = new TCanvas("c1a");
  c1a->SetLogy();
  h1a->Draw();

  TCanvas *c1b = new TCanvas("c1b");
  h1b->Draw(); 
 
  TCanvas *c1c = new TCanvas("c1c");
  c1c->SetLogy();
  h1a->SetMinimum(10);
  h1a->Draw();
  h1c->Draw("same");
  h1d->Draw("same");

  TLegend *leg1 = new TLegend(0.4,0.6,0.75,0.85);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(h1a,"All FCS ECal tower hits");
  leg1->AddEntry(h1c,"Hits matched to good track");
  leg1->AddEntry(h1d,"Hits matched to good track + isolation");
  leg1->Draw();

  TCanvas *c1d = new TCanvas("c1d");
  c1d->SetLogy();
  h1a->Draw();
  h1e->Draw("same");

  TLegend *leg1a = new TLegend(0.4,0.6,0.75,0.85);
  leg1a->SetBorderSize(0);
  leg1a->SetTextSize(0.03);
  leg1a->AddEntry(h1a,"All FCS ECal tower hits");
  leg1a->AddEntry(h1e,"Hits w/ isolation requirement only");
  leg1a->Draw();

  TCanvas *c2a = new TCanvas("c2a");
  c2a->SetLogy();
  h2a->Draw();

  TCanvas *c2b = new TCanvas("c2b");
  h2b->Draw();

  TCanvas *c2c = new TCanvas("c2c");
  c2c->SetLogy();
  h2a->SetMinimum(10);
  h2a->Draw();
  h2c->Draw("same");
  h2d->Draw("same");

  TLegend *leg2 = new TLegend(0.4,0.6,0.75,0.85);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(h2a,"All FCS HCal tower hits");
  leg2->AddEntry(h2c,"Hits matched to good track");
  leg2->AddEntry(h2d,"Hits matched to good track + isolation");
  leg2->Draw();

  TCanvas *c2d = new TCanvas("c2d");
  c2d->SetLogy();
  h2a->Draw();
  h2e->Draw("same");

  TLegend *leg2a = new TLegend(0.4,0.6,0.75,0.85);
  leg2a->SetBorderSize(0);
  leg2a->SetTextSize(0.03);
  leg2a->AddEntry(h1a,"All FCS HCal tower hits");
  leg2a->AddEntry(h1e,"Hits w/ isolation requirement only");
  leg2a->Draw();

  // Print plots to file
  c1a->Print("plots/mip_ana.pdf[");
  c1a->Print("plots/mip_ana.pdf");
  c1b->Print("plots/mip_ana.pdf");
  c1c->Print("plots/mip_ana.pdf");
  c1d->Print("plots/mip_ana.pdf");
  c2a->Print("plots/mip_ana.pdf");
  c2b->Print("plots/mip_ana.pdf");
  c2c->Print("plots/mip_ana.pdf");
  c2d->Print("plots/mip_ana.pdf");
  c2d->Print("plots/mip_ana.pdf]");

}

