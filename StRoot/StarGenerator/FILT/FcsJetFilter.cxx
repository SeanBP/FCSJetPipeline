#include "FcsJetFilter.h"

#include "StarGenerator/EVENT/StarGenParticle.h"
#include "StarGenerator/EVENT/StarGenEvent.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

static const float ZFCS    =  710.16 + 13.9 + 15.0;
static const float XFCSMin =  16.69;
static const float XFCSMax =  16.69 + 22.0*5.542;
static const float YFCS    =  34.0/2.0 * 5.542;
static const float FVCUT   =  20.0; //cm

static int NTOT=0;
static int NGOOD=0;

//_______________________________________________________________
FcsJetFilter::FcsJetFilter():StarFilterMaker("fcsJetFilter"), mEthr(50.0), mConePt(0.0), mConeR(0.7)
{
}
//_______________________________________________________________
void FcsJetFilter::SetEnergyThreshold(Float_t ethr)
{
    mEthr = ethr;
    cout<<"FCS JET filter is used!!! Energy threshold = "<<mEthr<<" GeV"<<endl;
}
//_______________________________________________________________
void FcsJetFilter::SetConePtThreshold(Float_t ptmin, Float_t radius)
{
    mConePt = ptmin;
    mConeR = radius;
    if(mConePt > 0.0)
	cout<<"FCS JET filter in cone mode!!! pT threshold = "<<mConePt<<" GeV in R = "<<mConeR<<" cones (energy threshold not used)"<<endl;
}
//_______________________________________________________________
Int_t FcsJetFilter::Filter( StarGenEvent *mEvent){
    NTOT++;
    // Get a reference to the current event
    StarGenEvent& event = *mEvent;

    //event.Print();
    if(event.GetNumberOfParticles() <= 0) {return kError;}

    TIter Iterator = event.IterAll();
    StarGenParticle *p = 0;

    float etot[2]={0.0,0.0};
    std::vector<float> cEta, cPhi, cPt; std::vector<int> cHalf; // cone mode: particles hitting the box
    while( ( p = (StarGenParticle*)Iterator.Next() ) ){
	if(p->GetStatus() != 1)continue;
	if(p->GetPz() < 0.0) continue;
	//simple box cut
	float x = p->GetVx()/10.0 + p->GetPx() / p->GetPz() * (ZFCS - p->GetVz()/10.0);
	float absx=fabs(x);
	if(absx<XFCSMin+FVCUT || absx>XFCSMax-FVCUT) continue;
	float y = fabs (p->GetVy()/10.0 + p->GetPy() / p->GetPz() * (ZFCS - p->GetVz()/10.0));
	if(y>YFCS-FVCUT) continue;
	//it is hitting FCS
	if(x>0.0)       {etot[1]+=p->GetEnergy();} //south
	else if(x<=0.0) {etot[0]+=p->GetEnergy();} //north
	if(mConePt > 0.0){
	    float pt = sqrt(p->GetPx()*p->GetPx() + p->GetPy()*p->GetPy());
	    if(pt <= 0.0) continue;
	    cEta.push_back(asinh(p->GetPz()/pt)); cPhi.push_back(atan2(p->GetPy(), p->GetPx()));
	    cPt.push_back(pt); cHalf.push_back(x>0.0 ? 1 : 0);
	}
    }
    bool accept = false;
    if(mConePt > 0.0){
	// any R cone around a box particle, same half, with summed pT >= mConePt
	for(size_t i=0; i<cPt.size() && !accept; i++){
	    float sum = 0.0;
	    for(size_t k=0; k<cPt.size(); k++){
		if(cHalf[k] != cHalf[i]) continue;
		float deta = cEta[k]-cEta[i];
		float dphi = remainder(cPhi[k]-cPhi[i], 2.0*M_PI);
		if(deta*deta + dphi*dphi < mConeR*mConeR) sum += cPt[k];
	    }
	    if(sum >= mConePt) accept = true;
	}
    }else{
	accept = (etot[0] > mEthr || etot[1] > mEthr);
    }
    if(accept){
	NGOOD++;
	cout << Form("FcsJetFilter : N_Genearted=%6d  N_Accepted=%6d  R=%6.4f",
		     NTOT,NGOOD,float(NGOOD)/float(NTOT) ) << endl;
	return StarGenEvent::kAccept;
    }
    return StarGenEvent::kReject;
}
