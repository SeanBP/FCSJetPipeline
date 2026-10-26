/*! \class FcsJetFilter

  This class is just used to filter non Jet events
  The important bits are:
    *change it so that it inherits from StarFilterMaker, not StMCFilter
    *replace RejectGE() with Filter()
    *fix the filter function so that it takes a StarGenEvent rather than a StarGenParticleMaster

  Local copy of the central STAR FcsJetFilter with the flux threshold
  (ETHR, GeV, default 50) settable per submission (submit_sim_to_jet.sh -e).

  Cone mode (SetConePtThreshold(), submit_sim_to_jet.sh -q): instead of the
  energy summed over an FCS half, accept an event if some R = 0.7 (eta, phi)
  cone around a final-state particle hitting the FCS box holds a summed pT
  of at least the threshold, counting particles of the same half only. The
  half-energy filter biases ECAL-only jets below its threshold (a 2-4 GeV pT
  jet at eta ~3 carries only 20-40 GeV, so accepted events carry extra
  energy in that half); a cone threshold below the jets' own cone pT does
  not (UnderlyingEventPlayground/energy_flow/cone_filter_check.cpp).
*/

#ifndef STAR_FcsJetFilter
#define STAR_FcsJetFilter

#include <vector>
#include <string>
#include "StarGenerator/FILT/StarFilterMaker.h"
#include "StarGenerator/EVENT/StarGenEvent.h"

class StarGenParticleMaster;
class StarGenParticle;
class StarGenEvent;

class FcsJetFilter : public StarFilterMaker
{
 public:
  FcsJetFilter(); ///constructor
  virtual ~FcsJetFilter(){;};///destructor
  Int_t Filter( StarGenEvent *mEvent );

  // Logs here, not the constructor: the threshold is always set right
  // after construction (see starsim_pythia8/6_filter.C), and a
  // constructor-time print would show the stale default instead of
  // whatever value actually ends up governing Filter().
  void SetEnergyThreshold(Float_t ethr);
  Float_t GetEnergyThreshold() const { return mEthr; }

  // Cone mode: > 0 replaces the half-energy test (pT in GeV, cone radius R)
  void SetConePtThreshold(Float_t ptmin, Float_t radius = 0.7);
  Float_t GetConePtThreshold() const { return mConePt; }

 private:
  Float_t mEthr; //! flux threshold (GeV); default 50.0, see header comment
  Float_t mConePt; //! cone-mode pT threshold (GeV); 0 = off (half-energy mode)
  Float_t mConeR;  //! cone radius in (eta, phi)

  ClassDef(FcsJetFilter,1);
};

#endif
