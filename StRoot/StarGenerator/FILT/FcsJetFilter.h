/*! \class FcsJetFilter

  This class is just used to filter non Jet events
  The important bits are:
    *change it so that it inherits from StarFilterMaker, not StMCFilter
    *replace RejectGE() with Filter()
    *fix the filter function so that it takes a StarGenEvent rather than a StarGenParticleMaster

  Local copy of the central STAR FcsJetFilter with the flux threshold
  (ETHR, GeV, default 50) settable per submission (submit_sim_to_jet.sh -e).
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

 private:
  Float_t mEthr; //! flux threshold (GeV); default 50.0, see header comment

  ClassDef(FcsJetFilter,1);
};

#endif
