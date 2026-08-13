// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim.C
//
// Herwig6 + FCS generator-level flux filter. Based on runHerwig6.C, with
// the pieces that macro was missing added: declared herwig6 (was an
// undeclared bare global there), the SetPtRange/SetEtaRange/SetPhiRange/
// SetVertex/SetSigma acceptance cuts every other generator macro here has,
// pT-hat threading via HWHARD's PTMIN, and the FCS filter. Herwig6 has no
// Pythia-style numbered tune -- this uses unmodified Herwig 6 defaults,
// only PTMIN is set.

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarFilterMaker;
StarFilterMaker *filter = 0;

class StarHerwig6;
StarHerwig6 *herwig6 = 0;

// ----------------------------------------------------------------------------
void geometry( TString tag, Bool_t agml=true )
{
  TString cmd = "DETP GEOM "; cmd += tag;
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
}
// ----------------------------------------------------------------------------
void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
void trig( Int_t n=1 )
{
  chain->EventLoop(n);
  _primary->event()->Print();
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Herwig6( TString mode="pp", Double_t ckin3=0.0 )
{
  gSystem->Load( "libHerwig6_5_20.so");

  herwig6 = new StarHerwig6("herwig6");
  if ( mode == "pp" )
  {
    Double_t pblue[]={0.,0.,320.0};
    Double_t pyell[]={0.,0.,-320.0};
    herwig6->SetFrame("3MOM", pblue, pyell );
    herwig6->SetBlue("proton");
    herwig6->SetYell("proton");
    herwig6->SetProcess(1000);   // QCD 2->2 hard subprocess

    // pT-hat minimum, same role as Pythia's ptHatMin. Must be set before
    // _primary->Init() -- StarHerwig6::Init() calls HWInit(), which runs
    // the actual Fortran HWIGIN initialization that reads HWHARD.
    herwig6->hwhard().PTMIN = ckin3;
  }

  _primary->AddGenerator(herwig6);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=10, Int_t rngSeed=1234, Double_t ckin3=10.0 )
{

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "y2023 geant gstar usexgeom agml ";
    bfc(0, simple );
  }

  gSystem->Load( "libVMC.so");

  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );

  gSystem->Load( "libMathMore.so"   );
  gSystem->Load( "xgeometry.so"     );

  gSystem->Load("$OPTSTAR/lib/libfastjet.so");
  gSystem->Load( "StarGeneratorFilt.so" );
  gSystem->Load( "FastJetFilter.so" );
  gSystem->Load( "FCSJetFilter.so" );

  // Setup RNG seed and capture ROOT TRandom
  StarRandom::seed(rngSeed);
  StarRandom::capture();

  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "herwig6.starsim.root");
    chain -> AddBefore( "geant", _primary );
  }

  //
  // Setup an event generator
  //
  Herwig6( "pp", ckin3 );

  geometry("y2023");

  //
  // Setup the generator filter (FCS forward energy flux, generator-agnostic
  // -- same filter used for the Pythia6/Pythia8 tunes)
  //
  filter = new FcsJetFilter();
  _primary -> AddFilter( filter );

  //
  // Setup cuts on which particles get passed to geant for simulation.
  //                    ptmin  ptmax
  _primary->SetPtRange  (0.0,  -1.0);         // GeV
  //                    etamin etamax
  _primary->SetEtaRange ( +1.0, +5.0 );
  //                    phimin phimax
  _primary->SetPhiRange ( 0., TMath::TwoPi() );

  //
  // Setup a realistic z-vertex distribution:
  //   x = 0 gauss width = 1mm
  //   y = 0 gauss width = 1mm
  //   z = 0 gauss width = 30cm
  //
  _primary->SetVertex( 0., 0., 0. );
  _primary->SetSigma( 0.1, 0.1, 30.0 );

  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();

  command("gkine -4 0");
  command("gfile o herwig6.starsim.fzd");

  //
  // Trigger on nevents
  //
  trig( nevents );

  //
  // Finish the chain
  //
  chain->Finish();

  command("call agexit");  // Make sure that STARSIM exits properly
}
// ----------------------------------------------------------------------------
