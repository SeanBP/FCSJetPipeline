// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim.C
//
// Pythia6 + FCS generator-level flux filter, mirroring
// starsim_pythia8_filter.C's structure (filter, ptHatMin threading) but
// using Pythia6/PyTune instead of Pythia8's Set() interface.

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarFilterMaker;
StarFilterMaker *filter = 0;

// bschmookler's dir has the MRST2007lomod grid but no pdfsets.index (and
// isn't writable by us to add one); the system LHAPDF share dir has the
// index but not this grid. Merged copy of both, writable, lives here:
TString LHAPDF_DATA_PATH="/star/u/seanp/FCSJetPipeline/sim_to_jet/lhapdf_data";

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
void Pythia6( TString mode="pp:minbias", Double_t ckin3=0.0, Int_t tune=325 )
{
  gSystem->Setenv("LHAPDF_DATA_PATH", LHAPDF_DATA_PATH.Data() );
  gSystem->Load( "/opt/star/$STAR_HOST_SYS/lib/libLHAPDF.so");
  gSystem->Load( "libPythia6_4_28.so");

  StarPythia6 *pythia6 = new StarPythia6("pythia6");
  if ( mode == "pp:minbias" )
  {
    pythia6->SetFrame("CMS", 500.0 );
    pythia6->SetBlue("proton");
    pythia6->SetYell("proton");

    // StarPythia6::Init() redundantly re-applies FRAME/Ecms/BLUE/YELL from
    // attributes ("if (SAttr('FRAME')) SetFrame(...)" etc.) -- if those
    // attributes are unset, SAttr() apparently still evaluates truthy here
    // (cause unclear), so the redundant calls fire with empty/zero values:
    // SetFrame(_, 0) crashes the mDirect=value/mRootS assertion (0/0=NaN),
    // and SetBlue("")/SetYell("") corrupt the beam species, which PyInit
    // then reports as "unrecognized beam particle ''...". Set all four
    // explicitly so the redundant calls are harmless no-ops.
    pythia6->SetAttr("FRAME", "CMS");
    pythia6->SetAttr("Ecms", 500.0);
    pythia6->SetAttr("BLUE", "proton");
    pythia6->SetAttr("YELL", "proton");

    if ( tune ) pythia6->PyTune( tune );

    // Setup phase space cut (pT-hat minimum), same role as Pythia8's
    // PhaseSpace:ptHatMin in starsim_pythia8_filter.C
    PySubs_t &pysubs = pythia6->pysubs();
    pysubs.ckin(3) = ckin3;
  }

  _primary->AddGenerator(pythia6);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=10, Int_t rngSeed=1234, Double_t ckin3=10.0, Int_t tune=325 )
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

  // Setup RNG seed and map all ROOT TRandom here
  StarRandom::seed( rngSeed );
  StarRandom::capture();

  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "pythia6.starsim.root");
    chain -> AddBefore( "geant", _primary );
  }

  //
  // Setup an event generator
  //
  Pythia6( "pp:minbias", ckin3, tune );
  command("call gstar_part");

  geometry("y2023");

  //
  // Setup the generator filter (FCS forward energy flux -- see
  // FcsJetFilter.cxx; operates on the generator-agnostic StarGenEvent, so
  // this is the same filter used for the Pythia8 tunes)
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
  command("gfile o pythia6.starsim.fzd");

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
