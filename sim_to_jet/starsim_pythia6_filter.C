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

class StarPythia6;
StarPythia6 *_pythia6 = 0;

// "FMS tune": the generator settings of the STAR FMS EM-jet analysis
// (BrightSTAR emjet-gpc, starSim/pythia.kumac -- see
// UnderlyingEventPlayground/fms_emjet_reference.txt), bundled so they
// switch on/off together via starsim()'s fmsTune argument:
//   - PyTune 370 (Perugia 2012; its PDF is CTEQ6L1 through LHAPDF, set 10042)
//   - PARP(90) = 0.213 (energy scaling of the MPI pT cutoff, i.e. the amount
//     of underlying event), overriding the tune's own value
//   - PARP(91) = 1 (primordial kT width)
//   - MSEL = 1 (hard QCD)
// Applied AFTER PyTune, which would otherwise overwrite the two PARPs.
// sqrt(s), the ptHat cuts and particle decays are deliberately NOT part of
// the bundle: they stay as this pipeline already sets them.
const Int_t    kFmsTunePyTune = 370;
const Double_t kFmsTuneParp90 = 0.213;
const Double_t kFmsTuneParp91 = 1.0;
const Int_t    kFmsTuneMsel   = 1;

// bschmookler's dir has the MRST2007lomod grid but no pdfsets.index (and
// isn't writable by us to add one); the system LHAPDF share dir has the
// index but not this grid. Merged copy of both, writable, lives here
// (plus the system dir's cteq6l1 grid, which PyTune 370 needs):
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
void Pythia6( TString mode="pp:minbias", Double_t ckin3=0.0, Int_t tune=325, Double_t ckin4=-1.0, Int_t fmsTune=0 )
{
  gSystem->Setenv("LHAPDF_DATA_PATH", LHAPDF_DATA_PATH.Data() );
  gSystem->Load( "/opt/star/$STAR_HOST_SYS/lib/libLHAPDF.so");
  gSystem->Load( "libPythia6_4_28.so");

  StarPythia6 *pythia6 = new StarPythia6("pythia6");
  _pythia6 = pythia6;
  if ( mode == "pp:minbias" )
  {
    // Real Run 22 energy: sqrt(s)~508.4 GeV (254.2 GeV/nucleon), not 500.
    pythia6->SetFrame("CMS", 508.4 );
    pythia6->SetBlue("proton");
    pythia6->SetYell("proton");

    // StarPythia6::Init() redundantly re-applies these from attributes;
    // if unset, that crashes (SetFrame(_,0)) or corrupts the beam species
    // (SetBlue("")), so set all four explicitly to make it a no-op.
    pythia6->SetAttr("FRAME", "CMS");
    pythia6->SetAttr("Ecms", 508.4);
    pythia6->SetAttr("BLUE", "proton");
    pythia6->SetAttr("YELL", "proton");

    if ( fmsTune ) tune = kFmsTunePyTune;
    if ( tune ) pythia6->PyTune( tune );

    if ( fmsTune ) {
      PyPars_t &pypars = pythia6->pypars();
      pypars.parp(90) = kFmsTuneParp90;
      pypars.parp(91) = kFmsTuneParp91;
      pythia6->pysubs().msel = kFmsTuneMsel;
    }

    // pT-hat minimum/maximum, same role as Pythia8's
    // PhaseSpace:ptHatMin/ptHatMax. PYTHIA6's own default is ckin(4)=-1,
    // which it treats as "no upper cut" (CKIN(4) < CKIN(3)) -- matches our
    // -1 sentinel convention exactly, so leaving ckin4 at its default here
    // reproduces the previous (never-set) behavior.
    PySubs_t &pysubs = pythia6->pysubs();
    pysubs.ckin(3) = ckin3;
    pysubs.ckin(4) = ckin4;
  }

  _primary->AddGenerator(pythia6);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=10, Int_t rngSeed=1234, Double_t ckin3=10.0, Int_t tune=325, Double_t fcsFilterEthr=50.0, Double_t ckin4=-1.0, Int_t fmsTune=0 )
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

  StarRandom::seed( rngSeed );
  StarRandom::capture();

  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "pythia6.starsim.root");
    chain -> AddBefore( "geant", _primary );
  }

  Pythia6( "pp:minbias", ckin3, tune, ckin4, fmsTune );
  command("call gstar_part");

  geometry("y2023");

  // FCS forward energy flux filter (generator-agnostic, same as pythia8)
  filter = new FcsJetFilter();
  ((FcsJetFilter*)filter) -> SetEnergyThreshold( fcsFilterEthr );
  _primary -> AddFilter( filter );

  // Particle acceptance cuts passed to geant (ptmin, ptmax / etamin, etamax / phimin, phimax)
  _primary->SetPtRange  (0.0,  -1.0);         // GeV
  _primary->SetEtaRange ( +1.0, +5.0 );
  _primary->SetPhiRange ( 0., TMath::TwoPi() );

  // z-vertex: gauss widths x=1mm, y=1mm, z=30cm
  _primary->SetVertex( 0., 0., 0. );
  _primary->SetSigma( 0.1, 0.1, 30.0 );

  _primary -> Init();

  command("gkine -4 0");
  command("gfile o pythia6.starsim.fzd");

  trig( nevents );

  // What PYTHIA actually ran with, read back after the event loop so it is
  // obvious from the log whether the FMS tune (or anything else) took effect
  // -- nothing between PyTune and here may have silently reset it.
  cout << "Tune check (fmsTune=" << fmsTune << "): MSEL=" << _pythia6->pysubs().msel
       << " PARP(90)=" << _pythia6->pypars().parp(90)
       << " PARP(91)=" << _pythia6->pypars().parp(91)
       << " PARP(82)=" << _pythia6->pypars().parp(82)
       << " MSTP(5)=" << _pythia6->pypars().mstp(5)
       << " MSTP(51)=" << _pythia6->pypars().mstp(51)
       << " MSTP(52)=" << _pythia6->pypars().mstp(52) << endl;

  chain->Finish();

  command("call agexit");  // ensure STARSIM exits properly
}
// ----------------------------------------------------------------------------
