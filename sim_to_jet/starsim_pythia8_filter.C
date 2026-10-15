// starsim with Pythia8 and the FCS generator-level flux filter.
// Run as root4star starsim.C(...)

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarFilterMaker;
StarFilterMaker *filter = 0;

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
// trig(n) generates n+1 events; last event in the FZD file will be corrupt.
void trig( Int_t n=1 )
{
  chain->EventLoop(n);
  _primary->event()->Print();
  
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Pythia8( TString config="pp:W", Double_t ckin3=0.0, Double_t ckin4=-1.0, Int_t pdfSet=-1, Int_t ppTune=-1 )
{

  StarPythia8 *pythia8 = new StarPythia8();
  if ( config=="pp:W" )
    {
      pythia8->SetFrame("CMS", 510.0);
      pythia8->SetBlue("proton");
      pythia8->SetYell("proton");
      
      pythia8->Set("WeakSingleBoson:all=off");
      pythia8->Set("WeakSingleBoson:ffbar2W=on");
      pythia8->Set("24:onMode=0");              // switch off all W+/- decays
      pythia8->Set("24:onIfAny 11 -11");        // switch on for decays to e+/-
      
    }
  if ( config=="pp:minbias" )
    {
      // Run 22 sqrt(s) = 508.4 GeV
      pythia8->SetFrame("CMS", 508.4);
      pythia8->SetBlue("proton");
      pythia8->SetYell("proton");    

      pythia8->Set("HardQCD:all = on");

    }
  if ( config=="pp:heavyflavor:D0jets" ) 
    {
      pythia8->Set("HardQCD:gg2ccbar = on");
      pythia8->Set("HardQCD:qqbar2ccbar = on");
      pythia8->Set("Charmonium:all = on");

      pythia8->Set("421:mayDecay = 0");

      pythia8->SetFrame("CMS", 200.0);
      pythia8->SetBlue("proton");
      pythia8->SetYell("proton");

    }

  // Setup phase space cuts
  pythia8 -> Set(Form("PhaseSpace:ptHatMin=%f", ckin3 ));
  pythia8 -> Set(Form("PhaseSpace:ptHatMax=%f", ckin4 ));

  // Pythia8 Tune:pp preset (e.g. 14=Monash, 33=Detroit, 3=Tune 2C, 5=Tune
  // 4C -- see pythia8.web.cern.ch/manuals/pythia8317/Tunes.html). Applied
  // before the PDF override below so an explicit pdfSet still wins.
  // StarPythia8::Init() re-sets Detroit's five MPI/ColourReconnection values
  // afterwards, so presets other than 33 keep Detroit's for those.
  if ( ppTune > 0 )
    pythia8 -> Set(Form("Tune:pp = %d", ppTune));

  // PDF tune override (see Pythia8 PDF:pSet options). pdfSet<=0 leaves
  // Pythia8's internal default PDF in place.
  if ( pdfSet > 0 )
    pythia8 -> Set(Form("PDF:pSet = %d", pdfSet));

  _primary -> AddGenerator( pythia8 );
  
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=10, Int_t rngSeed=1234, Double_t ckin3=10.0, Int_t pdfSet=-1, Double_t fcsFilterEthr=50.0, Double_t ckin4=-1.0, Int_t ppTune=-1 )
{

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "y2023 geant gstar usexgeom agml ";
    bfc(0, simple );
  }

  gSystem->Load( "libVMC.so");

  gSystem->Load( "StarGeneratorUtil.so");
  gSystem->Load( "StarGeneratorEvent.so");
  gSystem->Load( "StarGeneratorBase.so" );

  gSystem->Load( "Pythia8_3_03.so"  );

  gSystem->Load( "libMathMore.so"   );  

  // Force loading of xgeometry
  gSystem->Load( "xgeometry.so"     );

  gSystem->Load("$OPTSTAR/lib/libfastjet.so");
  gSystem->Load( "StarGeneratorFilt.so" );
  gSystem->Load( "FastJetFilter.so" );
  gSystem->Load( "FCSJetFilter.so" );

  StarRandom::seed( rngSeed );
  StarRandom::capture();

  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "pythia8.starsim.root");
    chain -> AddBefore( "geant", _primary );
  }

  Pythia8("pp:minbias", ckin3, ckin4, pdfSet, ppTune );
  command("call gstar_part");
  geometry("y2023");

#if 1
  filter = new FcsJetFilter();
  ((FcsJetFilter*)filter) -> SetEnergyThreshold( fcsFilterEthr );
  _primary -> AddFilter( filter );

#endif

  // Particle acceptance cuts passed to geant (ptmin, ptmax / etamin, etamax / phimin, phimax)
  _primary->SetPtRange  (0.0,  -1.0);         // GeV
  // Covers the EM-only truth particle band (JetParameters.h): the ECal
  // fiducial eta range +- 1.5R (1.29-4.62), 0.29 to spare below, 0.38 above
  _primary->SetEtaRange ( +1.0, +5.0 );
  _primary->SetPhiRange ( 0., TMath::TwoPi() );

  // z-vertex: gauss widths x=1mm, y=1mm, z=30cm
  _primary->SetVertex( 0., 0., 0. );
  _primary->SetSigma( 0.1, 0.1, 30.0 );

  _primary -> Init();

  command("gkine -4 0");
  command("gfile o pythia8.starsim.fzd");
  

  trig( nevents );
  chain->Finish();
  command("call agexit");  // ensure STARSIM exits properly

}
// ----------------------------------------------------------------------------

