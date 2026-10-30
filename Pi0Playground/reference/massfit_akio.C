// massfit.C -- pi0 invariant-mass fit macro, from Akio Ogawa's FCS pages.
// Stored here for reference (2026-10-07), not part of FCSJetPipeline's own
// calibration code. Fit model: a Gaussian pi0 peak plus a
// A/tau^2 * (x-x0)^p * exp(-(x-x0)/tau) background shape, background
// fitted alone first (via "xexp"), then both together (via "gxexp"),
// x0 fixed at 0, initial sigma >= 0.04. Companion analysis code (Liang's
// pi0 reconstruction) is at
// github.com/xilinliang/FCS-ECal-pi0-reconstruction. Original hardcoded
// input path (liangxl's StFcsPi0Result*_All.root) left as-is below --
// this macro was not run or adapted, just captured verbatim.

const float peakmin=0.1;
const float peakmax=0.2;
const float tailmax=0.4;

double xexp(double* x, double* p, int flag=0) {
  if(x[0]<p[1]) return 0.0;
  return p[0]/p[2]/p[2]*pow(x[0]-p[1],p[3])*exp(-(x[0]-p[1])/p[2]);
}

double gaus(double* x, double* p, int flag=0) {
  return p[0]*exp(-0.5*pow((x[0]-p[1])/p[2],2));
}

double gxexp(double* x, double* p, int flag=0) {
  return xexp(x,p) + gaus(x,&p[4]);
}

void massfit(){
  massfit(1);
  massfit(20);
  massfit(21);
  massfit(22);
  massfit(44);
  massfit(66);
  massfit(88);
  massfit(110);
  massfit(132);
  massfit(154);
  massfit(330);
  massfit(331);
  massfit(350);
  massfit(352);
  massfit(374);
  massfit(343);
  massfit(706);
  massfit(717);
}

void massfit(int id){
  TCanvas *c1 = new TCanvas("MassFit","MassFit",800,800);
  gStyle->SetOptStat(10);
  gStyle->SetTitleH(0.05);
  gStyle->SetStatH(0.1);

  //TFile *F = new TFile("/star/data01/pwg/liangxl/mywork/fcs_OO_Ecal/fcs/rootfile/StFcsPi0Result1002_All.root","old");
  TFile *F = new TFile("/star/data01/pwg/liangxl/mywork/fcs_OO_Ecal/fcs/rootfile/StFcsPi0Result1023_All.root","old");
  TH1F* h = F->Get(Form("mass_by_Ntower_%d",id));
  cout << h->GetTitle() << endl;
  TString t(h->GetTitle());
  t.ReplaceAll("invariant mass sorted by ","Id=");
  t.ReplaceAll(" tower","");
  cout << t.Data() << endl;
  h->SetTitle(t.Data());

  int ent = h->GetEntries();
  double max = h->GetMaximum();
  double peak = h->GetBinCenter(h->GetMaximumBin());
  if(peak<0.05 || peak>0.20) peak=0.135;
  printf("max=%f peak=%f\n",max,peak);

  TF1* ftx = new TF1("xexp",xexp,0,tailmax,4);
  ftx->SetLineColor(kGreen);
  ftx->SetParName(0,"A");     ftx->SetParameter(0,max/2);
  ftx->SetParName(1,"x0");    ftx->FixParameter(1,0.0);
  ftx->SetParName(2,"tau");   ftx->SetParameter(2,0.2);
  ftx->SetParName(3,"p");     ftx->SetParameter(3,1.0);
  TFitResultPtr rest = h->Fit("xexp","b","",0,tailmax);

  TF1* f = new TF1("gxexp",gxexp,0,tailmax,7);
  f->SetLineColor(kRed);
  f->SetParName(0,"A");     f->FixParameter(0,ftx->GetParameter(0));
  f->SetParName(1,"x0");    f->FixParameter(1,ftx->GetParameter(1));
  f->SetParName(2,"tau");   f->FixParameter(2,ftx->GetParameter(2));
  f->SetParName(3,"p");     f->FixParameter(3,ftx->GetParameter(3));
  f->SetParName(4,"peak");  f->FixParameter(4,max/2); f->SetParLimits(4,0.0,max);
  f->SetParName(5,"mean");  f->FixParameter(5,peak);
  f->SetParName(6,"sigma"); f->SetParameter(6,0.04); f->SetParLimits(6,0.01,0.1);
  TFitResultPtr resg = h->Fit("gxexp","","",0,tailmax);
  f->SetParameter(0,ftx->GetParameter(0));
  f->SetParameter(1,ftx->GetParameter(1));
  f->SetParameter(2,ftx->GetParameter(2));
  f->SetParameter(3,ftx->GetParameter(3));
  f->SetParameter(4,max/2); f->SetParLimits(4,0.0,max);
  f->SetParameter(5,peak);
  f->SetParameter(6,f->GetParameter(6)); f->SetParLimits(6,0.005,0.1);
  TFitResultPtr resg = h->Fit("gxexp","","",0,tailmax);

  TF1* bg = new TF1("bg",xexp,0,tailmax,4);
  bg->SetLineColor(kBlue);
  bg->SetParameter(0,f->GetParameter(0));
  bg->SetParameter(1,f->GetParameter(1));
  bg->SetParameter(2,f->GetParameter(2));
  bg->SetParameter(3,f->GetParameter(3));

  TF1* g = new TF1("gaus",gaus,0,tailmax,3);
  g->SetLineColor(kMagenta);
  g->SetParameter(0,f->GetParameter(4));
  g->SetParameter(1,f->GetParameter(5));
  g->SetParameter(2,f->GetParameter(6));

  h->Draw("e");
  ftx->Draw("same");
  f->Draw("same");
  bg->Draw("same");
  g->Draw("same");
  c1->SaveAs(Form("massfit_id%d.png",id));
}
