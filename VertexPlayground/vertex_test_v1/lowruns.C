void lowruns(){ TChain c("data"); c.Add("out/vtx_*.root"); c.Add("out2/vtx_*.root");
 c.Draw("event_run>>h(110000,23000000,23110000)","","goff"); c.Draw("event_run>>hb(110000,23000000,23110000)","bbc_vz_tac>-9000","goff");
 TH1*h=(TH1*)gROOT->FindObject("h"),*hb=(TH1*)gROOT->FindObject("hb");
 for(int i=1;i<=h->GetNbinsX();i++){double n=h->GetBinContent(i); if(n<1)continue; int r=(int)h->GetBinLowEdge(i);
  double f=hb->GetBinContent(i)/n; double fa=c.GetEntries(Form("event_run==%d&&abs(bbc_vz_tac)<150",r))/n;
  if(f<0.9||fa>0.84||fa<0.74) printf("%d n %.0f BBC %.2f BBC150 %.2f EPD %.2f\n",r,n,f,fa,c.GetEntries(Form("event_run==%d&&epd_vz>-9000",r))/n);} }
