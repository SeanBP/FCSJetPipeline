void r106(){ TChain c("data"); c.Add("out/vtx_*.root"); c.Add("out2/vtx_*.root");
 int rs[2]={23106010,23106038};
 for (int i=0;i<2;i++){ int r=rs[i]; double n=c.GetEntries(Form("event_run==%d",r)); printf("%d n %.0f BBC %.2f BBC150 %.2f EPD %.2f\n", r, n, c.GetEntries(Form("event_run==%d&&bbc_vz_tac>-9000",r))/n, c.GetEntries(Form("event_run==%d&&abs(bbc_vz_tac)<150",r))/n, c.GetEntries(Form("event_run==%d&&epd_vz>-9000",r))/n); } }
