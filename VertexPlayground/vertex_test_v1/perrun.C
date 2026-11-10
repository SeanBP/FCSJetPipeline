#include "TChain.h"
#include <cstdio>
void perrun(){ TChain c("data"); c.Add("out/vtx_*.root");
 const char *rs[]={"23011061","23023035","23031038","23039034","23050004","23058042","23068007","23077001","23083041","23091030","23097060","23104042"};
 for(auto r:rs){double n=c.GetEntries(Form("event_run==%s",r)); printf("%s n %5.0f EPD %3.0f%% EPD|z|<150 %3.0f%% nEPDhits>0(E,W) %3.0f%% %3.0f%%\n",r,n,100*c.GetEntries(Form("event_run==%s&&epd_vz>-9000",r))/n,100*c.GetEntries(Form("event_run==%s&&abs(epd_vz)<150",r))/n,100*c.GetEntries(Form("event_run==%s&&epd_ntac_e>0",r))/n,100*c.GetEntries(Form("event_run==%s&&epd_ntac_w>0",r))/n);} }
