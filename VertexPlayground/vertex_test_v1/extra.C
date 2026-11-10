#include "TChain.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include <cstdio>
void extra()
{
    TChain c("data"); c.Add("out/vtx_*.root");
    c.Draw("bbc_vz_det:epd_vz>>a(60,-300,300,60,-300,300)", "epd_vz>-9000&&abs(bbc_vz_det)<9000", "goff");
    TH2 *a = (TH2*)gROOT->FindObject("a"); printf("old zVertex vs EPD: n %.0f corr %.2f, old zVertex mean %.1f rms %.1f\n", a->GetEntries(), a->GetCorrelationFactor(), a->GetMean(2), a->GetRMS(2));
    // BBC TAC vs EPD: profile slope both ways
    c.Draw("bbc_vz_tac:epd_vz>>b(30,-150,150,120,-600,600)", "epd_vz>-9000&&bbc_vz_tac>-9000", "goff");
    TH2 *b = (TH2*)gROOT->FindObject("b");
    TProfile *pb = b->ProfileX(); pb->Fit("pol1", "Q0", "", -100, 100);
    TF1 *f = pb->GetFunction("pol1"); printf("<BBC TAC z> vs EPD z: slope %.2f +- %.2f, offset %.1f\n", f->GetParameter(1), f->GetParError(1), f->GetParameter(0));
    c.Draw("epd_vz:bbc_vz_tac>>e(30,-300,300,60,-300,300)", "epd_vz>-9000&&bbc_vz_tac>-9000", "goff");
    TProfile *pe = ((TH2*)gROOT->FindObject("e"))->ProfileX(); pe->Fit("pol1", "Q0", "", -200, 200);
    f = pe->GetFunction("pol1"); printf("<EPD z> vs BBC TAC z: slope %.3f +- %.3f\n", f->GetParameter(1), f->GetParError(1));
    // core width of BBC-EPD
    c.Draw("bbc_vz_tac-epd_vz>>d(120,-400,400)", "epd_vz>-9000&&bbc_vz_tac>-9000", "goff");
    TH1 *d = (TH1*)gROOT->FindObject("d"); d->Fit("gaus", "Q0", "", -150, 150);
    printf("BBC-EPD all: rms %.1f, gauss core sigma %.1f\n", d->GetRMS(), d->GetFunction("gaus")->GetParameter(2));
    // EPD availability by fcs trigger: EM vs JP via tcu? just by run
    for (const char *r : {"23011061","23039034","23068007","23104042"}) {
        double n = c.GetEntries(Form("event_run==%s", r)), e = c.GetEntries(Form("event_run==%s&&epd_vz>-9000", r)), bb = c.GetEntries(Form("event_run==%s&&abs(bbc_vz_tac)<150", r));
        c.Draw("bbc_vz_tac>>r1(100,-500,500)", Form("event_run==%s&&bbc_vz_tac>-9000", r), "goff");
        c.Draw("epd_vz>>r2(100,-500,500)", Form("event_run==%s&&epd_vz>-9000", r), "goff");
        printf("run %s: n %.0f EPD %.0f%% BBC|z|<150 %.0f%%  BBC rms %.0f EPD rms %.0f\n", r, n, 100*e/n, 100*bb/n, ((TH1*)gROOT->FindObject("r1"))->GetRMS(), ((TH1*)gROOT->FindObject("r2"))->GetRMS());
    }
}
