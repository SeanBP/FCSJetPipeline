// Trigger turn-on in data vs simulation, one panel per trigger: the fraction
// of jets with the (hardware / emulated) bit set vs the jet's trigger value /
// nominal threshold, from two trg_eff_scan outputs (trg_turnon_periods.root),
// with the cuts 1.15 (dashed) and 1.3 (solid). Data: the gain periods listed
// in data_periods (0 before 23027047, 1 23027047-23066055, 2 from 23066056);
// simulation: its one period sim_period (the period of its calibration run).
// Usage: root -b -q 'plot_turnon_data_sim.C("data/trg_turnon_periods.root", "sim/trg_turnon_periods.root", 2, "12", "turnon_data_sim.png")'

void plot_turnon_data_sim(const char *data_root, const char *sim_root, int sim_period, const char *data_periods, const char *out_png)
{
    gStyle->SetOptStat(0);
    TFile *fd = TFile::Open(data_root), *fs = TFile::Open(sim_root);
    if (!fd || !fs) { printf("cannot open inputs\n"); return; }
    const char *trg[6] = {"EM0", "EM1", "EM2", "JP0", "JP1", "JP2"};
    const double thr[6] = {2.0, 3.5, 5.0, 4.0, 6.0, 8.0};
    const char *plab[3] = {"data before 23027047", "data 23027047-23066055", "data from 23066056"};
    const int pcol[3] = {kAzure + 2, kOrange + 7, kBlack}, pmk[3] = {20, 21, 22};
    TCanvas *c = new TCanvas("c_ds", "", 1800, 1050);
    c->Divide(3, 2, 0.002, 0.002);
    for (int t = 0; t < 6; t++)
    {
        c->cd(t + 1);
        gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.12);
        TH1D *fr = new TH1D(Form("fds%d", t), "", 100, 0.4, 2.0);
        fr->SetMinimum(0); fr->SetMaximum(1.4);
        fr->GetXaxis()->SetTitle("jet trigger value / nominal threshold");
        fr->GetYaxis()->SetTitle("fraction with the trigger bit set");
        fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05);
        fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045);
        fr->GetYaxis()->SetTitleOffset(1.1);
        fr->Draw("axis");
        TLine *l1 = new TLine(1.15, 0, 1.15, 1.06); l1->SetLineColor(kRed + 1); l1->SetLineWidth(2); l1->SetLineStyle(7); l1->Draw();
        TLine *l2 = new TLine(1.30, 0, 1.30, 1.06); l2->SetLineColor(kRed + 1); l2->SetLineWidth(2); l2->Draw();
        TLegend *lg = new TLegend(0.13, 0.72, 0.96, 0.905);
        lg->SetBorderSize(0); lg->SetFillColor(kWhite); lg->SetFillStyle(1001); lg->SetTextSize(0.038);
        lg->SetHeader("cuts: 1.15 (red dashed), 1.3 (red solid)");
        for (int k = 0; data_periods[k]; k++)
        {
            const int p = data_periods[k] - '0';
            if (p < 0 || p > 2) continue;
            TGraphAsymmErrors *g = (TGraphAsymmErrors *)fd->Get(Form("turnon_%s_p%d", trg[t], p));
            if (!g) continue;
            g = (TGraphAsymmErrors *)g->Clone();
            for (int i = 0; i < g->GetN(); i++) g->SetPoint(i, g->GetX()[i] - 0.006, g->GetY()[i]);
            g->SetMarkerStyle(pmk[p]); g->SetMarkerSize(0.9); g->SetMarkerColor(pcol[p]); g->SetLineColor(pcol[p]);
            g->Draw("p same");
            lg->AddEntry(g, plab[p], "pe");
        }
        TGraphAsymmErrors *gsim = (TGraphAsymmErrors *)fs->Get(Form("turnon_%s_p%d", trg[t], sim_period));
        if (gsim)
        {
            gsim = (TGraphAsymmErrors *)gsim->Clone();
            for (int i = 0; i < gsim->GetN(); i++) gsim->SetPoint(i, gsim->GetX()[i] + 0.006, gsim->GetY()[i]);
            gsim->SetMarkerStyle(24); gsim->SetMarkerSize(1.0); gsim->SetMarkerColor(kGreen + 3); gsim->SetLineColor(kGreen + 3);
            gsim->Draw("p same");
            lg->AddEntry(gsim, "simulation (emulated trigger)", "pe");
        }
        lg->Draw();
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.06);
        tl.DrawLatex(0.12, 0.93, Form("#bf{%s}  (nominal %g GeV)", trg[t], thr[t]));
    }
    c->SaveAs(out_png);
    TString pdf = out_png; pdf.ReplaceAll(".png", ".pdf"); c->SaveAs(pdf);
}
