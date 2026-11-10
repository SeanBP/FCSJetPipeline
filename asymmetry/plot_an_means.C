// The per-bin means that go into A_N = eps / (<P> <|cos phi|>), from
// AnAsymmetry's table (an.txt): <|cos phi|> (left) and the jet-weighted beam
// polarization <P> (right) vs x_F, for one trigger level and the photon
// classes all, np1, np2, np12 and np3plus. Blue beam at x_F > 0, yellow at
// x_F < 0, as in the A_N plots; the x_F axis is the A_N plots' (x_range).
// Usage: root -b -q 'plot_an_means.C("an.txt", "an_means.png"[, "all"[, "sample"[, x_range]]])'

void plot_an_means(const char *table, const char *out_png, const char *level = "all", const char *sample = "", double x_range = 0.6)
{
    gStyle->SetOptStat(0);
    const int nc = 5;
    const char *cls[nc] = {"all", "np1", "np2", "np12", "np3plus"};
    const char *clab[nc] = {"all jets", "1 photon", "2 photons", "1 or 2 photons", "#geq 3 photons"};
    const int col[nc] = {kBlack, kAzure + 2, kOrange + 7, kGreen + 2, kMagenta + 1};
    const int mk[nc] = {20, 21, 22, 33, 23};
    const double msz[nc] = {1.1, 1.0, 1.2, 1.6, 1.2};
    TGraph *gc[nc], *gp[nc];
    for (int c = 0; c < nc; c++) { gc[c] = new TGraph(); gp[c] = new TGraph(); }
    ifstream in(table);
    if (!in) { printf("cannot open %s\n", table); return; }
    TString line;
    while (line.ReadLine(in))
    {
        if (line.BeginsWith("#")) continue;
        TObjArray *t = line.Tokenize(" ");
        if (t->GetEntries() < 11) { delete t; continue; }
        TString lv = ((TObjString *)t->At(0))->String(), cl = ((TObjString *)t->At(1))->String();
        const double x = ((TObjString *)t->At(5))->String().Atof();
        const double jets = ((TObjString *)t->At(8))->String().Atof();
        const double cs = ((TObjString *)t->At(9))->String().Atof(), p = ((TObjString *)t->At(10))->String().Atof();
        delete t;
        if (lv != level || jets < 100 || fabs(x) > x_range) continue;
        for (int c = 0; c < nc; c++)
            if (cl == cls[c]) { gc[c]->SetPoint(gc[c]->GetN(), x, cs); gp[c]->SetPoint(gp[c]->GetN(), x, p); }
    }
    // the trigger selection in words, as in plot_an.C
    TString lname = level;
    if (lname == "all") lname = TString(table).Contains("union") ? "EM + JP union" : "EM0 + EM1 + EM2";
    else if (lname == "allJP") lname = "JP0 + JP1 + JP2";
    TCanvas *cv = new TCanvas("c_means", "", 1600, 760);
    // header across both panels: what is plotted and the sample
    TLatex hd; hd.SetNDC(); hd.SetTextFont(42);
    hd.SetTextSize(0.034);
    hd.DrawLatex(0.015, 0.955, Form("Run 22 FCS %s jets: per-bin inputs to A_{N} = #varepsilon / (#LTP#GT #LT|cos #phi|#GT)", lname.Data()));
    if (sample[0]) { hd.SetTextSize(0.024); hd.DrawLatex(0.015, 0.915, sample); }
    const char *ptit[2] = {"Mean |cos #phi| of the jets", "Mean beam polarization, weighted by jets"};
    const char *ytit[2] = {"#LT|cos #phi|#GT", "#LTP#GT"};
    const double ylo[2] = {0.80, 0.53}, yhi[2] = {1.00, 0.62};  // top ~25% left for the legend
    for (int pnl = 0; pnl < 2; pnl++)
    {
        cv->cd();
        TPad *pad = new TPad(Form("pm%d", pnl), "", pnl * 0.5, 0.0, pnl * 0.5 + 0.5, 0.895);
        pad->Draw(); pad->cd();
        gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.12);
        TH1D *fr = new TH1D(Form("fr_means%d", pnl), Form(";x_{F}  (blue beam polarized at x_{F} > 0, yellow at x_{F} < 0);%s", ytit[pnl]), 100, -x_range, x_range);
        fr->SetMinimum(ylo[pnl]); fr->SetMaximum(yhi[pnl]);
        fr->GetYaxis()->SetTitleOffset(1.4); fr->GetXaxis()->SetTitleOffset(1.1);
        fr->GetXaxis()->SetTitleSize(0.04); fr->GetYaxis()->SetTitleSize(0.045);
        fr->Draw("axis");
        TLegend *lg = new TLegend(0.15, 0.755, 0.96, 0.905);
        lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
        for (int c = 0; c < nc; c++)
        {
            TGraph *g = pnl ? gp[c] : gc[c];
            g->Sort();
            g->SetMarkerStyle(mk[c]); g->SetMarkerSize(msz[c]); g->SetMarkerColor(col[c]); g->SetLineColor(col[c]); g->SetLineWidth(2);
            if (g->GetN()) g->Draw("P SAME");
            lg->AddEntry(g, clab[c], "p");
        }
        lg->Draw();
        TLatex tl; tl.SetNDC(); tl.SetTextFont(62); tl.SetTextSize(0.042);
        tl.DrawLatex(0.13, 0.94, ptit[pnl]);
    }
    cv->SaveAs(out_png);
    TString pdf = out_png; pdf.ReplaceAll(".png", ".pdf"); cv->SaveAs(pdf);
}
