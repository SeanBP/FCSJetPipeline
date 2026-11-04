// A_N vs x_F from AnAsymmetry's output: the blue beam at x_F > 0 and the
// yellow beam at x_F < 0. Each point sits at the mean corrected x_F of the
// jets in its bin (as the FMS EM-jet analysis; the falling x_F spectrum puts
// it below the bin centre), its horizontal bar spans the bin edges (so it is
// asymmetric), and its vertical bar is the statistical error (no systematic or
// polarization-scale uncertainty). The plot states the jet selection and the
// sample, not the error-bar conventions.
// Usage: root -b -q 'plot_an.C("an.root", "an_vs_xf.png"[, "all"[, "all"[, err_scale[, at_zero[, zoom[, "sample"[, y_lo, y_hi[, x_range]]]]]]]])'
//   third argument: EM level (EM0, EM1, EM2 or all = their sum, the final result)
//   fourth argument: photon class (np0 ... np5plus, np3plus = 3 or more,
//     np12 = 1 or 2, or all)
//   err_scale: multiplies the errors, for a projection to more data
//     (project_errors.py prints it); 1 = the measured errors
//   at_zero: 1 draws the points at A_N = 0, to show the projected errors alone
//   zoom: the axes cover the points whose error is at most zoom x the smallest
//     error (default 10); the others are still drawn, clipped at the frame.
//     0 = cover every point
//   y_lo, y_hi, x_range: fixed axes (A_N from y_lo to y_hi, x_F from
//     -x_range to x_range), to draw several plots alike; 0 = automatic. Keep
//     the top ~30% of the A_N range free for the text
//   sample: a line saying what data the plot shows, e.g. "560 of 2892 good
//     runs, 2.3% of the Run 22 good-run events" (smoke_report.sh fills it in)
// Writes the .png and a .pdf next to it.

// The photon-class requirement in words (AnConfig.h: reco_nphotons, FCS
// points above 1 GeV within R = 0.7 of the jet axis)
TString ClassText(const char *cls)
{
    TString c = cls;
    if (c == "all") return "all jets";
    if (c == "np5plus") return "jets with #geq 5 FCS points (E > 1 GeV) within R = 0.7 of the axis";
    if (c == "np3plus") return "jets with #geq 3 FCS points (E > 1 GeV) within R = 0.7 of the axis";
    if (c == "np12") return "jets with 1 or 2 FCS points (E > 1 GeV) within R = 0.7 of the axis";
    if (c.BeginsWith("np"))
    {
        TString n = c(2, c.Length() - 2);
        return Form("jets with exactly %s FCS point%s (E > 1 GeV) within R = 0.7 of the axis", n.Data(), n == "1" ? "" : "s");
    }
    return c;
}

void plot_an(const char *in_root, const char *out_png, const char *level = "all", const char *cls = "all",
             double err_scale = 1, int at_zero = 0, double zoom = 10, const char *sample = "",
             double y_lo = 0, double y_hi = 0, double x_range = 0)
{
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open(in_root);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", in_root); return; }
    TGraphAsymmErrors *g[2];
    const char *beam[2] = {"blue", "yellow"};
    const int col[2] = {kBlue + 1, kOrange + 1}, mk[2] = {20, 21};
    double ymin = 0, ymax = 0, xmax = 0, emin = 1e30;
    for (int b = 0; b < 2; b++)
    {
        TGraphAsymmErrors *g0 = (TGraphAsymmErrors *)f->Get(Form("gr_AN_%s_%s_%s", beam[b], level, cls));
        if (g0) for (int i = 0; i < g0->GetN(); i++) if (g0->GetErrorYhigh(i) > 0) emin = TMath::Min(emin, g0->GetErrorYhigh(i));
    }
    double ecut = 1e30;
    if (zoom > 0) ecut = zoom * emin * err_scale;
    for (int b = 0; b < 2; b++)
    {
        g[b] = (TGraphAsymmErrors *)f->Get(Form("gr_AN_%s_%s_%s", beam[b], level, cls));
        if (!g[b]) { printf("no gr_AN_%s_%s_%s in %s\n", beam[b], level, cls, in_root); return; }
        g[b] = (TGraphAsymmErrors *)g[b]->Clone(Form("g_%s", beam[b]));
        g[b]->SetMarkerStyle(mk[b]); g[b]->SetMarkerSize(1.0);
        g[b]->SetMarkerColor(col[b]); g[b]->SetLineColor(col[b]); g[b]->SetLineWidth(2);
        for (int i = 0; i < g[b]->GetN(); i++)
        {
            double y = g[b]->GetY()[i];
            if (at_zero) y = 0;
            g[b]->SetPoint(i, g[b]->GetX()[i], y);
            g[b]->SetPointEYlow(i, err_scale * g[b]->GetErrorYlow(i));
            g[b]->SetPointEYhigh(i, err_scale * g[b]->GetErrorYhigh(i));
            if (g[b]->GetErrorYhigh(i) > ecut) continue;
            ymin = TMath::Min(ymin, y - g[b]->GetErrorYlow(i));
            ymax = TMath::Max(ymax, y + g[b]->GetErrorYhigh(i));
            xmax = TMath::Max(xmax, TMath::Abs(g[b]->GetX()[i]) + TMath::Max(g[b]->GetErrorXlow(i), g[b]->GetErrorXhigh(i)));
        }
    }
    TCanvas *c = new TCanvas("c_an", "c_an", 1000, 760);
    c->SetGridy();
    c->SetLeftMargin(0.13);
    c->SetTopMargin(0.08);
    // x range: the covered points plus a margin, symmetric, in steps of 0.05
    double xr = 1.0;
    if (xmax > 0) xr = 0.05 * ceil((xmax + 0.02) / 0.05);
    if (xr > 1.0) xr = 1.0;
    if (x_range > 0) xr = x_range;
    const char *lname = level;
    if (TString(level) == "all") lname = "EM0 + EM1 + EM2";
    TH1D *frame = new TH1D("frame_an", Form("Run 22 FCS EM-jet A_{N}, %s;x_{F};A_{N}", lname), 100, -xr, xr);
    // the top 45% of the frame is left for the text
    const double span = ymax - ymin + 1e-6;
    frame->SetMinimum(ymin - 0.08 * span);
    frame->SetMaximum(ymax + (err_scale != 1 ? 1.05 : 0.85) * span);
    if (y_hi > y_lo) { frame->SetMinimum(y_lo); frame->SetMaximum(y_hi); }
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->Draw();
    TLine *zero = new TLine(-xr, 0, xr, 0);
    zero->SetLineStyle(2); zero->Draw();
    for (int b = 0; b < 2; b++) g[b]->Draw("P SAME");
    gPad->RedrawAxis();

    // legend and notes on a white box, so error bars running into the top
    // of the frame do not cross the text
    TPave *box = new TPave(0.135, err_scale != 1 ? 0.635 : 0.67, 0.895, 0.915, 0, "NDC");
    box->SetFillColor(kWhite); box->SetFillStyle(1001);
    box->Draw();
    TLegend *lg = new TLegend(0.15, 0.82, 0.55, 0.91);
    lg->SetTextSize(0.028);
    lg->AddEntry(g[0], "blue beam polarized (x_{F} > 0)", "pe");
    lg->AddEntry(g[1], "yellow beam polarized (x_{F} < 0)", "pe");
    lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->Draw();
    // the A_N formula (AnAsymmetry.cpp: cross ratio, other beam's spin held
    // fixed and the two states averaged), right of the legend
    TLatex eq; eq.SetNDC(); eq.SetTextFont(42); eq.SetTextSize(0.031); eq.SetTextAlign(12);
    eq.DrawLatex(0.535, 0.855, "A_{N} = #frac{1}{#LTP#GT #LT|cos #phi|#GT} #frac{#sqrt{N_{L}^{#uparrow}N_{R}^{#downarrow}} - #sqrt{N_{L}^{#downarrow}N_{R}^{#uparrow}}}{#sqrt{N_{L}^{#uparrow}N_{R}^{#downarrow}} + #sqrt{N_{L}^{#downarrow}N_{R}^{#uparrow}}}");
    TLatex t; t.SetNDC(); t.SetTextSize(0.024); t.SetTextFont(42);
    double yl = 0.79;
    t.DrawLatex(0.155, yl, ClassText(cls)); yl -= 0.035;
    if (TString(sample) != "") { t.DrawLatex(0.155, yl, sample); yl -= 0.035; }
    if (err_scale != 1)
        t.DrawLatex(0.155, yl, Form("statistical errors projected to all good runs (test errors #times %.4g)%s", err_scale,
                                    at_zero ? ", drawn at A_{N} = 0" : ""));
    t.DrawLatex(0.155, yl - (err_scale != 1 ? 0.035 : 0), "N_{L(R)}^{#uparrow(#downarrow)}: jets left (right) of the polarized beam, its spin up (down); per spin of the other beam, averaged");
    gPad->RedrawAxis();
    c->SaveAs(out_png);
    TString pdf = out_png;
    pdf.ReplaceAll(".png", ".pdf");
    c->SaveAs(pdf);
}
