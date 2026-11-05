// A_N vs x_F of two-photon jets inside vs outside the pi0 mass window, from
// two AnAsymmetry -m outputs (in:<peak> and out:<peak>; class np2, level
// "all": with -J the EM + JP union). Blue beam at x_F > 0, yellow at x_F < 0;
// filled = inside, open = outside (shifted slightly in x to stay readable).
// Usage: root -b -q 'plot_an_mass.C("an_mass_in.root", "an_mass_out.root", "an_vs_xf_mass.png", 0.115[, "sample"])'

void plot_an_mass(const char *in_root, const char *out_root, const char *out_png, double peak, const char *sample = "")
{
    gStyle->SetOptStat(0);
    TFile *f[2] = {TFile::Open(in_root), TFile::Open(out_root)};
    const char *beam[2] = {"blue", "yellow"};
    const int col[2] = {kBlue + 1, kOrange + 1};
    TGraphAsymmErrors *g[2][2];
    double ymin = 0, ymax = 0, xmax = 0;
    for (int w = 0; w < 2; w++)
        for (int b = 0; b < 2; b++)
        {
            TGraphAsymmErrors *g0 = f[w] ? (TGraphAsymmErrors *)f[w]->Get(Form("gr_AN_%s_all_np2", beam[b])) : 0;
            if (!g0) { printf("no gr_AN_%s_all_np2 in %s\n", beam[b], w ? out_root : in_root); return; }
            TGraphAsymmErrors *gg = (TGraphAsymmErrors *)g0->Clone(Form("g_%d_%d", w, b));
            for (int i = 0; i < gg->GetN(); i++)
            {
                const double x = gg->GetX()[i], y = gg->GetY()[i];
                gg->SetPoint(i, x + (w ? 0.006 : -0.006), y);
                // only points with a usable error set the axes
                if (gg->GetErrorYhigh(i) > 0.2) continue;
                const double lo = y - gg->GetErrorYlow(i), hi = y + gg->GetErrorYhigh(i);
                const double ax = (x < 0 ? -x : x) + gg->GetErrorXhigh(i);
                if (lo < ymin) ymin = lo;
                if (hi > ymax) ymax = hi;
                if (ax > xmax) xmax = ax;
            }
            gg->SetMarkerStyle(w ? 24 : 20); gg->SetMarkerSize(1.0);
            gg->SetMarkerColor(col[b]); gg->SetLineColor(col[b]); gg->SetLineWidth(2);
            g[w][b] = gg;
        }
    TCanvas *c = new TCanvas("c_anm", "", 1000, 760);
    c->SetGridy(); c->SetLeftMargin(0.13); c->SetTopMargin(0.08);
    double xr = 1.0;
    if (xmax > 0) xr = 0.05 * ceil((xmax + 0.02) / 0.05);
    if (xr > 1.0) xr = 1.0;
    TH1D *frame = new TH1D("frame_anm", "Run 22 FCS two-photon jets: A_{N} in vs out of the #pi^{0} mass window (EM + JP union);x_{F};A_{N}", 100, -xr, xr);
    const double span = ymax - ymin + 1e-6;
    frame->SetMinimum(ymin - 0.08 * span);
    frame->SetMaximum(ymax + 0.75 * span);
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->Draw();
    TLine *zero = new TLine(-xr, 0, xr, 0); zero->SetLineStyle(2); zero->Draw();
    for (int w = 0; w < 2; w++) for (int b = 0; b < 2; b++) g[w][b]->Draw("P SAME");
    TPave *box = new TPave(0.135, 0.70, 0.895, 0.915, 0, "NDC");
    box->SetFillColor(kWhite); box->SetFillStyle(1001); box->Draw();
    TLegend *lg = new TLegend(0.15, 0.77, 0.89, 0.91);
    lg->SetNColumns(2); lg->SetTextSize(0.027); lg->SetBorderSize(0); lg->SetFillStyle(0);
    // the window AnAsymmetry -m used: 25 MeV bins (AnConfig.h) whose centre is in 0.5-1.5 x peak
    double wlo = 1e9, whi = -1e9;
    for (int k = 0; k < 40; k++)
    {
        const double cen = (k + 0.5) * 0.025;
        if (cen >= 0.5 * peak && cen <= 1.5 * peak) { if (cen - 0.0125 < wlo) wlo = cen - 0.0125; if (cen + 0.0125 > whi) whi = cen + 0.0125; }
    }
    lg->AddEntry(g[0][0], Form("blue, M_{#gamma#gamma} in %.3f-%.3f GeV", wlo, whi), "pe");
    lg->AddEntry(g[1][0], "blue, M_{#gamma#gamma} outside", "pe");
    lg->AddEntry(g[0][1], "yellow (x_{F} < 0), in", "pe");
    lg->AddEntry(g[1][1], "yellow, outside", "pe");
    lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.024); t.SetTextFont(42);
    t.DrawLatex(0.155, 0.735, Form("jets with exactly 2 FCS points; window: 25 MeV mass bins with centre in 0.5-1.5 #times the fitted peak %.3f GeV%s%s", peak,
                                   TString(sample) != "" ? "; " : "", sample));
    gPad->RedrawAxis();
    c->SaveAs(out_png);
    TString pdf = out_png; pdf.ReplaceAll(".png", ".pdf"); c->SaveAs(pdf);
}
