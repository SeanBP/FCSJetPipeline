// A_N vs x_F of two-photon jets inside vs outside the pi0 mass window, from
// two AnAsymmetry -m outputs (in:<peak>:<sigma> and out:<peak>:<sigma>; class np2, level
// "all": with -J the EM + JP union). Blue beam at x_F > 0, yellow at x_F < 0;
// filled = inside, open = outside (shifted slightly in x to stay readable).
// Usage: root -b -q 'plot_an_mass.C("an_mass_in.root", "an_mass_out.root", "an_vs_xf_mass.png", 0.115[, "sample"[, y_lo, y_hi[, x_range[, err_max[, sigma]]]]])'
//   sigma: the fitted peak width (mgg_peak.C); > 0 = the window is peak +- 2
//     sigma (AnAsymmetry -m in:<peak>:<sigma>), 0 = the older 0.5-1.5 x peak
//   err_max: bins with a larger statistical error are left out, as in plot_an.C
//   y_lo, y_hi, x_range: fixed axes as in plot_an.C (0 = automatic), so the
//   plot matches the other A_N plots

// Width of a TLatex string in NDC at text size ts (font 42), measured on
// the current canvas (as plot_an.C)
double LatexWidth(const TString &str, double ts)
{
    TLatex l(0, 0, str.Data());
    l.SetNDC(); l.SetTextFont(42); l.SetTextSize(ts);
    UInt_t w = 0, h = 0;
    l.GetBoundingBox(w, h);
    return double(w) / gPad->GetWw();
}

// Splits s at spaces into lines at most wmax wide (NDC, text size ts), into
// out[0..]; returns the number of lines (as plot_an.C)
int WrapText(const TString &s, double wmax, double ts, TString *out, int nmax)
{
    int nout = 0;
    TObjArray *words = s.Tokenize(" ");
    TString line = "";
    for (int i = 0; i < words->GetEntriesFast(); i++)
    {
        TString w = ((TObjString *)words->At(i))->GetString();
        TString cand = line == "" ? w : line + " " + w;
        if (line != "" && LatexWidth(cand, ts) > wmax) { if (nout < nmax) out[nout++] = line; line = w; }
        else line = cand;
    }
    if (line != "" && nout < nmax) out[nout++] = line;
    delete words;
    return nout;
}

void plot_an_mass(const char *in_root, const char *out_root, const char *out_png, double peak, const char *sample = "",
                  double y_lo = 0, double y_hi = 0, double x_range = 0, double err_max = 0, double sigma = 0)
{
    gStyle->SetOptStat(0);
    TFile *f[2] = {TFile::Open(in_root), TFile::Open(out_root)};
    const char *beam[2] = {"blue", "yellow"};
    const int col[2] = {kBlue + 1, kOrange + 1};
    const int kMk[2][2] = {{20, 21}, {24, 25}};  // [inside, outside][blue, yellow]
    TGraphAsymmErrors *g[2][2];
    double ymin = 0, ymax = 0, xmax = 0;
    for (int w = 0; w < 2; w++)
        for (int b = 0; b < 2; b++)
        {
            TGraphAsymmErrors *g0 = f[w] ? (TGraphAsymmErrors *)f[w]->Get(Form("gr_AN_%s_all_np2", beam[b])) : 0;
            if (!g0) { printf("no gr_AN_%s_all_np2 in %s\n", beam[b], w ? out_root : in_root); return; }
            TGraphAsymmErrors *gg = (TGraphAsymmErrors *)g0->Clone(Form("g_%d_%d", w, b));
            // the innermost x_F bin (edge at 0) is not drawn, as in plot_an.C
            for (int i = gg->GetN() - 1; i >= 0; i--)
            {
                const double x0 = gg->GetX()[i];
                const double inner = x0 > 0 ? x0 - gg->GetErrorXlow(i) : -(x0 + gg->GetErrorXhigh(i));
                if (inner < 1e-6) { gg->RemovePoint(i); continue; }
                // too few jets (err_max as in plot_an.C; 0 = keep all)
                if (err_max > 0 && gg->GetErrorYhigh(i) > err_max) gg->RemovePoint(i);
            }
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
            // as plot_an.C: blue circle, yellow square; inside filled, outside open
            gg->SetMarkerStyle(kMk[w][b]); gg->SetMarkerSize(1.0);
            gg->SetMarkerColor(col[b]); gg->SetLineColor(col[b]); gg->SetLineWidth(2);
            g[w][b] = gg;
        }
    TCanvas *c = new TCanvas("c_anm", "", 1000, 760);
    c->SetGridy(); c->SetLeftMargin(0.13); c->SetTopMargin(0.05);
    double xr = 1.0;
    if (xmax > 0) xr = 0.05 * ceil((xmax + 0.02) / 0.05);
    if (xr > 1.0) xr = 1.0;
    if (x_range > 0) xr = x_range;
    // no title, as plot_an.C (FMS paper style); the data set heads the legend
    TH1D *frame = new TH1D("frame_anm", ";x_{F};A_{N}", 100, -xr, xr);
    const double span = ymax - ymin + 1e-6;
    frame->SetMinimum(ymin - 0.08 * span);
    frame->SetMaximum(ymax + 0.75 * span);
    if (y_hi > y_lo) { frame->SetMinimum(y_lo); frame->SetMaximum(y_hi); }
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->Draw();
    TLine *zero = new TLine(-xr, 0, xr, 0); zero->SetLineStyle(2); zero->Draw();
    // the window AnAsymmetry -m used: its TNamed "mass_window" ("lo hi width",
    // from 2026-11-11), else 25 MeV bins whose centre is in peak +- 2 sigma
    // (or, sigma = 0, in 0.5-1.5 x peak)
    const double mlo = sigma > 0 ? peak - 2 * sigma : 0.5 * peak, mhi = sigma > 0 ? peak + 2 * sigma : 1.5 * peak;
    double wlo = 1e9, whi = -1e9, bw = 0.025;
    TNamed *mw = f[0] ? (TNamed *)f[0]->Get("mass_window") : 0;
    if (mw) sscanf(mw->GetTitle(), "%lf %lf %lf", &wlo, &whi, &bw);
    else
        for (int k = 0; k < 40; k++)
        {
            const double cen = (k + 0.5) * 0.025;
            if (cen >= mlo && cen <= mhi) { if (cen - 0.0125 < wlo) wlo = cen - 0.0125; if (cen + 0.0125 > whi) whi = cen + 0.0125; }
        }

    // Layout as plot_an.C (Sean 2026-11-11): the notes in the upper left over
    // x_F < 0, wrapped to its width, on a white pane; points drawn over it
    const double xl = 0.16, xw = 0.505, dy = 0.030, ts = 0.022, yhead = 0.915;
    const char *ltext[4] = {Form("blue (x_{F} > 0), M_{#gamma#gamma} in %.3f-%.3f GeV", wlo, whi), "blue, M_{#gamma#gamma} outside",
                            "yellow (x_{F} < 0), in", "yellow, outside"};
    TString ltxt[4];
    for (int k = 0; k < 4; k++) ltxt[k] = ltext[k];
    TString notes[3];
    int nn = 0;
    if (sigma > 0)
        notes[nn++] = Form("jets with 2 photons; window: fitted #pi^{0} peak %.3f#pm2#sigma (#sigma = %.3f GeV)", peak, sigma);
    else
        notes[nn++] = Form("jets with 2 photons; window: %.0f MeV mass bins with centre in 0.5-1.5 #times the fitted peak %.3f GeV", 1000 * bw, peak);
    if (TString(sample) != "") notes[nn++] = sample;
    const int kMaxLines = 20;
    TString lines[kMaxLines], wrapped[10];
    int nl = 0;
    for (int k = 0; k < nn; k++)
    {
        const int nw = WrapText(notes[k], xw - xl - 0.01, ts, wrapped, 10);
        for (int j = 0; j < nw && nl < kMaxLines; j++) lines[nl++] = wrapped[j];
    }
    TString hl[2] = {"Run 22 FCS two-photon jets", "p^{#uparrow}+p #sqrt{s} = 508 GeV"};
    const int nhl = 2;
    const double ytop = yhead - 0.034 * nhl - 0.004;           // first marker row
    const double ynotes = ytop - 4 * 0.032 - 0.012;
    const double ybot = ynotes - dy * (nl - 0.5);
    TPave *box = new TPave(0.135, ybot - 0.012, xw, 0.945, 0, "NDC");
    box->SetFillColor(kWhite); box->SetFillStyle(1001); box->Draw();
    for (int ww = 0; ww < 2; ww++) for (int b = 0; b < 2; b++) g[ww][b]->Draw("P SAME");
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextAlign(12);
    t.SetTextSize(0.026);
    for (int k = 0; k < nhl; k++) t.DrawLatex(xl, yhead - 0.034 * k, hl[k]);
    // the four entries: in = filled, out = open; blue then yellow
    const int lw[4] = {0, 1, 0, 1}, lb[4] = {0, 0, 1, 1};
    for (int k = 0; k < 4; k++)
    {
        const double y = ytop - 0.032 * k;
        TMarker *m = new TMarker(0.0, 0.0, kMk[lw[k]][lb[k]]);
        m->SetNDC(); m->SetX(xl + 0.008); m->SetY(y); m->SetMarkerColor(col[lb[k]]); m->SetMarkerSize(1.0); m->Draw();
        t.SetTextSize(0.024); t.DrawLatex(xl + 0.025, y, ltxt[k]);
    }
    t.SetTextSize(ts);
    for (int k = 0; k < nl; k++) t.DrawLatex(xl, ynotes - dy * k, lines[k]);
    gPad->RedrawAxis();
    c->SaveAs(out_png);
    TString pdf = out_png; pdf.ReplaceAll(".png", ".pdf"); c->SaveAs(pdf);
}
