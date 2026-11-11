// A_N vs x_F from AnAsymmetry's output: the blue beam at x_F > 0 and the
// yellow beam at x_F < 0. Each point sits at the mean corrected x_F of the
// jets in its bin (as the FMS EM-jet analysis; the falling x_F spectrum puts
// it below the bin centre), its horizontal bar spans the bin edges (so it is
// asymmetric), and its vertical bar is the statistical error (no systematic or
// polarization-scale uncertainty). The plot states the jet selection and the
// sample, not the error-bar conventions.
// Usage: root -b -q 'plot_an.C("an.root", "an_vs_xf.png"[, "all"[, "all"[, err_scale[, at_zero[, zoom[, "sample"[, y_lo, y_hi[, x_range]]]]]]]])'
//   third argument: trigger level (EM0-2, JP0-2, all = the EM levels, or the
//   EM + JP union in a file from AnAsymmetry -J (name containing "union"),
//   allJP = the JP levels)
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
//   err_max: bins whose drawn statistical error exceeds it are left out
//     (too few jets); 0 = draw all. replot_report.sh passes AN_ERRMAX, 0.022
//     by default = the largest stat. error of the published FMS 500 GeV A_N
//     points (PRD 103, 092009, Figs. 4, 7, 9; Sean 2026-11-11)
//   show_formula: 1 draws the A_N formula (top right) and the N_L/N_R note;
//     default 0 (Sean 2026-11-11: no formula on the plots, it goes on the
//     slides)
//   ref: a reference file to overlay (fms_reference/*.txt): "set <label>"
//     starts a set, then lines "x_F A_N stat"; "#" starts a comment. Each set
//     is drawn with open markers and gets its own legend row. "" = none
// Writes the .png and a .pdf next to it.

// Removes the innermost x_F bin (the one with an edge at x_F = 0): very few
// jets, and the FMS results leave it out too (Sean 2026-11-11)
void DropInnermostBin(TGraphAsymmErrors *g)
{
    for (int i = g->GetN() - 1; i >= 0; i--)
    {
        const double x = g->GetX()[i];
        const double inner = x > 0 ? x - g->GetErrorXlow(i) : -(x + g->GetErrorXhigh(i));
        if (inner < 1e-6) g->RemovePoint(i);
    }
}

// Width of a TLatex string in NDC at text size ts (font 42), measured on
// the current canvas
double LatexWidth(const TString &str, double ts)
{
    TLatex l(0, 0, str.Data());
    l.SetNDC(); l.SetTextFont(42); l.SetTextSize(ts);
    UInt_t w = 0, h = 0;
    l.GetBoundingBox(w, h);
    return double(w) / gPad->GetWw();
}

// Splits s at spaces into lines at most wmax wide (NDC, text size ts), into
// out[0..]; returns the number of lines (at most nmax)
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

// Reference points (e.g. the FMS results) from a fms_reference/*.txt file
int ReadRef(const char *path, TGraphErrors **gr, TString *lab, double *shift, int nmax)
{
    for (int k = 0; k < nmax; k++) shift[k] = 0;
    std::ifstream in(path);
    if (!in) { printf("cannot open reference %s\n", path); return 0; }
    int n = -1;
    std::string line;
    while (std::getline(in, line))
    {
        TString l = line.c_str();
        if (l.Index("#") == 0) continue;
        if (l.BeginsWith("set "))
        {
            if (n + 1 >= nmax) break;
            n++;
            lab[n] = l(4, l.Length() - 4);
            gr[n] = new TGraphErrors();
            continue;
        }
        if (n < 0) continue;
        // "shift <dx>": draw this set's points dx to the right (for overlaps)
        if (l.BeginsWith("shift "))
        {
            double dx = atof(l(6, l.Length() - 6).Data());
            for (int i = 0; i < gr[n]->GetN(); i++) gr[n]->SetPoint(i, gr[n]->GetX()[i] + dx, gr[n]->GetY()[i]);
            shift[n] += dx;
            continue;
        }
        double x, y, e;
        if (sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &e) != 3) continue;
        int i = gr[n]->GetN();
        gr[n]->SetPoint(i, x + shift[n], y);
        gr[n]->SetPointError(i, 0, e);
    }
    return n + 1;
}

// The photon-class requirement in words (Sean 2026-11-06: "jets with n
// photons"; the classes count the jet's photon points, AnConfig.h)
TString ClassText(const char *cls)
{
    TString c = cls;
    if (c == "all") return "all jets";
    if (c == "np5plus") return "jets with #geq 5 photons";
    if (c == "np3plus") return "jets with > 2 photons";
    if (c == "np12") return "jets with 1 or 2 photons";
    if (c.BeginsWith("np"))
    {
        TString n = c(2, c.Length() - 2);
        return Form("jets with %s photon%s", n.Data(), n == "1" ? "" : "s");
    }
    return c;
}

void plot_an(const char *in_root, const char *out_png, const char *level = "all", const char *cls = "all",
             double err_scale = 1, int at_zero = 0, double zoom = 10, const char *sample = "",
             double y_lo = 0, double y_hi = 0, double x_range = 0, const char *ref = "", double err_max = 0, int show_formula = 0)
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
        DropInnermostBin(g[b]);
        // bins with too few jets: statistical error (as drawn, i.e. after
        // err_scale) above err_max are not drawn; 0 = keep all
        if (err_max > 0)
            for (int i = g[b]->GetN() - 1; i >= 0; i--)
                if (err_scale * g[b]->GetErrorYhigh(i) > err_max) g[b]->RemovePoint(i);
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
    c->SetTopMargin(0.05);
    // x range: the covered points plus a margin, symmetric, in steps of 0.05
    double xr = 1.0;
    if (xmax > 0) xr = 0.05 * ceil((xmax + 0.02) / 0.05);
    if (xr > 1.0) xr = 1.0;
    if (x_range > 0) xr = x_range;
    const char *lname = level;
    // EM + JP (the union) is the default and is not named (Sean 2026-11-11);
    // a title names the trigger only when the plot isolates some
    // (the pi0 mass-window files an_mass_in/out are -J union output too)
    const bool is_union = TString(in_root).Contains("union") || TString(in_root).Contains("mass_");
    if (TString(level) == "all") lname = is_union ? "" : "EM0 + EM1 + EM2";
    if (TString(level) == "allJP") lname = "JP0 + JP1 + JP2";
    // no title (Sean 2026-11-11: as the FMS paper's figures); the data set
    // is the first legend line instead
    TH1D *frame = new TH1D("frame_an", ";x_{F};A_{N}", 100, -xr, xr);
    const TString header = Form("Run 22 FCS EM-jets, p^{#uparrow}+p #sqrt{s} = 508 GeV%s%s", TString(lname) != "" ? ", " : "", lname);
    // automatic axes leave the top 45% of the frame for the text
    const double span = ymax - ymin + 1e-6;
    frame->SetMinimum(ymin - 0.08 * span);
    frame->SetMaximum(ymax + (err_scale != 1 ? 1.05 : 0.85) * span);
    if (y_hi > y_lo) { frame->SetMinimum(y_lo); frame->SetMaximum(y_hi); }
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->Draw();
    TLine *zero = new TLine(-xr, 0, xr, 0);
    zero->SetLineStyle(2); zero->Draw();
    // reference sets (open markers, drawn under our points)
    TGraphErrors *gref[4];
    TString lref[4];
    double sref[4];
    int nref = 0;
    if (TString(ref) != "") nref = ReadRef(ref, gref, lref, sref, 4);
    // open markers, larger than ours and in colours apart from blue/orange
    const int rcol[4] = {kBlack, kRed + 1, kGreen + 3, kMagenta + 2}, rmk[4] = {24, 25, 26, 27};
    for (int r = 0; r < nref; r++)
    {
        gref[r]->SetMarkerStyle(rmk[r]); gref[r]->SetMarkerSize(1.5);
        gref[r]->SetMarkerColor(rcol[r]); gref[r]->SetLineColor(rcol[r]); gref[r]->SetLineWidth(2);
        if (sref[r] != 0) lref[r] += Form(", shifted by %+.3g in x_{F}", sref[r]);
        lref[r] += " (stat. errors)";
    }

    // Layout (Sean 2026-11-11): the notes sit in the upper left, over the
    // x_F < 0 side, wrapped to its width, so the rising x_F > 0 points stay
    // visible up to the formula at the top right. White panes behind the
    // text hide the grid; the points are drawn over them, the text last.
    const double xl = 0.16, xw = 0.505, dy = 0.030, ts = 0.022;
    // the note lines, with an optional marker (style, colour) in front
    const int kMaxLines = 40;
    TString lines[kMaxLines];
    int lmk[kMaxLines], lcol[kMaxLines], nl = 0;
    TString ctext = ClassText(cls);
    // AnAsymmetry -m in/out output (np2 rows = two-photon jets in/out of the
    // pi0 mass window, the fitted peak +- 2 sigma): say so
    if (TString(in_root).Contains("mass_in")) ctext += ", pair mass within 2#sigma of the fitted #pi^{0} peak";
    if (TString(in_root).Contains("mass_out")) ctext += ", pair mass outside the #pi^{0} window";
    TString notes[4];
    int nn = 0;
    notes[nn++] = ctext;
    if (TString(sample) != "") notes[nn++] = sample;
    if (err_scale != 1)
        notes[nn++] = (Form("statistical errors projected to all good runs (test errors #times %.4g)%s", err_scale,
                             at_zero ? ", drawn at A_{N} = 0" : ""));
    // the formula's N_L/N_R, only with the formula
    if (show_formula) notes[nn++] = ("N_{L(R)}^{#uparrow(#downarrow)}: jets left (right) of the polarized beam, its spin up (down); the other beam's spin averaged");
    TString w[10];
    for (int k = 0; k < nn; k++)
    {
        const int nw = WrapText(notes[k], xw - xl - 0.01, ts, w, 10);
        for (int j = 0; j < nw && nl < kMaxLines; j++) { lines[nl] = w[j]; lmk[nl] = 0; lcol[nl] = 0; nl++; }
    }
    for (int r = 0; r < nref; r++)
    {
        const int nw = WrapText(lref[r], xw - xl - 0.035, ts, w, 10);
        for (int j = 0; j < nw && nl < kMaxLines; j++) { lines[nl] = w[j]; lmk[nl] = j == 0 ? rmk[r] : -1; lcol[nl] = rcol[r]; nl++; }
    }
    // the header (wrapped), then the two beam rows, then the notes
    TString hl[3];
    const int nhl = WrapText(header, xw - xl, 0.026, hl, 3);
    const double yhead = 0.915;
    const double ytop = yhead - 0.034 * nhl - 0.004;          // first beam row
    const double ynotes = ytop - 2 * 0.034 - 0.012;          // below the two beam rows
    const double ybot = ynotes - dy * (nl - 0.5);
    TPave *box = new TPave(0.135, ybot - 0.012, xw, 0.945, 0, "NDC");
    box->SetFillColor(kWhite); box->SetFillStyle(1001);
    box->Draw();
    if (show_formula)
    {
        TPave *eqbox = new TPave(0.525, 0.835, 0.895, 0.945, 0, "NDC");
        eqbox->SetFillColor(kWhite); eqbox->SetFillStyle(1001);
        eqbox->Draw();
    }
    for (int r = 0; r < nref; r++) gref[r]->Draw("P SAME");
    for (int b = 0; b < 2; b++) g[b]->Draw("P SAME");

    // the beams
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextAlign(12);
    t.SetTextSize(0.026);
    for (int k = 0; k < nhl; k++) t.DrawLatex(xl, yhead - 0.034 * k, hl[k]);
    const char *btext[2] = {"blue beam polarized (x_{F} > 0)", "yellow beam polarized (x_{F} < 0)"};
    for (int b = 0; b < 2; b++)
    {
        const double y = ytop - 0.034 * b;
        TMarker *m = new TMarker(0.0, 0.0, mk[b]);
        m->SetNDC(); m->SetX(xl + 0.008); m->SetY(y); m->SetMarkerColor(col[b]); m->SetMarkerSize(1.0); m->Draw();
        t.SetTextSize(0.026); t.DrawLatex(xl + 0.025, y, btext[b]);
    }
    // the notes and the reference sets
    t.SetTextSize(ts);
    for (int k = 0; k < nl; k++)
    {
        const double y = ynotes - dy * k;
        double x = xl;
        if (lmk[k] != 0)
        {
            x = xl + 0.025;
            if (lmk[k] > 0)
            {
                TMarker *m = new TMarker(0.0, 0.0, lmk[k]);
                m->SetNDC(); m->SetX(xl + 0.008); m->SetY(y); m->SetMarkerColor(lcol[k]); m->SetMarkerSize(1.3); m->Draw();
            }
        }
        t.DrawLatex(x, y, lines[k]);
    }
    // the A_N formula (AnAsymmetry.cpp: cross ratio, other beam's spin held
    // fixed and the two states averaged), top right
    TLatex eq; eq.SetNDC(); eq.SetTextFont(42); eq.SetTextSize(0.031); eq.SetTextAlign(12);
    if (show_formula) eq.DrawLatex(0.535, 0.89, "A_{N} = #frac{1}{#LTP#GT #LT|cos #phi|#GT} #frac{#sqrt{N_{L}^{#uparrow}N_{R}^{#downarrow}} - #sqrt{N_{L}^{#downarrow}N_{R}^{#uparrow}}}{#sqrt{N_{L}^{#uparrow}N_{R}^{#downarrow}} + #sqrt{N_{L}^{#downarrow}N_{R}^{#uparrow}}}");
    gPad->RedrawAxis();
    c->SaveAs(out_png);
    TString pdf = out_png;
    pdf.ReplaceAll(".png", ".pdf");
    c->SaveAs(pdf);
}
