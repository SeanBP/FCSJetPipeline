// The pi0 peak of the two-photon jets' pair mass, for AnAsymmetry -m: sums
// AnHistograms' x_F vs mass yields of the triggering jets in every yield file
// (EM levels at the default factor plus the JP jets not firing EM,
// h_xFm_<level>_<side>_<spin>[_jponly]: each jet once, E <= sqrt(s)/2),
// projects them on the mass and fits a Gaussian plus a linear background
// around the highest bin in 0.05-0.30 GeV. Prints "peak=<GeV>" (the Gaussian
// mean), "sigma=" (AnAsymmetry -m in:<peak>:<sigma> takes peak +- 2 sigma),
// and writes <out_png> with the fit and that window shaded.
// Usage: root -b -q 'mgg_peak.C("yields.list", "mgg_peak.png")'

void mgg_peak(const char *yield_list, const char *out_png)
{
    gStyle->SetOptStat(0);
    const char *lv[6] = {"EM0", "EM1", "EM2", "JP0", "JP1", "JP2"};
    const char *sd[2] = {"xneg", "xpos"};
    const char *sp[4] = {"BuYu", "BuYd", "BdYu", "BdYd"};
    TH1D *m = 0;
    ifstream in(yield_list);
    TString fn;
    int nf = 0;
    while (fn.ReadLine(in))
    {
        fn = fn.Strip(TString::kBoth);
        if (fn.Length() == 0) continue;
        TFile f(fn.Data());
        if (f.IsZombie()) continue;
        nf++;
        for (int l = 0; l < 6; l++)
            for (int s = 0; s < 2; s++)
                for (int p = 0; p < 4; p++)
                {
                    TH2D *h = (TH2D *)f.Get(Form("h_xFm_%s_%s_%s%s", lv[l], sd[s], sp[p], l >= 3 ? "_jponly" : ""));
                    if (!h) continue;
                    TH1D *py = h->ProjectionY("py_tmp");
                    if (!m) { m = (TH1D *)py->Clone("mgg_sum"); m->SetDirectory(0); m->Reset(); }
                    m->Add(py);
                    delete py;
                }
    }
    if (!m || m->Integral() <= 0) { printf("peak=0\nno h_xFm_* yields in %d files\n", nf); return; }
    int bmax = 0; double vmax = -1;
    for (int b = 1; b <= m->GetNbinsX(); b++)
    {
        const double x = m->GetBinCenter(b);
        if (x < 0.05 || x > 0.30) continue;
        if (m->GetBinContent(b) > vmax) { vmax = m->GetBinContent(b); bmax = b; }
    }
    const double x0 = m->GetBinCenter(bmax);
    TF1 *g = new TF1("g", "gaus(0)+pol1(3)", TMath::Max(0.0, x0 - 0.1), x0 + 0.15);
    g->SetParameters(vmax, x0, 0.03, 0, 0);
    g->SetParLimits(2, 0.005, 0.1);
    m->Fit(g, "RQ0");
    const double peak = g->GetParameter(1), sig = g->GetParameter(2);
    printf("peak=%.4f\nsigma=%.4f\nfiles=%d jets=%.0f\n", peak, sig, nf, m->Integral());
    TCanvas *c = new TCanvas("c_mgg", "", 900, 650);
    m->SetTitle(";M_{#gamma#gamma} of two-photon jets (GeV);jets");
    m->SetLineColor(kBlack); m->SetLineWidth(2);
    m->GetXaxis()->SetRangeUser(0, 0.6);
    m->SetMaximum(m->GetMaximum() * 1.15);
    m->Draw("hist");
    // the window AnAsymmetry -m uses: whole bins whose centre is in peak +- 2 sigma
    double wlo = 1e9, whi = -1e9;
    for (int b = 1; b <= m->GetNbinsX(); b++)
        if (m->GetBinCenter(b) >= peak - 2 * sig && m->GetBinCenter(b) <= peak + 2 * sig)
        {
            if (m->GetBinLowEdge(b) < wlo) wlo = m->GetBinLowEdge(b);
            if (m->GetBinLowEdge(b + 1) > whi) whi = m->GetBinLowEdge(b + 1);
        }
    TBox *in_w = new TBox(wlo, 0, whi, m->GetMaximum());
    in_w->SetFillColor(kAzure - 9); in_w->Draw();
    m->Draw("hist same"); g->SetLineColor(kRed + 1); g->Draw("same");
    gPad->RedrawAxis();
    TLatex t; t.SetNDC(); t.SetTextSize(0.035);
    t.DrawLatex(0.45, 0.85, Form("peak %.3f GeV, #sigma %.3f GeV", peak, sig));
    t.DrawLatex(0.45, 0.80, Form("shaded: %.3f-%.3f GeV", wlo, whi));
    t.DrawLatex(0.45, 0.75, "(bins with centre in peak #pm 2#sigma)");
    c->SaveAs(out_png);
}
