// A_N vs x_F from AnAsymmetry's output: the blue beam at x_F > 0 and the
// yellow beam at x_F < 0, each point at its bin's mean x_F with the bin as the
// horizontal bar and the statistical error as the vertical bar (no
// systematic or polarization-scale uncertainty).
// Usage: root -b -q 'plot_an.C("an.root", "an_vs_xf.png"[, "all"[, "all"]])'
//   third argument: EM level (EM0, EM1, EM2 or all = their sum, the final result)
//   fourth argument: photon class (np0 ... np5plus, or all)
// Writes the .png and a .pdf next to it.

void plot_an(const char *in_root, const char *out_png, const char *level = "all", const char *cls = "all")
{
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open(in_root);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", in_root); return; }
    TGraphAsymmErrors *g[2];
    const char *beam[2] = {"blue", "yellow"};
    const int col[2] = {kBlue + 1, kOrange + 1}, mk[2] = {20, 21};
    double ymin = 0, ymax = 0;
    for (int b = 0; b < 2; b++)
    {
        g[b] = (TGraphAsymmErrors *)f->Get(Form("gr_AN_%s_%s_%s", beam[b], level, cls));
        if (!g[b]) { printf("no gr_AN_%s_%s_%s in %s\n", beam[b], level, cls, in_root); return; }
        g[b]->SetMarkerStyle(mk[b]); g[b]->SetMarkerSize(1.0);
        g[b]->SetMarkerColor(col[b]); g[b]->SetLineColor(col[b]); g[b]->SetLineWidth(2);
        for (int i = 0; i < g[b]->GetN(); i++)
        {
            ymin = TMath::Min(ymin, g[b]->GetY()[i] - g[b]->GetErrorYlow(i));
            ymax = TMath::Max(ymax, g[b]->GetY()[i] + g[b]->GetErrorYhigh(i));
        }
    }
    const double pad = 0.15 * (ymax - ymin + 1e-3);
    TCanvas *c = new TCanvas("c_an", "c_an", 900, 650);
    c->SetGridy();
    c->SetLeftMargin(0.12);
    TH1D *frame = new TH1D("frame_an", Form("Run 22 FCS EM-jet A_{N}, %s levels, %s photon classes;x_{F};A_{N}",
                                            level, cls), 100, -1.0, 1.0);
    frame->SetMinimum(ymin - pad); frame->SetMaximum(ymax + pad);
    frame->Draw();
    TLine *zero = new TLine(-1.0, 0, 1.0, 0);
    zero->SetLineStyle(2); zero->Draw();
    for (int b = 0; b < 2; b++) g[b]->Draw("P SAME");
    TLegend *lg = new TLegend(0.13, 0.72, 0.55, 0.88);
    lg->AddEntry(g[0], "blue beam polarized (x_{F} > 0)", "pe");
    lg->AddEntry(g[1], "yellow beam polarized (x_{F} < 0)", "pe");
    lg->AddEntry((TObject *)0, "statistical errors only", "");
    lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->Draw();
    c->SaveAs(out_png);
    TString pdf = out_png;
    pdf.ReplaceAll(".png", ".pdf");
    c->SaveAs(pdf);
}
