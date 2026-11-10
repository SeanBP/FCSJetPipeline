// The two factors A_N is divided by, A_N = eps / (<P> <|cos phi|>), from
// AnAsymmetry's table (<out>.txt: level class beam xF_lo xF_hi <xF> A_N stat
// jets <|cos|> <P>), and the per-run polarization behind <P>:
//   left    <P> per x_F bin (jet-weighted mean of the runs' P), blue and
//           yellow vs |x_F|, for EM0 + EM1 + EM2 and each level
//   middle  <|cos phi|> per x_F bin (phi = jet azimuth from the beam's left),
//           the same for both beams, for all jets and each photon class
//   right   P_run of every run in the polarization table (grey) and of the
//           runs in the yields (coloured), blue and yellow, vs run index
// Bins with fewer than min_jets jets are not drawn.
// Usage (compiled, CINT loops on the stream parsing): root -b -q 'plot_pol_cos.C+("an.txt", "run22_polarization.txt", "yields.list", "pol_cos.png"[, min_jets = 100])'
//   yields.list: the yield files AnAsymmetry -r read (an_hists_<run>_<n>.root)

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>

#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>

TGraph *Pick(const char *table, const char *level, const char *cls, const char *beam, int col, int mk, bool pol, double min_jets)
{
    TGraph *g = new TGraph();
    std::ifstream in(table);
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string l, c, b;
        double lo, hi, mxf, an, st, jets, cs, P;
        if (!(ss >> l >> c >> b >> lo >> hi >> mxf >> an >> st >> jets >> cs >> P)) continue;
        if (l != level || c != cls || b != beam || jets < min_jets) continue;
        g->SetPoint(g->GetN(), TMath::Abs(mxf), pol ? P : cs);
    }
    g->Sort();
    g->SetMarkerStyle(mk); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerSize(1.0); g->SetLineWidth(2);
    return g;
}

void plot_pol_cos(const char *table, const char *poltab, const char *yields, const char *out_png, double min_jets = 100)
{
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c_pc", "", 1800, 600);
    c->Divide(3, 1);

    // left: <P>
    c->cd(1);
    gPad->SetGridy(); gPad->SetLeftMargin(0.13);
    TH1D *f1 = new TH1D("f1", "mean polarization of the jets per x_{F} bin;|x_{F}| (bin's mean);#LTP#GT", 10, 0, 0.65);
    f1->SetMinimum(0.50); f1->SetMaximum(0.60); f1->GetYaxis()->SetTitleOffset(1.5); f1->Draw();
    const char *lv[4] = {"all", "EM0", "EM1", "EM2"};
    const int mkb[4] = {20, 24, 25, 26}, mky[4] = {21, 24, 25, 26};
    TLegend *l1 = new TLegend(0.15, 0.12, 0.89, 0.34);
    l1->SetBorderSize(0); l1->SetFillStyle(0); l1->SetTextSize(0.030);
    for (int i = 0; i < 4; i++)
    {
        TGraph *gb = Pick(table, lv[i], "all", "blue", kBlue + 1, mkb[i], true, min_jets);
        TGraph *gy = Pick(table, lv[i], "all", "yellow", kOrange + 1, mky[i], true, min_jets);
        if (i) { gb->SetLineStyle(2); gy->SetLineStyle(2); gb->SetLineWidth(1); gy->SetLineWidth(1); }
        gb->Draw("PL SAME"); gy->Draw("PL SAME");
        const char *nm = i ? lv[i] : "EM0+1+2";
        l1->AddEntry(gb, Form("blue, %s", nm), "pl");
        l1->AddEntry(gy, Form("yellow, %s", nm), "pl");
    }
    l1->SetNColumns(2);
    l1->Draw();

    // middle: <|cos phi|>
    c->cd(2);
    gPad->SetGridy(); gPad->SetLeftMargin(0.13);
    TH1D *f2 = new TH1D("f2", "mean |cos #phi| of the jets per x_{F} bin;|x_{F}| (bin's mean);#LT|cos #phi|#GT", 10, 0, 0.65);
    f2->SetMinimum(0.80); f2->SetMaximum(1.0); f2->GetYaxis()->SetTitleOffset(1.5); f2->Draw();
    const char *cl[5] = {"all", "np1", "np2", "np3", "np5plus"};
    const char *cn[5] = {"all jets", "1 point", "2 points", "3 points", "#geq 5 points"};
    const int cc[5] = {kBlack, kGreen + 2, kAzure + 1, kMagenta + 1, kRed + 1}, cm[5] = {20, 24, 25, 26, 32};
    TLegend *l2 = new TLegend(0.15, 0.12, 0.60, 0.36);
    l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.033);
    for (int i = 0; i < 5; i++)
    {
        TGraph *g = Pick(table, "all", cl[i], "blue", cc[i], cm[i], false, min_jets);
        if (i) { g->SetLineStyle(2); g->SetLineWidth(1); }
        g->Draw("PL SAME");
        l2->AddEntry(g, Form("%s (EM0 + EM1 + EM2)", cn[i]), "pl");
    }
    l2->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.032); t.SetTextFont(42);
    t.DrawLatex(0.16, 0.84, "#phi: jet azimuth from the beam's left (#pm x),");
    t.DrawLatex(0.16, 0.79, "the same for both beams");

    // right: per-run P
    c->cd(3);
    gPad->SetLeftMargin(0.13); gPad->SetGridy();
    std::set<int> used;
    {
        std::ifstream in(yields);
        std::string f;
        while (in >> f)
        {
            size_t p = f.find("an_hists_");
            if (p != std::string::npos) used.insert(atoi(f.c_str() + p + 9));
        }
    }
    TGraph *ab = new TGraph(), *ay = new TGraph(), *ub = new TGraph(), *uy = new TGraph();
    std::ifstream in(poltab);
    std::string line;
    int idx = 0, first = 0, last = 0;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int run, fill; double tm, pb, dpb, py, dpy;
        if (!(ss >> run >> fill >> tm >> pb >> dpb >> py >> dpy)) continue;
        if (!first) first = run;
        last = run;
        if (pb > 0) ab->SetPoint(ab->GetN(), idx, pb);
        if (py > 0) ay->SetPoint(ay->GetN(), idx, py);
        if (used.count(run)) { if (pb > 0) ub->SetPoint(ub->GetN(), idx, pb); if (py > 0) uy->SetPoint(uy->GetN(), idx, py); }
        idx++;
    }
    TH1D *f3 = new TH1D("f3", Form("P per run (fill's P at the run's middle);run index (%d to %d);P", first, last), 10, 0, idx);
    f3->SetMinimum(0.2); f3->SetMaximum(0.9); f3->GetYaxis()->SetTitleOffset(1.5); f3->Draw();
    ab->SetMarkerStyle(7); ab->SetMarkerColor(kGray + 1); ab->Draw("P SAME");
    ay->SetMarkerStyle(7); ay->SetMarkerColor(kGray); ay->Draw("P SAME");
    ub->SetMarkerStyle(20); ub->SetMarkerSize(0.5); ub->SetMarkerColor(kBlue + 1); ub->Draw("P SAME");
    uy->SetMarkerStyle(21); uy->SetMarkerSize(0.5); uy->SetMarkerColor(kOrange + 1); uy->Draw("P SAME");
    TLegend *l3 = new TLegend(0.14, 0.12, 0.86, 0.26);
    l3->SetBorderSize(0); l3->SetFillStyle(1001); l3->SetFillColor(kWhite); l3->SetTextSize(0.029);
    l3->AddEntry(ub, Form("blue, runs in the yields (%d): mean %.3f", ub->GetN(), ub->GetMean(2)), "p");
    l3->AddEntry(uy, Form("yellow, runs in the yields (%d): mean %.3f", uy->GetN(), uy->GetMean(2)), "p");
    l3->AddEntry(ab, Form("every run in the polarization table (%d)", idx), "p");
    l3->Draw();

    c->SaveAs(out_png);
    TString pdf = out_png;
    pdf.ReplaceAll(".png", ".pdf");
    c->SaveAs(pdf);
}
