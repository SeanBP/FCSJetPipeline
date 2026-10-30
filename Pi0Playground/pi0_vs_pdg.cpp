// Reconstructed pi0 mass against the PDG mass (134.977 MeV), data and
// simulation through the same analysis: the ECal gain calibration's pi0
// selection (pi0_pairs' h3_liang_*, default clustering; "-f": the variant
// with the common 45 MeV tower floor, h3_liangf_*), Gaussian + quadratic fit
// as pi0_scale_plots. Left: the mass spectrum over all pair energies, which
// is what the per-tower calibration fits; right: the peak against pair
// energy. Peaks are printed as a table. "-s2 <sim2.root>" overlays a second
// simulation (e.g. k = 0.93), "-l1/-l2 <label>" name the two simulations.
//
// Build and run from this folder:
//   g++ -O2 pi0_vs_pdg.cpp -o pi0_vs_pdg `root-config --cflags --libs` -std=c++11
//   ./pi0_vs_pdg <data_all.root> <sim.root> <out.png> [-f] [-s2 <sim2.root>] [-l1 <label>] [-l2 <label>]

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPad.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

static const double kPdg = 0.134977;

struct Peak { double m, em, s; bool ok; };
// Gaussian + quadratic, as pi0_scale_plots' FitPi0
static Peak FitPi0(TH1D *h, TF1 *keep = 0)
{
    const double lo = 0.05, hi = 0.30, mlo = 0.09, mhi = 0.19, slo = 0.006, shi = 0.045;
    Peak p = {0, 0, 0, false};
    if (!h || h->GetEffectiveEntries() < 150) return p;
    double best = -1, at = 0.135;
    for (int i = h->FindBin(0.08); i <= h->FindBin(0.20); i++) if (h->GetBinContent(i) > best) { best = h->GetBinContent(i); at = h->GetBinCenter(i); }
    TF1 f("f", "[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])^2)", lo, hi);
    f.SetParameters(0.3 * best, 0, 0, 0.5 * best, at, 0.5 * (slo + shi));
    f.SetParLimits(3, 0, 3 * best); f.SetParLimits(4, mlo, mhi); f.SetParLimits(5, slo, shi);
    if (h->Fit(&f, "QRN0") != 0) return p;
    p.m = f.GetParameter(4); p.em = f.GetParError(4); p.s = f.GetParameter(5);
    p.ok = p.em > 0 && p.em < 0.02 && f.GetParameter(3) > 3 * f.GetParError(3) && p.m > mlo + 0.002 && p.m < mhi - 0.002 && p.s > slo * 1.05 && p.s < shi * 0.98;
    if (keep) { keep->SetRange(lo, hi); for (int i = 0; i < 6; i++) keep->SetParameter(i, f.GetParameter(i)); }
    return p;
}
// mass projection of a (mass, pair E, ring) histogram, all rings
static TH1D *Proj(TFile &f, const char *name, double e0, double e1)
{
    TH3D *h = (TH3D *)f.Get(name);
    if (!h) return 0;
    const int b0 = e0 < 0 ? 1 : h->GetYaxis()->FindBin(e0 + 1e-6), b1 = e1 < 0 ? h->GetNbinsY() : h->GetYaxis()->FindBin(e1 - 1e-6);
    return h->ProjectionX(Form("%s_%p_%g_%g", name, (void *)&f, e0, e1), b0, b1, 1, h->GetNbinsZ(), "e");
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 4) { std::fprintf(stderr, "usage: ./pi0_vs_pdg <data_all.root> <sim.root> <out.png> [-f] [-s2 <sim2.root>] [-l1 <label>] [-l2 <label>]\n"); return 1; }
    bool floor = false;
    std::string sim2, l1 = "simulation", l2 = "simulation 2";
    for (int a = 4; a < argc; a++)
    {
        if (!std::strcmp(argv[a], "-f")) floor = true;
        else if (!std::strcmp(argv[a], "-s2") && a + 1 < argc) sim2 = argv[++a];
        else if (!std::strcmp(argv[a], "-l1") && a + 1 < argc) l1 = argv[++a];
        else if (!std::strcmp(argv[a], "-l2") && a + 1 < argc) l2 = argv[++a];
    }
    const char *hn = floor ? "h3_liangf_m_E_ring_all" : "h3_liang_m_E_ring_all";
    TFile fd(argv[1]), fs(argv[2]);

    // inclusive
    TH1D *hd = Proj(fd, hn, -1, -1), *hs = Proj(fs, hn, -1, -1);
    if (!hd || !hs) { std::fprintf(stderr, "ERROR: %s missing\n", hn); return 1; }
    TF1 *fdk = new TF1("fdk", "[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])^2)", 0.05, 0.3), *fsk = new TF1("fsk", "[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])^2)", 0.05, 0.3);
    const Peak pd = FitPi0(hd, fdk), ps = FitPi0(hs, fsk);
    TFile *f2 = sim2.empty() ? 0 : new TFile(sim2.c_str());
    TH1D *h2 = f2 ? Proj(*f2, hn, -1, -1) : 0;
    TF1 *f2k = new TF1("f2k", "[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])^2)", 0.05, 0.3);
    const Peak p2 = h2 ? FitPi0(h2, f2k) : Peak{0, 0, 0, false};
    if (h2) std::printf(" %s: %.1f +- %.1f MeV (%+.1f%%)\n", l2.c_str(), 1e3 * p2.m, 1e3 * p2.em, 100 * (p2.m / kPdg - 1));
    std::printf("selection %s; PDG %.1f MeV\n all pair energies: data %.1f +- %.1f MeV (%+.1f%%) | sim %.1f +- %.1f MeV (%+.1f%%)\n", hn, 1e3 * kPdg,
                1e3 * pd.m, 1e3 * pd.em, 100 * (pd.m / kPdg - 1), 1e3 * ps.m, 1e3 * ps.em, 100 * (ps.m / kPdg - 1));

    // against pair energy
    const double eb[] = {4, 6, 8, 10, 12, 14, 16, 20, 30, 40};
    const int ne = sizeof(eb) / sizeof(eb[0]) - 1;
    TGraphErrors *gd = new TGraphErrors(), *gs = new TGraphErrors(), *g2 = new TGraphErrors();
    std::printf(" pair E (GeV)   data (MeV)        sim (MeV)\n");
    for (int i = 0; i < ne; i++)
    {
        const Peak a = FitPi0(Proj(fd, hn, eb[i], eb[i + 1])), b = FitPi0(Proj(fs, hn, eb[i], eb[i + 1]));
        const double x = 0.5 * (eb[i] + eb[i + 1]);
        if (a.ok) { const int k = gd->GetN(); gd->SetPoint(k, x - 0.2, 1e3 * a.m); gd->SetPointError(k, 0.5 * (eb[i + 1] - eb[i]), 1e3 * a.em); }
        if (b.ok) { const int k = gs->GetN(); gs->SetPoint(k, x + 0.2, 1e3 * b.m); gs->SetPointError(k, 0.5 * (eb[i + 1] - eb[i]), 1e3 * b.em); }
        if (f2)
        {
            const Peak c = FitPi0(Proj(*f2, hn, eb[i], eb[i + 1]));
            if (c.ok) { const int k = g2->GetN(); g2->SetPoint(k, x + 0.6, 1e3 * c.m); g2->SetPointError(k, 0.5 * (eb[i + 1] - eb[i]), 1e3 * c.em); }
            std::printf("            %s: %s\n", l2.c_str(), c.ok ? Form("%6.1f +- %.1f", 1e3 * c.m, 1e3 * c.em) : "-");
        }
        std::printf("  %4.0f-%-4.0f     %s          %s\n", eb[i], eb[i + 1], a.ok ? Form("%6.1f +- %.1f", 1e3 * a.m, 1e3 * a.em) : "     -     ",
                    b.ok ? Form("%6.1f +- %.1f", 1e3 * b.m, 1e3 * b.em) : "     -     ");
    }

    gStyle->SetOptStat(0);
    const int cd = TColor::GetColor("#222222"), cs = TColor::GetColor("#d9480f"), cp = TColor::GetColor("#2f9e44");
    // the two panels below a strip with the title and one shared legend
    TCanvas cv("cv", "", 1500, 720);
    TPad *pads[2] = {new TPad("pa", "", 0, 0, 0.5, 0.82), new TPad("pb", "", 0.5, 0, 1, 0.82)};
    pads[0]->Draw(); pads[1]->Draw();

    pads[0]->cd();
    gPad->SetMargin(0.12, 0.03, 0.12, 0.07);
    // normalise; the fits were made on the raw histograms, so scale their amplitude terms alike
    const double sd = 1.0 / hd->Integral(), ss = 1.0 / hs->Integral();
    hd->Scale(sd); hs->Scale(ss);
    const double s2 = h2 ? 1.0 / h2->Integral() : 1;
    if (h2) { h2->Scale(s2); for (int i = 0; i < 4; i++) f2k->SetParameter(i, f2k->GetParameter(i) * s2); h2->GetXaxis()->SetRangeUser(0, 0.32); }
    for (int i = 0; i < 4; i++) { fdk->SetParameter(i, fdk->GetParameter(i) * sd); fsk->SetParameter(i, fsk->GetParameter(i) * ss); }
    hd->GetXaxis()->SetRangeUser(0, 0.32); hs->GetXaxis()->SetRangeUser(0, 0.32);
    hd->SetTitle(";M_{#gamma#gamma} (GeV);pairs (normalised)");
    hd->SetMarkerStyle(20); hd->SetMarkerSize(0.8); hd->SetMarkerColor(cd); hd->SetLineColor(cd);
    hs->SetLineColor(cs); hs->SetLineWidth(2);
    hd->SetMaximum(1.1 * std::max(hd->GetMaximum(), hs->GetMaximum())); hd->SetMinimum(0);
    hd->GetYaxis()->SetTitleOffset(1.4);
    const int c2 = TColor::GetColor("#1c6bb0");
    if (h2) hd->SetMaximum(std::max(hd->GetMaximum(), 1.1 * h2->GetMaximum()));
    hd->Draw("E"); hs->Draw("HIST SAME");
    if (h2) { h2->SetLineColor(c2); h2->SetLineWidth(2); h2->Draw("HIST SAME"); f2k->SetLineColor(c2); f2k->SetLineStyle(2); f2k->Draw("SAME"); }
    fdk->SetLineColor(cd); fdk->SetLineStyle(2); fdk->Draw("SAME");
    fsk->SetLineColor(cs); fsk->SetLineStyle(2); fsk->Draw("SAME");
    TLine *lp = new TLine(kPdg, 0, kPdg, hd->GetMaximum()); lp->SetLineColor(cp); lp->SetLineWidth(2); lp->Draw();
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.04);
    t.DrawLatex(0.12, 0.95, "Mass spectrum, all pair energies (dashed: the peak fits)");

    pads[1]->cd();
    gPad->SetMargin(0.12, 0.03, 0.12, 0.07);
    TH1D *fr = new TH1D("fr", ";pair energy (GeV);#pi^{0} peak (MeV)", 1, 0, 40);
    fr->SetMinimum(110); fr->SetMaximum(185); fr->GetYaxis()->SetTitleOffset(1.4);
    fr->Draw();
    TLine *l2line = new TLine(0, 1e3 * kPdg, 40, 1e3 * kPdg); l2line->SetLineColor(cp); l2line->SetLineWidth(2); l2line->Draw();
    gd->SetMarkerStyle(20); gd->SetMarkerColor(cd); gd->SetLineColor(cd); gd->SetMarkerSize(1.2);
    gs->SetMarkerStyle(21); gs->SetMarkerColor(cs); gs->SetLineColor(cs); gs->SetMarkerSize(1.2);
    gd->Draw("P SAME"); gs->Draw("P SAME");
    if (f2) { g2->SetMarkerStyle(22); g2->SetMarkerColor(c2); g2->SetLineColor(c2); g2->SetMarkerSize(1.3); g2->Draw("P SAME"); }
    t.DrawLatex(0.12, 0.95, "Peak against pair energy");
    // shared legend: one line per sample, its text in the sample's colour, with the peak over all pair energies
    cv.cd();
    TLegend *lg = new TLegend(0.04, 0.83, 0.99, 0.94); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.028); lg->SetMargin(0.08); lg->SetNColumns(2);
    lg->AddEntry(gd, Form("data (run 23007048): peak %.1f MeV", 1e3 * pd.m), "pe")->SetTextColor(cd);
    lg->AddEntry(gs, Form("%s: peak %.1f MeV", l1.c_str(), 1e3 * ps.m), "pe")->SetTextColor(cs);
    if (f2) lg->AddEntry(g2, Form("%s: peak %.1f MeV", l2.c_str(), 1e3 * p2.m), "pe")->SetTextColor(c2);
    lg->AddEntry(lp, "PDG m_{#pi^{0}} = 135.0 MeV", "l")->SetTextColor(cp);
    lg->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(42); tt.SetTextSize(0.03);
    tt.DrawLatex(0.04, 0.955, Form("Calibration's #pi^{0} selection%s, data and simulation through the same analysis (peaks over all pair energies in the legend)", floor ? " (45 MeV floor)" : ""));
    cv.SaveAs(argv[3]);
    return 0;
}
