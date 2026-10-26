// The separation test on FCS data (run_sep_test.sh output): the UE step's
// estimate (reco_ue_E, UE energy x jet area) for the same jets three ways,
//   sep_veto  the pipeline method: separation 2.5R and the jet veto
//   sep       separation, no veto
//   nosep     no separation (and so no veto)
// Jets enter if all three have an estimate. Pages:
//   <base>_grid     per-jet distributions in jet-pT ranges (log y), with each
//                   variant's most probable value (Gaussian over the peak, the
//                   spike at 0 ignored; the mean where the fit fails, dotted)
//   <base>_profiles the MPV (open: mean where the fit fails) against jet pT,
//                   eta and energy, each over nosep (no separation, no veto) underneath
//   <base>_means    the same with the means
// The printout adds the proxy events per jet.
//
// Build and run from this folder:
//   g++ -m64 -O2 sep_test_compare.cpp ../../shared/VectorDict.cxx -o sep_test_compare -I../../shared `root-config --cflags --libs` -std=c++11
//   ./sep_test_compare <outdir of submit_sep_test.sh> [output_base]

#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

static const int kNV = 3;
static const char *kDir[kNV] = {"sep_veto", "sep", "nosep"};
static const char *kName[kNV] = {"separation + veto (pipeline)", "separation, no veto", "no separation, no veto"};
static const int kRef = 2; // ratios are to nosep
static const char *kColor[kNV] = {"#2a78d6", "#8a5cd1", "#1baf7a"};
static const int kVarMarker[kNV] = {20, 22, 33}, kVarOpen[kNV] = {24, 26, 27};
static const double kXMax = 6.0;
static const int kNBins = 120;

struct Fit { double mpv = -1, err = 0; bool ok = false; };
// Most probable value, ignoring the first bin: Gaussian over the FWHM, then over the fitted mean +- 1 sigma
static Fit FitPeak(TH1D *h, const char *name)
{
    Fit r;
    int bmax = 2;
    for (int b = 2; b <= h->GetNbinsX(); b++) if (h->GetBinContent(b) > h->GetBinContent(bmax)) bmax = b;
    const double half = 0.5 * h->GetBinContent(bmax);
    if (half <= 0) return r;
    int bl = bmax, br = bmax;
    while (bl > 2 && h->GetBinContent(bl - 1) > half) bl--;
    while (br < h->GetNbinsX() && h->GetBinContent(br + 1) > half) br++;
    double lo = h->GetBinLowEdge(bl), hi = h->GetBinLowEdge(br + 1);
    if (br - bl < 3) { lo = h->GetBinCenter(bmax) - 3 * h->GetBinWidth(1); hi = h->GetBinCenter(bmax) + 3 * h->GetBinWidth(1); }
    TF1 *g = new TF1(name, "gaus", lo, hi);
    g->SetParameters(h->GetBinContent(bmax), h->GetBinCenter(bmax), 0.5 * (hi - lo));
    if (h->Fit(g, "QNR0") != 0) return r;
    const double m = g->GetParameter(1), s = std::fabs(g->GetParameter(2));
    g->SetRange(std::max(h->GetBinLowEdge(2), m - s), m + s);
    if (h->Fit(g, "QNR0") != 0) return r;
    r.mpv = g->GetParameter(1); r.err = g->GetParError(1);
    const double sg = std::fabs(g->GetParameter(2));
    r.ok = r.mpv > h->GetBinLowEdge(2) && r.mpv < h->GetXaxis()->GetXmax() && sg > 0 && r.err < sg;
    return r;
}

struct Acc
{
    std::vector<TH1D *> h[kNV];
    std::vector<double> s[kNV], s2[kNV], np[kNV];
    std::vector<long> n;
    void Init(int nb, const char *tag)
    {
        n.assign(nb, 0);
        for (int v = 0; v < kNV; v++)
        {
            s[v].assign(nb, 0); s2[v].assign(nb, 0); np[v].assign(nb, 0); h[v].resize(nb);
            for (int b = 0; b < nb; b++) h[v][b] = new TH1D(Form("h_%s_%d_%d", tag, v, b), ";UE energy #times jet area (GeV);jets", kNBins, 0, kXMax);
        }
    }
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2) { std::cerr << "usage: ./sep_test_compare <outdir> [output_base]" << std::endl; return 1; }
    const std::string dir = argv[1], out = argc > 2 ? argv[2] : "sep_test";

    // File names from the first variant
    std::vector<std::string> names;
    TSystemDirectory d("d", (dir + "/" + kDir[0]).c_str());
    TList *files = d.GetListOfFiles();
    for (int i = 0; files && i < files->GetSize(); i++)
    {
        const std::string nm = files->At(i)->GetName();
        if (nm.find("jet_output_") == 0 && nm.find(".root") != std::string::npos) names.push_back(nm);
    }
    std::sort(names.begin(), names.end());

    const int kNVar = 3;
    const char *vtitle[kNVar] = {"jet p_{T} (GeV)", "jet #eta", "jet E (GeV)"};
    std::vector<double> edges[kNVar];
    edges[0] = {1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8};
    for (int i = 0; i <= 14; i++) edges[1].push_back(2.2 + 0.1 * i);
    edges[2] = {10, 15, 20, 25, 30, 40, 50, 60, 80, 120};
    Acc acc[kNVar];
    for (int k = 0; k < kNVar; k++) acc[k].Init((int)edges[k].size() - 1, Form("v%d", k));

    long n_jets = 0, n_all = 0;
    for (size_t i = 0; i < names.size(); i++)
    {
        TFile *tf[kNV]; TTree *t[kNV];
        bool ok = true;
        for (int v = 0; v < kNV; v++)
        {
            tf[v] = TFile::Open((dir + "/" + kDir[v] + "/" + names[i]).c_str());
            t[v] = tf[v] ? (TTree *)tf[v]->Get("jetTree") : 0;
            if (!t[v]) ok = false;
        }
        for (int v = 1; v < kNV && ok; v++) ok = t[v]->GetEntries() == t[0]->GetEntries();
        if (!ok) { std::cerr << "skip " << names[i] << std::endl; for (int v = 0; v < kNV; v++) delete tf[v]; continue; }
        std::vector<float> *pt = 0, *eta = 0, *E = 0, *ue[kNV] = {0};
        std::vector<int> *np[kNV] = {0};
        t[0]->SetBranchAddress("reco_pt", &pt); t[0]->SetBranchAddress("reco_eta", &eta); t[0]->SetBranchAddress("reco_E", &E);
        for (int v = 0; v < kNV; v++) { t[v]->SetBranchAddress("reco_ue_E", &ue[v]); t[v]->SetBranchAddress("reco_ue_nproxy", &np[v]); }
        for (Long64_t e = 0; e < t[0]->GetEntries(); e++)
        {
            for (int v = 0; v < kNV; v++) t[v]->GetEntry(e);
            for (size_t j = 0; j < pt->size(); j++)
            {
                n_all++;
                bool good = true;
                for (int v = 0; v < kNV; v++) good = good && j < ue[v]->size() && np[v]->at(j) > 0 && ue[v]->at(j) >= 0;
                if (!good) continue;
                n_jets++;
                const double x[kNVar] = {pt->at(j), eta->at(j), E->at(j)};
                for (int k = 0; k < kNVar; k++)
                {
                    const std::vector<double> &ed = edges[k];
                    if (x[k] < ed.front() || x[k] >= ed.back()) continue;
                    const int b = (int)(std::upper_bound(ed.begin(), ed.end(), x[k]) - ed.begin()) - 1;
                    acc[k].n[b]++;
                    for (int v = 0; v < kNV; v++)
                    {
                        const double u = ue[v]->at(j);
                        acc[k].h[v][b]->Fill(u); acc[k].s[v][b] += u; acc[k].s2[v][b] += u * u; acc[k].np[v][b] += np[v]->at(j);
                    }
                }
            }
        }
        for (int v = 0; v < kNV; v++) delete tf[v];
    }
    printf("sep_test_compare: %zu files, %ld reco jets, %ld with all three estimates\n", names.size(), n_all, n_jets);

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(7);
    const int c_ink = TColor::GetColor("#3a3a36");
    std::vector<TObject *> keep;
    // Values per variable, bin and variant
    std::vector<double> val[kNVar][kNV], err[kNVar][kNV], mean[kNVar][kNV], merr[kNVar][kNV];
    std::vector<int> fitted[kNVar][kNV];
    for (int k = 0; k < kNVar; k++)
    {
        const int nb = (int)edges[k].size() - 1;
        printf("  %s\n    bin          jets   ", vtitle[k]);
        for (int v = 0; v < kNV; v++) printf("  MPV %-9.9s     ", kDir[v]);
        for (int v = 0; v < kNV; v++) printf(" mean %-8.8s", kDir[v]);
        printf("  proxies per jet:");
        for (int v = 0; v < kNV; v++) printf(" %s", kDir[v]);
        printf("\n");
        for (int v = 0; v < kNV; v++) { val[k][v].assign(nb, -1); err[k][v].assign(nb, 0); mean[k][v].assign(nb, -1); merr[k][v].assign(nb, 0); fitted[k][v].assign(nb, 0); }
        for (int b = 0; b < nb; b++)
        {
            const long N = acc[k].n[b];
            if (N < 50) continue;
            const int rebin = N < 20000 ? 4 : (N < 200000 ? 2 : 1);
            for (int v = 0; v < kNV; v++)
            {
                mean[k][v][b] = acc[k].s[v][b] / N;
                merr[k][v][b] = std::sqrt(std::max(0.0, acc[k].s2[v][b] / N - mean[k][v][b] * mean[k][v][b]) / (N - 1));
                if (rebin > 1) acc[k].h[v][b]->Rebin(rebin);
                const Fit f = FitPeak(acc[k].h[v][b], Form("g%d_%d_%d", k, v, b));
                if (f.ok) { val[k][v][b] = f.mpv; err[k][v][b] = f.err; fitted[k][v][b] = 1; }
                else { val[k][v][b] = mean[k][v][b]; err[k][v][b] = merr[k][v][b]; }
            }
            printf("    %5.2f-%-6.2f %7ld", edges[k][b], edges[k][b + 1], N);
            for (int v = 0; v < kNV; v++) printf("  %s", fitted[k][v][b] ? Form("%.3f+-%.3f    ", val[k][v][b], err[k][v][b]) : Form("mean %.3f      ", val[k][v][b]));
            for (int v = 0; v < kNV; v++) printf("  %.3f    ", mean[k][v][b]);
            for (int v = 0; v < kNV; v++) printf(" %.0f", acc[k].np[v][b] / N);
            printf("\n");
        }
    }

    // Profiles (MPV, then means)
    for (int page = 0; page < 2; page++)
    {
        TCanvas *c = new TCanvas(Form("c%d", page), "", 1900, 1000);
        keep.push_back(c);
        c->Divide(kNVar, 1, 0.003, 0.003);
        for (int k = 0; k < kNVar; k++)
        {
            const int nb = (int)edges[k].size() - 1;
            c->cd(k + 1);
            TPad *top = new TPad(Form("t%d%d", page, k), "", 0, 0.34, 1, 1), *bot = new TPad(Form("b%d%d", page, k), "", 0, 0, 1, 0.34);
            top->SetBottomMargin(0.02); top->SetLeftMargin(0.15); top->SetRightMargin(0.03); top->SetTopMargin(0.08);
            bot->SetTopMargin(0.03); bot->SetBottomMargin(0.3); bot->SetLeftMargin(0.15); bot->SetRightMargin(0.03);
            top->Draw(); bot->Draw(); keep.push_back(top); keep.push_back(bot);
            TGraphErrors *g[kNV], *gf[kNV], *r[kNV], *rf[kNV];
            double ymax = 0;
            for (int v = 0; v < kNV; v++)
            {
                g[v] = new TGraphErrors(); gf[v] = new TGraphErrors(); r[v] = new TGraphErrors(); rf[v] = new TGraphErrors();
                for (int b = 0; b < nb; b++)
                {
                    const double y = page == 0 ? val[k][v][b] : mean[k][v][b], ey = page == 0 ? err[k][v][b] : merr[k][v][b];
                    const double y0 = page == 0 ? val[k][kRef][b] : mean[k][kRef][b], ey0 = page == 0 ? err[k][kRef][b] : merr[k][kRef][b];
                    if (!(y > 0)) continue;
                    const bool open = page == 0 && !fitted[k][v][b];
                    const double x = 0.5 * (edges[k][b] + edges[k][b + 1]), ex = 0.5 * (edges[k][b + 1] - edges[k][b]);
                    TGraphErrors *gg = open ? gf[v] : g[v];
                    int n = gg->GetN(); gg->SetPoint(n, x, y); gg->SetPointError(n, ex, ey);
                    ymax = std::max(ymax, y + ey);
                    if (v == kRef || !(y0 > 0)) continue;
                    TGraphErrors *rr = open ? rf[v] : r[v];
                    const double q = y / y0;
                    n = rr->GetN(); rr->SetPoint(n, x, q); rr->SetPointError(n, 0, q * std::sqrt(std::pow(ey / y, 2) + std::pow(ey0 / y0, 2)));
                }
            }
            top->cd();
            TH1D *fr = new TH1D(Form("fr%d%d", page, k), "", 1, edges[k].front(), edges[k].back());
            fr->SetMinimum(0); fr->SetMaximum(1.6 * ymax);
            fr->GetYaxis()->SetTitle(page == 0 ? "UE energy #times jet area, MPV (GeV)" : "UE energy #times jet area, mean (GeV)");
            fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetLabelSize(0.05); fr->GetXaxis()->SetLabelSize(0);
            fr->Draw("AXIS"); keep.push_back(fr);
            for (int v = 0; v < kNV; v++)
            {
                const int col = TColor::GetColor(kColor[v]);
                TGraphErrors *gs[2] = {g[v], gf[v]};
                for (int o = 0; o < 2; o++)
                {
                    gs[o]->SetMarkerStyle(o ? kVarOpen[v] : kVarMarker[v]); gs[o]->SetMarkerSize(kVarMarker[v] == 33 ? 1.6 : 1.2);
                    gs[o]->SetMarkerColor(col); gs[o]->SetLineColor(col); gs[o]->SetLineWidth(2);
                    if (gs[o]->GetN() > 0) gs[o]->Draw("P same");
                    keep.push_back(gs[o]);
                }
            }
            if (k == 0)
            {
                TLegend *leg = new TLegend(0.17, 0.68, 0.97, 0.91);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045); leg->SetTextColor(c_ink);
                for (int v = 0; v < kNV; v++) leg->AddEntry(g[v], kName[v], "lep");
                if (page == 0) leg->AddEntry((TObject *)0, "open: mean (peak fit failed)", "");
                leg->Draw(); keep.push_back(leg);
            }
            bot->cd();
            TH1D *fb = new TH1D(Form("fb%d%d", page, k), "", 1, edges[k].front(), edges[k].back());
            fb->SetMinimum(0.6); fb->SetMaximum(1.4);
            fb->GetXaxis()->SetTitle(vtitle[k]); fb->GetYaxis()->SetTitle("ratio to no sep., no veto");
            fb->GetXaxis()->SetTitleSize(0.11); fb->GetXaxis()->SetLabelSize(0.1); fb->GetYaxis()->SetTitleSize(0.085);
            fb->GetYaxis()->SetLabelSize(0.09); fb->GetYaxis()->SetTitleOffset(0.8); fb->GetYaxis()->SetNdivisions(505);
            fb->Draw("AXIS"); keep.push_back(fb);
            TLine *one = new TLine(edges[k].front(), 1, edges[k].back(), 1); one->SetLineStyle(2); one->SetLineColor(c_ink); one->Draw(); keep.push_back(one);
            for (int v = 0; v < kNV; v++)
            {
                if (v == kRef) continue;
                const int col = TColor::GetColor(kColor[v]);
                TGraphErrors *rs[2] = {r[v], rf[v]};
                for (int o = 0; o < 2; o++)
                {
                    rs[o]->SetMarkerStyle(o ? kVarOpen[v] : kVarMarker[v]); rs[o]->SetMarkerSize(kVarMarker[v] == 33 ? 1.6 : 1.2);
                    rs[o]->SetMarkerColor(col); rs[o]->SetLineColor(col); rs[o]->SetLineWidth(2);
                    if (rs[o]->GetN() > 0) rs[o]->Draw("P same");
                    keep.push_back(rs[o]);
                }
            }
        }
        c->SaveAs((out + (page == 0 ? "_profiles" : "_means") + ".png").c_str());
        c->SaveAs((out + (page == 0 ? "_profiles" : "_means") + ".pdf").c_str());
    }

    // Distribution grid in jet-pT ranges
    const int np_ = (int)edges[0].size() - 1;
    TCanvas cg("cg", "", 2400, 2100);
    cg.Divide(3, 3, 0.002, 0.002);
    for (int b = 0; b < np_ && b < 9; b++)
    {
        cg.cd(b + 1);
        gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09);
        double ymax = 0;
        for (int v = 0; v < kNV; v++) ymax = std::max(ymax, acc[0].h[v][b]->GetMaximum());
        for (int v = 0; v < kNV; v++)
        {
            TH1D *h = acc[0].h[v][b];
            h->SetLineColor(TColor::GetColor(kColor[v])); h->SetLineWidth(2);
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05); h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
            if (v == 0) { h->SetMaximum(40 * std::max(ymax, 1.0)); h->SetMinimum(0.5); h->Draw("hist"); }
            else h->Draw("hist same");
            if (val[0][v][b] > 0)
            {
                TLine *l = new TLine(val[0][v][b], 0.5, val[0][v][b], 3 * ymax);
                l->SetLineColor(TColor::GetColor(kColor[v])); l->SetLineWidth(2); l->SetLineStyle(fitted[0][v][b] ? 1 : 3);
                l->Draw(); keep.push_back(l);
            }
        }
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.055); tl.SetTextColor(c_ink); tl.SetTextFont(62);
        tl.DrawLatex(0.15, 0.93, Form("%.1f < jet p_{T} < %.1f GeV  (%ld jets)", edges[0][b], edges[0][b + 1], acc[0].n[b]));
        TLegend *leg = new TLegend(0.40, 0.70, 0.97, 0.89);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04); leg->SetTextColor(c_ink);
        for (int v = 0; v < kNV; v++)
            leg->AddEntry(acc[0].h[v][b], Form("%s: %s %.2f", kName[v], fitted[0][v][b] ? "MPV" : "mean", val[0][v][b]), "l");
        leg->Draw(); keep.push_back(leg);
    }
    cg.SaveAs((out + "_grid.png").c_str());
    cg.SaveAs((out + "_grid.pdf").c_str());
    return 0;
}
