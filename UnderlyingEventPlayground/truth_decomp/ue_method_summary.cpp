// Summary of the UE method choice (separation + veto, separation without
// veto, no separation), one page:
//   overall bias    each method's truth cone energy minus the same jets'
//                   truth UE (own-event +-90 deg cones), sigma-weighted over
//                   every sampled jet (ue_truth_decomposition output)
//   uniformity      the same bias against jet eta and against the method's
//                   band half-width, with chi2/ndf against no bias (bins with
//                   n_eff >= 20)
//   data agreement  the reco cone ratios nosep / sep and sep_veto / sep
//                   (mean proxy cone, reco_ue_cone_E), same jets, against jet
//                   pT: simulation (the decomposition's reco_cone) and data
//                   (run_sep_test.sh output, jets with all three estimates)
// Errors: linearized over the jets (ratio of sums) or the weighted spread
// (mean difference).
//
// Build and run from this folder:
//   g++ -m64 -O2 ue_method_summary.cpp ../../shared/VectorDict.cxx -o ue_method_summary -I../../shared \
//       `root-config --cflags --libs` -std=c++11
//   ./ue_method_summary <decomposition file list> <sep test outdir> <output base>

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "VectorDict.h"

static const int kNV = 3; // 0 sep_veto, 1 sep, 2 nosep
static const char *kName[kNV] = {"separation + veto", "separation, no veto (pipeline)", "no separation, no veto"};
static const char *kShort[kNV] = {"sep_veto", "sep", "nosep"};
static const char *kColor[kNV] = {"#2a78d6", "#8a5cd1", "#1baf7a"};
static const int kVMarker[kNV] = {20, 22, 33};

// Weighted sums for a vs b over jets: mean difference and ratio of sums
struct Pair
{
    double sw = 0, sw2 = 0, sd = 0, sd2 = 0, sa = 0, sb = 0, swa2 = 0, swb2 = 0, swab = 0;
    void Add(double w, double a, double b)
    {
        sw += w; sw2 += w * w; sd += w * (a - b); sd2 += w * (a - b) * (a - b);
        sa += w * a; sb += w * b; swa2 += w * w * a * a; swb2 += w * w * b * b; swab += w * w * a * b;
    }
    double Neff() const { return sw2 > 0 ? sw * sw / sw2 : 0; }
    double Diff() const { return sd / sw; }
    double DiffErr() const { const double m = Diff(); return std::sqrt(std::max(0.0, sd2 / sw - m * m) / Neff()); }
    double Ratio() const { return sa / sb; }
    double RatioErr() const { const double r = Ratio(); return std::sqrt(std::max(0.0, swa2 - 2 * r * swab + r * r * swb2)) / sb; }
};

static int Bin(const std::vector<double> &ed, double x)
{
    if (!(x >= ed.front() && x < ed.back())) return -1;
    return (int)(std::upper_bound(ed.begin(), ed.end(), x) - ed.begin()) - 1;
}

int main(int argc, char **argv)
{
    if (argc < 4) { std::cerr << "usage: ./ue_method_summary <decomposition file list> <sep test outdir> <output base>" << std::endl; return 1; }
    const std::vector<double> ed_eta = {2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.6};
    const std::vector<double> ed_w = {0, 10, 15, 20, 25, 30, 35, 40, 46};
    const std::vector<double> ed_pt = {2, 2.5, 3, 3.5, 4, 5, 6, 8, 12};

    // ---- simulation: bias vs truth UE, and reco cone ratios
    Pair all[kNV];
    std::vector<Pair> by_eta[kNV], by_w[kNV], sim_ratio[2]; // sim_ratio: 0 nosep/sep, 1 sep_veto/sep
    for (int v = 0; v < kNV; v++) { by_eta[v].resize(ed_eta.size() - 1); by_w[v].resize(ed_w.size() - 1); }
    for (int r = 0; r < 2; r++) sim_ratio[r].resize(ed_pt.size() - 1);
    {
        std::ifstream in(argv[1]);
        std::string fn;
        // reco cone per jet and variant, to form the ratios jet by jet: the
        // two builds sample the same jets in the same order (jet_index)
        while (in >> fn)
        {
            TFile f(fn.c_str());
            TTree *ti = (TTree *)f.Get("info"), *tj = (TTree *)f.Get("jets");
            if (!ti || !tj) { std::cerr << "skip " << fn << std::endl; continue; }
            Double_t sigma, trials; Long64_t n_est, n_samp;
            ti->SetBranchAddress("sigma_mb", &sigma); ti->SetBranchAddress("trials", &trials);
            ti->SetBranchAddress("n_estimable", &n_est); ti->SetBranchAddress("n_sampled", &n_samp);
            ti->GetEntry(0);
            if (!(sigma > 0) || !(trials > 0) || n_samp <= 0) continue;
            const double w = sigma / trials * (double)n_est / n_samp;
            Int_t variant; Float_t pt, eta, hw, ta, oa, rc;
            tj->SetBranchAddress("variant", &variant); tj->SetBranchAddress("jet_pt", &pt); tj->SetBranchAddress("jet_eta", &eta);
            tj->SetBranchAddress("half_width", &hw); tj->SetBranchAddress("t_all", &ta); tj->SetBranchAddress("own_all", &oa);
            tj->SetBranchAddress("reco_cone", &rc);
            for (Long64_t e = 0; e < tj->GetEntries(); e++)
            {
                tj->GetEntry(e);
                if (variant < 0 || variant >= kNV) continue;
                all[variant].Add(w, ta, oa);
                int b;
                if ((b = Bin(ed_eta, eta)) >= 0) by_eta[variant][b].Add(w, ta, oa);
                if ((b = Bin(ed_w, hw)) >= 0) by_w[variant][b].Add(w, ta, oa);
            }
        }
    }
    // sim reco cone ratios need the variants of one jet together: the sep
    // build writes sep_veto and sep rows of a jet back to back, the nosep
    // build its own file; pair them by file and row order
    {
        std::ifstream in(argv[1]);
        std::string fn;
        std::vector<std::string> sep_files, nosep_files;
        while (in >> fn) (fn.find("td_nosep_") != std::string::npos ? nosep_files : sep_files).push_back(fn);
        for (size_t i = 0; i < sep_files.size(); i++)
        {
            std::string partner = sep_files[i];
            partner.replace(partner.find("td_sep_"), 7, "td_nosep_");
            if (std::find(nosep_files.begin(), nosep_files.end(), partner) == nosep_files.end()) continue;
            TFile fs(sep_files[i].c_str()), fn2(partner.c_str());
            TTree *ts = (TTree *)fs.Get("jets"), *tn = (TTree *)fn2.Get("jets"), *ti = (TTree *)fs.Get("info");
            if (!ts || !tn || !ti) continue;
            Double_t sigma, trials; Long64_t n_est, n_samp;
            ti->SetBranchAddress("sigma_mb", &sigma); ti->SetBranchAddress("trials", &trials);
            ti->SetBranchAddress("n_estimable", &n_est); ti->SetBranchAddress("n_sampled", &n_samp);
            ti->GetEntry(0);
            const double w = sigma / trials * (double)n_est / n_samp;
            Int_t vs, vn; Float_t pts, ptn, rcs, rcn;
            ts->SetBranchAddress("variant", &vs); ts->SetBranchAddress("jet_pt", &pts); ts->SetBranchAddress("reco_cone", &rcs);
            tn->SetBranchAddress("variant", &vn); tn->SetBranchAddress("jet_pt", &ptn); tn->SetBranchAddress("reco_cone", &rcn);
            if (ts->GetEntries() != 2 * tn->GetEntries()) { std::cerr << "row counts differ: " << sep_files[i] << std::endl; continue; }
            for (Long64_t j = 0; j < tn->GetEntries(); j++)
            {
                ts->GetEntry(2 * j); const double c_veto = rcs, pt_s = pts; const int v0 = vs;
                ts->GetEntry(2 * j + 1); const double c_sep = rcs; const int v1 = vs;
                tn->GetEntry(j);
                if (v0 != 0 || v1 != 1 || std::fabs(pt_s - ptn) > 1e-4) { std::cerr << "rows not paired in " << sep_files[i] << std::endl; break; }
                const int b = Bin(ed_pt, ptn);
                if (b < 0) continue;
                sim_ratio[0][b].Add(w, rcn, c_sep);
                sim_ratio[1][b].Add(w, c_veto, c_sep);
            }
        }
    }

    // ---- data: reco cone ratios, same jets
    std::vector<Pair> data_ratio[2];
    for (int r = 0; r < 2; r++) data_ratio[r].resize(ed_pt.size() - 1);
    {
        const std::string dir = argv[2];
        std::vector<std::string> names;
        TSystemDirectory d("d", (dir + "/sep_veto").c_str());
        TIter next(d.GetListOfFiles());
        TSystemFile *sf;
        while ((sf = (TSystemFile *)next()))
        {
            const std::string n = sf->GetName();
            if (n.find(".root") != std::string::npos) names.push_back(n);
        }
        std::sort(names.begin(), names.end());
        for (size_t i = 0; i < names.size(); i++)
        {
            TFile *tf[kNV]; TTree *t[kNV];
            bool ok = true;
            for (int v = 0; v < kNV; v++)
            {
                tf[v] = TFile::Open((dir + "/" + kShort[v] + "/" + names[i]).c_str());
                t[v] = tf[v] ? (TTree *)tf[v]->Get("jetTree") : 0;
                ok = ok && t[v] && t[v]->GetEntries() == (t[0] ? t[0]->GetEntries() : -1);
            }
            if (!ok) { for (int v = 0; v < kNV; v++) delete tf[v]; continue; }
            std::vector<float> *pt = 0, *ue[kNV] = {0, 0, 0};
            std::vector<int> *np[kNV] = {0, 0, 0};
            t[0]->SetBranchAddress("reco_pt", &pt);
            for (int v = 0; v < kNV; v++) { t[v]->SetBranchAddress("reco_ue_cone_E", &ue[v]); t[v]->SetBranchAddress("reco_ue_nproxy", &np[v]); }
            for (Long64_t e = 0; e < t[0]->GetEntries(); e++)
            {
                for (int v = 0; v < kNV; v++) t[v]->GetEntry(e);
                for (size_t j = 0; j < pt->size(); j++)
                {
                    bool good = true;
                    for (int v = 0; v < kNV; v++) good = good && j < ue[v]->size() && np[v]->at(j) > 0 && ue[v]->at(j) >= 0;
                    const int b = Bin(ed_pt, pt->at(j));
                    if (!good || b < 0) continue;
                    data_ratio[0][b].Add(1, ue[2]->at(j), ue[1]->at(j));
                    data_ratio[1][b].Add(1, ue[0]->at(j), ue[1]->at(j));
                }
            }
            for (int v = 0; v < kNV; v++) delete tf[v];
        }
    }

    // ---- printout
    printf("ue_method_summary\n  overall bias vs truth UE (GeV) | ratio:\n");
    for (int v = 0; v < kNV; v++)
        printf("    %-9s %+.3f +- %.3f | %.3f +- %.3f   (n_eff %.0f)\n", kShort[v], all[v].Diff(), all[v].DiffErr(), all[v].Ratio(), all[v].RatioErr(), all[v].Neff());
    double chi2[2][kNV]; int ndf[2][kNV];
    for (int v = 0; v < kNV; v++)
        for (int k = 0; k < 2; k++)
        {
            const std::vector<Pair> &bins = k == 0 ? by_eta[v] : by_w[v];
            chi2[k][v] = 0; ndf[k][v] = 0;
            for (size_t b = 0; b < bins.size(); b++)
                if (bins[b].Neff() >= 20) { chi2[k][v] += std::pow(bins[b].Diff() / bins[b].DiffErr(), 2); ndf[k][v]++; }
        }
    printf("  chi2/ndf vs no bias:   eta | band half-width\n");
    for (int v = 0; v < kNV; v++) printf("    %-9s %5.1f/%d | %5.1f/%d\n", kShort[v], chi2[0][v], ndf[0][v], chi2[1][v], ndf[1][v]);
    printf("  reco cone ratios vs jet pT: data | sim\n    pT          nosep/sep               sep_veto/sep\n");
    for (size_t b = 0; b + 1 < ed_pt.size(); b++)
        printf("    %4.1f-%-4.1f  %.3f+-%.3f | %.3f+-%.3f   %.3f+-%.3f | %.3f+-%.3f\n", ed_pt[b], ed_pt[b + 1],
               data_ratio[0][b].Ratio(), data_ratio[0][b].RatioErr(), sim_ratio[0][b].Ratio(), sim_ratio[0][b].RatioErr(),
               data_ratio[1][b].Ratio(), data_ratio[1][b].RatioErr(), sim_ratio[1][b].Ratio(), sim_ratio[1][b].RatioErr());

    // ---- page
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(5);
    const int ink = TColor::GetColor("#3a3a36"), muted = TColor::GetColor("#8a8a84");
    std::vector<TObject *> keep;
    TCanvas c("c", "", 2100, 1150);
    c.Divide(3, 2, 0.004, 0.004);
    auto pad = [&](int i) {
        c.cd(i);
        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.09); gPad->SetRightMargin(0.03);
    };
    auto frame = [&](const char *name, double x0, double x1, double y0, double y1, const char *xt, const char *yt) {
        TH1D *fr = new TH1D(name, "", 1, x0, x1);
        keep.push_back(fr);
        fr->SetMinimum(y0); fr->SetMaximum(y1);
        fr->GetXaxis()->SetTitle(xt); fr->GetYaxis()->SetTitle(yt);
        fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleOffset(1.4);
        fr->Draw("AXIS");
        return fr;
    };
    auto hline = [&](double x0, double x1, double y) {
        TLine *l = new TLine(x0, y, x1, y);
        keep.push_back(l);
        l->SetLineColor(ink); l->SetLineStyle(2); l->Draw();
    };
    auto title = [&](const char *t) {
        TLatex *tl = new TLatex();
        keep.push_back(tl);
        tl->SetNDC(); tl->SetTextSize(0.05); tl->SetTextColor(ink); tl->DrawLatex(0.15, 0.94, t);
    };
    auto style = [&](TGraphErrors *g, int v, bool open) {
        const int col = v < 0 ? ink : TColor::GetColor(kColor[v]);
        g->SetMarkerStyle(v < 0 ? 20 : (open ? (kVMarker[v] == 33 ? 27 : kVMarker[v] + 4) : kVMarker[v]));
        g->SetMarkerSize(v >= 0 && kVMarker[v] == 33 ? 1.8 : 1.3);
        g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2);
    };

    // (1) overall bias
    pad(1);
    TH1D *f1 = frame("f1", 0.5, 3.5, -0.35, 0.2, "", "method #minus truth UE (GeV)");
    f1->GetXaxis()->SetLabelSize(0); f1->GetXaxis()->SetTickLength(0);
    f1->Draw("AXIS");
    hline(0.5, 3.5, 0);
    title("overall bias (all jets, #sigma-weighted)");
    for (int v = 0; v < kNV; v++)
    {
        TGraphErrors *g = new TGraphErrors();
        keep.push_back(g);
        g->SetPoint(0, v + 1, all[v].Diff()); g->SetPointError(0, 0, all[v].DiffErr());
        style(g, v, false);
        g->Draw("P same");
        TLatex *tl = new TLatex();
        keep.push_back(tl);
        tl->SetTextSize(0.042); tl->SetTextColor(ink); tl->SetTextAlign(21);
        tl->DrawLatex(v + 1, -0.33, Form("%s", kShort[v]));
        tl->DrawLatex(v + 1, all[v].Diff() + all[v].DiffErr() + 0.02, Form("%.2f #pm %.2f", all[v].Ratio(), all[v].RatioErr()));
    }
    {
        TLatex *tl = new TLatex();
        keep.push_back(tl);
        tl->SetNDC(); tl->SetTextSize(0.038); tl->SetTextColor(muted);
        tl->DrawLatex(0.18, 0.84, "labels: method / truth UE");
    }

    // (2), (3) uniformity: bias vs eta and vs band half-width
    for (int k = 0; k < 2; k++)
    {
        pad(2 + k);
        const std::vector<double> &ed = k == 0 ? ed_eta : ed_w;
        frame(Form("fu%d", k), ed.front(), ed.back(), -0.8, 0.8, k == 0 ? "jet #eta" : "band half-width of the method at the jet (deg)", "method #minus truth UE (GeV)");
        hline(ed.front(), ed.back(), 0);
        title(k == 0 ? "uniformity: bias vs jet #eta" : "uniformity: bias vs band width");
        TLegend *leg = new TLegend(0.18, 0.7, 0.97, 0.89);
        keep.push_back(leg);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04); leg->SetTextColor(ink);
        for (int v = 0; v < kNV; v++)
        {
            const std::vector<Pair> &bins = k == 0 ? by_eta[v] : by_w[v];
            TGraphErrors *g = new TGraphErrors();
            keep.push_back(g);
            for (size_t b = 0; b < bins.size(); b++)
            {
                // too few jets for a stable mean (the narrow edge bins of the separation variants)
                if (bins[b].Neff() < 20) continue;
                const int n = g->GetN();
                g->SetPoint(n, 0.5 * (ed[b] + ed[b + 1]) + (v - 1) * 0.012 * (ed.back() - ed.front()), bins[b].Diff());
                g->SetPointError(n, 0, bins[b].DiffErr());
            }
            style(g, v, false);
            g->Draw("PL same");
            leg->AddEntry(g, Form("%s   #chi^{2}/ndf %.1f/%d", kShort[v], chi2[k][v], ndf[k][v]), "lp");
        }
        leg->Draw();
    }

    // (4), (5) agreement with data: reco cone ratios vs jet pT
    for (int r = 0; r < 2; r++)
    {
        pad(4 + r);
        frame(Form("fd%d", r), ed_pt.front(), ed_pt.back(), 0.6, 1.5, "jet p_{T} (GeV)", r == 0 ? "reco cone: nosep / sep" : "reco cone: sep_veto / sep");
        hline(ed_pt.front(), ed_pt.back(), 1);
        title(r == 0 ? "data check: no separation vs separation" : "data check: veto vs no veto");
        TLegend *leg = new TLegend(0.18, 0.72, 0.97, 0.89);
        keep.push_back(leg);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042); leg->SetTextColor(ink);
        for (int s = 0; s < 2; s++)
        {
            const std::vector<Pair> &bins = s == 0 ? data_ratio[r] : sim_ratio[r];
            TGraphErrors *g = new TGraphErrors();
            keep.push_back(g);
            for (size_t b = 0; b < bins.size(); b++)
            {
                if (bins[b].Neff() < 20) continue;
                const int n = g->GetN();
                g->SetPoint(n, 0.5 * (ed_pt[b] + ed_pt[b + 1]) + (s - 0.5) * 0.15, bins[b].Ratio());
                g->SetPointError(n, 0, bins[b].RatioErr());
            }
            style(g, s == 0 ? -1 : (r == 0 ? 2 : 0), s == 1);
            g->Draw("PL same");
            leg->AddEntry(g, s == 0 ? "data, run 23007048 (reclustered)" : "simulation (unfiltered, #sigma-weighted)", "lp");
        }
        leg->Draw();
    }

    // (6) legend and reading
    c.cd(6);
    TLegend *leg = new TLegend(0.03, 0.55, 0.97, 0.9);
    keep.push_back(leg);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.05); leg->SetTextColor(ink);
    leg->SetHeader("proxy UE methods (truth UE: own event, #pm90#circ cones)");
    for (int v = 0; v < kNV; v++)
    {
        TGraphErrors *d = new TGraphErrors();
        keep.push_back(d);
        style(d, v, false);
        leg->AddEntry(d, kName[v], "lp");
    }
    leg->Draw();
    TLatex *tl = new TLatex();
    keep.push_back(tl);
    tl->SetNDC(); tl->SetTextSize(0.042); tl->SetTextColor(ink);
    tl->DrawLatex(0.05, 0.45, "sep and nosep tie on the overall bias (#pm3%, opposite signs);");
    tl->DrawLatex(0.05, 0.38, "sep is more uniform in band width, where nosep grows;");
    tl->DrawLatex(0.05, 0.31, "data vs sim: veto/no-veto agrees within ~0.05 at 2.5-5 GeV;");
    tl->DrawLatex(0.05, 0.24, "nosep/sep is 0.1-0.3 higher in sim than in data.");
    tl->DrawLatex(0.05, 0.17, "sep + veto: ~10% low everywhere.");

    const std::string out = argv[3];
    c.SaveAs((out + ".png").c_str());
    c.SaveAs((out + ".pdf").c_str());
    return 0;
}
