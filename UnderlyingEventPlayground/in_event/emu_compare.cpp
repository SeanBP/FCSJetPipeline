// FCS EM trigger levels free of prescales: data and simulation both put
// through the same trigger emulator at run 23007048's settings
// (Pi0Playground/trigemu/reemulate.C; friend tree "trgemu" per SimpleTree).
// Data: events the hardware recorded through fcsEM0 (prescale 221, so a
// random sample of events meeting the EM0 condition), their emulated level;
// simulation: events whose emulated EM0 fired. Per event, against the leading
// fiducial reco jet's E and pT:
//   P(EM1 emulated | EM0)  and  P(EM2 emulated | EM0)
// Also prints, for data, the hardware-bit ratios corrected by the run's
// prescales (EM1 17.9, EM2 3.32: README of Pi0Playground) as a cross-check.
//
// Build and run from this folder:
//   g++ -m64 -O2 emu_compare.cpp ../../shared/VectorDict.cxx -o emu_compare -I../../shared -I../../underlying_event -I../../jet_scale `root-config --cflags --libs` -std=c++11
//   ./emu_compare "<data pair lists>" <data emu dir> "<sim pair lists>" <sim emu dir> <out.png>

#include <TROOT.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "in_event_util.h"
#include "JetXSecWeights.h"

static const double kPs1 = 17.9, kPs2 = 3.32;  // run 23007048 prescales (RunLog.l0TriggerSet)

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}
struct Frac
{
    double w = 0, w2 = 0, p = 0;
    void Add(bool pass, double ww) { w += ww; w2 += ww * ww; if (pass) p += ww; }
    double Val() const { return w > 0 ? p / w : NAN; }
    double Err() const { const double f = Val(); return w > 0 ? std::sqrt(std::max(f * (1 - f), 1.0 / (w * w / w2)) * w2) / w : NAN; }
    double Neff() const { return w2 > 0 ? w * w / w2 : 0; }
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 6) { std::fprintf(stderr, "usage: ./emu_compare \"<data pair lists>\" <data emu dir> \"<sim pair lists>\" <sim emu dir> <out.png>\n"); return 1; }
    const std::vector<double> ed[2] = {{10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 120}, {2, 2.5, 3, 3.5, 4, 5, 6, 8, 12}};
    // [var][sample][level: EM1, EM2][bin]; hw: data hardware ratios (EM1|EM0 bits, EM2|EM1 bits)
    std::vector<Frac> fr[2][2][2], hw[2][2];
    for (int v = 0; v < 2; v++) for (int s = 0; s < 2; s++) for (int l = 0; l < 2; l++) fr[v][s][l].assign(ed[v].size() - 1, Frac());
    for (int v = 0; v < 2; v++) for (int l = 0; l < 2; l++) hw[v][l].assign(ed[v].size() - 1, Frac());

    for (int s = 0; s < 2; s++)
    {
        const std::string emudir = argv[2 + 2 * s];
        std::vector<std::pair<std::string, std::string> > files;
        const std::vector<std::string> lists = Split(argv[1 + 2 * s]);
        for (size_t l = 0; l < lists.size(); l++)
        {
            std::ifstream in(lists[l].c_str()); std::string a, b;
            while (in >> a >> b) files.push_back(std::make_pair(a, emudir + "/" + b.substr(b.find_last_of('/') + 1)));
        }
        std::map<int, double> sN, nN;
        std::vector<int> grp(files.size(), 0);
        if (s == 1)
            for (size_t f = 0; f < files.size(); f++)
            {
                const std::string base = files[f].first.substr(files[f].first.find_last_of('/') + 1);
                int a = 0, b = 0; const size_t p = base.rfind("_pt");
                if (p != std::string::npos && std::sscanf(base.c_str() + p, "_pt%d_%d.root", &a, &b) == 2) grp[f] = a;
                TFile x(files[f].first.c_str()); TTree *t = (TTree *)x.Get("jetTree");
                if (!t || t->GetEntries() == 0 || !t->GetBranch("sigma_job")) continue;
                Double_t sj; Int_t nj; t->SetBranchAddress("sigma_job", &sj); t->SetBranchAddress("N_job", &nj); t->GetEntry(0);
                const double sg = DecodePythia6Sigma(sj);
                if (sg > 0 && nj > 0) { sN[grp[f]] += sg * nj; nN[grp[f]] += nj; }
            }
        long nev = 0, nmiss = 0;
        for (size_t f = 0; f < files.size(); f++)
        {
            double w = 1;
            if (s == 1) { if (!(nN[grp[f]] > 0)) continue; w = sN[grp[f]] / nN[grp[f]] / nN[grp[f]]; }
            TFile x(files[f].first.c_str()), xe(files[f].second.c_str());
            TTree *t = (TTree *)x.Get("jetTree"), *te = xe.IsZombie() ? 0 : (TTree *)xe.Get("trgemu");
            if (!t || !te) { nmiss++; continue; }
            std::vector<float> *pt = 0, *eta = 0, *xx = 0, *yy = 0;
            Int_t e0, e1, e2, m0, m1, m2;
            Long64_t st = -1;
            const bool has_entry = t->GetBranch("simpletree_entry");
            t->SetBranchStatus("*", 0);
            const char *br[] = {"reco_pt", "reco_eta", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2"};
            for (int k = 0; k < 7; k++) t->SetBranchStatus(br[k], 1);
            t->SetBranchAddress("reco_pt", &pt); t->SetBranchAddress("reco_eta", &eta); t->SetBranchAddress("reco_x", &xx); t->SetBranchAddress("reco_y", &yy);
            t->SetBranchAddress("fcsEM0", &e0); t->SetBranchAddress("fcsEM1", &e1); t->SetBranchAddress("fcsEM2", &e2);
            if (has_entry) { t->SetBranchStatus("simpletree_entry", 1); t->SetBranchAddress("simpletree_entry", &st); }
            te->SetBranchAddress("fcsEM0", &m0); te->SetBranchAddress("fcsEM1", &m1); te->SetBranchAddress("fcsEM2", &m2);
            for (Long64_t i = 0; i < t->GetEntries(); i++)
            {
                t->GetEntry(i);
                if (te->GetEntry(has_entry ? st : i) <= 0) continue;
                if (!pt || pt->empty()) continue;
                int lead = -1;
                for (size_t j = 0; j < pt->size(); j++)
                    if (pt->at(j) >= reco_jet_pt_min_em && ProxyInFiducial(xx->at(j), yy->at(j)) && (lead < 0 || pt->at(j) > pt->at(lead))) lead = (int)j;
                if (lead < 0) continue;
                const double x_[2] = {pt->at(lead) * std::cosh(eta->at(lead)), pt->at(lead)};
                for (int v = 0; v < 2; v++)
                {
                    int b = -1;
                    for (size_t q = 0; q + 1 < ed[v].size(); q++) if (x_[v] >= ed[v][q] && x_[v] < ed[v][q + 1]) b = (int)q;
                    if (b < 0) continue;
                    // data: recorded through EM0 (hardware bit); simulation: emulated EM0
                    const bool base = s == 0 ? e0 != 0 : m0 != 0;
                    if (base) { fr[v][s][0][b].Add(m1 != 0, w); fr[v][s][1][b].Add(m2 != 0, w); }
                    if (s == 0)
                    {
                        if (e0) hw[v][0][b].Add(e1 != 0, 1);
                        if (e1) hw[v][1][b].Add(e2 != 0, 1);
                    }
                }
                nev++;
            }
        }
        std::printf("%s: %zu files (%ld without an emulation file), %ld events with a fiducial leading jet\n", s ? "simulation" : "data", files.size(), nmiss, nev);
    }

    const char *ln[2] = {"P(EM1 | EM0)", "P(EM2 | EM0)"};
    for (int v = 0; v < 2; v++)
    {
        std::printf("\nemulated level at run 23007048 settings, against leading-jet %s\n  bin           %s data | sim          %s data | sim        data hardware x prescale: EM1|EM0 x17.9  EM2|EM1 x3.32\n", v ? "pT" : "E", ln[0], ln[1]);
        for (size_t b = 0; b + 1 < ed[v].size(); b++)
            std::printf("  %5.1f-%-5.1f  %.3f+-%.3f | %.3f+-%.3f    %.3f+-%.3f | %.3f+-%.3f      %.3f  %.3f   (n data %.0f, sim n_eff %.0f)\n", ed[v][b], ed[v][b + 1],
                        fr[v][0][0][b].Val(), fr[v][0][0][b].Err(), fr[v][1][0][b].Val(), fr[v][1][0][b].Err(),
                        fr[v][0][1][b].Val(), fr[v][0][1][b].Err(), fr[v][1][1][b].Val(), fr[v][1][1][b].Err(),
                        hw[v][0][b].Val() * kPs1, hw[v][1][b].Val() * kPs2, fr[v][0][0][b].Neff(), fr[v][1][0][b].Neff());
    }

    gStyle->SetOptStat(0);
    const int cd = TColor::GetColor("#222222"), cs = TColor::GetColor("#d9480f");
    TCanvas cv("cv", "", 1500, 1000);
    cv.Divide(2, 2, 0.002, 0.002);
    std::vector<TObject *> keep;
    for (int l = 0; l < 2; l++)
        for (int v = 0; v < 2; v++)
        {
            cv.cd(l * 2 + v + 1);
            gPad->SetMargin(0.12, 0.03, 0.12, 0.09);
            if (v == 0) gPad->SetLogx();
            TH1D *h = new TH1D(Form("f%d%d", l, v), "", 1, ed[v].front(), ed[v].back()); keep.push_back(h);
            h->SetTitle(Form(";leading jet %s;%s fired (emulated)", v ? "p_{T} (GeV)" : "E (GeV)", l ? "EM2" : "EM1"));
            h->SetMinimum(0); h->SetMaximum(1.15);
            if (v == 0) h->GetXaxis()->SetMoreLogLabels();
            h->Draw();
            TLegend *lg = new TLegend(0.15, 0.76, 0.75, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04); keep.push_back(lg);
            for (int s = 0; s < 2; s++)
            {
                TGraphErrors *g = new TGraphErrors(); keep.push_back(g);
                for (size_t b = 0; b + 1 < ed[v].size(); b++)
                {
                    const Frac &f = fr[v][s][l][b];
                    if (f.Neff() < 5) continue;
                    const double xc = v ? 0.5 * (ed[v][b] + ed[v][b + 1]) : std::sqrt(ed[v][b] * ed[v][b + 1]);
                    const int i = g->GetN();
                    g->SetPoint(i, xc, f.Val()); g->SetPointError(i, 0, f.Err());
                }
                g->SetMarkerStyle(s ? 21 : 20); g->SetMarkerSize(1.3); g->SetMarkerColor(s ? cs : cd); g->SetLineColor(s ? cs : cd);
                g->Draw("PL SAME");
                lg->AddEntry(g, s ? "simulation, emulated EM0" : "data, recorded through EM0", "pl");
            }
            lg->Draw();
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
            t.DrawLatex(0.12, 0.93, Form("%s, same emulator (run 23007048 settings)", l ? "EM2" : "EM1"));
        }
    cv.SaveAs(argv[5]);
    return 0;
}
