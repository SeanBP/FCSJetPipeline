// How differently the FCS EM trigger levels respond in data and simulation:
// per event, against the leading reco jet's E and pT (pT >= 2 GeV, inside the
// reco fiducial region),
//   P(EM2 | EM1 fired)  and  P(EM1 | EM0 fired)
// data: fired trigger bits; simulation: StFcsTriggerSimMaker's emulated bits,
// weighted per ptHat bin (sigma / N). P(EM2 | EM1) is free of prescales if EM2
// is not prescaled (an event that recorded EM1 shows EM2 whenever its
// condition held); P(EM1 | EM0) in data carries EM1's prescale as an overall
// factor, so only its shape is comparable.
//
// Input: in_event_ue pair lists ("<jet_output> <SimpleTree>"; only the jet
// files are read).
// Build and run from this folder:
//   g++ -m64 -O2 trig_compare.cpp ../../shared/VectorDict.cxx -o trig_compare -I../../shared -I../../underlying_event -I../../jet_scale `root-config --cflags --libs` -std=c++11
//   ./trig_compare "<data pair lists>" "<sim pair lists>" <out.png>

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

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}
// weighted fraction: pass / all, error from the effective count
struct Frac
{
    double w = 0, w2 = 0, p = 0, p2 = 0;
    void Add(bool pass, double ww) { w += ww; w2 += ww * ww; if (pass) { p += ww; p2 += ww * ww; } }
    double Val() const { return w > 0 ? p / w : NAN; }
    double Err() const { const double f = Val(); return w > 0 ? std::sqrt(std::max(f * (1 - f), 1.0 / (w * w / w2)) * w2) / w : NAN; }
    double Neff() const { return w2 > 0 ? w * w / w2 : 0; }
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 4) { std::fprintf(stderr, "usage: ./trig_compare \"<data pair lists>\" \"<sim pair lists>\" <out.png>\n"); return 1; }
    const std::vector<double> ed[2] = {{10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 120}, {2, 2.5, 3, 3.5, 4, 5, 6, 8, 12}};
    // [var: E, pT][sample][ratio: EM2|EM1, EM1|EM0][bin]
    std::vector<Frac> fr[2][2][2];
    for (int v = 0; v < 2; v++) for (int s = 0; s < 2; s++) for (int r = 0; r < 2; r++) fr[v][s][r].assign(ed[v].size() - 1, Frac());

    for (int s = 0; s < 2; s++)
    {
        // jet files, grouped by ptHat bin in simulation for the weights
        std::vector<std::string> files;
        const std::vector<std::string> lists = Split(argv[1 + s]);
        for (size_t l = 0; l < lists.size(); l++)
        {
            std::ifstream in(lists[l].c_str()); std::string a, b;
            while (in >> a >> b) files.push_back(a);
        }
        // simulation: sigma / N per ptHat bin (first pass over the files)
        std::map<int, double> sN, nN;
        std::vector<int> grp(files.size(), 0);
        if (s == 1)
            for (size_t f = 0; f < files.size(); f++)
            {
                const std::string base = files[f].substr(files[f].find_last_of('/') + 1);
                int a = 0, b = 0; const size_t p = base.rfind("_pt");
                if (p != std::string::npos && std::sscanf(base.c_str() + p, "_pt%d_%d.root", &a, &b) == 2) grp[f] = a;
                TFile x(files[f].c_str()); TTree *t = (TTree *)x.Get("jetTree");
                if (!t || t->GetEntries() == 0 || !t->GetBranch("sigma_job")) continue;
                Double_t sj; Int_t nj; t->SetBranchAddress("sigma_job", &sj); t->SetBranchAddress("N_job", &nj); t->GetEntry(0);
                const double sg = DecodePythia6Sigma(sj);
                if (sg > 0 && nj > 0) { sN[grp[f]] += sg * nj; nN[grp[f]] += nj; }
            }
        long nev = 0;
        for (size_t f = 0; f < files.size(); f++)
        {
            double w = 1;
            if (s == 1) { if (!(nN[grp[f]] > 0)) continue; w = sN[grp[f]] / nN[grp[f]] / nN[grp[f]]; }
            TFile x(files[f].c_str()); TTree *t = (TTree *)x.Get("jetTree");
            if (!t) continue;
            std::vector<float> *pt = 0, *eta = 0, *xx = 0, *yy = 0;
            Int_t e0, e1, e2;
            t->SetBranchStatus("*", 0);
            const char *br[] = {"reco_pt", "reco_eta", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2"};
            for (int k = 0; k < 7; k++) t->SetBranchStatus(br[k], 1);
            t->SetBranchAddress("reco_pt", &pt); t->SetBranchAddress("reco_eta", &eta); t->SetBranchAddress("reco_x", &xx); t->SetBranchAddress("reco_y", &yy);
            t->SetBranchAddress("fcsEM0", &e0); t->SetBranchAddress("fcsEM1", &e1); t->SetBranchAddress("fcsEM2", &e2);
            for (Long64_t i = 0; i < t->GetEntries(); i++)
            {
                t->GetEntry(i);
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
                    if (e1) fr[v][s][0][b].Add(e2 != 0, w);
                    if (e0) fr[v][s][1][b].Add(e1 != 0, w);
                }
                nev++;
            }
        }
        std::printf("%s: %zu files, %ld events with a fiducial leading jet\n", s ? "simulation" : "data", files.size(), nev);
    }

    const char *rn[2] = {"P(EM2 | EM1 fired)", "P(EM1 | EM0 fired)"};
    for (int v = 0; v < 2; v++)
        for (int r = 0; r < 2; r++)
        {
            std::printf("\n%s against leading-jet %s: data | sim\n", rn[r], v ? "pT" : "E");
            for (size_t b = 0; b + 1 < ed[v].size(); b++)
                std::printf("  %5.1f-%-5.1f  %.3f+-%.3f (n %6.0f) | %.3f+-%.3f (n_eff %5.0f)\n", ed[v][b], ed[v][b + 1], fr[v][0][r][b].Val(), fr[v][0][r][b].Err(),
                            fr[v][0][r][b].Neff(), fr[v][1][r][b].Val(), fr[v][1][r][b].Err(), fr[v][1][r][b].Neff());
        }

    gStyle->SetOptStat(0);
    const int cd = TColor::GetColor("#222222"), cs = TColor::GetColor("#d9480f");
    TCanvas cv("cv", "", 1500, 1000);
    cv.Divide(2, 2, 0.002, 0.002);
    std::vector<TObject *> keep;
    for (int r = 0; r < 2; r++)
        for (int v = 0; v < 2; v++)
        {
            cv.cd(r * 2 + v + 1);
            gPad->SetMargin(0.12, 0.03, 0.12, 0.09);
            if (v == 0) gPad->SetLogx();
            TH1D *h = new TH1D(Form("f%d%d", r, v), "", 1, ed[v].front(), ed[v].back()); keep.push_back(h);
            h->SetTitle(Form(";leading jet %s;%s", v ? "p_{T} (GeV)" : "E (GeV)", rn[r]));
            h->SetMinimum(0); h->SetMaximum(1.15);
            if (v == 0) h->GetXaxis()->SetMoreLogLabels();
            h->Draw();
            TLegend *lg = new TLegend(0.15, 0.76, 0.6, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042); keep.push_back(lg);
            for (int s = 0; s < 2; s++)
            {
                TGraphErrors *g = new TGraphErrors(); keep.push_back(g);
                for (size_t b = 0; b + 1 < ed[v].size(); b++)
                {
                    const Frac &f = fr[v][s][r][b];
                    if (f.Neff() < 5) continue;
                    const double xc = v ? 0.5 * (ed[v][b] + ed[v][b + 1]) : std::sqrt(ed[v][b] * ed[v][b + 1]);
                    const int i = g->GetN();
                    g->SetPoint(i, xc, f.Val()); g->SetPointError(i, 0, f.Err());
                }
                g->SetMarkerStyle(s ? 21 : 20); g->SetMarkerSize(1.3); g->SetMarkerColor(s ? cs : cd); g->SetLineColor(s ? cs : cd);
                g->Draw("PL SAME");
                lg->AddEntry(g, s ? "simulation (emulated trigger)" : "data (fired trigger)", "pl");
            }
            lg->Draw();
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
            t.DrawLatex(0.12, 0.93, r == 0 ? "EM2 among events that fired EM1" : "EM1 among events that fired EM0 (data: times EM1 prescale)");
        }
    cv.SaveAs(argv[3]);
    return 0;
}
