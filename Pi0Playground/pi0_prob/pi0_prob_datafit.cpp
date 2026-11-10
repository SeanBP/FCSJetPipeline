// Data refit of the pi0 peak in S(m | E_pair) for the per-jet pi0 probability
// (shared/Pi0Prob.h). The calibration's S is fitted on truth pi0 pairs of
// simulated two-point jets; the data peak sits lower (point energy scale,
// ~8-17% for the cluster-pair study, unknown for points) and may be wider
// (vertex spread: sim sigma_z 30 cm, data larger). Same fit, data and sim:
// every two-point jet's pair mass, per pair-energy bin, fitted with a Crystal
// Ball (tail alpha, n fixed to the calibration's truth S of that bin) plus a
// quadratic background on [0.03, 0.40] GeV. Then, per bin,
//   S_data_mean  = S_mean  * mu_data / mu_sim
//   S_data_sigma = S_sigma * sigma_data / sigma_sim
// so the fit's own biases cancel. A bin whose data or sim fit fails takes
// the ratios of the nearest good bin below it (above for the first).
// Sim: the training half of the dumps (cross-section weighted); data: kept
// jet files of a Run 22 pass (jetTree, reco_const_*), no LED, ECal point
// constituents only. Writes <out>.json (the calibration with the two data
// arrays added), <out>.pdf and a table.
//
// Build:
//   g++ -O2 pi0_prob_datafit.cpp -o pi0_prob_datafit -I../../shared -I../../jet_scale \
//       `root-config --cflags --libs` -std=c++11
// Run:  ./pi0_prob_datafit <calib.json> "<dump glob>" <SimpleTree list, or - for in-job dumps> <data jet file list> <out prefix>
#include "Pi0Prob.h"
#include "Pi0DumpWeights.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

static double CBq(double *x, double *p)
{
    // p0 norm, p1 mu, p2 sigma, p3 alpha (fixed), p4 n (fixed), p5-p7 quadratic
    const double t = (x[0] - p[1]) / p[2], a = fabs(p[3]), n = p[4];
    double s;
    if (t > -a) s = exp(-0.5 * t * t);
    else s = pow(n / a, n) * exp(-0.5 * a * a) * pow(n / a - a - t, -n);
    return p[0] * s + p[5] + p[6] * x[0] + p[7] * x[0] * x[0];
}

struct Fit { double mu, sig, emu, esig, nsig; bool ok; };

static Fit DoFit(TH1D *h, double alpha, double n, double mu0, double sig0, TF1 *keep)
{
    Fit r = {0, 0, 0, 0, 0, false};
    if (h->GetEffectiveEntries() < 200) return r;
    TF1 f("f", CBq, 0.03, 0.40, 8);
    const double pk = h->GetBinContent(h->FindBin(mu0));
    f.SetParameters(pk, mu0, sig0, alpha, n, 0.2 * pk, 0, 0);
    f.FixParameter(3, alpha); f.FixParameter(4, n);
    f.SetParLimits(0, 0, 10 * h->GetMaximum() + 1);
    f.SetParLimits(1, 0.07, 0.22); f.SetParLimits(2, 0.004, 0.08);
    if (h->Fit(&f, "QRN0") != 0) return r;
    r.mu = f.GetParameter(1); r.sig = f.GetParameter(2); r.emu = f.GetParError(1); r.esig = f.GetParError(2);
    r.nsig = f.GetParameter(0) * r.sig * sqrt(2 * M_PI) / h->GetBinWidth(1);
    r.ok = r.emu > 0 && r.emu < 0.01 && r.mu > 0.075 && r.mu < 0.215 && r.sig > 0.0045 && r.sig < 0.075 && f.GetParameter(0) > 3 * f.GetParError(0);
    if (keep) for (int i = 0; i < 8; i++) keep->SetParameter(i, f.GetParameter(i));
    return r;
}

int main(int argc, char **argv)
{
    if (argc < 6) { fprintf(stderr, "Usage: %s <calib.json> \"<dump glob>\" <SimpleTree list> <data jet list> <out prefix>\n", argv[0]); return 1; }
    gROOT->SetBatch(1);
    gStyle->SetOptStat(0);
    Pi0ProbCalib cal;
    string err;
    if (!cal.Load(argv[1], false, &err)) { fprintf(stderr, "%s\n", err.c_str()); return 1; }
    const int nE = (int)cal.eS.size() - 1;
    const string out = argv[5];
    vector<TH1D *> hs(nE), hd(nE);
    for (int e = 0; e < nE; e++)
    {
        hs[e] = new TH1D(Form("hsim_%d", e), Form("2-point jets, E_{pair} %g-%g GeV;m_{#gamma#gamma} (GeV)", cal.eS[e], cal.eS[e + 1]), 100, 0, 0.5);
        hd[e] = new TH1D(Form("hdat_%d", e), "", 100, 0, 0.5);
        hs[e]->Sumw2(); hd[e]->Sumw2();
    }

    // ---- sim: training half of the dumps, weighted ----
    {
        TChain ch("jets");
        ch.Add(argv[2]);
        map<int, double> wOf;
        if (!Pi0DumpXSecWeights(argv[2], argv[3], wOf)) return 1;
        Int_t file_idx, n; Float_t pE[200], px[200], py[200], pz[200];
        ch.SetBranchStatus("*", 0);
        for (auto b : {"file_idx", "n", "pt_E", "pt_x", "pt_y", "pt_z"}) ch.SetBranchStatus(b, 1);
        ch.SetBranchAddress("file_idx", &file_idx); ch.SetBranchAddress("n", &n);
        ch.SetBranchAddress("pt_E", pE); ch.SetBranchAddress("pt_x", px); ch.SetBranchAddress("pt_y", py); ch.SetBranchAddress("pt_z", pz);
        for (Long64_t i = 0; i < ch.GetEntries(); i++)
        {
            ch.GetEntry(i);
            if (n != 2 || file_idx % 2) continue;
            auto wit = wOf.find(file_idx);
            if (wit == wOf.end()) continue;
            const double w = wit->second;
            if (!(w > 0)) continue;
            double v[2][4];
            for (int a = 0; a < 2; a++)
            {
                const double r = sqrt(px[a] * px[a] + py[a] * py[a] + pz[a] * pz[a]);
                v[a][0] = pE[a]; v[a][1] = pE[a] * px[a] / r; v[a][2] = pE[a] * py[a] / r; v[a][3] = pE[a] * pz[a] / r;
            }
            const double m = Pi0PairMass(v[0][0], v[0][1], v[0][2], v[0][3], v[1][0], v[1][1], v[1][2], v[1][3]);
            hs[Pi0ProbCalib::Bin(cal.eS, pE[0] + pE[1])]->Fill(m, w);
        }
    }
    // ---- data: kept jet files ----
    long nd = 0;
    {
        TChain ch("jetTree");
        ifstream in(argv[4]); string l;
        while (getline(in, l)) if (!l.empty()) ch.Add(l.c_str());
        vector<float> *E = 0, *X = 0, *Y = 0, *Z = 0; vector<int> *ji = 0, *det = 0;
        Int_t led = 0;
        ch.SetBranchStatus("*", 0);
        for (auto b : {"reco_const_E", "reco_const_px", "reco_const_py", "reco_const_pz", "reco_const_jetIndex", "reco_const_detid", "fcs_led"}) ch.SetBranchStatus(b, 1);
        ch.SetBranchAddress("reco_const_E", &E); ch.SetBranchAddress("reco_const_px", &X); ch.SetBranchAddress("reco_const_py", &Y);
        ch.SetBranchAddress("reco_const_pz", &Z); ch.SetBranchAddress("reco_const_jetIndex", &ji); ch.SetBranchAddress("reco_const_detid", &det);
        ch.SetBranchAddress("fcs_led", &led);
        for (Long64_t i = 0; i < ch.GetEntries(); i++)
        {
            ch.GetEntry(i);
            if (led) continue;
            map<int, vector<int> > byJet;
            for (size_t c = 0; c < E->size(); c++) if ((*det)[c] == 0 || (*det)[c] == 1) byJet[(*ji)[c]].push_back((int)c);
            // a jet with exactly two constituents, both ECal points
            map<int, int> ncon;
            for (size_t c = 0; c < E->size(); c++) ncon[(*ji)[c]]++;
            for (auto &kv : byJet)
            {
                if (kv.second.size() != 2 || ncon[kv.first] != 2) continue;
                const int a = kv.second[0], b = kv.second[1];
                const double m = Pi0PairMass((*E)[a], (*X)[a], (*Y)[a], (*Z)[a], (*E)[b], (*X)[b], (*Y)[b], (*Z)[b]);
                hd[Pi0ProbCalib::Bin(cal.eS, (*E)[a] + (*E)[b])]->Fill(m);
                nd++;
            }
        }
    }
    printf("data two-point jets: %ld\n", nd);

    // ---- fits ----
    vector<Fit> fs(nE), fd(nE);
    TCanvas cv("cv", "", 1200, 900);
    cv.Divide(3, 3);
    printf("%-10s %8s %8s | %8s %8s | %8s %8s | %7s %7s | %8s %8s\n", "E_pair", "mu_sim", "sig_sim", "mu_dat", "sig_dat", "mu_d/s", "sig_d/s",
           "Ssim", "Sdat", "S_mean", "S_sigma");
    for (int e = 0; e < nE; e++)
    {
        TF1 *ks = new TF1(Form("ks_%d", e), CBq, 0.03, 0.40, 8), *kd = new TF1(Form("kd_%d", e), CBq, 0.03, 0.40, 8);
        fs[e] = DoFit(hs[e], cal.sAlpha[e], cal.sN[e], cal.sMean[e], cal.sSigma[e], ks);
        fd[e] = DoFit(hd[e], cal.sAlpha[e], cal.sN[e], 0.88 * cal.sMean[e], cal.sSigma[e], kd);
        cv.cd(e + 1);
        TH1D *a = (TH1D *)hs[e]->Clone(), *b = (TH1D *)hd[e]->Clone();
        const double ia = a->Integral(), ib = b->Integral();
        if (ia > 0) a->Scale(1 / ia);
        if (ib > 0) b->Scale(1 / ib);
        a->SetLineColor(kBlue + 1); b->SetLineColor(kBlack); b->SetMarkerStyle(20); b->SetMarkerSize(0.5);
        a->SetMaximum(1.3 * max(a->GetMaximum(), b->GetMaximum()));
        a->Draw("hist"); b->Draw("P same");
        if (fs[e].ok) { ks->SetParameter(0, ks->GetParameter(0) / ia); for (int q = 5; q < 8; q++) ks->SetParameter(q, ks->GetParameter(q) / ia); ks->SetLineColor(kBlue + 1); ks->SetLineStyle(2); ks->Draw("same"); }
        if (fd[e].ok) { kd->SetParameter(0, kd->GetParameter(0) / ib); for (int q = 5; q < 8; q++) kd->SetParameter(q, kd->GetParameter(q) / ib); kd->SetLineColor(kRed + 1); kd->Draw("same"); }
        TLegend *L = new TLegend(0.45, 0.68, 0.89, 0.89);
        L->SetBorderSize(0); L->SetFillStyle(0);
        L->AddEntry(a, fs[e].ok ? Form("sim: #mu %.3f #sigma %.3f", fs[e].mu, fs[e].sig) : "sim", "l");
        L->AddEntry(b, fd[e].ok ? Form("data: #mu %.3f #sigma %.3f", fd[e].mu, fd[e].sig) : "data", "p");
        L->Draw();
    }
    cv.Print((out + ".pdf").c_str());

    // ratios, with the nearest good bin's for failed bins
    vector<double> rmu(nE, 0), rsig(nE, 0);
    vector<bool> good(nE);
    for (int e = 0; e < nE; e++)
    {
        good[e] = fs[e].ok && fd[e].ok;
        if (good[e]) { rmu[e] = fd[e].mu / fs[e].mu; rsig[e] = fd[e].sig / fs[e].sig; }
    }
    for (int e = 0; e < nE; e++)
    {
        if (good[e]) continue;
        int src = -1;
        for (int q = e - 1; q >= 0 && src < 0; q--) if (good[q]) src = q;
        for (int q = e + 1; q < nE && src < 0; q++) if (good[q]) src = q;
        if (src >= 0) { rmu[e] = rmu[src]; rsig[e] = rsig[src]; }
        else { rmu[e] = 1; rsig[e] = 1; }
    }
    vector<double> dm(nE), ds(nE);
    for (int e = 0; e < nE; e++)
    {
        dm[e] = cal.sMean[e] * rmu[e]; ds[e] = cal.sSigma[e] * rsig[e];
        printf("%4g-%-5g %8.4f %8.4f | %8.4f %8.4f | %8.3f %8.3f | %7.0f %7.0f | %8.4f %8.4f%s\n", cal.eS[e], cal.eS[e + 1],
               fs[e].mu, fs[e].sig, fd[e].mu, fd[e].sig, rmu[e], rsig[e], fs[e].nsig, fd[e].nsig, dm[e], ds[e], good[e] ? "" : "  (ratios from a neighbour bin)");
    }

    // ---- JSON: the calibration plus the data arrays ----
    ifstream in(argv[1]);
    stringstream ss;
    ss << in.rdbuf();
    string js = ss.str();
    const size_t endb = js.rfind('}');
    size_t last = js.find_last_not_of(" \n\t\r", endb - 1);
    ostringstream add;
    add << ",\n  \"S_data_source\": \"two-point jets of " << argv[4] << " (" << nd << " jets); CB(alpha, n fixed)+pol2, data/sim ratios\",\n";
    add << "  \"S_data_mean\": [";
    for (int e = 0; e < nE; e++) add << (e ? ", " : "") << dm[e];
    add << "],\n  \"S_data_sigma\": [";
    for (int e = 0; e < nE; e++) add << (e ? ", " : "") << ds[e];
    add << "]\n";
    js = js.substr(0, last + 1) + add.str() + "}\n";
    ofstream o((out + ".json").c_str());
    o << js;
    printf("wrote %s.json, %s.pdf\n", out.c_str(), out.c_str());
    return 0;
}
