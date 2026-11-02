// The UE energy each method subtracts from a truth jet, per jet, on FCS
// simulation: the truth-level counterpart of the FMS data plot
// (FMSJetPipeline/compare/ue_subtracted.cpp -c -n).
//   FMS two cones  the two cones at the jet's eta and phi +- 90 deg: mean of
//                  their truth-photon energies, a full cone being pi R^2
//   pixel ring     truth_ue_E of the jet tree: the pipeline's truth ring
//                  (TruthPixelRing, |dphi| <= 150 deg, other jets kept)
// Photons as the truth ring: IsEmTruthParticle, InTruthParticleEtaRange, no
// energy cut. Full azimuth at truth, so the FMS plot's full and FCS-like
// rings are one here. Jets as truth_scan: truth pT >= reco_jet_pt_min_em, E
// 10-120 GeV, position inside the reco fiducial region. Weight sigma / N of the
// ptHat bin. Histograms normalised to unit weight; linear and log panels (the
// linear one scaled to the bins above zero; the last bin holds overflow).
// "-r": the subtracted UE as a fraction of the jet energy.
//
// Build and run from this folder:
//   g++ -m64 -O2 truth_ue_subtracted.cpp ../../shared/VectorDict.cxx -o truth_ue_subtracted -I../../shared -I../../underlying_event -I../../jet_scale `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//   ./truth_ue_subtracted <pairs list> <out.root>     (pairs "<jet_output> <SimpleTree>")
//   ./truth_ue_subtracted plot <out.root> <out.png> [-r]

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "in_event_util.h"
#include "JetXSecWeights.h"

static const int MAX_MC = 20000;

static int Plot(const char *in, const char *png, bool frac)
{
    TFile f(in);
    TTree *t = (TTree *)f.Get("jets");
    Float_t cone, ring, w, jet_E;
    t->SetBranchAddress("cone_dE", &cone); t->SetBranchAddress("ring_dE", &ring); t->SetBranchAddress("weight", &w); t->SetBranchAddress("jet_E", &jet_E);
    const int kNM = 2;
    const char *name[kNM] = {"FMS two cones #pm90#circ (per #piR^{2})", "pixel ring (#leq150#circ, jets kept)"};
    const int col[kNM] = {TColor::GetColor("#222222"), TColor::GetColor("#1c4f8a")};
    const double xmax = frac ? 0.25 : 8.0;
    const char *unit = frac ? "" : " GeV";
    TH1D *h[kNM];
    double sw[kNM] = {0}, sw2[kNM] = {0}, s[kNM] = {0}, ss[kNM] = {0}, nz[kNM] = {0};
    for (int k = 0; k < kNM; k++)
    {
        h[k] = new TH1D(Form("h%d", k), frac ? ";UE energy subtracted / jet energy;jets (normalised)" : ";UE energy subtracted from the jet (GeV);jets (normalised)", 80, 0, xmax);
        h[k]->Sumw2();
    }
    const Long64_t nj = t->GetEntries();
    for (Long64_t i = 0; i < nj; i++)
    {
        t->GetEntry(i);
        double v[kNM] = {cone, ring};
        if (frac) for (int k = 0; k < kNM; k++) if (v[k] >= 0) v[k] /= jet_E;
        for (int k = 0; k < kNM; k++)
        {
            if (v[k] < 0) continue;
            h[k]->Fill(std::min(v[k], xmax * 0.9999), w);
            sw[k] += w; sw2[k] += (double)w * w; s[k] += w * v[k]; ss[k] += w * v[k] * v[k];
            if (v[k] == 0) nz[k] += w;
        }
    }
    double m[kNM], r[kNM];
    std::printf("truth_ue_subtracted: %lld truth jets%s\n method                                              n_eff     mean          RMS         share = 0\n", (long long)nj, frac ? " (UE / jet energy)" : " (GeV)");
    for (int k = 0; k < kNM; k++)
    {
        m[k] = s[k] / sw[k]; r[k] = std::sqrt(std::max(0.0, ss[k] / sw[k] - m[k] * m[k]));
        const double neff = sw[k] * sw[k] / sw2[k];
        std::printf(" %-50s %7.0f   %.4f+-%.4f   %.4f      %.3f\n", name[k], neff, m[k], r[k] / std::sqrt(neff), r[k], nz[k] / sw[k]);
    }
    gStyle->SetOptStat(0);
    TCanvas cv("cv", "", 1500, 760);
    TPad *pads[2] = {new TPad("pa", "", 0, 0, 0.5, 0.82), new TPad("pb", "", 0.5, 0, 1, 0.82)};
    pads[0]->Draw(); pads[1]->Draw();
    TLegend *lg = new TLegend(0.04, 0.83, 0.99, 0.94); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.026); lg->SetNColumns(2); lg->SetMargin(0.12);
    for (int p = 0; p < 2; p++)
    {
        pads[p]->cd();
        gPad->SetMargin(0.12, 0.03, 0.12, 0.07);
        if (p == 1) gPad->SetLogy();
        double ymax = 0;
        for (int k = 0; k < kNM; k++)
        {
            h[k]->Scale(p == 0 ? 1.0 / h[k]->Integral() : 1.0);
            for (int b = p == 0 ? 2 : 1; b <= h[k]->GetNbinsX(); b++) ymax = std::max(ymax, h[k]->GetBinContent(b));
        }
        for (int k = 0; k < kNM; k++)
        {
            h[k]->SetLineColor(col[k]); h[k]->SetLineWidth(3); h[k]->SetMarkerColor(col[k]);
            h[k]->SetMaximum(p == 0 ? 1.15 * ymax : 3 * ymax); h[k]->SetMinimum(p == 0 ? 0 : 1e-5);
            h[k]->GetYaxis()->SetTitleOffset(1.3);
            h[k]->DrawCopy(k == 0 ? "HIST" : "HIST SAME");
            if (p == 1) continue;
            TLegendEntry *le = lg->AddEntry(h[k], frac ? Form("%s: mean %.3f, RMS %.3f, zero %.0f%%", name[k], m[k], r[k], 100 * nz[k] / sw[k])
                                                       : Form("%s: mean %.2f%s, RMS %.2f, zero %.0f%%", name[k], m[k], unit, r[k], 100 * nz[k] / sw[k]), "l");
            le->SetTextColor(col[k]);
        }
        TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.04);
        l.DrawLatex(0.12, 0.95, p ? "log scale" : "linear scale");
    }
    cv.cd();
    lg->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextFont(42); tt.SetTextSize(0.03);
    tt.DrawLatex(0.04, 0.955, Form("FCS simulation (Pythia 6, k = 0.947), %lld truth jets, E 10-120 GeV: UE energy subtracted per jet%s, truth level", (long long)nj, frac ? ", over the jet energy" : ""));
    cv.SaveAs(png);
    return 0;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc >= 4 && !std::strcmp(argv[1], "plot")) return Plot(argv[2], argv[3], argc > 4 && !std::strcmp(argv[4], "-r"));
    if (argc < 3) { std::fprintf(stderr, "usage: ./truth_ue_subtracted <pairs list> <out.root> | plot <out.root> <out.png> [-r]\n"); return 1; }

    std::vector<std::pair<std::string, std::string> > pairs;
    {
        std::ifstream in(argv[1]); std::string a, b;
        while (in >> a >> b) pairs.push_back(std::make_pair(a, b));
    }
    TFile fout(argv[2], "RECREATE");
    TTree *out = new TTree("jets", "UE subtracted from truth jets");
    Float_t o_E, o_w, o_cone, o_ring;
    out->Branch("jet_E", &o_E, "jet_E/F"); out->Branch("weight", &o_w, "weight/F");
    out->Branch("cone_dE", &o_cone, "cone_dE/F"); out->Branch("ring_dE", &o_ring, "ring_dE/F");

    // ptHat-bin weights: sigma / N over the files of each bin
    std::map<int, double> sN, nN;
    std::vector<int> grp(pairs.size(), 0);
    for (size_t f = 0; f < pairs.size(); f++)
    {
        const std::string base = pairs[f].first.substr(pairs[f].first.find_last_of('/') + 1);
        const size_t p = base.rfind("_pt");
        int a = 0;
        if (p != std::string::npos) std::sscanf(base.c_str() + p, "_pt%d", &a);
        grp[f] = a;
        TFile x(pairs[f].first.c_str()); TTree *t = (TTree *)x.Get("jetTree");
        if (!t || t->GetEntries() == 0) continue;
        Double_t sj; Int_t nj; t->SetBranchAddress("sigma_job", &sj); t->SetBranchAddress("N_job", &nj); t->GetEntry(0);
        const double sg = DecodePythia6Sigma(sj);
        if (sg > 0 && nj > 0) { sN[a] += sg * nj; nN[a] += nj; }
    }

    static Int_t mc_n, mc_geid[MAX_MC], mc_idVtx[MAX_MC];
    static Float_t mc_px[MAX_MC], mc_py[MAX_MC], mc_pz[MAX_MC], mc_E[MAX_MC], mc_vx[MAX_MC], mc_vy[MAX_MC], mc_vz[MAX_MC];
    long nev = 0, njet = 0, nnoring = 0;
    for (size_t f = 0; f < pairs.size(); f++)
    {
        if (!(nN[grp[f]] > 0)) continue;
        const double w = sN[grp[f]] / nN[grp[f]] / nN[grp[f]];
        TFile fj(pairs[f].first.c_str()), fs(pairs[f].second.c_str());
        TTree *jt = (TTree *)fj.Get("jetTree"), *st = (TTree *)fs.Get(fs.Get("SimpleTree") ? "SimpleTree" : "simpleTree");
        if (!st) { TIter it(fs.GetListOfKeys()); TObject *k; while ((k = it())) { TObject *o = fs.Get(k->GetName()); if (o && o->InheritsFrom("TTree")) { st = (TTree *)o; break; } } }
        if (!jt || !st) { std::fprintf(stderr, "skip %s\n", pairs[f].first.c_str()); continue; }
        std::vector<float> *tpt = 0, *teta = 0, *tphi = 0, *tE = 0, *tx = 0, *ty = 0, *tue = 0;
        Long64_t sent = -1;
        const bool has_entry = jt->GetBranch("simpletree_entry");
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"truth_pt", "truth_eta", "truth_phi", "truth_E", "truth_x", "truth_y", "truth_ue_E"};
        for (int i = 0; i < 7; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("truth_pt", &tpt); jt->SetBranchAddress("truth_eta", &teta); jt->SetBranchAddress("truth_phi", &tphi);
        jt->SetBranchAddress("truth_E", &tE); jt->SetBranchAddress("truth_x", &tx); jt->SetBranchAddress("truth_y", &ty);
        jt->SetBranchAddress("truth_ue_E", &tue);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &sent); }
        st->SetBranchStatus("*", 0);
        const char *mb[] = {"mcpart_num", "mcpart_geid", "mcpart_idVtx", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
        for (int i = 0; i < 10; i++) st->SetBranchStatus(mb[i], 1);
        st->SetBranchAddress("mcpart_num", &mc_n); st->SetBranchAddress("mcpart_geid", mc_geid); st->SetBranchAddress("mcpart_idVtx", mc_idVtx);
        st->SetBranchAddress("mcpart_px", mc_px); st->SetBranchAddress("mcpart_py", mc_py); st->SetBranchAddress("mcpart_pz", mc_pz);
        st->SetBranchAddress("mcpart_E", mc_E); st->SetBranchAddress("mcpart_Vtx_x", mc_vx); st->SetBranchAddress("mcpart_Vtx_y", mc_vy);
        st->SetBranchAddress("mcpart_Vtx_z", mc_vz);
        std::vector<double> pe, pp, pE;
        for (Long64_t ev = 0; ev < jt->GetEntries(); ev++)
        {
            if (jt->GetEntry(ev) <= 0 || !tpt || tpt->empty()) continue;
            std::vector<int> sel;
            for (size_t j = 0; j < tpt->size(); j++)
                if (tpt->at(j) >= reco_jet_pt_min_em && tE->at(j) >= 10 && tE->at(j) < 120 && ProxyInFiducial(tx->at(j), ty->at(j))) sel.push_back((int)j);
            if (sel.empty()) continue;
            if (st->GetEntry(has_entry ? sent : ev) <= 0) continue;
            float pvx, pvy, pvz;
            const int n = std::min((int)mc_n, MAX_MC);
            if (!FindEmTruthPrimaryVertex(n, mc_idVtx, mc_vx, mc_vy, mc_vz, pvx, pvy, pvz)) continue;
            pe.clear(); pp.clear(); pE.clear();
            for (int i = 0; i < n; i++)
            {
                if (!IsEmTruthParticle(mc_geid[i], mc_vx[i], mc_vy[i], mc_vz[i], pvx, pvy, pvz)) continue;
                const double ptp = std::hypot(mc_px[i], mc_py[i]);
                if (!(ptp > 0)) continue;
                const double e = std::asinh(mc_pz[i] / ptp);
                if (!InTruthParticleEtaRange(e)) continue;
                pe.push_back(e); pp.push_back(std::atan2(mc_py[i], mc_px[i])); pE.push_back(mc_E[i]);
            }
            nev++;
            for (size_t q = 0; q < sel.size(); q++)
            {
                const int j = sel[q];
                const double je = teta->at(j), jp = tphi->at(j);
                double c = 0;
                for (size_t i = 0; i < pE.size(); i++)
                    for (int s = -1; s <= 1; s += 2)
                        if (std::hypot(pe[i] - je, std::remainder(pp[i] - jp - s * M_PI / 2, 2 * M_PI)) < R) c += pE[i];
                o_E = tE->at(j); o_w = w; o_cone = 0.5 * c; o_ring = tue->at(j);
                if (o_ring < 0) nnoring++;
                out->Fill();
                njet++;
            }
        }
    }
    std::printf("truth_ue_subtracted: %zu pairs, %ld events with a selected truth jet, %ld jets, %ld without a ring estimate\n", pairs.size(), nev, njet, nnoring);
    fout.Write();
    return 0;
}
