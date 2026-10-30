// Truth-level version of the two azimuthal UE scans of the FMS test
// (FMSJetPipeline/compare/wedge_density.cpp), on FCS simulation: truth EM jets
// and truth photons (IsEmTruthParticle, E > 0.2 GeV, eta 1.5-5.5), full
// azimuth (no detector), other jets kept.
//   wedges  the band |eta - eta_jet| < R with no jet exclusion (the jet's own
//           photons included, to see where the jet gives way to the UE), in
//           wedges of |dphi| to the jet (kEdges); each wedge weighted on
//           its own so that its area-per-eta matches the jet cone's (eta bins
//           of kBin, bins the wedge covers under kMinRatio of the jet's carry
//           no weight), rho = sum E w(eta) / sum_b w_b A_wedge(b). With no
//           detector the weights depend only on eta - eta_jet.
//   cones   the two cones at the jet's eta and phi +- dphi (kOff, from 0: on
//           the jet), density
//           (E1 + E2) / (2 pi R^2)
// Jets: truth pT >= reco_jet_pt_min_em, E 10-120 GeV, truth position inside
// the reco fiducial region (as the reco jets of in_event_ue); every event, and
// the simulation's EM bits stored (em: highest level, -1 none) to select
// triggered events later.
// Output tree "jets": jet_pt, jet_eta, jet_E, em, group (-1 - ptHatMin),
// weight (sigma / N of the ptHat bin), w_rho[kNW], c_rho[kNOff]; histogram
// "h_wcov" holds each wedge's coverage of the jet's eta range.
//
// Build and run from this folder:
//   g++ -m64 -O2 truth_scan.cpp ../../shared/VectorDict.cxx -o truth_scan -I../../shared -I../../underlying_event -I../../jet_scale `root-config --cflags --libs` -std=c++11
//   ./truth_scan <pairs list> <out.root>     (pairs "<jet_output> <SimpleTree>")

#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "in_event_util.h"
#include "JetXSecWeights.h"

static const double kBin = 0.05, kMinRatio = 0.1;
static const int kNB = 28;   // 2R / kBin
static const int kNW = 18;
static const double kEdges[kNW + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180};
static const int kNOff = 19;
static const double kOff[kNOff] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180};
static const int MAX_MC = 20000;

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::fprintf(stderr, "usage: ./truth_scan <pairs list> <out.root>\n"); return 1; }

    // wedge weights per eta bin of (eta - eta_jet + R): jet-cone and wedge areas from fine eta steps
    double ratio[kNW][kNB], den[kNW], cov[kNW], area_jet = 0;
    {
        double hj[kNB] = {0}, hw[kNW][kNB] = {{0}};
        const int nsub = 50;
        for (int b = 0; b < kNB; b++)
            for (int q = 0; q < nsub; q++)
            {
                const double d = -R + (b + (q + 0.5) / nsub) * kBin, de = kBin / nsub;
                const double s = std::sqrt(std::max(0.0, R * R - d * d));
                hj[b] += 2 * s * de;
                for (int k = 0; k < kNW; k++)
                {
                    const double lo = kEdges[k] * M_PI / 180, hi = kEdges[k + 1] * M_PI / 180;
                    hw[k][b] += 2 * (hi - lo) * de;
                }
            }
        for (int b = 0; b < kNB; b++) area_jet += hj[b];
        for (int k = 0; k < kNW; k++)
        {
            den[k] = 0;
            for (int b = 0; b < kNB; b++)
            {
                ratio[k][b] = hw[k][b] >= kMinRatio * hj[b] && hw[k][b] > 0 ? hj[b] / hw[k][b] : 0;
                den[k] += ratio[k][b] * hw[k][b];
            }
            cov[k] = den[k] / area_jet;
        }
    }

    std::vector<std::pair<std::string, std::string> > pairs;
    {
        std::ifstream in(argv[1]); std::string a, b;
        while (in >> a >> b) pairs.push_back(std::make_pair(a, b));
    }
    TFile fout(argv[2], "RECREATE");
    TTree *out = new TTree("jets", "truth-level azimuthal UE scans");
    Float_t o_pt, o_eta, o_E, o_w, w_rho[kNW], c_rho[kNOff];
    Int_t o_em, o_group;
    out->Branch("jet_pt", &o_pt, "jet_pt/F"); out->Branch("jet_eta", &o_eta, "jet_eta/F"); out->Branch("jet_E", &o_E, "jet_E/F");
    out->Branch("em", &o_em, "em/I"); out->Branch("group", &o_group, "group/I"); out->Branch("weight", &o_w, "weight/F");
    out->Branch("w_rho", w_rho, Form("w_rho[%d]/F", kNW)); out->Branch("c_rho", c_rho, Form("c_rho[%d]/F", kNOff));
    TH1D *h_wcov = new TH1D("h_wcov", ";wedge;coverage of the jet's eta range", kNW, 0, kNW);
    for (int k = 0; k < kNW; k++) h_wcov->SetBinContent(k + 1, cov[k]);

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
    long nev = 0, njet = 0;
    for (size_t f = 0; f < pairs.size(); f++)
    {
        if (!(nN[grp[f]] > 0)) continue;
        const double w = sN[grp[f]] / nN[grp[f]] / nN[grp[f]];
        TFile fj(pairs[f].first.c_str()), fs(pairs[f].second.c_str());
        TTree *jt = (TTree *)fj.Get("jetTree"), *st = (TTree *)fs.Get(fs.Get("SimpleTree") ? "SimpleTree" : "simpleTree");
        if (!st) { TIter it(fs.GetListOfKeys()); TObject *k; while ((k = it())) { TObject *o = fs.Get(k->GetName()); if (o && o->InheritsFrom("TTree")) { st = (TTree *)o; break; } } }
        if (!jt || !st) { std::fprintf(stderr, "skip %s\n", pairs[f].first.c_str()); continue; }
        std::vector<float> *tpt = 0, *teta = 0, *tphi = 0, *tE = 0, *tx = 0, *ty = 0;
        Int_t em0, em1, em2; Long64_t sent = -1;
        const bool has_entry = jt->GetBranch("simpletree_entry");
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"truth_pt", "truth_eta", "truth_phi", "truth_E", "truth_x", "truth_y", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int i = 0; i < 9; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("truth_pt", &tpt); jt->SetBranchAddress("truth_eta", &teta); jt->SetBranchAddress("truth_phi", &tphi);
        jt->SetBranchAddress("truth_E", &tE); jt->SetBranchAddress("truth_x", &tx); jt->SetBranchAddress("truth_y", &ty);
        jt->SetBranchAddress("fcsEM0", &em0); jt->SetBranchAddress("fcsEM1", &em1); jt->SetBranchAddress("fcsEM2", &em2);
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
                if (mc_E[i] < 0.2 || !IsEmTruthParticle(mc_geid[i], mc_vx[i], mc_vy[i], mc_vz[i], pvx, pvy, pvz)) continue;
                const double ptp = std::hypot(mc_px[i], mc_py[i]);
                if (!(ptp > 0)) continue;
                const double e = std::asinh(mc_pz[i] / ptp);
                if (e < 1.5 || e > 5.5) continue;
                pe.push_back(e); pp.push_back(std::atan2(mc_py[i], mc_px[i])); pE.push_back(mc_E[i]);
            }
            nev++;
            o_em = em2 ? 2 : em1 ? 1 : em0 ? 0 : -1;
            o_group = -1 - grp[f]; o_w = w;
            for (size_t q = 0; q < sel.size(); q++)
            {
                const int j = sel[q];
                const double je = teta->at(j), jp = tphi->at(j);
                double num[kNW] = {0}, cE[kNOff] = {0};
                for (size_t i = 0; i < pE.size(); i++)
                {
                    const double de = pe[i] - je, adp = std::fabs(std::remainder(pp[i] - jp, 2 * M_PI));
                    if (std::fabs(de) < R)
                    {
                        const int b = (int)std::floor((de + R) / kBin), k = std::min(kNW - 1, (int)(adp * 180 / M_PI / 10));
                        if (b >= 0 && b < kNB) num[k] += pE[i] * ratio[k][b];
                    }
                    for (int c = 0; c < kNOff; c++)
                        for (int s = -1; s <= 1; s += 2)
                            if (std::hypot(de, std::remainder(pp[i] - jp - s * kOff[c] * M_PI / 180, 2 * M_PI)) < R) cE[c] += pE[i];
                }
                for (int k = 0; k < kNW; k++) w_rho[k] = den[k] > 0 ? num[k] / den[k] : -1;
                for (int c = 0; c < kNOff; c++) c_rho[c] = cE[c] / (2 * M_PI * R * R);
                o_pt = tpt->at(j); o_eta = je; o_E = tE->at(j);
                out->Fill();
                njet++;
            }
        }
    }
    std::printf("truth_scan: %zu pairs, %ld events with a selected truth jet, %ld jets\n", pairs.size(), nev, njet);
    fout.Write();
    return 0;
}
