// FCS counterpart of the FMS wedge scan (FMSJetPipeline/compare/wedge_density.cpp,
// FmsUnderlyingEvent wedge_rho): the UE density around each jet in 10 deg
// wedges of |dphi| to the jet, 0 to 180 deg, on real FCS data, to compare the
// recoil rise towards 180 deg with the FMS at 200 GeV.
//
// Each wedge is the pipeline's eta-weighted side wedges (WeightedWedgeUE.h:
// the jet's eta band |eta - eta_jet| < R on the ECal tiles, other jets kept,
// each eta bin weighted so the wedge's area-per-eta matches the jet cone's)
// restricted to |dphi| in [k*10, (k+1)*10) deg and with no exclusion around
// the jet (as the FMS scan: the first wedges hold the jet itself). Points as
// the reco jets' inputs (ECal points above the point energy cut).
//
// Jets: as the FMS FCS-like sample and ring_cap_check, reco pT >= 2 GeV,
// E 10-120 GeV, axis inside the ECal fiducial region with the 20 cm buffer;
// events as written by data_to_jet (EM-triggered). Also stored: whether the
// jet passes the A_N trigger cut (uncorrected pT > 1.15 x the threshold of a
// listed EM level, AnConfig.h values).
//
// Output tree "wedges", one row per jet: run, jet_E, jet_pt, jet_eta,
// jet_phi, an_trig, wedge_rho[18], wedge_cov[18] (density GeV per unit
// eta-phi area, -1 if the wedge reaches none of the jet's eta bins; share of
// the jet's live area in eta bins the wedge reaches).
//
// Build and run from this folder:
//   g++ -m64 -O2 wedge_scan.cpp ../../shared/VectorDict.cxx -o wedge_scan -I../../shared -I../../underlying_event `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//   ./wedge_scan <pairs list: "<jet_output> <SimpleTree>"> <out.root>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "VectorDict.h"
#include "WeightedWedgeUE.h"

static const int kNW = 18;
static const int MAX_PTS = 2000;
static const double kThr[3] = {2.0, 3.5, 5.0}, kTrigFactor = 1.15;

// RecoWeightedWedges for every wedge |dphi| in [k, k+1) x 10 deg at once, with
// no jet exclusion: one pass over the tiles fills the jet cone's and each
// wedge's area per eta bin, one pass over the points the weighted energies
static void AllWedges(const UeGeometry &g, double jeta, double jphi, const std::vector<float> &p_eta, const std::vector<float> &p_phi,
                      const std::vector<float> &p_E, const std::vector<int> &p_tile, float *rho, float *cov)
{
    const double eta0 = jeta - R - 0.1, wdeg = 10 * M_PI / 180;
    static double hr[kNW][kUeNB];
    double hj[kUeNB] = {0}, area = 0;
    for (int k = 0; k < kNW; k++) for (int b = 0; b < kUeNB; b++) hr[k][b] = 0;
    std::vector<signed char> wk(g.tiles.size(), -1);
    const int ns = kUeSub * kUeSub;
    for (size_t c = 0; c < g.tiles.size(); c++)
    {
        if (std::fabs(g.tiles[c].eta - jeta) > R + 0.2) continue;
        if (std::fabs(g.tiles[c].eta - jeta) < R)
            wk[c] = std::min(kNW - 1, (int)(std::fabs(std::remainder(g.tiles[c].phi - jphi, 2.0 * M_PI)) / wdeg));
        for (int q = 0; q < ns; q++)
        {
            const UeGeometry::Sub &u = g.sub[c * ns + q];
            const int b = (int)std::floor((u.eta - eta0) / kUeBin);
            if (b < 0 || b >= kUeNB) continue;
            if (UeDeltaR(u.eta, u.phi, jeta, jphi) < R) { hj[b] += u.area; area += u.area; }
            if (wk[c] >= 0) hr[wk[c]][b] += u.area;
        }
    }
    double ratio[kNW][kUeNB], den[kNW], sE[kNW];
    for (int k = 0; k < kNW; k++)
    {
        den[k] = 0; sE[k] = 0;
        for (int b = 0; b < kUeNB; b++)
        {
            ratio[k][b] = hj[b] > 0 && hr[k][b] >= kUeMinRatio * hj[b] ? hj[b] / hr[k][b] : 0;
            den[k] += ratio[k][b] * hr[k][b];
        }
    }
    for (size_t i = 0; i < p_E.size(); i++)
    {
        if (p_tile[i] < 0 || wk[p_tile[i]] < 0) continue;
        const int b = (int)std::floor((p_eta[i] - eta0) / kUeBin);
        if (b < 0 || b >= kUeNB) continue;
        const int k = wk[p_tile[i]];
        sE[k] += p_E[i] * ratio[k][b];
    }
    for (int k = 0; k < kNW; k++)
    {
        cov[k] = area > 0 ? den[k] / area : 0;
        rho[k] = den[k] > 0 ? sE[k] / den[k] : -1;
    }
}

int main(int argc, char **argv)
{
    if (argc < 3) { std::cerr << "Usage: " << argv[0] << " <pairs list> <out.root>" << std::endl; return 1; }
    UeGeometry geo;
    geo.Build();
    TFile out(argv[2], "RECREATE");
    TTree *t = new TTree("wedges", "FCS UE density in 10 deg wedges of |dphi| to the jet");
    Int_t o_run, o_trig;
    Float_t o_E, o_pt, o_eta, o_phi, o_rho[kNW], o_cov[kNW];
    t->Branch("run", &o_run, "run/I");
    t->Branch("jet_E", &o_E, "jet_E/F"); t->Branch("jet_pt", &o_pt, "jet_pt/F");
    t->Branch("jet_eta", &o_eta, "jet_eta/F"); t->Branch("jet_phi", &o_phi, "jet_phi/F");
    t->Branch("an_trig", &o_trig, "an_trig/I");
    t->Branch("wedge_rho", o_rho, Form("wedge_rho[%d]/F", kNW));
    t->Branch("wedge_cov", o_cov, Form("wedge_cov[%d]/F", kNW));

    std::ifstream in(argv[1]);
    std::string jf, sf;
    long njets = 0;
    while (in >> jf >> sf)
    {
        TFile fj(jf.c_str()), fs(sf.c_str());
        TTree *jt = (TTree *)fj.Get("jetTree"), *st = (TTree *)fs.Get("data");
        if (!jt || !st) { std::cerr << "skip " << jf << std::endl; continue; }
        int run = 0, fn = 0;
        std::sscanf(jf.substr(jf.find_last_of('/') + 1).c_str(), "jet_output_%d_%d.root", &run, &fn);
        std::vector<float> *eta = 0, *phi = 0, *E = 0, *pt = 0, *x = 0, *y = 0;
        Long64_t ste = -1;
        Double_t sqrt_s = 510;
        Int_t em[3] = {0}, emt[3] = {0};
        jt->SetBranchAddress("reco_eta", &eta); jt->SetBranchAddress("reco_phi", &phi);
        jt->SetBranchAddress("reco_E", &E); jt->SetBranchAddress("reco_pt", &pt);
        jt->SetBranchAddress("reco_x", &x); jt->SetBranchAddress("reco_y", &y);
        const bool has_entry = jt->GetBranch("simpletree_entry");
        if (has_entry) jt->SetBranchAddress("simpletree_entry", &ste);
        if (jt->GetBranch("event_sqrt_s_gev")) jt->SetBranchAddress("event_sqrt_s_gev", &sqrt_s);
        for (int l = 0; l < 3; l++)
        {
            if (jt->GetBranch(Form("fcsEM%d", l))) jt->SetBranchAddress(Form("fcsEM%d", l), &em[l]);
            if (jt->GetBranch(Form("fcsEM%d_tpc", l))) jt->SetBranchAddress(Form("fcsEM%d_tpc", l), &emt[l]);
        }
        Int_t npts = 0, detid[MAX_PTS];
        Float_t pe[MAX_PTS], px[MAX_PTS], py[MAX_PTS], pz[MAX_PTS];
        st->SetBranchAddress("Cal_npts", &npts); st->SetBranchAddress("Cal_pt_energy", pe);
        st->SetBranchAddress("Cal_pt_x", px); st->SetBranchAddress("Cal_pt_y", py); st->SetBranchAddress("Cal_pt_z", pz);
        st->SetBranchAddress("Cal_pt_detid", detid);
        std::vector<float> p_eta, p_phi, p_E;
        std::vector<int> p_tile;
        for (Long64_t ev = 0; ev < jt->GetEntries(); ev++)
        {
            jt->GetEntry(ev);
            std::vector<size_t> sel;
            for (size_t j = 0; j < eta->size(); j++)
                if ((*pt)[j] >= 2.0 && (*E)[j] >= 10 && (*E)[j] <= 120 &&
                    pass_fiducial_cut((*x)[j], (*y)[j], reco_fiducial_buffer_ecal, kFiducialRectEcal))
                    sel.push_back(j);
            if (sel.empty() || st->GetEntry(has_entry ? ste : ev) <= 0) continue;
            p_eta.clear(); p_phi.clear(); p_E.clear(); p_tile.clear();
            for (int i = 0; i < std::min((int)npts, MAX_PTS); i++)
            {
                if (detid[i] != 0 && detid[i] != 1) continue;
                if (!pass_point_energy_cut(pe[i], sqrt_s)) continue;
                const double rr = std::sqrt((double)px[i] * px[i] + (double)py[i] * py[i]);
                p_eta.push_back(std::asinh(pz[i] / rr)); p_phi.push_back(std::atan2(py[i], px[i]));
                p_E.push_back(pe[i]); p_tile.push_back(geo.TileOf(px[i], py[i]));
            }
            for (size_t s = 0; s < sel.size(); s++)
            {
                const size_t j = sel[s];
                o_run = run; o_E = (*E)[j]; o_pt = (*pt)[j]; o_eta = (*eta)[j]; o_phi = (*phi)[j];
                o_trig = 0;
                for (int l = 0; l < 3; l++) if ((em[l] || emt[l]) && (*pt)[j] > kTrigFactor * kThr[l]) o_trig = 1;
                AllWedges(geo, (*eta)[j], (*phi)[j], p_eta, p_phi, p_E, p_tile, o_rho, o_cov);
                t->Fill();
                njets++;
            }
        }
    }
    out.cd();
    t->Write();
    out.Close();
    std::printf("wedge_scan: %ld jets -> %s\n", njets, argv[2]);
    return 0;
}
