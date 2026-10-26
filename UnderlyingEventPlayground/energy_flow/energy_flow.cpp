// Energy flow around a reco jet, data and simulation alike: the ECal point
// energy density (GeV per unit eta-phi area) as a function of the distance
// dR and of the azimuthal angle |dphi| from the jet. The proxy-event method
// samples this density, in other events, at |dphi| ~ 90 deg (no separation)
// or further out (separation; ~115 deg on average), so it shows directly
// where data and simulation differ.
//
// Reads pairs "<jet_output> <SimpleTree>" (as underlying_event/
// UnderlyingEvent.cpp; same point selection and directions). For each event
// with an EM trigger (fcsEM0/1/2; emulated in simulation) and each fiducial
// reco jet (pT >= reco_jet_pt_min_em), every point of the event is binned in
// dR and, for dR >= proxy_min_dr (the separation, 2.5R), in |dphi|. The area
// over which a point could have been found in each bin is the fiducial area
// (ProxyFiducialCells) in that bin around the jet, summed per jet, so the
// ratio is a density independent of the detector shape. Points are not
// removed for belonging to other jets (as the no-veto proxy cones); only
// points whose direction is inside the fiducial region count.
//
// Simulation (files ..._pt<ptHat>_<job>.root with sigma_job, N_job) is
// weighted per ptHat bin by sigma / (events in the bin), so the bins mix as
// in nature; data events have weight 1.
//
// Output: histograms per jet pT bin (all EM levels, and per EM level)
//   h_E_dr_<pt>[_em<l>], h_A_dr_<pt>[...]       point energy, area (per jet sum)
//   h_E_dphi_<pt>[...],  h_A_dphi_<pt>[...]     the same in |dphi|, dR >= 2.5R
//   h_J_dr_<pt>[...]                            other fiducial reco jets in dR
//   h_Pt_ring_<pt>[...]                         point energies, 0.8 <= dR < 1.6
//   h_njet_<pt>[...]                            jets (weighted)
//   h_Pt_all_em<l>                              energies of all fiducial points
// and prints the density tables.
//
// Build and run from this folder:
//   g++ -m64 -O2 energy_flow.cpp ../../shared/VectorDict.cxx -o energy_flow -I../../shared -I../../underlying_event -I../../jet_scale \
//       `root-config --cflags --libs` -std=c++11
//   ./energy_flow <pairs.list> <out.root> [all]      all: also untriggered events

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "ProxySelection.h"
#include "JetXSecWeights.h"

static const int kNPt = 4;
static const double kPtEdges[kNPt + 1] = {2, 3, 4, 6, 10};
static const int kNLev = 4; // all, em0, em1, em2

struct Hists
{
    TH1D *E_dr, *A_dr, *E_dphi, *A_dphi, *njet, *J_dr, *Pt_ring;
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::cerr << "usage: ./energy_flow <pairs.list> <out.root> [all]" << std::endl; return 1; }
    const bool all_events = argc > 3 && std::string(argv[3]) == "all";

    std::vector<std::pair<std::string, std::string> > pairs;
    {
        std::ifstream in(argv[1]);
        std::string a, b;
        while (in >> a >> b) pairs.push_back(std::make_pair(a, b));
    }
    if (pairs.empty()) { std::cerr << "ERROR: no pairs in " << argv[1] << std::endl; return 1; }

    // Simulation: sigma / events per ptHat bin, from the files' sigma_job and N_job
    std::vector<double> file_w(pairs.size(), 1.0);
    {
        typedef std::pair<float, float> Key;
        std::map<Key, std::pair<double, double> > bins; // sum N, sum N sigma
        std::vector<Key> key(pairs.size(), Key(-1, -1));
        bool sim = false;
        for (size_t f = 0; f < pairs.size(); f++)
        {
            TFile fj(pairs[f].first.c_str());
            TTree *jt = (TTree *)fj.Get("jetTree");
            if (!jt || jt->GetEntries() == 0 || !jt->GetBranch("mc_pthatmin_gev")) continue;
            sim = true;
            Double_t sigma_job = -1; Int_t N_job = -1; Float_t ptmin = -1, ptmax = -1;
            jt->SetBranchStatus("*", 0);
            const char *b[] = {"sigma_job", "N_job", "mc_pthatmin_gev", "mc_pthatmax_gev"};
            for (int i = 0; i < 4; i++) jt->SetBranchStatus(b[i], 1);
            jt->SetBranchAddress("sigma_job", &sigma_job); jt->SetBranchAddress("N_job", &N_job);
            jt->SetBranchAddress("mc_pthatmin_gev", &ptmin); jt->SetBranchAddress("mc_pthatmax_gev", &ptmax);
            jt->GetEntry(0);
            const double sg = DecodePythia6Sigma(sigma_job);
            key[f] = Key(ptmin, ptmax);
            if (sg > 0 && N_job > 0) { bins[key[f]].first += N_job; bins[key[f]].second += N_job * sg; }
        }
        if (sim)
        {
            std::map<Key, double> w;
            double wmax = 0;
            for (std::map<Key, std::pair<double, double> >::iterator it = bins.begin(); it != bins.end(); ++it)
            {
                w[it->first] = it->second.first > 0 ? it->second.second / it->second.first / it->second.first : 0;
                wmax = std::max(wmax, w[it->first]);
            }
            for (size_t f = 0; f < pairs.size(); f++) file_w[f] = w.count(key[f]) && wmax > 0 ? w[key[f]] / wmax : 0;
            std::cout << "energy_flow: simulation, weight per ptHat bin (relative):";
            for (std::map<Key, double>::iterator it = w.begin(); it != w.end(); ++it)
                std::cout << " [" << it->first.first << "," << it->first.second << ") " << it->second / wmax;
            std::cout << std::endl;
        }
    }

    TFile out(argv[2], "RECREATE");
    Hists h[kNLev][kNPt];
    for (int l = 0; l < kNLev; l++)
        for (int p = 0; p < kNPt; p++)
        {
            const std::string s = Form("_%d%s", p, l == 0 ? "" : Form("_em%d", l - 1));
            h[l][p].E_dr = new TH1D(("h_E_dr" + s).c_str(), ";#DeltaR to the jet;point energy (GeV)", 40, 0, 4);
            h[l][p].A_dr = new TH1D(("h_A_dr" + s).c_str(), ";#DeltaR to the jet;fiducial area", 40, 0, 4);
            h[l][p].E_dphi = new TH1D(("h_E_dphi" + s).c_str(), ";|#Delta#phi| to the jet (deg);point energy (GeV)", 18, 0, 180);
            h[l][p].A_dphi = new TH1D(("h_A_dphi" + s).c_str(), ";|#Delta#phi| to the jet (deg);fiducial area", 18, 0, 180);
            h[l][p].J_dr = new TH1D(("h_J_dr" + s).c_str(), ";#DeltaR to the jet;other reco jets (p_{T} >= 2 GeV)", 40, 0, 4);
            h[l][p].Pt_ring = new TH1D(("h_Pt_ring" + s).c_str(), ";point energy (GeV), 0.8 <= #DeltaR < 1.6;points", 40, 0, 20);
            h[l][p].njet = new TH1D(("h_njet" + s).c_str(), ";;jets", 1, 0, 1);
        }
    const std::vector<ProxyCell> &cells = ProxyFiducialCells();
    // Every fiducial point of the triggered events with a jet, per EM level
    TH1D *h_pt_all[3];
    for (int l = 0; l < 3; l++) h_pt_all[l] = new TH1D(Form("h_Pt_all_em%d", l), ";point energy (GeV);points", 200, 0, 50);
    double n_ev_w[3] = {0, 0, 0};

    const int MAX_PTS = 200;
    long n_ev = 0, n_jets = 0;
    for (size_t f = 0; f < pairs.size(); f++)
    {
        if (file_w[f] <= 0) continue;
        TFile *fj = TFile::Open(pairs[f].first.c_str());
        TFile *fs = TFile::Open(pairs[f].second.c_str());
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        const bool has_entry = jt && jt->GetBranch("simpletree_entry");
        if (!jt || !st || !st->GetBranch("Cal_npts") || (!has_entry && jt->GetEntries() != st->GetEntries()))
        {
            std::cerr << "skipping " << pairs[f].first << std::endl;
            delete fj; delete fs; continue;
        }
        std::vector<float> *reco_pt = 0, *reco_eta = 0, *reco_phi = 0, *reco_x = 0, *reco_y = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        Long64_t simpletree_entry = -1;
        Double_t sqrt_s = kSqrtS_MainRun22_GeV;
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int i = 0; i < 8; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &reco_pt); jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi); jt->SetBranchAddress("reco_x", &reco_x);
        jt->SetBranchAddress("reco_y", &reco_y);
        jt->SetBranchAddress("fcsEM0", &em0); jt->SetBranchAddress("fcsEM1", &em1); jt->SetBranchAddress("fcsEM2", &em2);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &simpletree_entry); }
        if (jt->GetBranch("event_sqrt_s_gev")) { jt->SetBranchStatus("event_sqrt_s_gev", 1); jt->SetBranchAddress("event_sqrt_s_gev", &sqrt_s); }
        Int_t Cal_npts = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        Int_t Cal_pt_detid[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pb[i], 1);
        st->SetBranchAddress("Cal_npts", &Cal_npts); st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x); st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z); st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);

        const double w = file_w[f];
        std::vector<float> p_eta, p_phi, p_E;
        for (Long64_t ev = 0; ev < jt->GetEntries(); ev++)
        {
            if (jt->GetEntry(ev) <= 0 || reco_pt->empty()) continue;
            const int em = ProxyEmLevel(em0, em1, em2);
            if (em < 0 && !all_events) continue;
            if (st->GetEntry(has_entry ? simpletree_entry : ev) <= 0) continue;
            n_ev++;
            if (em >= 0) n_ev_w[em] += w;
            p_eta.clear(); p_phi.clear(); p_E.clear();
            for (int i = 0; i < std::min((int)Cal_npts, MAX_PTS); i++)
            {
                if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                if (!pass_point_energy_cut(Cal_pt_energy[i], sqrt_s)) continue;
                // only points whose direction falls in the fiducial region, the
                // area the densities are normalised to (positions at z_proj)
                if (!(Cal_pt_z[i] > 0) || !ProxyInFiducial(Cal_pt_x[i] * z_proj / Cal_pt_z[i], Cal_pt_y[i] * z_proj / Cal_pt_z[i])) continue;
                const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                p_eta.push_back(std::asinh(Cal_pt_z[i] / rho)); p_phi.push_back(std::atan2(Cal_pt_y[i], Cal_pt_x[i]));
                p_E.push_back(Cal_pt_energy[i]);
                if (em >= 0) h_pt_all[em]->Fill(Cal_pt_energy[i], w);
            }
            for (size_t j = 0; j < reco_pt->size(); j++)
            {
                const double pt = reco_pt->at(j);
                if (pt < reco_jet_pt_min_em || !ProxyInFiducial(reco_x->at(j), reco_y->at(j))) continue;
                int p = -1;
                for (int k = 0; k < kNPt; k++) if (pt >= kPtEdges[k] && pt < kPtEdges[k + 1]) p = k;
                if (p < 0) continue;
                n_jets++;
                const double eta = reco_eta->at(j), phi = reco_phi->at(j);
                const int lev[2] = {0, em + 1};
                for (int li = 0; li < (em >= 0 ? 2 : 1); li++)
                {
                    Hists &hh = h[lev[li]][p];
                    hh.njet->Fill(0.5, w);
                    for (size_t c = 0; c < cells.size(); c++)
                    {
                        const double dr = ProxyDeltaR(cells[c].eta, cells[c].phi, eta, phi);
                        hh.A_dr->Fill(dr, w * cells[c].area);
                        if (dr >= proxy_min_dr)
                            hh.A_dphi->Fill(std::fabs(std::remainder(cells[c].phi - phi, 2.0 * M_PI)) * 180.0 / M_PI, w * cells[c].area);
                    }
                    for (size_t o = 0; o < reco_pt->size(); o++)
                        if (o != j && reco_pt->at(o) >= reco_jet_pt_min_em && ProxyInFiducial(reco_x->at(o), reco_y->at(o)))
                            hh.J_dr->Fill(ProxyDeltaR(reco_eta->at(o), reco_phi->at(o), eta, phi), w);
                    for (size_t i = 0; i < p_E.size(); i++)
                    {
                        const double dr = ProxyDeltaR(p_eta[i], p_phi[i], eta, phi);
                        hh.E_dr->Fill(dr, w * p_E[i]);
                        if (dr >= 0.8 && dr < 1.6) hh.Pt_ring->Fill(p_E[i], w);
                        if (dr >= proxy_min_dr)
                            hh.E_dphi->Fill(std::fabs(std::remainder(p_phi[i] - phi, 2.0 * M_PI)) * 180.0 / M_PI, w * p_E[i]);
                    }
                }
            }
        }
        fj->Close(); delete fj;
        fs->Close(); delete fs;
    }
    std::printf("energy_flow: %zu file pairs, %ld events, %ld jets\n", pairs.size(), n_ev, n_jets);

    // Density tables, all EM levels
    for (int p = 0; p < kNPt; p++)
    {
        Hists &hh = h[0][p];
        std::printf("  jet pT %.0f-%.0f GeV (%.0f jets weighted): density GeV / unit area\n    dphi (deg, dR>=2.5R):", kPtEdges[p], kPtEdges[p + 1],
                    hh.njet->GetBinContent(1));
        for (int b = 1; b <= hh.E_dphi->GetNbinsX(); b++)
            if (hh.A_dphi->GetBinContent(b) > 0) std::printf(" %.0f:%.3f", hh.E_dphi->GetBinCenter(b), hh.E_dphi->GetBinContent(b) / hh.A_dphi->GetBinContent(b));
        std::printf("\n    dR:");
        for (int b = 1; b <= hh.E_dr->GetNbinsX(); b += 2)
        {
            const double e = hh.E_dr->GetBinContent(b) + hh.E_dr->GetBinContent(b + 1), a = hh.A_dr->GetBinContent(b) + hh.A_dr->GetBinContent(b + 1);
            if (a > 0) std::printf(" %.1f:%.3f", hh.E_dr->GetBinLowEdge(b + 1), e / a);
        }
        std::printf("\n    other jets per unit area, dR:");
        for (int b = 1; b <= hh.J_dr->GetNbinsX(); b += 2)
        {
            const double n = hh.J_dr->GetBinContent(b) + hh.J_dr->GetBinContent(b + 1), a = hh.A_dr->GetBinContent(b) + hh.A_dr->GetBinContent(b + 1);
            if (a > 0) std::printf(" %.1f:%.3f", hh.J_dr->GetBinLowEdge(b + 1), n / a);
        }
        std::printf("\n    points per jet in 0.8<=dR<1.6 by energy (GeV):");
        const double nj = hh.njet->GetBinContent(1);
        const double pe[] = {0, 0.5, 1, 2, 4, 20};
        for (int k = 0; k < 5; k++)
            std::printf(" %.1f-%.0f:%.3f", pe[k], pe[k + 1], nj > 0 ? hh.Pt_ring->Integral(hh.Pt_ring->FindBin(pe[k] + 1e-6), hh.Pt_ring->FindBin(pe[k + 1] - 1e-6)) / nj : 0.);
        std::printf("\n");
    }
    std::printf("  fiducial points per event, by point energy (GeV):\n");
    const double pe[] = {1, 1.5, 2, 3, 4, 6, 10, 20, 50};
    for (int l = 0; l < 3; l++)
    {
        std::printf("    EM%d:", l);
        for (int k = 0; k < 8; k++)
            std::printf(" %g-%g:%.3f", pe[k], pe[k + 1], n_ev_w[l] > 0 ? h_pt_all[l]->Integral(h_pt_all[l]->FindBin(pe[k] + 1e-6), h_pt_all[l]->FindBin(pe[k + 1] - 1e-6)) / n_ev_w[l] : 0.);
        std::printf("\n");
    }
    out.Write();
    out.Close();
    return 0;
}
