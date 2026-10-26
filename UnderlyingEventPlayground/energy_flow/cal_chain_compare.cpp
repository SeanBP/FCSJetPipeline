// ECal reconstruction chain, data against simulation: where the data's
// deficit of low-energy points (1-2 GeV; energy_flow.cpp) comes from.
//
// Reads SimpleTrees (data or simulation, one kind per list) and, for events
// with an EM trigger (fcsEM0/1/2; emulated in simulation), histograms the
// ECal (detid 0, 1) hits, clusters and points:
//   hits      energy, ADC sum, energy / ADC (the gain), towers above energy
//             thresholds per event
//   clusters  energy, towers per cluster, energy vs towers
//   points    energy; points per parent cluster; point energy split by the
//             number of points sharing its cluster (1, 2, 3+); point energy
//             per eta band
// Simulation (mc_pthatmin_gev present) is weighted per ptHat bin by sigma /
// (events in the bin); tables are per triggered event.
//
// Build and run from this folder:
//   g++ -m64 -O2 cal_chain_compare.cpp -o cal_chain_compare -I../../shared -I../../jet_scale \
//       `root-config --cflags --libs` `fastjet-config --cxxflags` -std=c++11
//   ./cal_chain_compare <SimpleTree list> <out.root> [seed_GeV=1.0]
//
// Seed emulation: P24ia (SL24a) seeds ECal clusters only at towers >= 1 GeV;
// the DEV library of the simulation at 0.05 GeV. Each cluster's highest
// tower is estimated as the highest ECal hit within 1.5 towers of the
// cluster centroid; clusters below seed_GeV, and their points (simulation:
// parentClusterId = cluster index within its detector), are dropped from the
// "seeded" point spectrum. Data clusters all pass by construction (check).

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "JetXSecWeights.h"

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::cerr << "usage: ./cal_chain_compare <SimpleTree list> <out.root> [seed_GeV]" << std::endl; return 1; }
    const double seed = argc > 3 ? std::atof(argv[3]) : 1.0;
    const double kTower = 5.542, kReach = 1.5 * kTower; // ECal tower size (cm), search half-width
    std::vector<std::string> files;
    {
        std::ifstream in(argv[1]);
        std::string a;
        while (in >> a) files.push_back(a);
    }

    // Simulation: sigma / events per ptHat bin
    std::vector<double> file_w(files.size(), 1.0);
    {
        typedef std::pair<float, float> Key;
        std::map<Key, std::pair<double, double> > bins;
        std::vector<Key> key(files.size(), Key(-1, -1));
        bool sim = false;
        for (size_t f = 0; f < files.size(); f++)
        {
            TFile fs(files[f].c_str());
            TTree *t = (TTree *)fs.Get("data");
            if (!t || t->GetEntries() == 0 || !t->GetBranch("mc_pthatmin_gev")) continue;
            Double_t sigma_job = -1; Int_t N_job = -1; Float_t ptmin = -1, ptmax = -1;
            t->SetBranchStatus("*", 0);
            const char *b[] = {"sigma_job", "N_job", "mc_pthatmin_gev", "mc_pthatmax_gev"};
            for (int i = 0; i < 4; i++) t->SetBranchStatus(b[i], 1);
            t->SetBranchAddress("sigma_job", &sigma_job); t->SetBranchAddress("N_job", &N_job);
            t->SetBranchAddress("mc_pthatmin_gev", &ptmin); t->SetBranchAddress("mc_pthatmax_gev", &ptmax);
            t->GetEntry(0);
            if (ptmax <= 0) continue; // data SimpleTrees carry the branch unset
            sim = true;
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
            for (size_t f = 0; f < files.size(); f++) file_w[f] = w.count(key[f]) && wmax > 0 ? w[key[f]] / wmax : 0;
            std::printf("cal_chain_compare: simulation, %zu ptHat bins\n", w.size());
        }
        else std::printf("cal_chain_compare: data\n");
    }

    TFile out(argv[2], "RECREATE");
    TH1D *h_hitE = new TH1D("h_hitE", ";ECal hit energy (GeV);hits", 500, 0, 5);
    TH1D *h_hitADC = new TH1D("h_hitADC", ";ECal hit ADC sum;hits", 400, 0, 400);
    TH1D *h_gain = new TH1D("h_gain", ";ECal hit energy / ADC sum (GeV);hits", 200, 0, 0.02);
    TH1D *h_clE = new TH1D("h_clE", ";ECal cluster energy (GeV);clusters", 200, 0, 50);
    TH1D *h_clNT = new TH1D("h_clNT", ";towers per cluster;clusters", 60, 0, 60);
    TH2D *h_clE_NT = new TH2D("h_clE_NT", ";cluster energy (GeV);towers", 50, 0, 10, 30, 0, 30);
    TH1D *h_ptE = new TH1D("h_ptE", ";ECal point energy (GeV);points", 200, 0, 50);
    TH1D *h_ptE_n[3];
    for (int k = 0; k < 3; k++) h_ptE_n[k] = new TH1D(Form("h_ptE_n%d", k + 1), ";ECal point energy (GeV);points", 200, 0, 50);
    TH1D *h_npc = new TH1D("h_npc", ";points per cluster;clusters", 10, 0, 10);
    const double eta_edges[4] = {2.4, 2.9, 3.3, 4.2};
    TH1D *h_ptE_eta[3];
    for (int k = 0; k < 3; k++) h_ptE_eta[k] = new TH1D(Form("h_ptE_eta%d", k), ";ECal point energy (GeV);points", 200, 0, 50);
    TH1D *h_ptE_seed = new TH1D("h_ptE_seed", ";ECal point energy (GeV), seed emulated;points", 200, 0, 50);
    TH1D *h_clE_seed = new TH1D("h_clE_seed", ";ECal cluster energy (GeV), seed emulated;clusters", 200, 0, 50);
    TH1D *h_clMax = new TH1D("h_clMax", ";highest hit near the cluster (GeV);clusters", 100, 0, 10);
    double n_ev = 0;
    long n_raw = 0;

    for (size_t f = 0; f < files.size(); f++)
    {
        if (file_w[f] <= 0) continue;
        TFile *fs = TFile::Open(files[f].c_str());
        TTree *t = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        if (!t) { delete fs; continue; }
        static Int_t nh, det[5000], adc[5000], nc, cdet[100], cnt[100], np, pdet[200], ppar[200];
        static Float_t hE[5000], hx[5000], hy[5000], cE[100], cx[100], cy[100], pE[200], px[200], py[200], pz[200];
        Int_t em0 = 0, em1 = 0, em2 = 0;
        t->SetBranchStatus("*", 0);
        const char *b[] = {"Cal_nhits", "Cal_detid", "Cal_adcsum", "Cal_hit_energy", "Cal_nclus", "Cal_clus_detid", "Cal_clus_ntowers",
                           "Cal_clus_energy", "Cal_npts", "Cal_pt_detid", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z",
                           "Cal_pt_parentClusterId", "fcsEM0", "fcsEM1", "fcsEM2", "Cal_hit_posx", "Cal_hit_posy", "Cal_clus_x", "Cal_clus_y"};
        for (int i = 0; i < 22; i++) t->SetBranchStatus(b[i], 1);
        t->SetBranchAddress("Cal_hit_posx", hx); t->SetBranchAddress("Cal_hit_posy", hy);
        t->SetBranchAddress("Cal_clus_x", cx); t->SetBranchAddress("Cal_clus_y", cy);
        t->SetBranchAddress("Cal_nhits", &nh); t->SetBranchAddress("Cal_detid", det); t->SetBranchAddress("Cal_adcsum", adc);
        t->SetBranchAddress("Cal_hit_energy", hE); t->SetBranchAddress("Cal_nclus", &nc); t->SetBranchAddress("Cal_clus_detid", cdet);
        t->SetBranchAddress("Cal_clus_ntowers", cnt); t->SetBranchAddress("Cal_clus_energy", cE); t->SetBranchAddress("Cal_npts", &np);
        t->SetBranchAddress("Cal_pt_detid", pdet); t->SetBranchAddress("Cal_pt_energy", pE); t->SetBranchAddress("Cal_pt_x", px);
        t->SetBranchAddress("Cal_pt_y", py); t->SetBranchAddress("Cal_pt_z", pz); t->SetBranchAddress("Cal_pt_parentClusterId", ppar);
        t->SetBranchAddress("fcsEM0", &em0); t->SetBranchAddress("fcsEM1", &em1); t->SetBranchAddress("fcsEM2", &em2);
        const double w = file_w[f];
        for (Long64_t e = 0; e < t->GetEntries(); e++)
        {
            if (t->GetEntry(e) <= 0 || !(em0 || em1 || em2)) continue;
            n_ev += w; n_raw++;
            for (int i = 0; i < std::min(nh, 5000); i++)
            {
                if (det[i] != 0 && det[i] != 1) continue;
                h_hitE->Fill(hE[i], w); h_hitADC->Fill(adc[i], w);
                if (adc[i] > 0) h_gain->Fill(hE[i] / adc[i], w);
            }
            for (int i = 0; i < std::min(nc, 100); i++)
            {
                if (cdet[i] != 0 && cdet[i] != 1) continue;
                h_clE->Fill(cE[i], w); h_clNT->Fill(cnt[i], w); h_clE_NT->Fill(cE[i], cnt[i], w);
            }
            // Seed emulation: clusters' highest nearby hit; index within detector
            std::map<int, bool> seeded;
            int idx_in_det[2] = {0, 0};
            for (int i = 0; i < std::min(nc, 100); i++)
            {
                if (cdet[i] != 0 && cdet[i] != 1) continue;
                double hmax = 0;
                for (int k = 0; k < std::min(nh, 5000); k++)
                    if (det[k] == cdet[i] && std::fabs(hx[k] - cx[i]) <= kReach && std::fabs(hy[k] - cy[i]) <= kReach) hmax = std::max(hmax, (double)hE[k]);
                h_clMax->Fill(hmax, w);
                const bool ok = hmax >= seed;
                seeded[cdet[i] * 1000 + idx_in_det[cdet[i]]++] = ok;
                if (ok) h_clE_seed->Fill(cE[i], w);
            }
            std::map<int, int> per_cluster;
            for (int i = 0; i < std::min(np, 200); i++)
                if (pdet[i] == 0 || pdet[i] == 1) per_cluster[pdet[i] * 1000 + ppar[i]]++;
            for (std::map<int, int>::iterator it = per_cluster.begin(); it != per_cluster.end(); ++it) h_npc->Fill(it->second, w);
            for (int i = 0; i < std::min(np, 200); i++)
            {
                if (pdet[i] != 0 && pdet[i] != 1) continue;
                h_ptE->Fill(pE[i], w);
                {
                    std::map<int, bool>::iterator it = seeded.find(pdet[i] * 1000 + ppar[i]);
                    if (it == seeded.end() || it->second) h_ptE_seed->Fill(pE[i], w); // unknown parent (data): kept
                }
                h_ptE_n[std::min(per_cluster[pdet[i] * 1000 + ppar[i]], 3) - 1]->Fill(pE[i], w);
                const double rho = std::sqrt((double)px[i] * px[i] + (double)py[i] * py[i]);
                const double eta = std::asinh(pz[i] / rho);
                for (int k = 0; k < 3; k++) if (eta >= eta_edges[k] && eta < eta_edges[k + 1]) h_ptE_eta[k]->Fill(pE[i], w);
            }
        }
        fs->Close(); delete fs;
    }
    std::printf("cal_chain_compare: %zu files, %ld triggered events\n", files.size(), n_raw);

    // Per triggered event, in energy ranges
    struct R { static double In(TH1 *h, double lo, double hi) { return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6)); } };
    const double he[] = {0, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 5};
    std::printf("  hits per event by energy (GeV):");
    for (int k = 0; k < 8; k++) std::printf(" %g-%g:%.2f", he[k], he[k + 1], R::In(h_hitE, he[k], he[k + 1]) / n_ev);
    std::printf("\n  lowest hit energy with entries: ");
    for (int b = 1; b <= h_hitE->GetNbinsX(); b++) if (h_hitE->GetBinContent(b) > 0) { std::printf("%.3f GeV", h_hitE->GetBinLowEdge(b)); break; }
    const double ae[] = {0, 5, 10, 15, 20, 30, 50, 100, 400};
    std::printf("\n  hits per event by ADC sum:");
    for (int k = 0; k < 8; k++) std::printf(" %g-%g:%.2f", ae[k], ae[k + 1], R::In(h_hitADC, ae[k], ae[k + 1]) / n_ev);
    std::printf("\n  hit gain E/ADC: mean %.5f rms %.5f\n", h_gain->GetMean(), h_gain->GetRMS());
    const double ce[] = {0, 0.5, 1, 1.5, 2, 3, 4, 6, 10, 20, 50};
    std::printf("  clusters per event by energy:");
    for (int k = 0; k < 10; k++) std::printf(" %g-%g:%.3f", ce[k], ce[k + 1], R::In(h_clE, ce[k], ce[k + 1]) / n_ev);
    std::printf("\n  towers per cluster: mean %.2f;", h_clNT->GetMean());
    for (int k = 0; k < 6; k++)
    {
        TH1D *pr = h_clE_NT->ProjectionY("pr", h_clE_NT->GetXaxis()->FindBin(ce[k + 1] + 1e-6), h_clE_NT->GetXaxis()->FindBin(ce[k + 2] - 1e-6));
        std::printf(" E %g-%g: %.2f", ce[k + 1], ce[k + 2], pr->GetMean());
        delete pr;
    }
    std::printf("\n  points per cluster: mean %.3f; fractions 1:%.3f 2:%.3f 3+:%.3f\n", h_npc->GetMean(),
                h_npc->GetBinContent(2) / h_npc->Integral(), h_npc->GetBinContent(3) / h_npc->Integral(), h_npc->Integral(4, 11) / h_npc->Integral());
    std::printf("  points per event by energy (all | cluster with 1 | 2 | 3+ points):\n");
    for (int k = 0; k < 10; k++)
        std::printf("    %4g-%-4g %.3f | %.3f %.3f %.3f\n", ce[k], ce[k + 1], R::In(h_ptE, ce[k], ce[k + 1]) / n_ev,
                    R::In(h_ptE_n[0], ce[k], ce[k + 1]) / n_ev, R::In(h_ptE_n[1], ce[k], ce[k + 1]) / n_ev, R::In(h_ptE_n[2], ce[k], ce[k + 1]) / n_ev);
    std::printf("  clusters with highest nearby hit below %g GeV: %.3f of all\n", seed, h_clMax->Integral(1, h_clMax->FindBin(seed - 1e-6)) / h_clMax->Integral());
    std::printf("  seed emulated (%g GeV), per event by energy: points | clusters\n", seed);
    for (int k = 0; k < 10; k++)
        std::printf("    %4g-%-4g %.3f | %.3f\n", ce[k], ce[k + 1], R::In(h_ptE_seed, ce[k], ce[k + 1]) / n_ev, R::In(h_clE_seed, ce[k], ce[k + 1]) / n_ev);
    std::printf("  points per event by energy in eta bands:\n");
    for (int j = 0; j < 3; j++)
    {
        std::printf("    eta %.1f-%.1f:", eta_edges[j], eta_edges[j + 1]);
        for (int k = 1; k < 8; k++) std::printf(" %g-%g:%.3f", ce[k], ce[k + 1], R::In(h_ptE_eta[j], ce[k], ce[k + 1]) / n_ev);
        std::printf("\n");
    }
    out.Write();
    out.Close();
    return 0;
}
