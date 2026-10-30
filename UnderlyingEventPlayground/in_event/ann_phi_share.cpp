// Is the low-area closure bias of the whole-ring annulus an away-side bias?
// For each jet, recompute the ring's weights from the geometry (no other-jet
// exclusion, as an_*[3]) and the share of the weighted live area at
// |dphi| > 120 deg; bin the simulation closure (an_tl[3] / an_tf[3]) and the
// data mean estimate in that share. 20-120 GeV jets.
//
// Build and run from this folder:
//   g++ -O2 ann_phi_share.cpp -o ann_phi_share -I../../shared -I../../underlying_event `root-config --cflags --libs` -std=c++11
//   ./ann_phi_share "<data files>" "<sim files>"

#include <TROOT.h>
#include <TChain.h>

#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "in_event_geo.h"

static const int kNAV = 6, kNS_ = 6;
static const double kShEdge[kNS_ + 1] = {0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.01};
static const double kMaxFill = 1 / 0.7;

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}

// share of the ring's weighted live area at |dphi| > 120 deg
static double AwayShare(const Geometry &g, double eta, double phi)
{
    std::vector<int> cl;
    g.Near(eta, phi, R, cl);
    double Ajs[kNS] = {0}, Tjs[kNS] = {0}, T = 0;
    for (size_t q = 0; q < cl.size(); q++)
    {
        const Cell &c = g.cells[cl[q]];
        const int s = StripOf(c.eta - eta);
        if (s >= 0) { Ajs[s] += c.area; Tjs[s] += c.area; T += c.area; }
    }
    static std::vector<double> rho(g.cells.size(), 1.0);
    static const std::vector<std::pair<float, float> > none;
    Annulus an;
    BuildAnnulus(g, rho, eta, phi, none, Ajs, T, Tjs, 0, an);
    double all = 0, away = 0;
    for (size_t c = 0; c < g.cells.size(); c++)
    {
        if (!(an.w[c] > 0)) continue;
        const double x = an.w[c] * g.cells[c].area;
        all += x;
        if (std::fabs(std::remainder(g.cells[c].phi - phi, 2 * M_PI)) > 2 * M_PI / 3) away += x;
    }
    return all > 0 ? away / all : -1;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::fprintf(stderr, "usage: ./ann_phi_share \"<data files>\" \"<sim files>\"\n"); return 1; }
    Geometry geo; geo.Build();
    double n[2][kNS_] = {}, ntot[2] = {}, cx[kNS_] = {}, cy[kNS_] = {}, e[kNS_] = {}, eall = 0, nall = 0;
    for (int sim = 0; sim < 2; sim++)
    {
        TChain ch("jets");
        std::vector<std::string> f = Split(argv[1 + sim]);
        for (size_t i = 0; i < f.size(); i++) ch.Add(f[i].c_str());
        Float_t E, w, eta, phi, an_e[kNAV], an_fill[kNAV], an_tl[kNAV], an_tf[kNAV];
        ch.SetBranchAddress("jet_E", &E); ch.SetBranchAddress("weight", &w); ch.SetBranchAddress("jet_eta", &eta); ch.SetBranchAddress("jet_phi", &phi);
        ch.SetBranchAddress("an_e", an_e); ch.SetBranchAddress("an_fill", an_fill); ch.SetBranchAddress("an_tl", an_tl); ch.SetBranchAddress("an_tf", an_tf);
        const Long64_t N = sim ? ch.GetEntries() : std::min(ch.GetEntries(), (Long64_t)30000);
        for (Long64_t i = 0; i < N; i++)
        {
            ch.GetEntry(i);
            if (E < 20 || E >= 120 || !(an_fill[3] > 0 && an_fill[3] <= kMaxFill)) continue;
            const double sh = AwayShare(geo, eta, phi);
            int b = -1;
            for (int q = 0; q < kNS_; q++) if (sh >= kShEdge[q] && sh < kShEdge[q + 1]) b = q;
            if (b < 0) continue;
            const double ww = sim ? w : 1.0;
            n[sim][b] += ww; ntot[sim] += ww;
            if (sim) { cx[b] += ww * an_tl[3]; cy[b] += ww * an_tf[3]; }
            else { e[b] += an_e[3]; eall += an_e[3]; nall += 1; }
        }
    }
    std::printf("whole ring, no other-jet exclusion, 20-120 GeV (data: first 30k jets)\n");
    std::printf(" away share   jets data  sim   sim closure   data mean / all\n");
    for (int q = 0; q < kNS_; q++)
        std::printf("  %.1f-%.1f     %5.3f %5.3f    %6.3f        %6.3f\n", kShEdge[q], std::min(kShEdge[q + 1], 1.0), n[0][q] / ntot[0], n[1][q] / ntot[1],
                    cy[q] > 0 ? cx[q] / cy[q] : NAN, n[0][q] > 0 ? (e[q] / n[0][q]) / (eall / nall) : NAN);
    return 0;
}
