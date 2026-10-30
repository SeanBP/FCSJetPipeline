// Smoke test of the pixel ring's point weighting (in_event_geo.h): the same
// rings, the same points, weighted per cell (pr.w[cell], the cell's
// area-averaged ratio) and per point (pr.Weight: the ratio of the point's own
// eta bin). The denominator, the ring's live area weighted per eta bin, is
// common, so empty ring regions count the same in both. Reco jets as in
// in_event_ue (pT >= 2 GeV, E 10-120 GeV, fiducial), every event; ring at 2R,
// no other-jet exclusion, uncut and with |dphi| <= 150 deg. Prints the mean
// densities, their ratio and the per-jet spread of the difference, overall
// and in jet eta bins (the beam-hole edge is where they can differ).
//
// Build and run from this folder:
//   g++ -m64 -O2 px_weight_check.cpp ../../shared/VectorDict.cxx -o px_weight_check -I../../shared -I../../underlying_event `root-config --cflags --libs` -std=c++11
//   ./px_weight_check <pairs list> [max files = 10]

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "in_event_geo.h"

static const int MAX_PTS = 2000;
static const int kNE = 4;
static const double kEtaEdges[kNE + 1] = {2.5, 2.9, 3.2, 3.5, 4.2};

struct Acc
{
    double n = 0, sc = 0, sp = 0, sd = 0, sdd = 0;
    void Add(double c, double p) { n++; sc += c; sp += p; sd += p - c; sdd += (p - c) * (p - c); }
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2) { std::fprintf(stderr, "usage: ./px_weight_check <pairs list> [max files]\n"); return 1; }
    const int max_files = argc > 2 ? std::atoi(argv[2]) : 10;
    Geometry geo; geo.Build();
    const std::vector<std::pair<float, float> > none;
    Acc acc[2][kNE + 1]; // [uncut, cut][eta bins, all]
    std::ifstream in(argv[1]); std::string jf, sf;
    int nf = 0;
    PixelRing pr;
    while (in >> jf >> sf && nf < max_files)
    {
        nf++;
        TFile fj(jf.c_str()), fs(sf.c_str());
        TTree *jt = (TTree *)fj.Get("jetTree"), *st = 0;
        TIter it(fs.GetListOfKeys()); TObject *k;
        while ((k = it())) { TObject *o = fs.Get(k->GetName()); if (o && o->InheritsFrom("TTree")) { st = (TTree *)o; break; } }
        if (!jt || !st) continue;
        std::vector<float> *pt = 0, *eta = 0, *phi = 0, *E = 0, *x = 0, *y = 0;
        Long64_t sent = -1;
        const bool has_entry = jt->GetBranch("simpletree_entry");
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_E", "reco_x", "reco_y"};
        for (int i = 0; i < 6; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &pt); jt->SetBranchAddress("reco_eta", &eta); jt->SetBranchAddress("reco_phi", &phi);
        jt->SetBranchAddress("reco_E", &E); jt->SetBranchAddress("reco_x", &x); jt->SetBranchAddress("reco_y", &y);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &sent); }
        Int_t np = 0, det[MAX_PTS];
        Float_t pE[MAX_PTS], px[MAX_PTS], py[MAX_PTS], pz[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pb[i], 1);
        st->SetBranchAddress("Cal_npts", &np); st->SetBranchAddress("Cal_pt_energy", pE); st->SetBranchAddress("Cal_pt_x", px);
        st->SetBranchAddress("Cal_pt_y", py); st->SetBranchAddress("Cal_pt_z", pz); st->SetBranchAddress("Cal_pt_detid", det);
        for (Long64_t ev = 0; ev < jt->GetEntries(); ev++)
        {
            if (jt->GetEntry(ev) <= 0 || !pt || pt->empty()) continue;
            std::vector<int> sel;
            for (size_t j = 0; j < pt->size(); j++)
                if (pt->at(j) >= reco_jet_pt_min_em && E->at(j) >= 10 && E->at(j) < 120 && ProxyInFiducial(x->at(j), y->at(j))) sel.push_back((int)j);
            if (sel.empty() || st->GetEntry(has_entry ? sent : ev) <= 0) continue;
            std::vector<float> te, tE; std::vector<int> tc;
            for (int i = 0; i < std::min((int)np, MAX_PTS); i++)
            {
                if (det[i] > 1 || !pass_point_energy_cut(pE[i], kSqrtS_MainRun22_GeV)) continue;
                const double r = std::hypot(px[i], py[i]);
                te.push_back(std::asinh(pz[i] / r)); tE.push_back(pE[i]); tc.push_back(geo.CellOf(px[i], py[i]));
            }
            for (size_t q = 0; q < sel.size(); q++)
            {
                const int j = sel[q];
                int eb = -1;
                for (int b = 0; b < kNE; b++) if (eta->at(j) >= kEtaEdges[b] && eta->at(j) < kEtaEdges[b + 1]) eb = b;
                for (int c = 0; c < 2; c++)
                {
                    BuildPixelRing(geo, eta->at(j), phi->at(j), none, kAnnJetExcl, pr, c ? 150 * M_PI / 180 + 1e-9 : M_PI);
                    if (!(pr.den > 0)) continue;
                    double ec = 0, ep = 0;
                    for (size_t i = 0; i < tE.size(); i++)
                    {
                        if (tc[i] < 0) continue;
                        ec += tE[i] * pr.w[tc[i]];
                        ep += tE[i] * pr.Weight(tc[i], te[i]);
                    }
                    acc[c][kNE].Add(ec / pr.den, ep / pr.den);
                    if (eb >= 0) acc[c][eb].Add(ec / pr.den, ep / pr.den);
                }
            }
        }
    }
    std::printf("px_weight_check: %d files\n", nf);
    for (int c = 0; c < 2; c++)
    {
        std::printf("%s\n  jet eta        jets     rho per cell   rho per point   point/cell   per-jet (point - cell): mean    RMS     RMS / mean rho\n", c ? "ring |dphi| <= 150 deg" : "ring uncut");
        for (int b = 0; b <= kNE; b++)
        {
            const Acc &a = acc[c][b];
            if (a.n < 1) continue;
            const double mc = a.sc / a.n, mp = a.sp / a.n, md = a.sd / a.n, rd = std::sqrt(std::max(0.0, a.sdd / a.n - md * md));
            if (b < kNE) std::printf("  %.1f-%.1f   ", kEtaEdges[b], kEtaEdges[b + 1]); else std::printf("  all       ");
            std::printf("%7.0f     %.4f         %.4f          %.4f              %+.4f   %.4f   %.3f\n", a.n, mc, mp, mp / mc, md, rd, rd / mc);
        }
    }
    return 0;
}
