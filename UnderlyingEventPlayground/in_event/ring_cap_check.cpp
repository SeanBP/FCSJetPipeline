// The FCS pixel ring (in_event_geo.h, per-point weights; 2R, no other-jet
// exclusion) with cells at |dphi| above a cap left out, caps kCaps, on real
// FCS data jets (reco pT >= 2 GeV, E 10-120 GeV, fiducial; every event). Per
// cap: the mean density, the per-jet spread, the eta coverage (share of the
// jet's live area in eta bins the ring reaches), the share of jets below full,
// 0.9 and 0.5 coverage and with no ring (fallback), the mean largest cell
// weight and live ring area, and the per-jet spread of (capped - uncut).
//
// Build and run from this folder:
//   g++ -m64 -O2 ring_cap_check.cpp ../../shared/VectorDict.cxx -o ring_cap_check -I../../shared -I../../underlying_event `root-config --cflags --libs` -std=c++11
//   ./ring_cap_check <pairs list> [max files = 120]

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
static const int kNC = 3;
static const double kCaps[kNC] = {180, 150, 120};
struct Acc
{
    double n = 0, sr = 0, srr = 0, scov = 0, nc1 = 0, nc09 = 0, nc05 = 0, nfb = 0, sw = 0, sa = 0, sd = 0, sdd = 0, nd = 0;
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2) { std::fprintf(stderr, "usage: ./ring_cap_check <pairs list> [max files]\n"); return 1; }
    const int max_files = argc > 2 ? std::atoi(argv[2]) : 120;
    Geometry geo; geo.Build();
    const std::vector<std::pair<float, float> > none;
    Acc acc[kNC];
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
                double rho0 = -1;
                for (int c = 0; c < kNC; c++)
                {
                    BuildPixelRing(geo, eta->at(j), phi->at(j), none, kAnnJetExcl, pr, kCaps[c] * M_PI / 180 + 1e-9);
                    Acc &a = acc[c];
                    a.n++;
                    const double cov = pr.area_jet > 0 ? pr.den / pr.area_jet : 0;
                    a.scov += cov; a.sw += pr.w_max; a.sa += pr.area_ring;
                    if (cov < 0.999) a.nc1++;
                    if (cov < 0.9) a.nc09++;
                    if (cov < 0.5) a.nc05++;
                    if (!(pr.den > 0)) { a.nfb++; continue; }
                    double e = 0;
                    for (size_t i = 0; i < tE.size(); i++) e += tE[i] * pr.Weight(tc[i], te[i]);
                    const double rho = e / pr.den;
                    a.sr += rho; a.srr += rho * rho;
                    if (c == 0) rho0 = rho;
                    else if (rho0 >= 0) { a.sd += rho - rho0; a.sdd += (rho - rho0) * (rho - rho0); a.nd++; }
                }
            }
        }
    }
    std::printf("ring_cap_check: %d files, %.0f jets (FCS pixel ring, 2R, other jets kept, per-point weights)\n", nf, acc[0].n);
    std::printf("  cut (|dphi| <=)  mean rho   per-jet RMS   eta cov  <1      <0.9    <0.5    no ring  mean w_max  mean ring area   capped - uncut: mean   RMS\n");
    for (int c = 0; c < kNC; c++)
    {
        const Acc &a = acc[c];
        const double nr = a.n - a.nfb, m = a.sr / nr, r = std::sqrt(std::max(0.0, a.srr / nr - m * m));
        std::printf("  %3.0f deg          %.4f     %.4f        %.4f   %.4f  %.4f  %.4f  %.4f   %.3f       %.3f", kCaps[c], m, r, a.scov / a.n, a.nc1 / a.n, a.nc09 / a.n, a.nc05 / a.n, a.nfb / a.n, a.sw / a.n, a.sa / a.n);
        if (a.nd > 0) { const double md = a.sd / a.nd; std::printf("            %+.4f  %.4f", md, std::sqrt(std::max(0.0, a.sdd / a.nd - md * md))); }
        std::printf("\n");
    }
    return 0;
}
