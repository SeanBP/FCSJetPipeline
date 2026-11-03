// The eta-weighted wedges (WeightedWedgeUE.h) when eta bins they can't sample give their
// area to nearby bins that they can, against the current scheme, where the whole
// band stands in for them. A source bin's area M goes to the wedges' other bins
// as r(b) += M K(b) / sum_b K(b) h_wedge(b), K a Gaussian in eta, so the
// oversampling falls off gradually from the source. Variants: current; sigma
// 0.1, 0.2, 0.4 with unreached bins as sources; sigma 0.2 with bins above r = 3
// also giving their excess (thin bins, not only empty ones).
//
// Reco jets, pT >= 2 GeV, E 10-120 GeV, fiducial. Per variant, for all jets and
// for jets with a bin unreached: mean density, shift from current, its per-jet
// RMS, mean largest weight. With truth in the files (-t), a closure on the same
// jet axes: truth photons on the live wedge tiles, each variant's weights,
// against a target from the full-azimuth truth wedges: its density per eta bin,
// averaged with the jet's live h_jet(b). Jets are unweighted.
//
// Build and run from this folder, in the container:
//   g++ -m64 -O2 ring_taper_check.cpp ../../shared/VectorDict.cxx -o ring_taper_check -I../../shared -I../../underlying_event `root-config --cflags --libs` -std=c++11
//   ./ring_taper_check <pairs list> [max files = 120] [-t]

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "WeightedWedgeUE.h"
#include "in_event_util.h"

static const int kNV = 5;
static const char *kName[kNV] = {"current", "taper 0.1", "taper 0.2", "taper 0.4", "taper 0.2 + cap 3"};
static const double kSigma[kNV] = {0, 0.1, 0.2, 0.4, 0.2};
static const double kCap[kNV] = {0, 0, 0, 0, 3};

// the jet cone's and the wedges' areas per eta bin, as RecoWeightedWedges
static double WedgeAreas(const UeGeometry &g, double jeta, double jphi, double *hj, double *hr, std::vector<char> &in)
{
    const double eta0 = jeta - R - 0.1;
    double area = 0;
    for (int b = 0; b < kUeNB; b++) hj[b] = hr[b] = 0;
    in.assign(g.tiles.size(), 0);
    const int ns = kUeSub * kUeSub;
    for (size_t c = 0; c < g.tiles.size(); c++)
    {
        if (std::fabs(g.tiles[c].eta - jeta) > R + 0.2) continue;
        in[c] = InUeWedges(g.tiles[c].eta, g.tiles[c].phi, jeta, jphi);
        for (int k = 0; k < ns; k++)
        {
            const UeGeometry::Sub &q = g.sub[c * ns + k];
            const int b = (int)std::floor((q.eta - eta0) / kUeBin);
            if (b < 0 || b >= kUeNB) continue;
            if (UeDeltaR(q.eta, q.phi, jeta, jphi) < R) { hj[b] += q.area; area += q.area; }
            if (in[c]) hr[b] += q.area;
        }
    }
    return area;
}

// the per-area weights of variant v; returns sum r(b) h_wedge(b)
static double Weights(int v, const double *hj, const double *hr, double *r)
{
    double src[kUeNB];
    for (int b = 0; b < kUeNB; b++)
    {
        const bool reached = hj[b] > 0 && hr[b] >= kUeMinRatio * hj[b];
        r[b] = reached ? hj[b] / hr[b] : 0;
        src[b] = hj[b] > 0 && !reached ? hj[b] : 0;
        if (kCap[v] > 0 && r[b] > kCap[v]) { src[b] = (r[b] - kCap[v]) * hr[b]; r[b] = kCap[v]; }
    }
    if (kSigma[v] > 0)
    {
        double add[kUeNB] = {0};
        for (int s = 0; s < kUeNB; s++)
        {
            if (!(src[s] > 0)) continue;
            double norm = 0, k[kUeNB];
            for (int b = 0; b < kUeNB; b++)
            {
                const double d = (b - s) * kUeBin / kSigma[v];
                // targets: bins the wedges reach, the jet's or not, other than the source
                k[b] = b != s && (hj[b] > 0 ? hr[b] >= kUeMinRatio * hj[b] : hr[b] > 0) ? std::exp(-0.5 * d * d) : 0;
                norm += k[b] * hr[b];
            }
            if (norm > 0) for (int b = 0; b < kUeNB; b++) add[b] += src[s] * k[b] / norm;
        }
        for (int b = 0; b < kUeNB; b++) r[b] += add[b];
    }
    double den = 0;
    for (int b = 0; b < kUeNB; b++) den += r[b] * hr[b];
    return den;
}

struct Acc
{
    double n, s, sd, sdd, sw, sc, scc;
    void Add(double rho, double rho0, double wmax, double closure)
    {
        n++; s += rho; sd += rho - rho0; sdd += (rho - rho0) * (rho - rho0); sw += wmax;
        sc += closure; scc += closure * closure;
    }
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2) { std::fprintf(stderr, "usage: ./ring_taper_check <pairs list> [max files] [-t]\n"); return 1; }
    const int max_files = argc > 2 ? std::atoi(argv[2]) : 120;
    const bool do_truth = argc > 3 && !std::strcmp(argv[3], "-t");
    UeGeometry g; g.Build();
    Acc all[kNV], low[kNV];
    std::memset(all, 0, sizeof(all)); std::memset(low, 0, sizeof(low));
    double n_low[2] = {0, 0}, s_cov = 0, m_tot = 0, m_hi = 0, s_target[2] = {0, 0}, n_target[2] = {0, 0};
    const int MAX_PTS = 200, MAX_MC = 10000;
    std::vector<Int_t> geid(MAX_MC), idv(MAX_MC);
    std::vector<Float_t> mpx(MAX_MC), mpy(MAX_MC), mpz(MAX_MC), mE(MAX_MC), mvx(MAX_MC), mvy(MAX_MC), mvz(MAX_MC);
    std::ifstream in_list(argv[1]); std::string jf, sf;
    int nf = 0;
    while (in_list >> jf >> sf && nf < max_files)
    {
        TFile fj(jf.c_str()), fs(sf.c_str());
        TTree *jt = (TTree *)fj.Get("jetTree"), *st = (TTree *)fs.Get("data");
        if (!jt || !st) continue;
        nf++;
        std::vector<float> *pt = 0, *eta = 0, *phi = 0, *E = 0, *x = 0, *y = 0;
        Long64_t sent = -1;
        const bool has_entry = jt->GetBranch("simpletree_entry");
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_E", "reco_x", "reco_y"};
        for (int i = 0; i < 6; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &pt); jt->SetBranchAddress("reco_eta", &eta); jt->SetBranchAddress("reco_phi", &phi);
        jt->SetBranchAddress("reco_E", &E); jt->SetBranchAddress("reco_x", &x); jt->SetBranchAddress("reco_y", &y);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &sent); }
        Int_t np = 0, det[MAX_PTS], nmc = 0;
        Float_t pE[MAX_PTS], px[MAX_PTS], py[MAX_PTS], pz[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pb[i], 1);
        st->SetBranchAddress("Cal_npts", &np); st->SetBranchAddress("Cal_pt_energy", pE); st->SetBranchAddress("Cal_pt_x", px);
        st->SetBranchAddress("Cal_pt_y", py); st->SetBranchAddress("Cal_pt_z", pz); st->SetBranchAddress("Cal_pt_detid", det);
        if (do_truth)
        {
            const char *mb[] = {"mcpart_num", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_geid", "mcpart_idVtx", "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
            for (int i = 0; i < 10; i++) st->SetBranchStatus(mb[i], 1);
            st->SetBranchAddress("mcpart_num", &nmc); st->SetBranchAddress("mcpart_px", &mpx[0]); st->SetBranchAddress("mcpart_py", &mpy[0]);
            st->SetBranchAddress("mcpart_pz", &mpz[0]); st->SetBranchAddress("mcpart_E", &mE[0]); st->SetBranchAddress("mcpart_geid", &geid[0]);
            st->SetBranchAddress("mcpart_idVtx", &idv[0]); st->SetBranchAddress("mcpart_Vtx_x", &mvx[0]);
            st->SetBranchAddress("mcpart_Vtx_y", &mvy[0]); st->SetBranchAddress("mcpart_Vtx_z", &mvz[0]);
        }
        std::vector<char> in;
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
                if ((det[i] != 0 && det[i] != 1) || !pass_point_energy_cut(pE[i], kSqrtS_MainRun22_GeV)) continue;
                te.push_back(std::asinh(pz[i] / std::hypot(px[i], py[i]))); tE.push_back(pE[i]); tc.push_back(g.TileOf(px[i], py[i]));
            }
            std::vector<float> ge, gp, gE; std::vector<int> gc;
            if (do_truth)
            {
                float vx = 0, vy = 0, vz = 0;
                const int n = std::min((int)nmc, MAX_MC);
                if (FindEmTruthPrimaryVertex(n, &idv[0], &mvx[0], &mvy[0], &mvz[0], vx, vy, vz))
                    for (int i = 0; i < n; i++)
                    {
                        if (!IsEmTruthParticle(geid[i], mvx[i], mvy[i], mvz[i], vx, vy, vz)) continue;
                        const double p = std::hypot((double)mpx[i], (double)mpy[i]);
                        if (!(p > 0) || !(mpz[i] > 0)) continue;
                        const double e = std::asinh(mpz[i] / p), f = std::atan2(mpy[i], mpx[i]), rr = z_proj / std::sinh(e);
                        ge.push_back(e); gp.push_back(f); gE.push_back(mE[i]); gc.push_back(g.TileOf(rr * std::cos(f), rr * std::sin(f)));
                    }
            }
            for (size_t q = 0; q < sel.size(); q++)
            {
                const int j = sel[q];
                const double jeta = eta->at(j), jphi = phi->at(j), eta0 = jeta - R - 0.1;
                double hj[kUeNB], hr[kUeNB], r[kUeNB];
                const double area = WedgeAreas(g, jeta, jphi, hj, hr, in);
                // the unreached area, and where it lies in the jet's band
                double miss = 0, miss_hi = 0;
                for (int b = 0; b < kUeNB; b++)
                    if (hj[b] > 0 && hr[b] < kUeMinRatio * hj[b]) { miss += hj[b]; if (eta0 + (b + 0.5) * kUeBin > jeta) miss_hi += hj[b]; }
                const bool is_low = miss > 0.001 * area;
                n_low[0]++; if (is_low) n_low[1]++;
                s_cov += area > 0 ? 1 - miss / area : 0; m_tot += miss; m_hi += miss_hi;
                double target = 0;
                if (do_truth)
                {
                    double Ef[kUeNB] = {0}, hf[kUeNB], w = 0;
                    for (size_t i = 0; i < gE.size(); i++)
                    {
                        if (!InUeWedges(ge[i], gp[i], jeta, jphi)) continue;
                        const int b = (int)std::floor((ge[i] - eta0) / kUeBin);
                        if (b >= 0 && b < kUeNB) Ef[b] += gE[i];
                    }
                    for (int b = 0; b < kUeNB; b++)
                    {
                        hf[b] = 0;
                        for (int s = 0; s < 50; s++)
                        {
                            const double d = eta0 + (b + (s + 0.5) / 50) * kUeBin - jeta;
                            if (std::fabs(d) < R) hf[b] += 2 * std::max(0.0, kUeDphiMax - std::sqrt(kUeJetExcl * kUeJetExcl - d * d)) * kUeBin / 50;
                        }
                        if (hf[b] > 0 && hj[b] > 0) { target += hj[b] * Ef[b] / hf[b]; w += hj[b]; }
                    }
                    target = w > 0 ? target / w : 0;
                    s_target[0] += target; n_target[0]++;
                    if (is_low) { s_target[1] += target; n_target[1]++; }
                }
                double rho0 = 0;
                for (int v = 0; v < kNV; v++)
                {
                    const double den = Weights(v, hj, hr, r);
                    if (!(den > 0)) continue;
                    double e = 0, wmax = 0, c = 0;
                    for (int b = 0; b < kUeNB; b++) wmax = std::max(wmax, r[b]);
                    for (size_t i = 0; i < tE.size(); i++)
                    {
                        if (tc[i] < 0 || !in[tc[i]]) continue;
                        const int b = (int)std::floor((te[i] - eta0) / kUeBin);
                        if (b >= 0 && b < kUeNB) e += tE[i] * r[b];
                    }
                    const double rho = e / den;
                    if (v == 0) rho0 = rho;
                    if (do_truth)
                    {
                        for (size_t i = 0; i < gE.size(); i++)
                        {
                            if (gc[i] < 0 || !in[gc[i]]) continue;
                            const int b = (int)std::floor((ge[i] - eta0) / kUeBin);
                            if (b >= 0 && b < kUeNB) c += gE[i] * r[b];
                        }
                        c = c / den - target;
                    }
                    all[v].Add(rho, rho0, wmax, c);
                    if (is_low) low[v].Add(rho, rho0, wmax, c);
                }
            }
        }
    }
    std::printf("ring_taper_check: %d files, %.0f jets; %.0f (%.2f%%) with unreached bins; mean eta coverage %.4f; share of the unreached area on the high-eta side %.2f\n",
                nf, n_low[0], n_low[1], 100 * n_low[1] / n_low[0], s_cov / n_low[0], m_tot > 0 ? m_hi / m_tot : 0);
    for (int set = 0; set < 2; set++)
    {
        const Acc *a = set ? low : all;
        std::printf("\n%s\n  %-18s %8s %9s %9s %9s %7s", set ? "jets with unreached bins" : "all jets", "variant", "jets", "mean rho", "- current", "RMS", "w_max");
        if (do_truth) std::printf(" %14s %11s", "truth bias", "truth RMS");
        std::printf("\n");
        if (do_truth && n_target[set] > 0) std::printf("  truth target mean density %.4f\n", s_target[set] / n_target[set]);
        for (int v = 0; v < kNV; v++)
        {
            if (!(a[v].n > 0)) continue;
            const double md = a[v].sd / a[v].n, rd = std::sqrt(std::max(0.0, a[v].sdd / a[v].n - md * md));
            std::printf("  %-18s %8.0f %9.4f %+9.4f %9.4f %7.3f", kName[v], a[v].n, a[v].s / a[v].n, md, rd, a[v].sw / a[v].n);
            if (do_truth)
            {
                const double mc = a[v].sc / a[v].n, rc = std::sqrt(std::max(0.0, a[v].scc / a[v].n - mc * mc));
                std::printf(" %+7.4f+-%.4f %11.4f", mc, rc / std::sqrt(a[v].n), rc);
            }
            std::printf("\n");
        }
    }
    return 0;
}
