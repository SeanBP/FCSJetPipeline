// Low-area edge cases of the whole-ring annulus (in_event_ue v4+ outputs),
// binned in the live annulus area (units of the cone area pi R^2), 20-120 GeV
// jets:
//   share of jets per area bin (data, simulation), and failing the validity
//   cut (an_fill = 0 or > 1/0.7)
//   simulation closure (truth estimate / truth target) of the ring without
//   the other-jet exclusion, binned in its own area
//   data per-jet RMS / mean of the estimate, and the bin's mean estimate over
//   the all-jet mean in the same energy bin (a bias or selection indicator)
//   mean jet eta and |x| per bin
//
// Build and run from this folder:
//   g++ -m64 -O2 ann_edge.cpp -o ann_edge `root-config --cflags --libs` -std=c++11
//   ./ann_edge "<data files>" "<sim files>"

#include <TROOT.h>
#include <TChain.h>

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

static const int kNA = 6, kNE = 3, kNBA = 7;
static const double kE[kNE + 1] = {20, 35, 60, 120};
static const double kAEdge[kNBA + 1] = {0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 99};
static const double kMaxFill = 1 / 0.7, kCone = M_PI * 0.49;

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}
struct Acc { double w = 0, wx = 0, wxx = 0; void Add(double x, double ww) { w += ww; wx += ww * x; wxx += ww * x * x; } double M() const { return w > 0 ? wx / w : NAN; } double Rms() const { const double m = M(); return w > 0 ? std::sqrt(std::max(0.0, wxx / w - m * m)) : NAN; } };

// ratio of weighted sums with the error from per-jet covariance
struct Ratio
{
    double w = 0, w2 = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    void Add(double x, double y, double ww) { w += ww; w2 += ww * ww; sx += ww * x; sy += ww * y; sxx += ww * x * x; syy += ww * y * y; sxy += ww * x * y; }
    double Val() const { return sy != 0 ? sx / sy : NAN; }
    double Err() const
    {
        if (!(w > 0) || sy == 0) return NAN;
        const double mx = sx / w, my = sy / w, r = mx / my;
        const double vx = sxx / w - mx * mx, vy = syy / w - my * my, cxy = sxy / w - mx * my;
        return std::fabs(r) * std::sqrt(std::max(0.0, vx / (mx * mx) + vy / (my * my) - 2 * cxy / (mx * my)) / (w * w / w2));
    }
};

static int ABin(double a) { for (int q = 0; q < kNBA; q++) if (a >= kAEdge[q] && a < kAEdge[q + 1]) return q; return -1; }
static int EBin(double E) { for (int q = 0; q < kNE; q++) if (E >= kE[q] && E < kE[q + 1]) return q; return -1; }

struct Res
{
    double n[kNBA] = {}, nfail[kNBA] = {}, ntot = 0;
    Ratio cl[kNBA];
    Acc est[kNBA][kNE], all[kNE], eta[kNBA], ax[kNBA];
};

static void Fill(const std::vector<std::string> &files, bool sim, Res &r)
{
    TChain ch("jets");
    for (size_t i = 0; i < files.size(); i++) ch.Add(files[i].c_str());
    Float_t E, w, jeta, jx, an_e[kNA], an_area[kNA], an_fill[kNA], an_tl[kNA], an_tf[kNA];
    ch.SetBranchAddress("jet_E", &E); ch.SetBranchAddress("weight", &w); ch.SetBranchAddress("jet_eta", &jeta); ch.SetBranchAddress("jet_x", &jx);
    ch.SetBranchAddress("an_e", an_e); ch.SetBranchAddress("an_area", an_area); ch.SetBranchAddress("an_fill", an_fill);
    ch.SetBranchAddress("an_tl", an_tl); ch.SetBranchAddress("an_tf", an_tf);
    for (Long64_t i = 0; i < ch.GetEntries(); i++)
    {
        ch.GetEntry(i);
        const int eb = EBin(E);
        if (eb < 0) continue;
        const double ww = sim ? w : 1.0;
        const int ab = ABin(an_area[0] / kCone);
        if (ab < 0) continue;
        const bool ok = an_fill[0] > 0 && an_fill[0] <= kMaxFill;
        r.ntot += ww; r.n[ab] += ww; if (!ok) r.nfail[ab] += ww;
        r.eta[ab].Add(jeta, ww); r.ax[ab].Add(std::fabs(jx), ww);
        if (ok) { r.est[ab][eb].Add(an_e[0], ww); r.all[eb].Add(an_e[0], ww); }
        const int ab3 = ABin(an_area[3] / kCone);
        if (sim && ab3 >= 0 && an_fill[3] > 0 && an_fill[3] <= kMaxFill) r.cl[ab3].Add(an_tl[3], an_tf[3], ww);
    }
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::fprintf(stderr, "usage: ./ann_edge \"<data files>\" \"<sim files>\"\n"); return 1; }
    Res D, S;
    Fill(Split(argv[1]), false, D);
    Fill(Split(argv[2]), true, S);
    std::printf("whole ring, other jets excluded, 20-120 GeV jets; area in units of pi R^2\n");
    std::printf(" area       share data  sim   fail data  sim   closure(sim, no excl)  data RMS/mean (20-35 35-60 60-120)  mean/all (20-35 35-60 60-120)   <eta>  <|x|>\n");
    // the no-exclusion ring is larger: its area bins hold different jets than the share columns
    for (int q = 0; q < kNBA; q++)
    {
        std::printf(" %4.2f-%-5.2f  %5.3f %5.3f   %5.3f %5.3f   %6.3f+-%.3f         ", kAEdge[q], kAEdge[q + 1] > 50 ? 9.99 : kAEdge[q + 1],
                    D.n[q] / D.ntot, S.n[q] / S.ntot, D.n[q] > 0 ? D.nfail[q] / D.n[q] : NAN, S.n[q] > 0 ? S.nfail[q] / S.n[q] : NAN,
                    S.cl[q].Val(), S.cl[q].Err());
        for (int e = 0; e < kNE; e++) std::printf(" %5.2f", D.est[q][e].Rms() / D.est[q][e].M());
        std::printf("                  ");
        for (int e = 0; e < kNE; e++) std::printf(" %5.2f", D.est[q][e].M() / D.all[e].M());
        std::printf("          %5.2f  %5.1f\n", D.eta[q].M(), D.ax[q].M());
    }
    return 0;
}
