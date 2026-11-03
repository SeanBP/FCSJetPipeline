// Toy closure of the A_N sign convention (AnAsymmetry.cpp): jets uniform over
// both halves of the buffered ECal fiducial region, x_F (blue frame) uniform,
// spin states with unequal luminosities, and 20% fewer jets on the x < 0 side
// (both should cancel in the cross ratio). Each jet is kept with probability
// proportional to
//   1 + sum over beams P A(x_F of that beam) (s x p_beam) . pT_hat
// with s = +-y the beam's spin and p_beam its direction (blue +z, yellow -z),
// i.e. the physical definition of A_N (left of the beam = spin x momentum),
// written without AnAsymmetry's side assignment. Injected: blue
// A_N = 0.10 x_F (x_F > 0), yellow A_N = -0.04 (x_F < 0).
// Writes h_xF_*, h_abscos_* and h_xfsum_* histograms and the "runs" entry as
// AnHistograms does (all jets in EM1, np2) and prints, per x_F bin of 0.05,
// the injected A_N averaged over the bin's jets, which AnAsymmetry should
// return. Toy files made with different polarizations and run numbers test
// AnAsymmetry's per-run polarization (-r).
//
// Usage: ./an_toy_closure <out_yields.root> [jets = 20000000] [P_blue = 0.55] [P_yellow = 0.50] [run = 1] [seed = 12345]
// Build: g++ -m64 -O2 an_toy_closure.cpp -o an_toy_closure -I../shared `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11

#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TRandom3.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "JetParameters.h"
#include "AnConfig.h"

static double InjectedAN(int beam, double xf_beam) { return beam == 0 ? 0.10 * xf_beam : -0.04; }

int main(int argc, char **argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: ./an_toy_closure <out.root> [jets] [P_blue] [P_yellow]\n"); return 1; }
    const long n_want = argc > 2 ? std::atol(argv[2]) : 20000000L;
    const double pol[2] = {argc > 3 ? std::atof(argv[3]) : 0.55, argc > 4 ? std::atof(argv[4]) : 0.50};
    const int run = argc > 5 ? std::atoi(argv[5]) : 1;
    const int seed = argc > 6 ? std::atoi(argv[6]) : 12345;
    const double lumi[kNSpinStates] = {0.30, 0.20, 0.25, 0.25};  // BuYu, BuYd, BdYu, BdYd
    const int blue_spin[kNSpinStates] = {+1, +1, -1, -1}, yellow_spin[kNSpinStates] = {+1, -1, +1, -1};
    const double beam_z[2] = {+1, -1};

    TFile out(argv[1], "RECREATE");
    TH1D *h[kNEmLevels][kNPhotonClasses][kNSides][kNSpinStates], *hc[kNEmLevels][kNPhotonClasses], *hx[kNEmLevels][kNPhotonClasses];
    for (int l = 0; l < kNEmLevels; l++)
        for (int k = 0; k < kNPhotonClasses; k++)
            for (int s = 0; s < kNSides; s++)
                for (int p = 0; p < kNSpinStates; p++)
                {
                    TString name = TString::Format("h_xF_%s_%s_%s_%s", kEmLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]);
                    h[l][k][s][p] = new TH1D(name, name + ";x_{F};jets", kNXfBins, kXfMin, kXfMax);
                }
    for (int l = 0; l < kNEmLevels; l++)
        for (int k = 0; k < kNPhotonClasses; k++)
        {
            TString name = TString::Format("h_abscos_%s_%s", kEmLevelName[l], PhotonClassName(k));
            hc[l][k] = new TH1D(name, name, kNXfBins, kXfMin, kXfMax);
            name = TString::Format("h_xfsum_%s_%s", kEmLevelName[l], PhotonClassName(k));
            hx[l][k] = new TH1D(name, name, kNXfBins, kXfMin, kXfMax);
        }
    const int nb = kNXfBins / 2;
    double sum_exp[2][kNXfBins / 2] = {{0}}, n_exp[kNXfBins / 2] = {0};

    TRandom3 rng(seed);
    const FiducialRect &r = kFiducialRectEcal;
    const double b = reco_fiducial_buffer_ecal;
    const double wmax = 1 + pol[0] * 0.10 + pol[1] * 0.04;
    long n_kept = 0;
    while (n_kept < n_want)
    {
        // position at the ECal plane, either half
        const double side_sign = rng.Rndm() < 0.5 ? -1 : 1;
        const double x = side_sign * rng.Uniform(r.x_inner + b, r.x_outer - b), y = rng.Uniform(r.y_min + b, r.y_max - b);
        if (side_sign < 0 && rng.Rndm() < 0.2) continue;  // side-dependent efficiency
        const double xf = rng.Uniform(0.025, 0.975);
        // spin state by luminosity
        double u = rng.Rndm(); int p = 0;
        while (p < kNSpinStates - 1 && u >= lumi[p]) u -= lumi[p++];
        const double pt_hat[3] = {x / std::hypot(x, y), y / std::hypot(x, y), 0};
        double w = 1, lr[2];
        for (int beam = 0; beam < 2; beam++)
        {
            // left = s x p_beam, s = (0, spin, 0), p_beam = (0, 0, beam_z)
            const double spin = beam == 0 ? blue_spin[p] : yellow_spin[p];
            const double left[3] = {spin * beam_z[beam], 0, 0};  // (s_y p_z - 0, 0 - 0, 0 - 0)
            lr[beam] = left[0] * pt_hat[0] + left[1] * pt_hat[1] + left[2] * pt_hat[2];
            w += pol[beam] * InjectedAN(beam, beam == 0 ? xf : -xf) * lr[beam];
        }
        if (rng.Rndm() * wmax > w) continue;
        n_kept++;
        h[1][2][x > 0 ? 1 : 0][p]->Fill(xf);
        hc[1][2]->Fill(xf, std::fabs(pt_hat[0]));
        hx[1][2]->Fill(xf, xf);
        // expected: the injected A_N averaged over the bin's jets
        const int ib = (int)(xf * nb);
        n_exp[ib]++;
        for (int beam = 0; beam < 2; beam++) sum_exp[beam][ib] += InjectedAN(beam, beam == 0 ? xf : -xf);
    }
    TNamed("runs", TString::Format("%d", run).Data()).Write();
    out.Write();
    out.Close();
    std::printf("an_toy_closure: %ld jets, P_blue %.2f, P_yellow %.2f\n  expected A_N (injected, mean over the bin's jets):\n  %-6s %13s %10s\n", n_kept, pol[0], pol[1], "beam", "x_F", "A_N");
    for (int beam = 0; beam < 2; beam++)
        for (int i = 0; i < nb; i++)
        {
            const double lo = (double)i / nb, hi = (double)(i + 1) / nb;
            std::printf("  %-6s %6.3f-%6.3f %+9.5f\n", beam == 0 ? "blue" : "yellow", beam == 0 ? lo : -hi, beam == 0 ? hi : -lo, sum_exp[beam][i] / n_exp[i]);
        }
    return 0;
}
