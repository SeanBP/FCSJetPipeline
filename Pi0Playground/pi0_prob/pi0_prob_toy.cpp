// Closure test of the pi0 probability code (shared/Pi0Prob.h) on toy jets
// drawn from the model itself, so the predicted P must be calibrated:
//   n uniform in 2..12; k resolved pi0 drawn from a known prior pi_k(n)
//   (Poisson, mean 0.25 n, truncated to k <= n/2); the k signal pairs a
//   uniformly random k-matching; signal pair masses Gaussian (0.135, 0.020)
//   GeV, every other pair's mass exponential (slope 4/GeV) on [0, 1.2] GeV,
//   all drawn independently (the model's assumption).
// Prints the calibration table (true fraction of k >= 1 against mean P in P
// bins, overall and by n; and of point 0 being in a signal pair against
// P_leadpair) and the largest deviation in units of its error.
//
// Build: g++ -O2 pi0_prob_toy.cpp -o pi0_prob_toy -I../../shared `root-config --cflags --libs` -std=c++11
// Run:   ./pi0_prob_toy [jets = 200000] [seed = 1]
#include "Pi0Prob.h"

#include <TMath.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

int main(int argc, char **argv)
{
    const long N = argc > 1 ? atol(argv[1]) : 200000;
    TRandom3 rng(argc > 2 ? atoi(argv[2]) : 1);
    const double mu = 0.135, sig = 0.020, slope = 4.0, mmax = 1.2;
    const int nMin = 2, nMax = 12, kmax = 6;

    Pi0ProbCalib c;
    c.eS = {0, 1e4}; c.sMean = {mu}; c.sSigma = {sig}; c.sAlpha = {20}; c.sN = {3};
    c.eB = {0, 1e4};
    const int nb = 1200;
    for (int i = 0; i <= nb; i++) c.mB.push_back(mmax * i / nb);
    c.bPdf.assign(1, vector<double>(nb));
    const double norm = slope / (1 - exp(-slope * mmax));
    for (int i = 0; i < nb; i++) c.bPdf[0][i] = norm * exp(-slope * 0.5 * (c.mB[i] + c.mB[i + 1]));
    c.eps = 1e-12;
    c.kmax = kmax;
    c.eJ = {0, 1e4};
    for (int n = nMin; n <= nMax; n++) c.nLo.push_back(n);
    vector<vector<double> > pri(nMax + 1, vector<double>(kmax + 1, 0));
    for (int n = nMin; n <= nMax; n++)
    {
        const double lam = 0.25 * n;
        double t = 0;
        for (int k = 0; k <= min(kmax, n / 2); k++) { pri[n][k] = TMath::Poisson(k, lam); t += pri[n][k]; }
        for (int k = 0; k <= kmax; k++) { pri[n][k] /= t; c.prior.push_back(pri[n][k]); }
    }
    c.Normalise();
    c.ok = true;

    const int NB = 10;
    // calibration sums: [n class][P bin] (class 0 all, 1: n 2-3, 2: 4-6, 3: 7+)
    double sw[4][NB] = {{0}}, sp[4][NB] = {{0}}, sy[4][NB] = {{0}};
    // P_leadpair (point 0 in a signal pair), all n
    double lw[NB] = {0}, lp[NB] = {0}, ly[NB] = {0};
    for (long it = 0; it < N; it++)
    {
        const int n = nMin + (int)rng.Integer(nMax - nMin + 1);
        double u = rng.Uniform(), cum = 0;
        int k = 0;
        for (int q = 0; q <= kmax; q++) { cum += pri[n][q]; if (u < cum) { k = q; break; } }
        vector<int> perm(n);
        for (int i = 0; i < n; i++) perm[i] = i;
        for (int i = n - 1; i > 0; i--) swap(perm[i], perm[rng.Integer(i + 1)]);
        vector<vector<char> > sigp(n, vector<char>(n, 0));
        for (int q = 0; q < k; q++) sigp[perm[2 * q]][perm[2 * q + 1]] = sigp[perm[2 * q + 1]][perm[2 * q]] = 1;
        vector<vector<double> > r(n, vector<double>(n, 0));
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                double m;
                if (sigp[i][j]) m = rng.Gaus(mu, sig);
                else { do m = -log(1 - rng.Uniform() * (1 - exp(-slope * mmax))) / slope; while (m >= mmax); }
                r[i][j] = r[j][i] = c.S(m, 1) / c.B(m, 1);
            }
        Pi0ProbResult res;
        res.npairs = n * (n - 1) / 2;
        res = Pi0ProbFromRatios(c, r, 1, res);
        const int b = min(NB - 1, (int)(res.prob * NB));
        const int cls = n <= 3 ? 1 : (n <= 6 ? 2 : 3);
        for (int q : {0, cls}) { sw[q][b] += 1; sp[q][b] += res.prob; sy[q][b] += (k >= 1); }
        bool leadSig = false;
        for (int j = 1; j < n; j++) if (sigp[0][j]) leadSig = true;
        const int bl = min(NB - 1, (int)(res.lead_pair * NB));
        lw[bl] += 1; lp[bl] += res.lead_pair; ly[bl] += leadSig;
    }
    const char *nm[4] = {"all", "n 2-3", "n 4-6", "n 7+"};
    double worst = 0;
    for (int q = 0; q < 4; q++)
    {
        printf("%-6s  <P> -> true fraction (pull) [jets]\n", nm[q]);
        for (int b = 0; b < NB; b++)
        {
            if (sw[q][b] < 1) continue;
            const double mp = sp[q][b] / sw[q][b], f = sy[q][b] / sw[q][b], e = sqrt(max(f * (1 - f), 1e-6) / sw[q][b]);
            const double pull = (f - mp) / e;
            if (sw[q][b] >= 100) worst = max(worst, fabs(pull));
            printf("   %.3f -> %.3f (%+.1f) [%.0f]\n", mp, f, pull, sw[q][b]);
        }
    }
    printf("leading point in a signal pair: <P_leadpair> -> true fraction (pull) [jets]\n");
    for (int b = 0; b < NB; b++)
    {
        if (lw[b] < 1) continue;
        const double mp = lp[b] / lw[b], f = ly[b] / lw[b], e = sqrt(max(f * (1 - f), 1e-6) / lw[b]);
        if (lw[b] >= 100) worst = max(worst, fabs((f - mp) / e));
        printf("   %.3f -> %.3f (%+.1f) [%.0f]\n", mp, f, (f - mp) / e, lw[b]);
    }
    printf("largest |pull| (bins with >= 100 jets): %.2f\n", worst);
    return 0;
}
