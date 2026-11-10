#ifndef PI0_PROB_H
#define PI0_PROB_H

// Per-jet probability that an EM jet contains at least one resolved pi0 (both
// decay photons reconstructed as separate points of the jet), given the
// invariant masses of all its point pairs (Bayesian matching model, 2026-11-09).
//
// For a jet with n points and M = n(n-1)/2 pairs of mass m_ij:
//   r_ij = S(m_ij | E_ij) / B(m_ij | E_ij)   (pair energy E_ij = E_i + E_j)
//   e_k  = sum over the k-matchings A (k disjoint pairs) of prod_{(ij) in A} r_ij
//   N_k  = n! / (k! 2^k (n-2k)!), the number of k-matchings
//   O    = sum_{k>=1} (pi_k(n)/pi_0(n)) e_k / N_k   (posterior odds)
//   P    = O / (1 + O)
// S: mass PDF of a true pi0 pair (Crystal Ball, low-side tail), B: PDF of
// every other pair, pi_k(n, E_jet): prior of k resolved pi0 among n points.
// The 1/N_k is the trials factor: with many points one pair near the pi0 mass
// is weak evidence. Masses are treated as independent given the matching
// (pairs share points, so P is approximate at large n; see the validation).
//
// All three come from simulation (Pi0Playground/pi0_prob/pi0_prob_derive.cpp)
// and are read from a versioned JSON (Pi0ProbCalib::Load). For data the JSON
// can carry S peak position and width refitted on data point pairs
// ("S_data_mean", "S_data_sigma"), used when Load(..., data = true).
//
// The pair mass uses the point 4-vectors as the jets do (energy along the
// direction from the nominal vertex at the origin), m^2 = (p_i + p_j)^2.
//
// Leading point (Sean 2026-11-09: most EM jets contain some pi0, so the jet
// classes are "leading point from a pi0" or not). With the points ordered by
// energy, the leading one first, the same sum restricted to the matchings
// that pair the leading point gives
//   P_leadpair = P(leading point in a resolved pi0 pair | masses)
//              = sum_k (pi_k/N_k) e_k^lead / sum_k (pi_k/N_k) e_k  (pi_0 e_0 = pi_0)
// (e_k^lead = e_k minus the sum over the k-matchings leaving it unpaired).
// A merged pi0 (both photons in one point) or one whose partner photon is no
// jet point leaves no mass to see; its share comes from simulation,
// q(n, E_jet) = fraction of jets whose leading point is from a pi0 among
// those whose leading point is in no resolved pi0 pair, so
//   P_lead = P_leadpair + (1 - P_leadpair) q(n, E_jet),
// defined for n = 1 too (P_lead = q(1, E)).
//
// Status: kPi0ProbOk; kPi0ProbFewPoints (n < 2, P = 0, P_lead = q); kPi0ProbTooMany
// (n > max_points, P = P_lead = logodds = -1, not approximated); kPi0ProbNoCalib.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

enum { kPi0ProbOk = 0, kPi0ProbFewPoints = 1, kPi0ProbTooMany = 2, kPi0ProbNoCalib = 3 };

struct Pi0ProbResult
{
    float prob = 0, logodds = -99;
    int npairs = 0, status = kPi0ProbNoCalib;
    float best_r = 0;
    int best_i = -1, best_j = -1;
    float lead_pair = 0, lead = -1;  // P_leadpair, P_lead (-1: not computed)
};

struct Pi0ProbCalib
{
    std::string version;
    double eps = 1e-6;          // floor of S and B, GeV^-1
    int max_points = 20;
    // S: per pair-energy bin (edges eS)
    std::vector<double> eS, sMean, sSigma, sAlpha, sN, sNorm;
    // B: per pair-energy bin (edges eB), density on mass bins (edges mB)
    std::vector<double> eB, mB;
    std::vector<std::vector<double> > bPdf;
    // prior: n bins (lower edges nLo, last open), jet-energy bins (edges eJ),
    // pi_k for k = 0..kmax
    std::vector<double> nLo, eJ;
    int kmax = 0;
    std::vector<double> prior; // [(in * nE + ie) * (kmax + 1) + k]
    // q(n, E_jet) for the leading point: n bins (lower edges qNLo, from 1),
    // the prior's jet-energy bins; absent from the JSON -> P_lead not computed
    std::vector<double> qNLo, q; // q[iq * nE + ie]
    bool ok = false;
    bool hasLead() const { return !q.empty(); }
    double Q(int n, double ejet) const
    {
        const int nE = (int)eJ.size() - 1;
        return q[LoBin(qNLo, n) * nE + Bin(eJ, ejet)];
    }

    static int Bin(const std::vector<double> &edges, double x)
    {
        // edges e0 < e1 < ... ; x below e0 -> 0, above the last -> last bin
        const int nb = (int)edges.size() - 1;
        if (nb < 1) return 0;
        for (int i = 0; i < nb; i++) if (x < edges[i + 1]) return i;
        return nb - 1;
    }
    static int LoBin(const std::vector<double> &lo, double x)
    {
        int b = 0;
        for (size_t i = 0; i < lo.size(); i++) if (x >= lo[i]) b = (int)i;
        return b;
    }

    double SRaw(int b, double m) const
    {
        const double t = (m - sMean[b]) / sSigma[b], a = std::fabs(sAlpha[b]), nn = sN[b];
        if (t > -a) return std::exp(-0.5 * t * t);
        const double A = std::pow(nn / a, nn) * std::exp(-0.5 * a * a), Bc = nn / a - a;
        return A * std::pow(Bc - t, -nn);
    }
    double S(double m, double epair) const
    {
        const int b = Bin(eS, epair);
        const double v = SRaw(b, m) / (sNorm[b] * sSigma[b]);
        return v > eps ? v : eps;
    }
    double B(double m, double epair) const
    {
        const int b = Bin(eB, epair);
        const std::vector<double> &d = bPdf[b];
        const int nm = (int)mB.size() - 1;
        double v;
        if (m <= 0.5 * (mB[0] + mB[1])) v = d[0];
        else if (m >= 0.5 * (mB[nm - 1] + mB[nm])) v = d[nm - 1];
        else
        {
            int i = 0;
            while (i < nm - 1 && m >= 0.5 * (mB[i + 1] + mB[i + 2])) i++;
            const double c0 = 0.5 * (mB[i] + mB[i + 1]), c1 = 0.5 * (mB[i + 1] + mB[i + 2]);
            v = d[i] + (d[i + 1] - d[i]) * (m - c0) / (c1 - c0);
        }
        return v > eps ? v : eps;
    }
    // pi_k(n, E_jet), k = 0..kmax (zero beyond)
    double Prior(int n, double ejet, int k) const
    {
        if (k > kmax) return 0;
        const int in = LoBin(nLo, n), nE = (int)eJ.size() - 1, ie = Bin(eJ, ejet);
        return prior[(in * nE + ie) * (kmax + 1) + k];
    }

    // normalisation of each S bin's Crystal Ball over [0, 2 GeV], in units of sigma
    void Normalise()
    {
        sNorm.assign(sMean.size(), 1);
        for (size_t b = 0; b < sMean.size(); b++)
        {
            double s = 0;
            const int N = 4000;
            for (int i = 0; i < N; i++) s += SRaw((int)b, (i + 0.5) * 2.0 / N) * 2.0 / N;
            sNorm[b] = s / sSigma[b];
        }
    }

    // --- minimal JSON reading (flat structure written by pi0_prob_derive) ---
    static bool FindKey(const std::string &js, const std::string &key, size_t &pos)
    {
        pos = js.find("\"" + key + "\"");
        if (pos == std::string::npos) return false;
        pos = js.find(':', pos);
        return pos != std::string::npos;
    }
    static bool Num(const std::string &js, const std::string &key, double &v)
    {
        size_t p;
        if (!FindKey(js, key, p)) return false;
        v = std::strtod(js.c_str() + p + 1, 0);
        return true;
    }
    static bool Str(const std::string &js, const std::string &key, std::string &v)
    {
        size_t p;
        if (!FindKey(js, key, p)) return false;
        const size_t a = js.find('"', p), b = js.find('"', a + 1);
        if (a == std::string::npos || b == std::string::npos) return false;
        v = js.substr(a + 1, b - a - 1);
        return true;
    }
    // a flat array of numbers (nested brackets flattened)
    static bool Arr(const std::string &js, const std::string &key, std::vector<double> &v)
    {
        size_t p;
        v.clear();
        if (!FindKey(js, key, p)) return false;
        p = js.find('[', p);
        if (p == std::string::npos) return false;
        int depth = 0;
        for (size_t i = p; i < js.size(); i++)
        {
            const char c = js[i];
            if (c == '[') depth++;
            else if (c == ']') { if (--depth == 0) return true; }
            else if (c == '-' || c == '.' || (c >= '0' && c <= '9'))
            {
                char *end = 0;
                v.push_back(std::strtod(js.c_str() + i, &end));
                i = end - js.c_str() - 1;
            }
        }
        return false;
    }

    bool Load(const std::string &path, bool data, std::string *err = 0)
    {
        ok = false;
        std::ifstream in(path.c_str());
        if (!in) { if (err) *err = "cannot open " + path; return false; }
        std::stringstream ss;
        ss << in.rdbuf();
        const std::string js = ss.str();
        double x;
        Str(js, "version", version);
        if (Num(js, "eps", x)) eps = x;
        if (Num(js, "max_points", x)) max_points = (int)x;
        if (Num(js, "kmax", x)) kmax = (int)x;
        std::vector<double> flat;
        bool good = Arr(js, "S_pairE_edges", eS) && Arr(js, "S_mean", sMean) && Arr(js, "S_sigma", sSigma) &&
                    Arr(js, "S_alpha", sAlpha) && Arr(js, "S_n", sN) && Arr(js, "B_pairE_edges", eB) &&
                    Arr(js, "B_mass_edges", mB) && Arr(js, "B_pdf", flat) && Arr(js, "prior_n_lo", nLo) &&
                    Arr(js, "prior_jetE_edges", eJ) && Arr(js, "prior", prior);
        if (good && data)
        {
            std::vector<double> dm, ds;
            if (Arr(js, "S_data_mean", dm) && dm.size() == sMean.size()) sMean = dm;
            if (Arr(js, "S_data_sigma", ds) && ds.size() == sSigma.size()) sSigma = ds;
        }
        const size_t nS = eS.size() - 1, nBe = eB.size() - 1, nm = mB.size() - 1;
        good = good && sMean.size() == nS && sSigma.size() == nS && sAlpha.size() == nS && sN.size() == nS &&
               flat.size() == nBe * nm && prior.size() == nLo.size() * (eJ.size() - 1) * (size_t)(kmax + 1);
        if (!good) { if (err) *err = "malformed pi0 probability calibration " + path; return false; }
        if (Arr(js, "lead_q_n_lo", qNLo) && Arr(js, "lead_q", q) && q.size() != qNLo.size() * (eJ.size() - 1))
        { if (err) *err = "malformed lead_q in " + path; return false; }
        bPdf.assign(nBe, std::vector<double>(nm));
        for (size_t b = 0; b < nBe; b++)
            for (size_t i = 0; i < nm; i++) bPdf[b][i] = flat[b * nm + i];
        Normalise();
        ok = true;
        return true;
    }
};

// m^2 = (p_i + p_j)^2 of two massless-ish 4-vectors
inline double Pi0PairMass(double e1, double x1, double y1, double z1, double e2, double x2, double y2, double z2)
{
    const double E = e1 + e2, px = x1 + x2, py = y1 + y2, pz = z1 + z2;
    const double m2 = E * E - px * px - py * py - pz * pz;
    return m2 > 0 ? std::sqrt(m2) : 0;
}

// e_k for k = 0..n/2: sums over the k-matchings of the product of r over
// their pairs. Exact; memoised over the subsets of points still unpaired
// (the lowest unpaired point is either left alone or paired with each other
// one). Memory ~ 2^n (n/2 + 1) doubles: 0.1 MB at n = 12, 90 MB at n = 20.
inline void Pi0MatchingSums(const std::vector<std::vector<double> > &r, std::vector<double> &e, std::vector<double> *e0free = 0)
{
    const int n = (int)r.size(), K = n / 2 + 1;
    const unsigned full = (n >= 32) ? 0u : ((1u << n) - 1u);
    std::vector<double> memo((size_t)(full + 1) * K, 0.0);
    std::vector<char> done(full + 1, 0);
    // iterate masks in increasing order: every submask is smaller, so done
    for (unsigned mask = 0; mask <= full; mask++)
    {
        double *out = &memo[(size_t)mask * K];
        if (mask == 0) { out[0] = 1; done[0] = 1; continue; }
        int i = 0;
        while (!((mask >> i) & 1u)) i++;
        const unsigned rest = mask & ~(1u << i);
        const double *a = &memo[(size_t)rest * K];
        for (int k = 0; k < K; k++) out[k] = a[k];
        for (int j = i + 1; j < n; j++)
        {
            if (!((rest >> j) & 1u)) continue;
            const double *b = &memo[(size_t)(rest & ~(1u << j)) * K];
            const double rij = r[i][j];
            for (int k = 0; k + 1 < K; k++) out[k + 1] += rij * b[k];
        }
        done[mask] = 1;
    }
    e.assign(memo.begin() + (size_t)full * K, memo.begin() + (size_t)full * K + K);
    // the matchings leaving point 0 unpaired: the sums over the other points
    if (e0free) e0free->assign(memo.begin() + (size_t)(full & ~1u) * K, memo.begin() + (size_t)(full & ~1u) * K + K);
}

// log of the number of k-matchings of n points, n! / (k! 2^k (n-2k)!)
inline double Pi0LogNk(int n, int k)
{
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - k * std::log(2.0) - std::lgamma(n - 2.0 * k + 1.0);
}

// The posterior from the pair likelihood ratios r[i][j] (symmetric, n x n,
// leading point at index 0); res carries npairs and the best pair already.
// Used by Pi0Prob and the toy test.
inline Pi0ProbResult Pi0ProbFromRatios(const Pi0ProbCalib &c, const std::vector<std::vector<double> > &r, double ejet, Pi0ProbResult res)
{
    const int n = (int)r.size();
    res.status = kPi0ProbOk;
    const double p0 = c.Prior(n, ejet, 0);
    if (!(p0 > 0)) { res.prob = 1; res.logodds = 99; res.lead_pair = -1; return res; } // (floored away in the JSON)
    std::vector<double> e, ef;
    Pi0MatchingSums(r, e, &ef);
    // odds relative to the empty matching; the part of each k term with point 0 paired
    double odds = 0, oddsLead = 0;
    for (int k = 1; k < (int)e.size(); k++)
    {
        const double pk = c.Prior(n, ejet, k);
        if (!(pk > 0) || !(e[k] > 0)) continue;
        const double w = pk / p0 * std::exp(-Pi0LogNk(n, k));
        odds += w * e[k];
        const double el = e[k] - (k < (int)ef.size() ? ef[k] : 0.0);
        if (el > 0) oddsLead += w * el;
    }
    res.logodds = odds > 0 ? (float)std::log(odds) : -99.f;
    res.prob = (float)(odds / (1 + odds));
    res.lead_pair = (float)(oddsLead / (1 + odds));
    if (c.hasLead()) res.lead = (float)(res.lead_pair + (1 - res.lead_pair) * c.Q(n, ejet));
    return res;
}

// The probabilities for one jet. E, px, py, pz: its points (ECal point
// constituents, any order; best_i/j index them), ejet: the jet energy the
// prior is binned in (raw reco_E).
inline Pi0ProbResult Pi0Prob(const Pi0ProbCalib &c, const std::vector<double> &E, const std::vector<double> &px,
                             const std::vector<double> &py, const std::vector<double> &pz, double ejet)
{
    Pi0ProbResult res;
    const int n = (int)E.size();
    res.npairs = n * (n - 1) / 2;
    if (!c.ok) { res.status = kPi0ProbNoCalib; res.prob = res.logodds = res.lead_pair = res.lead = -1; return res; }
    if (n < 2)
    {
        res.status = kPi0ProbFewPoints; res.prob = 0; res.lead_pair = 0;
        res.lead = (n == 1 && c.hasLead()) ? (float)c.Q(1, ejet) : -1.f;
        return res;
    }
    if (n > c.max_points) { res.status = kPi0ProbTooMany; res.prob = res.logodds = res.lead_pair = res.lead = -1; return res; }
    // points by decreasing energy: the leading point is index 0
    std::vector<int> o(n);
    for (int i = 0; i < n; i++) o[i] = i;
    for (int i = 1; i < n; i++)
        for (int j = i; j > 0 && E[o[j]] > E[o[j - 1]]; j--) std::swap(o[j], o[j - 1]);
    std::vector<std::vector<double> > r(n, std::vector<double>(n, 0.0));
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
        {
            const int i = o[a], j = o[b];
            const double m = Pi0PairMass(E[i], px[i], py[i], pz[i], E[j], px[j], py[j], pz[j]), ep = E[i] + E[j];
            r[a][b] = r[b][a] = c.S(m, ep) / c.B(m, ep);
            if (r[a][b] > res.best_r) { res.best_r = (float)r[a][b]; res.best_i = std::min(i, j); res.best_j = std::max(i, j); }
        }
    return Pi0ProbFromRatios(c, r, ejet, res);
}

#endif
