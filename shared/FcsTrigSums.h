#ifndef FCS_TRIG_SUMS_H
#define FCS_TRIG_SUMS_H

// Per-jet FCS trigger quantities: what the Run 22 trigger algorithm compared
// with its thresholds, rebuilt from offline hit energies. The SimpleTree
// (StSimpleReaderMaker) holds each event's trigger 2x2 E_T sums in the
// stage-2 layout (trg_e2x2[ns][16][10], trg_h2x2[ns][10][6]) and their
// centres (tree trg_geom); this header forms, as stage_2_202207.cxx
// (StRoot/RTS/src/TRG_FCS):
//   - stage-1 truncation: each 2x2 in counts of kTrgGeVPerCount (0.0316 GeV
//     E_T), rounded down, saturating at 255
//   - EM: overlapping ECal 4x4 windows (2x2 of 2x2s, [15][9] per side); the
//     HCal 4x4s ([9][5]) nearest (EtoHmap, the 202203 E+H sum) and around
//     (EtoH3map, the 202207 veto: largest <= ECal / 4)
//   - JP: the five fixed patches per side (A, B, C, D, E), ECal + HCal 2x2s
//     with the B/C cutouts near the beam
// Per jet: the largest EM window whose centre is within R of the jet axis
// (em_et) with that window's HCal (em_hmax, em_hnear); and, among the
// patches containing the jet axis (the ECal 2x2 nearest the axis lies in
// the patch's ECal range), the largest patch sum (jp_et), which patch
// (jp_patch = ns*5 + i, i = 0 A ... 4 E) and all containing patches
// (jp_mask bit ns*5 + i). JPA = patch 0, JPBC = 1 or 2, JPDE = 3 or 4. The
// *_raw values skip the truncation. GeV E_T throughout; -1 if not computed.

#include <cmath>

static const double kTrgGeVPerCount = 0.0316;

namespace FcsTrg {
  static const int kTrgEtoHmap[15][9][2] = {
    { { 0, 0},{ 0, 1},{ 0, 1},{ 0, 2},{ 0, 2},{ 0, 3},{ 0, 4},{ 0, 4},{ 0, 4}},
    { { 0, 0},{ 0, 1},{ 0, 1},{ 0, 2},{ 0, 2},{ 0, 3},{ 0, 4},{ 0, 4},{ 0, 4}},
    { { 1, 0},{ 1, 1},{ 1, 1},{ 1, 2},{ 1, 2},{ 1, 3},{ 1, 4},{ 1, 4},{ 1, 4}},
    { { 2, 0},{ 2, 1},{ 2, 1},{ 2, 2},{ 2, 2},{ 2, 3},{ 2, 4},{ 2, 4},{ 2, 4}},
    { { 2, 0},{ 2, 1},{ 2, 1},{ 2, 2},{ 2, 2},{ 2, 3},{ 2, 4},{ 2, 4},{ 2, 4}},
    { { 3, 0},{ 3, 1},{ 3, 1},{ 3, 2},{ 3, 2},{ 3, 3},{ 3, 4},{ 3, 4},{ 3, 4}},
    { { 3, 0},{ 3, 1},{ 3, 1},{ 3, 2},{ 3, 2},{ 3, 3},{ 3, 4},{ 3, 4},{ 3, 4}},
    { { 4, 0},{ 4, 1},{ 4, 1},{ 4, 2},{ 4, 2},{ 4, 3},{ 4, 4},{ 4, 4},{ 4, 4}},
    { { 5, 0},{ 5, 1},{ 5, 1},{ 5, 2},{ 5, 2},{ 5, 3},{ 5, 4},{ 5, 4},{ 5, 4}},
    { { 5, 0},{ 5, 1},{ 5, 1},{ 5, 2},{ 5, 2},{ 5, 3},{ 5, 4},{ 5, 4},{ 5, 4}},
    { { 6, 0},{ 6, 1},{ 6, 1},{ 6, 2},{ 6, 2},{ 6, 3},{ 6, 4},{ 6, 4},{ 6, 4}},
    { { 6, 0},{ 6, 1},{ 6, 1},{ 6, 2},{ 6, 2},{ 6, 3},{ 6, 4},{ 6, 4},{ 6, 4}},
    { { 7, 0},{ 7, 1},{ 7, 1},{ 7, 2},{ 7, 2},{ 7, 3},{ 7, 4},{ 7, 4},{ 7, 4}},
    { { 8, 0},{ 8, 1},{ 8, 1},{ 8, 2},{ 8, 2},{ 8, 3},{ 8, 4},{ 8, 4},{ 8, 4}},
    { { 8, 0},{ 8, 1},{ 8, 1},{ 8, 2},{ 8, 2},{ 8, 3},{ 8, 4},{ 8, 4},{ 8, 4}}
  } ;
  static const int kTrgEtoH3map[15][9][4] = {
    {{-1,-1,-1, 0},{-1,-1, 0, 1},{-1,-1, 1, 2},{-1,-1, 1, 2},{-1,-1, 2, 3},{-1,-1, 2, 3},{-1,-1, 3, 4},{-1,-1, 4,-1},{-1,-1, 4,-1}},
    {{-1, 0,-1, 5},{ 0, 1, 5, 6},{ 1, 2, 6, 7},{ 1, 2, 6, 7},{ 2, 3, 7, 8},{ 2, 3, 7, 8},{ 3, 4, 8, 9},{ 4,-1, 9,-1},{ 4,-1, 9,-1}},
    {{-1, 0,-1, 5},{ 0, 1, 5, 6},{ 1, 2, 6, 7},{ 1, 2, 6, 7},{ 2, 3, 7, 8},{ 2, 3, 7, 8},{ 3, 4, 8, 9},{ 4,-1, 9,-1},{ 4,-1, 9,-1}},
    {{-1, 5,-1,10},{ 5, 6,10,11},{ 6, 7,11,12},{ 6, 7,11,12},{ 7, 8,12,13},{ 7, 8,12,13},{ 8, 9,13,14},{ 9,-1,14,-1},{ 9,-1,14,-1}},
    {{-1,10,-1,15},{10,11,15,16},{11,12,16,17},{11,12,16,17},{12,13,17,18},{12,13,17,18},{13,14,18,19},{14,-1,19,-1},{14,-1,19,-1}},
    {{-1,10,-1,15},{10,11,15,16},{11,12,16,17},{11,12,16,17},{12,13,17,18},{12,13,17,18},{13,14,18,19},{14,-1,19,-1},{14,-1,19,-1}},
    {{-1,15,-1,20},{15,16,20,21},{16,17,21,22},{16,17,21,22},{17,18,22,23},{17,18,22,23},{18,19,23,24},{19,-1,24,-1},{19,-1,24,-1}},
    {{-1,20,-1,25},{20,21,25,26},{21,22,26,27},{21,22,26,27},{22,23,27,28},{22,23,27,28},{23,24,28,29},{24,-1,29,-1},{24,-1,29,-1}},
    {{-1,20,-1,25},{20,21,25,26},{21,22,26,27},{21,22,26,27},{22,23,27,28},{22,23,27,28},{23,24,28,29},{24,-1,29,-1},{24,-1,29,-1}},
    {{-1,25,-1,30},{25,26,30,31},{26,27,31,32},{26,27,31,32},{27,28,32,33},{27,28,32,33},{28,29,33,34},{29,-1,34,-1},{29,-1,34,-1}},
    {{-1,25,-1,30},{25,26,30,31},{26,27,31,32},{26,27,31,32},{27,28,32,33},{27,28,32,33},{28,29,33,34},{29,-1,34,-1},{29,-1,34,-1}},
    {{-1,30,-1,35},{30,31,35,36},{31,32,36,37},{31,32,36,37},{32,33,37,38},{32,33,37,38},{33,34,38,39},{34,-1,39,-1},{34,-1,39,-1}},
    {{-1,35,-1,40},{35,36,40,41},{36,37,41,42},{36,37,41,42},{37,38,42,43},{37,38,42,43},{38,39,43,44},{39,-1,44,-1},{39,-1,44,-1}},
    {{-1,35,-1,40},{35,36,40,41},{36,37,41,42},{36,37,41,42},{37,38,42,43},{37,38,42,43},{38,39,43,44},{39,-1,44,-1},{39,-1,44,-1}},
    {{-1,40,-1,-1},{40,41,-1,-1},{41,42,-1,-1},{41,42,-1,-1},{42,43,-1,-1},{42,43,-1,-1},{43,44,-1,-1},{44,-1,-1,-1},{44,-1,-1,-1}}
  };
  static const int kECol0[5] = {0, 0, 0, 3, 3}, kECol1[5] = {3, 6, 6, 9, 9};
  static const int kERow0[5] = {3, 0, 5, 0, 5}, kERow1[5] = {12, 10, 15, 10, 15};
  static const int kHCol0[5] = {0, 0, 0, 2, 2}, kHCol1[5] = {2, 3, 3, 5, 5};
  static const int kHRow0[5] = {2, 0, 3, 0, 3}, kHRow1[5] = {7, 6, 9, 6, 9};
}

struct FcsTrgGeom
{
    float ex[2][16][10], ey[2][16][10], ez[2][16][10];
    float hx[2][10][6], hy[2][10][6], hz[2][10][6];
    bool ok = false;
};

struct FcsTrgJet
{
    float em_et = -1, em_hmax = -1, em_hnear = -1, em_et_raw = -1;
    float jp_et = -1, jp_et_raw = -1;
    int jp_patch = -1, jp_mask = 0;
};

// counts of one 2x2 after stage 1
inline int FcsTrgCounts(float et)
{
    if (!(et > 0)) return 0;
    const double n = std::floor(et / kTrgGeVPerCount);
    return n > 255 ? 255 : (int)n;
}

inline void FcsTrgEtaPhi(double x, double y, double z, double &eta, double &phi)
{
    const double rho = std::sqrt(x * x + y * y);
    eta = rho > 0 ? std::asinh(z / rho) : 0;
    phi = std::atan2(y, x);
}

inline double FcsTrgDR(double eta1, double phi1, double eta2, double phi2)
{
    double dphi = std::fabs(phi1 - phi2);
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    return std::sqrt((eta1 - eta2) * (eta1 - eta2) + dphi * dphi);
}

// e2x2/h2x2: one event's sums (GeV E_T); jet axis eta, phi; cone radius R
inline FcsTrgJet FcsTrgForJet(const float e2x2[2][16][10], const float h2x2[2][10][6], const FcsTrgGeom &g,
                              double jeta, double jphi, double R)
{
    using namespace FcsTrg;
    FcsTrgJet out;
    if (!g.ok) return out;
    // counts per 2x2
    int ec[2][16][10], hc[2][10][6];
    for (int ns = 0; ns < 2; ns++)
    {
        for (int r = 0; r < 16; r++) for (int c = 0; c < 10; c++) ec[ns][r][c] = FcsTrgCounts(e2x2[ns][r][c]);
        for (int r = 0; r < 10; r++) for (int c = 0; c < 6; c++) hc[ns][r][c] = FcsTrgCounts(h2x2[ns][r][c]);
    }
    // EM windows with centres inside the cone
    int best = -1; double best_raw = -1;
    for (int ns = 0; ns < 2; ns++)
    {
        int hs[9][5];
        for (int r = 0; r < 9; r++) for (int c = 0; c < 5; c++)
            hs[r][c] = hc[ns][r][c] + hc[ns][r][c + 1] + hc[ns][r + 1][c] + hc[ns][r + 1][c + 1];
        for (int r = 0; r < 15; r++) for (int c = 0; c < 9; c++)
        {
            const double x = 0.25 * (g.ex[ns][r][c] + g.ex[ns][r][c + 1] + g.ex[ns][r + 1][c] + g.ex[ns][r + 1][c + 1]);
            const double y = 0.25 * (g.ey[ns][r][c] + g.ey[ns][r][c + 1] + g.ey[ns][r + 1][c] + g.ey[ns][r + 1][c + 1]);
            const double z = 0.25 * (g.ez[ns][r][c] + g.ez[ns][r][c + 1] + g.ez[ns][r + 1][c] + g.ez[ns][r + 1][c + 1]);
            if (!(z > 0)) continue;
            double weta, wphi;
            FcsTrgEtaPhi(x, y, z, weta, wphi);
            if (FcsTrgDR(weta, wphi, jeta, jphi) >= R) continue;
            const int esum = ec[ns][r][c] + ec[ns][r][c + 1] + ec[ns][r + 1][c] + ec[ns][r + 1][c + 1];
            const double eraw = e2x2[ns][r][c] + e2x2[ns][r][c + 1] + e2x2[ns][r + 1][c] + e2x2[ns][r + 1][c + 1];
            if (eraw > best_raw) best_raw = eraw;
            if (esum > best)
            {
                best = esum;
                int hmax = 0;
                for (int k = 0; k < 4; k++)
                {
                    const int id = kTrgEtoH3map[r][c][k];
                    if (id < 0) continue;
                    if (hs[id / 5][id % 5] > hmax) hmax = hs[id / 5][id % 5];
                }
                out.em_et = esum * kTrgGeVPerCount;
                out.em_hmax = hmax * kTrgGeVPerCount;
                out.em_hnear = hs[kTrgEtoHmap[r][c][0]][kTrgEtoHmap[r][c][1]] * kTrgGeVPerCount;
            }
        }
    }
    if (best_raw >= 0) out.em_et_raw = best_raw;
    // JP: the ECal 2x2 nearest the jet axis decides which patches contain the jet
    int nns = -1, nr = -1, ncl = -1; double dmin = 1e9;
    for (int ns = 0; ns < 2; ns++) for (int r = 0; r < 16; r++) for (int c = 0; c < 10; c++)
    {
        if (!(g.ez[ns][r][c] > 0)) continue;
        double ueta, uphi;
        FcsTrgEtaPhi(g.ex[ns][r][c], g.ey[ns][r][c], g.ez[ns][r][c], ueta, uphi);
        const double d = FcsTrgDR(ueta, uphi, jeta, jphi);
        if (d < dmin) { dmin = d; nns = ns; nr = r; ncl = c; }
    }
    if (nns < 0) return out;
    for (int i = 0; i < 5; i++)
    {
        if (nr < kERow0[i] || nr > kERow1[i] || ncl < kECol0[i] || ncl > kECol1[i]) continue;
        out.jp_mask |= 1 << (nns * 5 + i);
        int sum = 0; double raw = 0;
        for (int c = kECol0[i]; c <= kECol1[i]; c++) for (int r = kERow0[i]; r <= kERow1[i]; r++)
        {
            if ((i == 1 || i == 2) && (c == 0 || c == 1) && (r >= 5 && r <= 10)) continue;  // cutout
            sum += ec[nns][r][c]; raw += e2x2[nns][r][c];
        }
        for (int c = kHCol0[i]; c <= kHCol1[i]; c++) for (int r = kHRow0[i]; r <= kHRow1[i]; r++)
        {
            if ((i == 1 || i == 2) && c == 0 && (r >= 3 && r <= 6)) continue;  // cutout
            sum += hc[nns][r][c]; raw += h2x2[nns][r][c];
        }
        if (sum * kTrgGeVPerCount > out.jp_et)
        {
            out.jp_et = sum * kTrgGeVPerCount;
            out.jp_et_raw = raw;
            out.jp_patch = nns * 5 + i;
        }
    }
    return out;
}

#endif
