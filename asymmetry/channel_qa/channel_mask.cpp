// Hot / dead FCS channel masks per period (Sean 2026-11-24) from the
// channel_occupancy outputs of a sample (submit_an_job.sh -Q 1 on
// make_mask_sample.py's list).
//
// Per period, the occ_*.root of its sample runs are summed. For each channel
// of ECal N/S and HCal N/S:
//   rate   = hits above the tower threshold per event
//   meanE  = mean energy of those hits
// compared with its neighbours in the same detector (centres within 1.5 tower
// pitches: ECal 5.57 cm, HCal 9.99 cm; the median of the neighbours, so one bad
// neighbour does not move it). A channel is flagged
//   hot     rate / neighbour median > kHotRate (and > kSig sigma above it)
//   dead    no hits above 0 at all, or rate / neighbour median < kDeadRate
//   hotE    meanE / neighbour median > kHotE (gain too high or a noisy ADC)
//   coldE   meanE / neighbour median < kColdE (gain too low)
// Channels with fewer than kMinNb neighbours (corners) are compared with the
// median of their ring instead (same distance from the beam within half a
// pitch). The thresholds are first values to tune by eye on the plots.
//
// Outputs in <outdir>:
//   channel_mask.txt     one row per flagged channel and period:
//                        period first_run last_run det id x y flag rate_ratio e_ratio n_thr
//   channel_mask_summary.txt  flagged channels per period and detector
//   channel_ratio_<det>_p<k>.png  per period: rate and mean-energy ratio maps
//
//   g++ -m64 -O2 channel_mask.cpp -o channel_mask `root-config --cflags --libs` -std=c++11
//   ./channel_mask <periods.txt> <runs.txt> <occ dir> <outdir>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

static const int kNDet = 4;
static const int kNCh[kNDet] = {748, 748, 260, 260};
static const char *kDetName[kNDet] = {"ecal_n", "ecal_s", "hcal_n", "hcal_s"};
static const double kPitch[kNDet] = {5.572, 5.572, 9.99, 9.99};
static const double kHotRate = 3.0, kDeadRate = 0.2, kHotE = 2.0, kColdE = 0.5, kSig = 5.0;
static const int kMinNb = 3;

struct Period { int idx, first, last; };

static double Median(std::vector<double> v)
{
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    if (argc < 5) { std::fprintf(stderr, "usage: ./channel_mask <periods.txt> <runs.txt> <occ dir> <outdir>\n"); return 1; }
    std::vector<Period> periods;
    {
        std::ifstream f(argv[1]); std::string l;
        while (std::getline(f, l)) { if (l.empty() || l[0] == '#') continue; std::istringstream s(l); Period p; s >> p.idx >> p.first >> p.last; periods.push_back(p); }
    }
    const std::string occdir = argv[3], outdir = argv[4];
    gSystem->mkdir(outdir.c_str(), kTRUE);
    // occ files by run: the SimpleTree suffix holds the 8-digit run number
    std::map<int, std::vector<std::string> > files;
    {
        void *d = gSystem->OpenDirectory(occdir.c_str());
        const char *e;
        while ((e = gSystem->GetDirEntry(d)))
        {
            const std::string n = e;
            if (n.compare(0, 4, "occ_") != 0 || n.size() < 6 || n.compare(n.size() - 5, 5, ".root") != 0) continue;
            for (size_t i = 0; i + 8 <= n.size(); i++)
                if (n[i] == '2' && n[i + 1] == '3' && std::all_of(n.begin() + i, n.begin() + i + 8, ::isdigit) && (i + 8 == n.size() || !::isdigit(n[i + 8])))
                { files[std::atoi(n.substr(i, 8).c_str())].push_back(occdir + "/" + n); break; }
        }
        gSystem->FreeDirectory(d);
    }
    FILE *fm = std::fopen((outdir + "/channel_mask.txt").c_str(), "w");
    FILE *fs = std::fopen((outdir + "/channel_mask_summary.txt").c_str(), "w");
    std::fprintf(fm, "# hot / dead channel mask (channel_mask.cpp; hot rate > %.1f x neighbours and > %.0f sigma, dead rate < %.1f x or no hits, hotE > %.1f x, coldE < %.1f x)\n",
                 kHotRate, kSig, kDeadRate, kHotE, kColdE);
    std::fprintf(fm, "# period first_run last_run det id x y flag rate_ratio e_ratio n_thr\n");
    std::fprintf(fs, "# period first_run last_run files events  ecal_n ecal_s hcal_n hcal_s  (flagged channels: hot dead hotE coldE)\n");
    for (const Period &p : periods)
    {
        std::vector<std::string> pf;
        for (auto &kv : files) if (kv.first >= p.first && kv.first <= p.last) pf.insert(pf.end(), kv.second.begin(), kv.second.end());
        if (pf.empty()) { std::fprintf(fs, "%d %d %d 0 0  no sample files\n", p.idx, p.first, p.last); continue; }
        std::vector<double> any[kNDet], thr[kNDet], se[kNDet], px[kNDet], py[kNDet], np[kNDet];
        for (int d = 0; d < kNDet; d++) { any[d].assign(kNCh[d], 0); thr[d] = se[d] = px[d] = py[d] = np[d] = any[d]; }
        double nev = 0;
        for (const auto &fn : pf)
        {
            TFile f(fn.c_str());
            if (!f.IsOpen()) continue;
            TH1D *he = (TH1D *)f.Get("n_events");
            if (!he) continue;
            nev += he->GetBinContent(1);
            for (int d = 0; d < kNDet; d++)
            {
                TH1D *ha = (TH1D *)f.Get(Form("n_any_%s", kDetName[d])), *ht = (TH1D *)f.Get(Form("n_thr_%s", kDetName[d])),
                     *hs = (TH1D *)f.Get(Form("sumE_%s", kDetName[d])), *hx = (TH1D *)f.Get(Form("posx_%s", kDetName[d])),
                     *hy = (TH1D *)f.Get(Form("posy_%s", kDetName[d])), *hn = (TH1D *)f.Get(Form("npos_%s", kDetName[d]));
                if (!ha || !ht || !hs || !hx || !hy || !hn) continue;
                for (int c = 0; c < kNCh[d]; c++)
                {
                    any[d][c] += ha->GetBinContent(c + 1); thr[d][c] += ht->GetBinContent(c + 1); se[d][c] += hs->GetBinContent(c + 1);
                    px[d][c] += hx->GetBinContent(c + 1); py[d][c] += hy->GetBinContent(c + 1); np[d][c] += hn->GetBinContent(c + 1);
                }
            }
        }
        int counts[kNDet][4] = {{0}};
        for (int d = 0; d < kNDet; d++)
        {
            const int n = kNCh[d];
            std::vector<double> x(n, NAN), y(n, NAN), rate(n, 0), me(n, NAN);
            for (int c = 0; c < n; c++)
            {
                if (np[d][c] > 0) { x[c] = px[d][c] / np[d][c]; y[c] = py[d][c] / np[d][c]; }
                rate[c] = nev > 0 ? thr[d][c] / nev : 0;
                if (thr[d][c] > 0) me[c] = se[d][c] / thr[d][c];
            }
            // per-channel ratios to the neighbour medians (or the ring's, at corners)
            std::vector<double> rr(n, NAN), er(n, NAN);
            for (int c = 0; c < n; c++)
            {
                if (!(x[c] == x[c])) continue;
                std::vector<double> nr, ne;
                for (int o = 0; o < n; o++)
                    if (o != c && x[o] == x[o] && std::hypot(x[o] - x[c], y[o] - y[c]) < 1.5 * kPitch[d]) { nr.push_back(rate[o]); if (me[o] == me[o]) ne.push_back(me[o]); }
                if ((int)nr.size() < kMinNb)
                {
                    nr.clear(); ne.clear();
                    const double r0 = std::hypot(x[c], y[c]);
                    for (int o = 0; o < n; o++)
                        if (o != c && x[o] == x[o] && std::fabs(std::hypot(x[o], y[o]) - r0) < 0.5 * kPitch[d]) { nr.push_back(rate[o]); if (me[o] == me[o]) ne.push_back(me[o]); }
                }
                const double mr = Median(nr), mE = Median(ne);
                if (mr > 0) rr[c] = rate[c] / mr;
                if (mE > 0 && me[c] == me[c]) er[c] = me[c] / mE;
                std::string flag;
                const double expect = mr * nev, sig = expect > 0 ? (thr[d][c] - expect) / std::sqrt(expect) : 0;
                if (any[d][c] == 0 || (rr[c] == rr[c] && rr[c] < kDeadRate)) flag = "dead";
                else if (rr[c] > kHotRate && sig > kSig) flag = "hot";
                else if (er[c] == er[c] && er[c] > kHotE && thr[d][c] >= 20) flag = "hotE";
                else if (er[c] == er[c] && er[c] < kColdE && thr[d][c] >= 20) flag = "coldE";
                if (flag.empty()) continue;
                counts[d][flag == "hot" ? 0 : flag == "dead" ? 1 : flag == "hotE" ? 2 : 3]++;
                std::fprintf(fm, "%d %d %d %d %d %.2f %.2f %s %.3f %.3f %.0f\n", p.idx, p.first, p.last, d, c, x[c], y[c], flag.c_str(), rr[c], er[c], thr[d][c]);
            }
            // maps of the two ratios
            TH2D *hr = new TH2D(Form("rr_%d_%d", d, p.idx), Form("%s, period %d (%d-%d): hit rate / neighbour median;x (cm);y (cm)", kDetName[d], p.idx, p.first, p.last),
                                d < 2 ? 22 : 13, 0, 1, d < 2 ? 34 : 20, 0, 1);
            TH2D *he = (TH2D *)hr->Clone(Form("er_%d_%d", d, p.idx));
            he->SetTitle(Form("%s, period %d: mean hit energy / neighbour median;x (cm);y (cm)", kDetName[d], p.idx));
            double xl = 1e9, xh = -1e9, yl = 1e9, yh = -1e9;
            for (int c = 0; c < n; c++) if (x[c] == x[c]) { xl = std::min(xl, x[c]); xh = std::max(xh, x[c]); yl = std::min(yl, y[c]); yh = std::max(yh, y[c]); }
            const int nx = d < 2 ? 22 : 13, ny = d < 2 ? 34 : 20;
            const double dx = (xh - xl) / (nx - 1), dy = (yh - yl) / (ny - 1);
            hr->SetBins(nx, xl - dx / 2, xh + dx / 2, ny, yl - dy / 2, yh + dy / 2);
            he->SetBins(nx, xl - dx / 2, xh + dx / 2, ny, yl - dy / 2, yh + dy / 2);
            for (int c = 0; c < n; c++) if (x[c] == x[c]) { if (rr[c] == rr[c]) hr->Fill(x[c], y[c], std::max(rr[c], 1e-3)); if (er[c] == er[c]) he->Fill(x[c], y[c], er[c]); }
            TCanvas cv("cv", "", 1300, 620);
            cv.Divide(2, 1);
            cv.cd(1); gPad->SetRightMargin(0.14); gPad->SetLogz(); hr->SetMinimum(0.1); hr->SetMaximum(10); hr->Draw("colz");
            cv.cd(2); gPad->SetRightMargin(0.14); he->SetMinimum(0.3); he->SetMaximum(3); he->Draw("colz");
            cv.SaveAs(Form("%s/channel_ratio_%s_p%02d.png", outdir.c_str(), kDetName[d], p.idx));
        }
        std::fprintf(fs, "%d %d %d %zu %.0f ", p.idx, p.first, p.last, pf.size(), nev);
        for (int d = 0; d < kNDet; d++) std::fprintf(fs, " %d/%d/%d/%d", counts[d][0], counts[d][1], counts[d][2], counts[d][3]);
        std::fprintf(fs, "\n");
    }
    std::fclose(fm);
    std::fclose(fs);
    std::printf("channel_mask: %zu periods -> %s/channel_mask.txt, channel_mask_summary.txt, channel_ratio_*.png\n", periods.size(), outdir.c_str());
    return 0;
}
