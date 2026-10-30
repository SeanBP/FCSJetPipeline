// Run QA of the EM2 trigger: real fcsEM2 events with no hard offline ECal
// cluster ("empty": highest cluster E_T < 1.5 GeV) against the rest, per run.
// Usage: ./empty_em2_events <Liang depch_id.txt> <SimpleTree list> [...]
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
// Real fcsEM2 events with no hard offline ECal cluster (highest cluster E_T < 1.5 GeV) against the
// rest: fraction per file list, the tower and DEP board of the event's largest ECal hit, and
// ECal / HCal energy totals. depch: Liang's depch_id.txt ("ehp ns dep ch id").
int main(int argc, char **argv)
{
    std::map<std::pair<int, int>, std::pair<int, int>> id2dep; // (ns,id) -> (dep,ch)
    { std::ifstream d(argv[1]); int e, ns, dep, ch, id; while (d >> e >> ns >> dep >> ch >> id) id2dep[{ns, id}] = {dep, ch}; }
    for (int a = 2; a < argc; a++)
    {
        std::ifstream in(argv[a]); std::string f;
        long nem2 = 0, nempty = 0;
        std::map<std::pair<int, int>, long> twr, twrAll, dep, depAll;
        double eE[2] = {0}, eH[2] = {0}, maxh[2] = {0}, nhit[2] = {0}, n[2] = {0};
        while (in >> f)
        {
            TFile x(f.c_str()); TTree *t = (TTree *)x.Get("data");
            static int nh, hd[5000], hid[5000], nc, cd[200], em2;
            static float he[5000], ce[200], cx[200], cy[200], cz[200];
            t->SetBranchStatus("*", 0);
            for (const char *b : {"Cal_nhits", "Cal_detid", "Cal_hitid", "Cal_hit_energy", "Cal_nclus", "Cal_clus_detid", "Cal_clus_energy", "Cal_clus_x", "Cal_clus_y", "Cal_clus_z", "fcsEM2"}) t->SetBranchStatus(b, 1);
            t->SetBranchAddress("Cal_nhits", &nh); t->SetBranchAddress("Cal_detid", hd); t->SetBranchAddress("Cal_hitid", hid); t->SetBranchAddress("Cal_hit_energy", he);
            t->SetBranchAddress("Cal_nclus", &nc); t->SetBranchAddress("Cal_clus_detid", cd); t->SetBranchAddress("Cal_clus_energy", ce);
            t->SetBranchAddress("Cal_clus_x", cx); t->SetBranchAddress("Cal_clus_y", cy); t->SetBranchAddress("Cal_clus_z", cz);
            t->SetBranchAddress("fcsEM2", &em2);
            for (Long64_t i = 0; i < t->GetEntries(); i++)
            {
                t->GetEntry(i);
                if (!em2) continue;
                nem2++;
                double m = 0;
                for (int k = 0; k < nc && k < 200; k++) if (cd[k] <= 1) m = std::max(m, (double)ce[k] * std::hypot(cx[k], cy[k]) / std::sqrt(cx[k] * cx[k] + cy[k] * cy[k] + cz[k] * cz[k]));
                const int cls = m < 1.5 ? 0 : 1;
                if (cls == 0) nempty++;
                int kmax = -1; double emax = 0, se = 0, sh = 0, nn = 0;
                for (int k = 0; k < nh && k < 5000; k++)
                {
                    if (hd[k] <= 1) { se += he[k]; nn++; if (he[k] > emax) { emax = he[k]; kmax = k; } }
                    else if (hd[k] <= 3) sh += he[k];
                }
                eE[cls] += se; eH[cls] += sh; maxh[cls] += emax; nhit[cls] += nn; n[cls]++;
                if (kmax >= 0)
                {
                    const std::pair<int, int> key(hd[kmax], hid[kmax]);
                    (cls == 0 ? twr : twrAll)[key]++;
                    auto it = id2dep.find(key);
                    if (it != id2dep.end()) (cls == 0 ? dep : depAll)[{hd[kmax], it->second.first}]++;
                }
            }
        }
        std::printf("%s\n  real EM2 events %ld, with highest ECal cluster E_T < 1.5 GeV: %ld (%.1f%%)\n", argv[a], nem2, nempty, 100.0 * nempty / std::max(1L, nem2));
        const char *nm[2] = {"empty", "normal"};
        for (int c = 0; c < 2; c++)
            if (n[c] > 0) std::printf("  %-6s: mean ECal E %.2f GeV, HCal E %.2f GeV, largest ECal hit %.2f GeV, ECal hits %.1f\n", nm[c], eE[c] / n[c], eH[c] / n[c], maxh[c] / n[c], nhit[c] / n[c]);
        auto top = [](std::map<std::pair<int, int>, long> &mp, long tot, const char *what) {
            std::vector<std::pair<long, std::pair<int, int>>> v;
            for (auto &kv : mp) v.push_back({kv.second, kv.first});
            std::sort(v.rbegin(), v.rend());
            std::printf("  %s:", what);
            for (size_t i = 0; i < v.size() && i < 6; i++) std::printf("  (ns %d, %d) %.1f%%", v[i].second.first, v[i].second.second, 100.0 * v[i].first / std::max(1L, tot));
            std::printf("  [%zu distinct]\n", v.size());
        };
        top(twr, nempty, "empty events, tower of the largest hit (ns, id)");
        top(twrAll, nem2 - nempty, "normal events, tower of the largest hit");
        top(dep, nempty, "empty events, DEP board of the largest hit (ns, dep)");
        top(depAll, nem2 - nempty, "normal events, DEP board of the largest hit");
    }
}
