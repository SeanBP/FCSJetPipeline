// Hot-tower scan for run QA: for each SimpleTree list (one run), the fraction
// of events in which each ECal tower is the event's largest ECal hit. Towers
// above a threshold fraction are flagged, with their position and the mean
// energy of the hit when it is the largest. Hit positions are those stored
// in the SimpleTree.
//
// Usage: ./hot_tower_scan <threshold fraction, e.g. 0.02> <label>:<list> [...]
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

int main(int argc, char **argv)
{
    if (argc < 3) { std::printf("usage: ./hot_tower_scan <threshold fraction> <label>:<list> ...\n"); return 1; }
    const double thr = std::atof(argv[1]);
    for (int a = 2; a < argc; a++)
    {
        const std::string arg = argv[a];
        const size_t c = arg.find(':');
        std::ifstream in(arg.substr(c + 1).c_str());
        std::string f;
        std::map<std::pair<int, int>, double> count, esum, px, py;
        long nev = 0;
        while (in >> f)
        {
            TFile x(f.c_str());
            TTree *t = (TTree *)x.Get("data");
            static int nh, hd[5000], hid[5000];
            static float he[5000], hx[5000], hy[5000];
            t->SetBranchStatus("*", 0);
            for (const char *b : {"Cal_nhits", "Cal_detid", "Cal_hitid", "Cal_hit_energy", "Cal_hit_posx", "Cal_hit_posy"}) t->SetBranchStatus(b, 1);
            t->SetBranchAddress("Cal_nhits", &nh); t->SetBranchAddress("Cal_detid", hd); t->SetBranchAddress("Cal_hitid", hid);
            t->SetBranchAddress("Cal_hit_energy", he); t->SetBranchAddress("Cal_hit_posx", hx); t->SetBranchAddress("Cal_hit_posy", hy);
            for (Long64_t i = 0; i < t->GetEntries(); i++)
            {
                t->GetEntry(i);
                int km = -1; double em = 0;
                for (int k = 0; k < nh && k < 5000; k++) if (hd[k] <= 1 && he[k] > em) { em = he[k]; km = k; }
                if (km < 0) continue;
                nev++;
                const std::pair<int, int> key(hd[km], hid[km]);
                count[key]++; esum[key] += em; px[key] = hx[km]; py[key] = hy[km];
            }
        }
        std::vector<std::pair<double, std::pair<int, int>>> v;
        for (auto &kv : count) v.push_back({kv.second, kv.first});
        std::sort(v.rbegin(), v.rend());
        std::printf("%s: %ld events, %zu towers ever the largest hit; flagged (> %.1f%% of events):\n", arg.substr(0, c).c_str(), nev, v.size(), 100 * thr);
        int nflag = 0;
        for (auto &e : v)
        {
            if (e.first / nev <= thr) break;
            const auto &k = e.second;
            std::printf("    det %d id %3d  (x %6.1f, y %6.1f cm, r %5.1f)  %5.2f%% of events, mean largest-hit E %.2f GeV\n",
                        k.first, k.second, px[k], py[k], std::hypot(px[k], py[k]), 100 * e.first / nev, esum[k] / e.first);
            nflag++;
        }
        if (!nflag) std::printf("    none; busiest tower %.2f%%\n", v.empty() ? 0 : 100 * v[0].first / nev);
    }
    return 0;
}
