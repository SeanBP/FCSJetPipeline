// Run QA of the EM2 trigger: radius of the largest ECal hit in "empty" and
// normal real fcsEM2 events (first 40 files of the list).
// Usage: ./empty_em2_radius <SimpleTree list>
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
// Radius (and |x|, y) of the largest ECal hit in real EM2 events, empty (highest cluster E_T < 1.5 GeV) vs normal
int main(int argc, char **argv)
{
    std::ifstream in(argv[1]); std::string f; int nf = 0;
    const double rb[8] = {0, 30, 40, 50, 66, 88, 110, 200};
    double cnt[2][7] = {{0}};
    while (in >> f && nf++ < 40)
    {
        TFile x(f.c_str()); TTree *t = (TTree *)x.Get("data");
        static int nh, hd[5000], nc, cd[200], em2; static float he[5000], hx[5000], hy[5000], ce[200], cx[200], cy[200], cz[200];
        t->SetBranchStatus("*", 0);
        for (const char *b : {"Cal_nhits", "Cal_detid", "Cal_hit_energy", "Cal_hit_posx", "Cal_hit_posy", "Cal_nclus", "Cal_clus_detid", "Cal_clus_energy", "Cal_clus_x", "Cal_clus_y", "Cal_clus_z", "fcsEM2"}) t->SetBranchStatus(b, 1);
        t->SetBranchAddress("Cal_nhits", &nh); t->SetBranchAddress("Cal_detid", hd); t->SetBranchAddress("Cal_hit_energy", he);
        t->SetBranchAddress("Cal_hit_posx", hx); t->SetBranchAddress("Cal_hit_posy", hy);
        t->SetBranchAddress("Cal_nclus", &nc); t->SetBranchAddress("Cal_clus_detid", cd); t->SetBranchAddress("Cal_clus_energy", ce);
        t->SetBranchAddress("Cal_clus_x", cx); t->SetBranchAddress("Cal_clus_y", cy); t->SetBranchAddress("Cal_clus_z", cz);
        t->SetBranchAddress("fcsEM2", &em2);
        for (Long64_t i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);
            if (!em2) continue;
            double m = 0;
            for (int k = 0; k < nc && k < 200; k++) if (cd[k] <= 1) m = std::max(m, (double)ce[k] * std::hypot(cx[k], cy[k]) / std::sqrt(cx[k] * cx[k] + cy[k] * cy[k] + cz[k] * cz[k]));
            int kmax = -1; double emax = 0;
            for (int k = 0; k < nh && k < 5000; k++) if (hd[k] <= 1 && he[k] > emax) { emax = he[k]; kmax = k; }
            if (kmax < 0) continue;
            const double r = std::hypot(hx[kmax], hy[kmax]);
            for (int b = 0; b < 7; b++) if (r >= rb[b] && r < rb[b + 1]) cnt[m < 1.5 ? 0 : 1][b]++;
        }
    }
    double s[2] = {0, 0};
    for (int c = 0; c < 2; c++) for (int b = 0; b < 7; b++) s[c] += cnt[c][b];
    std::printf("radius of the largest ECal hit (cm): share of empty | normal real-EM2 events\n");
    for (int b = 0; b < 7; b++) std::printf("  %3.0f-%-3.0f  %5.1f%% | %5.1f%%\n", rb[b], rb[b + 1], 100 * cnt[0][b] / s[0], 100 * cnt[1][b] / s[1]);
}
