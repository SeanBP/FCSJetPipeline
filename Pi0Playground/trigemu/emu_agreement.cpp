// Agreement of the re-emulated EM bits (reemulate.C friend files) with the
// recorded ones, overall and split by where the event's highest-E_T ECal
// cluster sits: inside the jet fiducial region (kFiducialRectEcal shrunk by
// reco_fiducial_buffer_ecal), or outside it. Also by the radius of the
// event's largest ECal hit (edge: r < 50 cm or r > 110 cm), and, inside the
// fiducial region, by the top cluster's E_T.
//
// Usage: ./emu_agreement <SimpleTree list> <trgemu dir>
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include "JetParameters.h"

int main(int argc, char **argv)
{
    if (argc < 3) { std::printf("usage: ./emu_agreement <SimpleTree list> <trgemu dir>\n"); return 1; }
    std::ifstream in(argv[1]);
    const std::string dir = argv[2];
    // [region][level]: real, emulated, both
    const int kNG = 9;
    const char *rname[kNG] = {"all events", "top cluster in jet fiducial", "top cluster outside", "largest hit 50-110 cm", "largest hit at an edge",
                              "fiducial, top E_T < 3 GeV", "fiducial, top E_T 3-5 GeV", "fiducial, top E_T 5-7 GeV", "fiducial, top E_T > 7 GeV"};
    double real[kNG][3] = {{0}}, emu[kNG][3] = {{0}}, both[kNG][3] = {{0}}, nev[kNG] = {0};
    std::string f;
    while (in >> f)
    {
        TFile x(f.c_str());
        TTree *t = (TTree *)x.Get("data");
        t->AddFriend("trgemu", (dir + "/" + f.substr(f.rfind('/') + 1)).c_str());
        static int nh, hd[5000], nc, cd[200], r[3], e[3];
        static float he[5000], hx[5000], hy[5000], ce[200], cx[200], cy[200], cz[200];
        t->SetBranchAddress("Cal_nhits", &nh); t->SetBranchAddress("Cal_detid", hd); t->SetBranchAddress("Cal_hit_energy", he);
        t->SetBranchAddress("Cal_hit_posx", hx); t->SetBranchAddress("Cal_hit_posy", hy);
        t->SetBranchAddress("Cal_nclus", &nc); t->SetBranchAddress("Cal_clus_detid", cd); t->SetBranchAddress("Cal_clus_energy", ce);
        t->SetBranchAddress("Cal_clus_x", cx); t->SetBranchAddress("Cal_clus_y", cy); t->SetBranchAddress("Cal_clus_z", cz);
        for (int k = 0; k < 3; k++)
        {
            t->SetBranchAddress(Form("fcsEM%d", k), &r[k]);
            t->SetBranchAddress(Form("trgemu.fcsEM%d", k), &e[k]);
        }
        for (Long64_t i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);
            double best = -1; int kb = -1;
            for (int k = 0; k < nc && k < 200; k++)
            {
                if (cd[k] > 1) continue;
                const double et = ce[k] * std::hypot(cx[k], cy[k]) / std::sqrt(cx[k] * cx[k] + cy[k] * cy[k] + cz[k] * cz[k]);
                if (et > best) { best = et; kb = k; }
            }
            int kh = -1; double emax = 0;
            for (int k = 0; k < nh && k < 5000; k++) if (hd[k] <= 1 && he[k] > emax) { emax = he[k]; kh = k; }
            const bool fid = kb >= 0 && pass_fiducial_cut(cx[kb], cy[kb], reco_fiducial_buffer_ecal, kFiducialRectEcal);
            const double rh = kh >= 0 ? std::hypot(hx[kh], hy[kh]) : -1;
            const bool edge = rh >= 0 && (rh < 50 || rh > 110);
            const bool in[kNG] = {true, fid, !fid, rh >= 50 && rh <= 110, edge,
                                  fid && best < 3, fid && best >= 3 && best < 5, fid && best >= 5 && best < 7, fid && best >= 7};
            for (int g = 0; g < kNG; g++)
            {
                if (!in[g]) continue;
                nev[g]++;
                for (int k = 0; k < 3; k++) { real[g][k] += r[k]; emu[g][k] += e[k]; both[g][k] += r[k] && e[k]; }
            }
        }
    }
    std::printf("%-30s %8s", "region", "events");
    for (int k = 0; k < 3; k++) std::printf("   EM%d: P(real|emu) P(emu|real)", k);
    std::printf("\n");
    for (int g = 0; g < kNG; g++)
    {
        std::printf("%-30s %8.0f", rname[g], nev[g]);
        for (int k = 0; k < 3; k++) std::printf("        %5.3f        %5.3f", emu[g][k] > 0 ? both[g][k] / emu[g][k] : -1, real[g][k] > 0 ? both[g][k] / real[g][k] : -1);
        std::printf("\n");
    }
    return 0;
}
