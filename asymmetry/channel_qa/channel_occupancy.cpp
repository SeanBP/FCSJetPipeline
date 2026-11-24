// Per-channel FCS hit statistics for the hot / dead channel mask (Sean
// 2026-11-24: the run divided into periods, a sample of events per period, a
// mask per period). Reads SimpleTrees (runMudst.C output: the A_N jobs' EM/JP
// trigger filter) and, for every channel of ECal north / south (det 0, 1) and
// HCal north / south (det 2, 3), counts:
//   n_any   hits with energy > 0 (the channel read out above its pedestal cut)
//   n_thr   hits above the reclustering tower threshold (ECal 0.4 GeV, HCal
//           0.1 GeV; what can enter clusters and points)
//   sumE    summed energy of the hits above that threshold, and sumE2
// plus the number of events read, and each channel's position (x, y). Output
// histograms are indexed by channel id (hitid) and add with hadd, so a period
// is the hadd of its files' outputs; channel_mask.py turns a period into a mask.
//
//   g++ -m64 -O2 channel_occupancy.cpp -o channel_occupancy `root-config --cflags --libs` -std=c++11
//   ./channel_occupancy <SimpleTree.root or list> <out.root>
#include <TFile.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TTree.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

static const int kNDet = 4;
static const int kNCh[kNDet] = {748, 748, 260, 260};       // ECal 22 x 34, HCal 13 x 20 per side
static const float kThr[kNDet] = {0.4f, 0.4f, 0.1f, 0.1f}; // reclustering tower thresholds, GeV
static const char *kDetName[kNDet] = {"ecal_n", "ecal_s", "hcal_n", "hcal_s"};

int main(int argc, char **argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: ./channel_occupancy <SimpleTree.root | file list> <out.root>\n"); return 1; }
    std::vector<std::string> files;
    const std::string in = argv[1];
    if (in.size() > 5 && in.compare(in.size() - 5, 5, ".root") == 0) files.push_back(in);
    else { std::ifstream f(in.c_str()); std::string l; while (std::getline(f, l)) if (!l.empty()) files.push_back(l); }

    TH1D *h_any[kNDet], *h_thr[kNDet], *h_e[kNDet], *h_e2[kNDet], *h_x[kNDet], *h_y[kNDet];
    for (int d = 0; d < kNDet; d++)
    {
        h_any[d] = new TH1D(Form("n_any_%s", kDetName[d]), ";channel id;hits with E > 0", kNCh[d], -0.5, kNCh[d] - 0.5);
        h_thr[d] = new TH1D(Form("n_thr_%s", kDetName[d]), Form(";channel id;hits with E > %.1f GeV", kThr[d]), kNCh[d], -0.5, kNCh[d] - 0.5);
        h_e[d] = new TH1D(Form("sumE_%s", kDetName[d]), ";channel id;summed E above threshold (GeV)", kNCh[d], -0.5, kNCh[d] - 0.5);
        h_e2[d] = new TH1D(Form("sumE2_%s", kDetName[d]), ";channel id;summed E^{2} above threshold", kNCh[d], -0.5, kNCh[d] - 0.5);
        // positions: filled once per channel with the channel's x / y (bin content), counted in h_npos
        h_x[d] = new TH1D(Form("posx_%s", kDetName[d]), ";channel id;x (cm) x n_pos", kNCh[d], -0.5, kNCh[d] - 0.5);
        h_y[d] = new TH1D(Form("posy_%s", kDetName[d]), ";channel id;y (cm) x n_pos", kNCh[d], -0.5, kNCh[d] - 0.5);
    }
    TH1D *h_npos[kNDet];
    for (int d = 0; d < kNDet; d++) h_npos[d] = new TH1D(Form("npos_%s", kDetName[d]), ";channel id;files with the position", kNCh[d], -0.5, kNCh[d] - 0.5);

    static int detid[5000], hitid[5000];
    static float energy[5000], posx[5000], posy[5000];
    int nhits = 0;
    double n_events = 0;
    for (const auto &fn : files)
    {
        TFile f(fn.c_str());
        TTree *t = f.IsOpen() ? (TTree *)f.Get("data") : nullptr;
        if (!t) { std::fprintf(stderr, "channel_occupancy: no tree in %s, skipped\n", fn.c_str()); continue; }
        t->SetBranchStatus("*", 0);
        for (const char *b : {"Cal_nhits", "Cal_detid", "Cal_hitid", "Cal_hit_energy", "Cal_hit_posx", "Cal_hit_posy"}) t->SetBranchStatus(b, 1);
        t->SetBranchAddress("Cal_nhits", &nhits);
        t->SetBranchAddress("Cal_detid", detid);
        t->SetBranchAddress("Cal_hitid", hitid);
        t->SetBranchAddress("Cal_hit_energy", energy);
        t->SetBranchAddress("Cal_hit_posx", posx);
        t->SetBranchAddress("Cal_hit_posy", posy);
        std::vector<char> seen[kNDet];
        for (int d = 0; d < kNDet; d++) seen[d].assign(kNCh[d], 0);
        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; i++)
        {
            t->GetEntry(i);
            for (int h = 0; h < nhits && h < 5000; h++)
            {
                const int d = detid[h], c = hitid[h];
                if (d < 0 || d >= kNDet || c < 0 || c >= kNCh[d]) continue;
                const float e = energy[h];
                if (!seen[d][c]) { h_x[d]->Fill(c, posx[h]); h_y[d]->Fill(c, posy[h]); h_npos[d]->Fill(c); seen[d][c] = 1; }
                if (e > 0) h_any[d]->Fill(c);
                if (e > kThr[d]) { h_thr[d]->Fill(c); h_e[d]->Fill(c, e); h_e2[d]->Fill(c, (double)e * e); }
            }
        }
        n_events += n;
    }
    TFile out(argv[2], "RECREATE");
    for (int d = 0; d < kNDet; d++) { h_any[d]->Write(); h_thr[d]->Write(); h_e[d]->Write(); h_e2[d]->Write(); h_x[d]->Write(); h_y[d]->Write(); h_npos[d]->Write(); }
    TH1D h_ev("n_events", ";;events read", 1, 0, 1);
    h_ev.SetBinContent(1, n_events);
    h_ev.Write();
    out.Close();
    std::printf("channel_occupancy: %zu file(s), %.0f events -> %s\n", files.size(), n_events, argv[2]);
    return 0;
}
