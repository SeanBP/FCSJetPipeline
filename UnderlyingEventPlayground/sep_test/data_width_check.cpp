// The separation test on FCS data, against the band half-width each jet's
// position gives: does the no-separation estimate grow with the window
// width, as the simulation's truth decomposition (truth_decomp/) finds?
//
// Matches the jets of the three variants (run_sep_test.sh output) jet by
// jet, as sep_test_compare.cpp, and bins them by the no-separation
// half-width (the band of the 5% area with no exclusion) and, separately,
// by the pipeline half-width (the 5% area outside 2.5R), both from the jet
// position on a 2 cm grid. Prints per bin the mean reco cone energy
// (reco_ue_cone_E, the mean proxy cone, not area scaled) of each variant and
// the ratios nosep / sep and sep_veto / sep.
//
// Build and run from this folder:
//   g++ -m64 -O2 data_width_check.cpp ../../shared/VectorDict.cxx -o data_width_check -I../../shared -I../../underlying_event \
//       `root-config --cflags --libs` `fastjet-config --cxxflags` -std=c++11
//   ./data_width_check <outdir of submit_sep_test.sh>

#include <TFile.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "VectorDict.h"
#include "ProxySelection.h"

static double HalfWidthNoExclusion(double ue_phi, double target)
{
    const std::vector<ProxyCell> &cells = ProxyFiducialCells();
    std::vector<std::pair<double, double> > dev_area;
    for (size_t i = 0; i < cells.size(); i++)
    {
        const double dev = ProxyOffAxisDeviationDeg(ue_phi, cells[i].phi);
        if (dev <= proxy_half_width_max_deg) dev_area.push_back(std::make_pair(dev, cells[i].area));
    }
    std::sort(dev_area.begin(), dev_area.end());
    double area = 0.0;
    for (size_t i = 0; i < dev_area.size(); i++)
        if ((area += dev_area[i].second) >= target) return dev_area[i].first;
    return proxy_half_width_max_deg;
}

// Both half-widths on a grid over one half of the fiducial region
struct Tables
{
    double x_lo, y_lo, step; int nx, ny;
    std::vector<float> hw_nosep, hw_sep;
    void Build(double grid)
    {
        const FiducialRect &r = kFiducialRectEcal;
        const double b = reco_fiducial_buffer_ecal, target = ProxyTargetArea();
        step = grid; x_lo = r.x_inner + b; y_lo = r.y_min + b;
        nx = (int)std::ceil((r.x_outer - b - x_lo) / step); ny = (int)std::ceil((r.y_max - b - y_lo) / step);
        hw_nosep.assign(nx * ny, 0.f); hw_sep.assign(nx * ny, 0.f);
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
            {
                double eta, phi;
                computeEtaPhi(x_lo + (i + 0.5) * step, y_lo + (j + 0.5) * step, eta, phi);
                hw_nosep[i * ny + j] = HalfWidthNoExclusion(phi, target);
                hw_sep[i * ny + j] = ProxyHalfWidthDeg(eta, phi);
            }
    }
    int Index(double x, double y) const
    {
        const int i = std::min(nx - 1, std::max(0, (int)((std::fabs(x) - x_lo) / step)));
        const int j = std::min(ny - 1, std::max(0, (int)((y - y_lo) / step)));
        return i * ny + j;
    }
};

int main(int argc, char **argv)
{
    if (argc < 2) { printf("usage: ./data_width_check <outdir>\n"); return 1; }
    const std::string dir = argv[1];
    const int NV = 3;
    const char *vdir[NV] = {"sep_veto", "sep", "nosep"};
    Tables tab;
    tab.Build(2.0);

    std::vector<std::string> names;
    {
        TSystemDirectory d("d", (dir + "/sep_veto").c_str());
        TList *l = d.GetListOfFiles();
        TIter next(l);
        TSystemFile *f;
        while ((f = (TSystemFile *)next()))
        {
            const std::string n = f->GetName();
            if (n.find(".root") != std::string::npos) names.push_back(n);
        }
        std::sort(names.begin(), names.end());
    }
    const double edges[] = {0, 10, 15, 20, 25, 30, 35, 40, 46};
    const int nb = 8;
    double s[2][nb][NV], n[2][nb];
    for (int k = 0; k < 2; k++) for (int b = 0; b < nb; b++) { n[k][b] = 0; for (int v = 0; v < NV; v++) s[k][b][v] = 0; }
    long n_jets = 0;
    for (size_t i = 0; i < names.size(); i++)
    {
        TFile *tf[NV]; TTree *t[NV];
        bool ok = true;
        for (int v = 0; v < NV; v++)
        {
            tf[v] = TFile::Open((dir + "/" + vdir[v] + "/" + names[i]).c_str());
            t[v] = tf[v] ? (TTree *)tf[v]->Get("jetTree") : 0;
            ok = ok && t[v] && t[v]->GetEntries() == (t[0] ? t[0]->GetEntries() : -1);
        }
        if (!ok) { for (int v = 0; v < NV; v++) delete tf[v]; continue; }
        std::vector<float> *x = 0, *y = 0, *ue[NV] = {0, 0, 0};
        std::vector<int> *np[NV] = {0, 0, 0};
        t[0]->SetBranchAddress("reco_x", &x); t[0]->SetBranchAddress("reco_y", &y);
        for (int v = 0; v < NV; v++) { t[v]->SetBranchAddress("reco_ue_cone_E", &ue[v]); t[v]->SetBranchAddress("reco_ue_nproxy", &np[v]); }
        for (Long64_t e = 0; e < t[0]->GetEntries(); e++)
        {
            for (int v = 0; v < NV; v++) t[v]->GetEntry(e);
            for (size_t j = 0; j < x->size(); j++)
            {
                bool good = true;
                for (int v = 0; v < NV; v++) good = good && j < ue[v]->size() && np[v]->at(j) > 0 && ue[v]->at(j) >= 0;
                if (!good) continue;
                n_jets++;
                const int idx = tab.Index(x->at(j), y->at(j));
                const double w[2] = {tab.hw_nosep[idx], tab.hw_sep[idx]};
                for (int k = 0; k < 2; k++)
                    for (int b = 0; b < nb; b++)
                        if (w[k] >= edges[b] && w[k] < edges[b + 1])
                        {
                            n[k][b]++;
                            for (int v = 0; v < NV; v++) s[k][b][v] += ue[v]->at(j);
                        }
            }
        }
        for (int v = 0; v < NV; v++) delete tf[v];
    }
    printf("data_width_check: %zu files, %ld jets with all three estimates\n", names.size(), n_jets);
    const char *kt[2] = {"no-separation half-width (deg)", "pipeline half-width (deg)"};
    for (int k = 0; k < 2; k++)
    {
        printf("  binned by the %s: mean reco cone energy (GeV)\n    bin       jets   sep_veto   sep     nosep    nosep/sep  sep_veto/sep\n", kt[k]);
        for (int b = 0; b < nb; b++)
        {
            if (n[k][b] < 50) continue;
            const double m[NV] = {s[k][b][0] / n[k][b], s[k][b][1] / n[k][b], s[k][b][2] / n[k][b]};
            printf("    %2.0f-%-3.0f %7.0f   %6.3f   %6.3f   %6.3f     %5.3f      %5.3f\n", edges[b], edges[b + 1], n[k][b], m[0], m[1], m[2], m[2] / m[1], m[0] / m[1]);
        }
    }
    return 0;
}
