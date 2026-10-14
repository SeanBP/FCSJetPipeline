// Is the scatter of the FMS-style (TSVDUnfold, kreg = number of bins)
// multiplicity unfolding in fms_multiplicity_unfold.cpp statistical?
//
// Same closure test (odd files = response, even files = Asimov pseudo-data
// with A(n) = A0 / n_particle), repeated on disjoint sub-samples holding
// 1/nsub of the files, for nsub = 1, 2, 4, 8. Both the response and the
// pseudo-data come from the same physics, so with infinite statistics the
// matrix inversion closes exactly; any non-closure is sample fluctuations
// amplified by the inversion. If so, the RMS relative deviation from the
// injected asymmetry grows like sqrt(nsub), and its size at nsub = 1 tells
// how much more simulation is needed. Also prints the Kish effective number
// of jets, (sum w)^2 / sum w^2, since the cross-section weights shrink it.
//
// Usage: ./fms_unfold_stat_scan <file list> [A0=0.10]
// Build: g++ -m64 -O2 fms_unfold_stat_scan.cpp ../shared/VectorDict.cxx -o fms_unfold_stat_scan
//            -I../shared -I../jet_scale `root-config --cflags --libs` -std=c++11

#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSVDUnfold.h>
#include <TError.h>

#include <algorithm>
#include <cmath>
#include "UnfoldingInput.h"

const int n_nph = 5;
inline int NphBin(int n) { return std::min(std::max(n, 1), n_nph); }
const int n_ebin = 6;
const double e_edges[n_ebin + 1] = {0, 20, 40, 60, 80, 100, 150};
const int n_sets = n_ebin + 1;
const char *set_name[n_sets] = {"E 0-20", "E 20-40", "E 40-60", "E 60-80", "E 80-100", "E 100-150", "all E"};

struct Result { double rms[2]; double neff_resp, neff_data; bool ok; };

// One closure test on the pairs with (file / 2) % nsub == j
Result RunOne(const std::vector<JetPair> &pairs, int s_sel, int nsub, int j, double A0)
{
    static int uid = 0;
    uid++;
    // heap, never freed: TSVDUnfold and the projections keep pointers to these
    TH2D &resp = *new TH2D(Form("r%d", uid), "", n_nph, 0.5, n_nph + 0.5, n_nph, 0.5, n_nph + 0.5);
    TH1D &D = *new TH1D(Form("D%d", uid), "", n_nph, 0.5, n_nph + 0.5), &S = *new TH1D(Form("S%d", uid), "", n_nph, 0.5, n_nph + 0.5);
    TH1D &Dt = *new TH1D(Form("Dt%d", uid), "", n_nph, 0.5, n_nph + 0.5), &St = *new TH1D(Form("St%d", uid), "", n_nph, 0.5, n_nph + 0.5);
    resp.Sumw2(); D.Sumw2(); S.Sumw2();
    double sw[2] = {0, 0}, sw2[2] = {0, 0};
    for (size_t i = 0; i < pairs.size(); i++)
    {
        const JetPair &p = pairs[i];
        if (!p.has_ue || p.em_level < 0 || p.rptsub() < 2.0) continue;
        if ((p.file / 2) % nsub != j) continue;
        const int ie = (int)(std::upper_bound(e_edges, e_edges + n_ebin + 1, p.rEsub()) - e_edges) - 1;
        if (s_sel < n_ebin && ie != s_sel) continue;
        if (s_sel == n_ebin && (ie < 0 || ie >= n_ebin)) continue;
        const int bd = NphBin(p.rnph), bt = NphBin(p.tnph);
        const int h = p.file % 2;
        sw[h] += p.w; sw2[h] += p.w * p.w;
        if (h == 1) resp.Fill(bd, bt, p.w);
        else
        {
            const double c2 = std::cos(p.rphi) * std::cos(p.rphi), a = A0 / p.tnph;
            D.Fill(bd, p.w * a * c2); S.Fill(bd, p.w * c2);
            Dt.Fill(bt, p.w * a * c2); St.Fill(bt, p.w * c2);
        }
    }
    Result r;
    r.neff_resp = sw2[1] > 0 ? sw[1] * sw[1] / sw2[1] : 0;
    r.neff_data = sw2[0] > 0 ? sw[0] * sw[0] / sw2[0] : 0;
    r.ok = resp.Integral() > 0 && S.Integral() > 0;
    if (!r.ok) return r;
    TH1D *bini = resp.ProjectionX(Form("bi%d", uid));
    TH1D *xini = resp.ProjectionY(Form("xi%d", uid));
    const int kregs[2] = {n_nph, 2};
    for (int k = 0; k < 2; k++)
    {
        TH1D *unf[2];
        for (int m = 0; m < 2; m++)
        {
            TH1D *meas = m == 0 ? &D : &S;
            TH2D *cov = new TH2D(Form("c%d_%d_%d", uid, k, m), "", n_nph, 0.5, n_nph + 0.5, n_nph, 0.5, n_nph + 0.5);
            for (int b = 1; b <= n_nph; b++) cov->SetBinContent(b, b, std::max(meas->GetBinError(b) * meas->GetBinError(b), 1e-12));
            TSVDUnfold *svd = new TSVDUnfold(meas, cov, bini, xini, &resp);
            unf[m] = svd->Unfold(kregs[k]);
        }
        double sum = 0;
        for (int b = 1; b <= n_nph; b++)
        {
            const double at = St.GetBinContent(b) > 0 ? Dt.GetBinContent(b) / St.GetBinContent(b) : 0;
            const double au = unf[1]->GetBinContent(b) != 0 ? unf[0]->GetBinContent(b) / unf[1]->GetBinContent(b) : 0;
            const double d = at > 0 ? (au - at) / at : 0;
            sum += d * d;
        }
        r.rms[k] = std::sqrt(sum / n_nph);
    }
    return r;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning + 1;
    if (argc < 2) { std::cerr << "usage: ./fms_unfold_stat_scan <file list> [A0]" << std::endl; return 1; }
    const double A0 = argc > 2 ? std::atof(argv[2]) : 0.10;
    std::vector<JetPair> pairs;
    if (!LoadJetPairs(ReadFileList(argv[1]), kMatchPipeline, pairs)) return 1;

    const int nsubs[4] = {1, 2, 4, 8};
    std::printf("RMS relative deviation of unfolded A_N from injected, over the 5 multiplicity bins,\n"
                "averaged over the nsub disjoint sub-samples (ratio column: rms(nsub)/rms(1), sqrt(nsub) if statistical)\n");
    for (int s = 0; s < n_sets; s++)
    {
        std::printf("\n== %s GeV ==\n  nsub  Neff(resp)  Neff(data)   k=5 rms   ratio  sqrt(nsub)   k=2 rms\n", set_name[s]);
        double base = 0;
        for (int q = 0; q < 4; q++)
        {
            double m5 = 0, m2 = 0, nr = 0, nd = 0;
            int n = 0;
            for (int j = 0; j < nsubs[q]; j++)
            {
                Result r = RunOne(pairs, s, nsubs[q], j, A0);
                if (!r.ok) continue;
                m5 += r.rms[0]; m2 += r.rms[1]; nr += r.neff_resp; nd += r.neff_data; n++;
            }
            if (n == 0) continue;
            m5 /= n; m2 /= n; nr /= n; nd /= n;
            if (q == 0) base = m5;
            std::printf("  %4d  %10.0f  %10.0f  %8.3f  %6.2f  %10.2f  %8.3f\n", nsubs[q], nr, nd, m5, base > 0 ? m5 / base : 0,
                        std::sqrt((double)nsubs[q]), m2);
        }
    }
    return 0;
}
