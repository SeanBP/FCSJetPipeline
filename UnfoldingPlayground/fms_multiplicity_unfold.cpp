// Closure test of the STAR FMS EM-jet analysis's photon-multiplicity
// correction of A_N, on our FCS simulation (BrightSTAR
// emJetAnalysis/src/EjUnfolding.C, EjRooUnfolding.C and
// emJetAnalysis/scripts/unfolding/file/anaUnfolding.C).
//
// FMS measured A_N in bins of detector-level photon multiplicity and used a
// simulated migration matrix (detector vs particle multiplicity, per energy
// and pT bin) to correct it. Its code holds three approaches:
//   1. TSVDUnfold (ROOT; RooUnfoldSvd in EjRooUnfolding.C) with the
//      regularisation kreg equal to the number of bins (5), i.e. effectively
//      no regularisation: a matrix inversion
//   2. a weighted sum, A(n_part) = sum over n_det of A_det(n_det) P(n_det | n_part)
//      (anaUnfolding.C ph2, "Making weighted adding")
//   3. a direct linear solve of the matrix equation (anaUnfolding.C linSol)
// Here (1) is applied to the spin-sorted yields, the correct quantity to
// unfold, and (2) as FMS wrote it. (3) is (1) without regularisation and
// without the toy-based errors, so it is not repeated.
//
// Test: the sample is split by file, odd files for the response, even files
// as pseudo-data. Pseudo-data jets get an injected asymmetry that depends
// only on the particle-level multiplicity, A(n) = A0 / n. Without statistical
// fluctuations ("Asimov" yields: each jet adds w (1 +/- A cos phi) / 2 to the
// up/down yields), the transverse-spin moments are
//   D(n) = sum w A cos^2(phi),   S(n) = sum w cos^2(phi),   A(n) = D(n) / S(n)
// and D and S are each unfolded with the same response, so the unfolded
// asymmetry is D_unf / S_unf. Any difference from A0 / n is then the bias of
// the method, not noise. Matched pairs only (as FMS): no misses or fakes.
//
// Usage: ./fms_multiplicity_unfold <file list> [output_base] [A0=0.10] [pipeline|fms matching]
// Build: g++ -m64 -O2 fms_multiplicity_unfold.cpp ../shared/VectorDict.cxx -o fms_multiplicity_unfold
//            -I../shared -I../jet_scale `root-config --cflags --libs` -std=c++11

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSVDUnfold.h>

#include <algorithm>
#include "UnfoldingInput.h"

const int n_nph = 5;
inline int NphBin(int n) { return std::min(std::max(n, 1), n_nph); }

// detector-level UE-subtracted energy bins, as in fms_response.cpp, plus "all"
const int n_ebin = 6;
const double e_edges[n_ebin + 1] = {0, 20, 40, 60, 80, 100, 150};

struct Moments
{
    TH1D *D, *S;
    Moments(const char *name)
    {
        D = new TH1D(Form("%s_D", name), "", n_nph, 0.5, n_nph + 0.5);
        S = new TH1D(Form("%s_S", name), "", n_nph, 0.5, n_nph + 0.5);
        D->Sumw2(); S->Sumw2();
    }
    double A(int b) const { return S->GetBinContent(b) > 0 ? D->GetBinContent(b) / S->GetBinContent(b) : 0.0; }
};

TH2D *DiagCov(const TH1D *h, const char *name)
{
    TH2D *c = new TH2D(name, "", n_nph, 0.5, n_nph + 0.5, n_nph, 0.5, n_nph + 0.5);
    for (int b = 1; b <= n_nph; b++) c->SetBinContent(b, b, std::max(h->GetBinError(b) * h->GetBinError(b), 1e-12));
    return c;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    if (argc < 2)
    {
        std::cerr << "usage: ./fms_multiplicity_unfold <file list> [output_base] [A0] [pipeline|fms]" << std::endl;
        return 1;
    }
    const std::string out_base = argc > 2 ? argv[2] : "fms_multiplicity_unfold";
    const double A0 = argc > 3 ? std::atof(argv[3]) : 0.10;
    const MatchMode mode = (argc > 4 && std::string(argv[4]) == "fms") ? kMatchFms : kMatchPipeline;

    std::vector<JetPair> pairs;
    if (!LoadJetPairs(ReadFileList(argv[1]), mode, pairs)) return 1;

    const int n_sets = n_ebin + 1;   // the energy bins, then all energies
    std::vector<TH2D *> resp(n_sets);
    std::vector<TH1D *> resp_det(n_sets), resp_part(n_sets);
    std::vector<Moments *> det(n_sets), part(n_sets);
    for (int s = 0; s < n_sets; s++)
    {
        resp[s] = new TH2D(Form("resp_%d", s), ";detector multiplicity;particle multiplicity", n_nph, 0.5, n_nph + 0.5, n_nph, 0.5, n_nph + 0.5);
        resp[s]->Sumw2();
        det[s] = new Moments(Form("det_%d", s));
        part[s] = new Moments(Form("part_%d", s));
    }

    long n_resp = 0, n_data = 0;
    for (size_t i = 0; i < pairs.size(); i++)
    {
        const JetPair &p = pairs[i];
        if (!p.has_ue || p.em_level < 0 || p.rptsub() < 2.0) continue;   // analysis selection, as for data
        const int ie = (int)(std::upper_bound(e_edges, e_edges + n_ebin + 1, p.rEsub()) - e_edges) - 1;
        const int bd = NphBin(p.rnph), bt = NphBin(p.tnph);
        const int sets[2] = {ie >= 0 && ie < n_ebin ? ie : -1, n_ebin};
        for (int k = 0; k < 2; k++)
        {
            const int s = sets[k];
            if (s < 0) continue;
            if (p.file % 2 == 1)
            {
                resp[s]->Fill(bd, bt, p.w);
                if (k == 1) n_resp++;
            }
            else
            {
                const double c2 = std::cos(p.rphi) * std::cos(p.rphi);
                const double a = A0 / std::max(p.tnph, 1);   // injected: depends only on the particle multiplicity (0 counts as 1)
                det[s]->D->Fill(bd, p.w * a * c2); det[s]->S->Fill(bd, p.w * c2);
                part[s]->D->Fill(bt, p.w * a * c2); part[s]->S->Fill(bt, p.w * c2);
                if (k == 1) n_data++;
            }
        }
    }
    std::cout << n_resp << " response pairs (odd files), " << n_data << " pseudo-data pairs (even files); injected A(n) = " << A0
              << " / n_particle" << std::endl;

    const char *set_name[n_sets] = {"E 0-20", "E 20-40", "E 40-60", "E 60-80", "E 80-100", "E 100-150", "all E"};
    TCanvas c("c", "", 1600, 900);
    c.Divide(4, 2, 0.001, 0.001);
    std::vector<TObject *> keep;
    TFile fout((out_base + ".root").c_str(), "RECREATE");

    for (int s = 0; s < n_sets; s++)
    {
        if (resp[s]->Integral() <= 0 || det[s]->S->Integral() <= 0) continue;
        TH1D *bini = resp[s]->ProjectionX(Form("bini_%d", s));
        TH1D *xini = resp[s]->ProjectionY(Form("xini_%d", s));

        // (1) TSVDUnfold of D and S, kreg = number of bins (FMS) and 2 (regularised)
        double A_svd[2][n_nph + 1];
        const int kregs[2] = {n_nph, 2};
        for (int r = 0; r < 2; r++)
        {
            TH1D *unf[2];
            for (int m = 0; m < 2; m++)
            {
                TH1D *meas = m == 0 ? det[s]->D : det[s]->S;
                TH2D *cov = DiagCov(meas, Form("cov_%d_%d_%d", s, r, m));
                TSVDUnfold svd(meas, cov, bini, xini, resp[s]);
                unf[m] = svd.Unfold(kregs[r]);
                unf[m]->SetName(Form("unf_%d_k%d_%s", s, kregs[r], m == 0 ? "D" : "S"));
            }
            for (int b = 1; b <= n_nph; b++)
                A_svd[r][b] = unf[1]->GetBinContent(b) != 0 ? unf[0]->GetBinContent(b) / unf[1]->GetBinContent(b) : 0.0;
        }

        // (2) FMS weighted sum: A(n_part) = sum_d A_det(d) P(d | n_part)
        double A_wsum[n_nph + 1];
        for (int t = 1; t <= n_nph; t++)
        {
            double row = 0, v = 0;
            for (int d = 1; d <= n_nph; d++) row += resp[s]->GetBinContent(d, t);
            for (int d = 1; d <= n_nph; d++) v += det[s]->A(d) * (row > 0 ? resp[s]->GetBinContent(d, t) / row : 0);
            A_wsum[t] = v;
        }

        std::cout << "\n== " << set_name[s] << " GeV ==\n  n     true    detector  SVD k=" << n_nph << "  SVD k=2   weighted sum" << std::endl;
        for (int b = 1; b <= n_nph; b++)
            std::printf("  %s  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n", b < n_nph ? Form("%d ", b) : ">=5", part[s]->A(b), det[s]->A(b),
                        A_svd[0][b], A_svd[1][b], A_wsum[b]);

        c.cd(s + 1);
        gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08);
        TH1F *frame = gPad->DrawFrame(0.5, 0.0, n_nph + 0.5, 1.25 * A0);
        frame->GetXaxis()->SetTitle("photon multiplicity (5 = 5 or more)");
        frame->GetYaxis()->SetTitle("A_{N}");
        frame->GetXaxis()->SetNdivisions(5);
        TGraph *g[5];
        const char *gname[5] = {"injected (particle level)", "detector level", Form("TSVDUnfold k=%d (FMS)", n_nph), "TSVDUnfold k=2",
                                "FMS weighted sum"};
        int col[5] = {kBlack, kGray + 1, kRed + 1, kAzure + 2, kGreen + 2};
        int mk[5] = {20, 24, 21, 25, 22};
        for (int q = 0; q < 5; q++)
        {
            g[q] = new TGraph();
            for (int b = 1; b <= n_nph; b++)
            {
                const double v = q == 0 ? part[s]->A(b) : q == 1 ? det[s]->A(b) : q == 2 ? A_svd[0][b] : q == 3 ? A_svd[1][b] : A_wsum[b];
                g[q]->SetPoint(b - 1, b + 0.06 * (q - 2), v);
            }
            g[q]->SetMarkerStyle(mk[q]); g[q]->SetMarkerColor(col[q]); g[q]->SetLineColor(col[q]); g[q]->SetMarkerSize(1.1);
            g[q]->Draw(q == 0 ? "PL same" : "P same");
            keep.push_back(g[q]);
        }
        TLatex *tl = new TLatex();
        tl->SetNDC(); tl->SetTextSize(0.055); tl->SetTextFont(62);
        tl->DrawLatex(0.14, 0.935, Form("%s GeV", set_name[s]));
        keep.push_back(tl);
        if (s == 0)
        {
            c.cd(8);
            TLegend *leg = new TLegend(0.05, 0.3, 0.95, 0.85);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.06);
            for (int q = 0; q < 5; q++) leg->AddEntry(g[q], gname[q], "p");
            leg->Draw();
            keep.push_back(leg);
        }
        fout.cd();
        resp[s]->Write();
    }
    c.SaveAs((out_base + ".png").c_str());
    c.SaveAs((out_base + ".pdf").c_str());
    fout.Close();
    return 0;
}
