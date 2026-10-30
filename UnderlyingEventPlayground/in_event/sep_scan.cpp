// How close to the jet can the eta-weighted annulus (whole ring, other jets
// excluded) start before jet energy leaks into it? From in_event_ue v5:
//   top     shell estimate / reference against the shell's dR from the jet
//   bottom  ring estimate with dR >= d_min / reference against d_min
// The reference is the same jet's ring with dR >= 2.5R (1.75); both must
// pass the validity cut (covered strips >= 70% of the jet template). Ratios
// of the weighted sums over jets; data reco, simulation reco and simulation
// truth photons on live cells. Also prints the tables.
//
// Build and run from this folder:
//   g++ -m64 -O2 sep_scan.cpp -o sep_scan `root-config --cflags --libs` -std=c++11
//   ./sep_scan "<data files>" "<sim files>" <out.png>

#include <TROOT.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

static const int kNB = 4, kNSh = 15, kNDm = 8;
static const double kE[kNB + 1] = {10, 20, 35, 60, 120};
static const double kShEdge[kNSh + 1] = {0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.3};
static const double kDmin[kNDm] = {0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.75};
static const int kRef = kNDm - 1;
static const double kMaxFill = 1 / 0.7, R = 0.7;

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}
// ratio of weighted sums with the error from per-jet covariance
struct Ratio
{
    double w, w2, sx, sy, sxx, syy, sxy;
    Ratio() : w(0), w2(0), sx(0), sy(0), sxx(0), syy(0), sxy(0) {}
    void Add(double x, double y, double ww) { w += ww; w2 += ww * ww; sx += ww * x; sy += ww * y; sxx += ww * x * x; syy += ww * y * y; sxy += ww * x * y; }
    double Val() const { return sy != 0 ? sx / sy : NAN; }
    double Err() const
    {
        if (!(w > 0) || sy == 0) return NAN;
        const double mx = sx / w, my = sy / w, r = mx / my;
        const double vx = sxx / w - mx * mx, vy = syy / w - my * my, cxy = sxy / w - mx * my;
        return std::fabs(r) * std::sqrt(std::max(0.0, vx / (mx * mx) + vy / (my * my) - 2 * cxy / (mx * my)) / (w * w / w2));
    }
};
// [source: 0 data reco, 1 sim reco, 2 sim truth][E bin][point]
struct Res { Ratio sh[3][kNB][kNSh], dm[3][kNB][kNDm]; };

static bool Ok(float f) { return f > 0 && f <= kMaxFill; }

static void Fill(const std::vector<std::string> &files, bool sim, Res &r)
{
    TChain ch("jets");
    for (size_t i = 0; i < files.size(); i++) ch.Add(files[i].c_str());
    Float_t E, w, sh_e[kNSh], sh_fill[kNSh], sh_tl[kNSh], dm_e[kNDm], dm_fill[kNDm], dm_tl[kNDm];
    ch.SetBranchAddress("jet_E", &E); ch.SetBranchAddress("weight", &w);
    ch.SetBranchAddress("sh_e", sh_e); ch.SetBranchAddress("sh_fill", sh_fill); ch.SetBranchAddress("sh_tl", sh_tl);
    ch.SetBranchAddress("dm_e", dm_e); ch.SetBranchAddress("dm_fill", dm_fill); ch.SetBranchAddress("dm_tl", dm_tl);
    for (Long64_t i = 0; i < ch.GetEntries(); i++)
    {
        ch.GetEntry(i);
        int b = -1;
        for (int q = 0; q < kNB; q++) if (E >= kE[q] && E < kE[q + 1]) b = q;
        if (b < 0 || !Ok(dm_fill[kRef])) continue;
        const double ww = sim ? w : 1.0;
        for (int q = 0; q < kNSh; q++)
            if (Ok(sh_fill[q]))
            {
                r.sh[sim ? 1 : 0][b][q].Add(sh_e[q], dm_e[kRef], ww);
                if (sim) r.sh[2][b][q].Add(sh_tl[q], dm_tl[kRef], ww);
            }
        for (int q = 0; q < kNDm; q++)
            if (Ok(dm_fill[q]))
            {
                r.dm[sim ? 1 : 0][b][q].Add(dm_e[q], dm_e[kRef], ww);
                if (sim) r.dm[2][b][q].Add(dm_tl[q], dm_tl[kRef], ww);
            }
    }
    std::printf("%s: %lld jets\n", sim ? "simulation" : "data", ch.GetEntries());
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 4) { std::fprintf(stderr, "usage: ./sep_scan \"<data files>\" \"<sim files>\" <out.png>\n"); return 1; }
    Res r;
    Fill(Split(argv[1]), false, r);
    Fill(Split(argv[2]), true, r);

    const char *src[3] = {"data reco", "sim reco", "sim truth"};
    for (int b = 0; b < kNB; b++)
    {
        std::printf("\njet E %.0f-%.0f GeV: shell / ring(dR >= 1.75)      data reco | sim reco | sim truth\n", kE[b], kE[b + 1]);
        for (int q = 0; q < kNSh; q++)
            std::printf("  shell %.1f-%.1f   %.3f+-%.3f | %.3f+-%.3f | %.3f+-%.3f\n", kShEdge[q], kShEdge[q + 1],
                        r.sh[0][b][q].Val(), r.sh[0][b][q].Err(), r.sh[1][b][q].Val(), r.sh[1][b][q].Err(), r.sh[2][b][q].Val(), r.sh[2][b][q].Err());
        for (int q = 0; q < kNDm; q++)
            std::printf("  dR >= %.2f      %.3f+-%.3f | %.3f+-%.3f | %.3f+-%.3f\n", kDmin[q],
                        r.dm[0][b][q].Val(), r.dm[0][b][q].Err(), r.dm[1][b][q].Val(), r.dm[1][b][q].Err(), r.dm[2][b][q].Val(), r.dm[2][b][q].Err());
    }

    gStyle->SetOptStat(0);
    const int col[3] = {TColor::GetColor("#222222"), TColor::GetColor("#d9480f"), TColor::GetColor("#1c6bb0")};
    const int mk[3] = {20, 21, 22};
    TCanvas cv("cv", "", 1900, 950);
    cv.Divide(kNB, 2, 0.002, 0.002);
    std::vector<TObject *> keep;
    for (int row = 0; row < 2; row++)
        for (int b = 0; b < kNB; b++)
        {
            cv.cd(row * kNB + b + 1);
            gPad->SetMargin(0.14, 0.03, 0.13, 0.09);
            TH1D *fr = new TH1D(Form("fr%d%d", row, b), "", 1, 0.6, row == 0 ? 3.0 : 1.85);
            fr->SetTitle(row == 0 ? ";shell #DeltaR from the jet axis;shell estimate / ring (#DeltaR #geq 2.5R)"
                                  : ";inner edge d_{min} of the ring;ring (#DeltaR #geq d_{min}) / ring (#DeltaR #geq 2.5R)");
            fr->SetMinimum(row == 0 ? 0.6 : 0.85); fr->SetMaximum(row == 0 ? 2.6 : 1.45);
            fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045);
            fr->GetYaxis()->SetTitleOffset(1.3);
            fr->Draw(); keep.push_back(fr);
            const double xhi = row == 0 ? 3.0 : 1.85;
            TBox *band = new TBox(0.6, 0.97, xhi, 1.03); band->SetFillColor(TColor::GetColor("#eeeeee")); band->Draw(); keep.push_back(band);
            TLine *l1 = new TLine(0.6, 1, xhi, 1); l1->SetLineStyle(2); l1->Draw(); keep.push_back(l1);
            for (int k = 0; k < 2; k++)
            {
                const double x = k == 0 ? 1.0 * R : 2 * R;
                TLine *l = new TLine(x, fr->GetMinimum(), x, fr->GetMaximum()); l->SetLineColor(TColor::GetColor("#868e96")); l->SetLineStyle(7); l->Draw(); keep.push_back(l);
            }
            for (int s = 0; s < 3; s++)
            {
                TGraphErrors *g = new TGraphErrors();
                const int n = row == 0 ? kNSh : kNDm;
                for (int q = 0; q < n; q++)
                {
                    const Ratio &rr = row == 0 ? r.sh[s][b][q] : r.dm[s][b][q];
                    if (!(rr.w > 0)) continue;
                    const double x = (row == 0 ? 0.5 * (kShEdge[q] + kShEdge[q + 1]) : kDmin[q]) + 0.012 * (s - 1);
                    const int p = g->GetN();
                    g->SetPoint(p, x, rr.Val()); g->SetPointError(p, 0, rr.Err());
                }
                g->SetMarkerColor(col[s]); g->SetLineColor(col[s]); g->SetMarkerStyle(mk[s]); g->SetMarkerSize(1.2);
                g->Draw("PL SAME"); keep.push_back(g);
                if (row == 0 && b == 0)
                {
                    static TLegend lg(0.45, 0.62, 0.96, 0.86);
                    lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.045);
                    lg.AddEntry(g, src[s], "pl");
                    if (s == 2) lg.Draw();
                }
            }
            gPad->RedrawAxis();
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.05);
            t.DrawLatex(0.14, 0.93, Form("jets %.0f-%.0f GeV", kE[b], kE[b + 1]));
            t.SetTextSize(0.038); t.SetTextColor(TColor::GetColor("#868e96"));
            t.DrawLatex(0.15, 0.17, "dashed: #DeltaR = R, 2R; band: #pm3%");
        }
    cv.SaveAs(argv[3]);
    return 0;
}
