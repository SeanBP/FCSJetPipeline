// Diagram of the masked FCS channels per period (Sean 2026-11-24): ECal (left)
// and HCal (right), north and south, at their STAR x, y; every channel drawn as
// an outlined box, masked ones filled by reason (hot red, dead dark blue, hotE
// orange, coldE green if present). One PNG per period, <out>_p<k>.png, and all
// periods in <out>.pdf.
// Channel positions come from any channel_occupancy output (posx / posy / npos).
// Usage: root -b -q 'plot_channel_mask.C+("channel_mask_v1.txt","occ_23007037_0.root","mask_v1","mask_v1_periods.txt")'
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

struct MaskRow { int period, first, last, det, id; double x, y; std::string flag; };

void plot_channel_mask(const char *mask, const char *occfile, const char *out = "channel_mask", const char *periods = "")
{
    gStyle->SetOptStat(0);
    const int nch[4] = {748, 748, 260, 260};
    const char *dn[4] = {"ecal_n", "ecal_s", "hcal_n", "hcal_s"};
    const double pitch[4] = {5.572, 5.572, 9.99, 9.99};
    std::vector<double> x[4], y[4];
    TFile f(occfile);
    if (!f.IsOpen()) { printf("plot_channel_mask: cannot open %s\n", occfile); return; }
    for (int d = 0; d < 4; d++)
    {
        TH1D *hx = (TH1D *)f.Get(Form("posx_%s", dn[d])), *hy = (TH1D *)f.Get(Form("posy_%s", dn[d])), *hn = (TH1D *)f.Get(Form("npos_%s", dn[d]));
        x[d].assign(nch[d], 1e9); y[d].assign(nch[d], 1e9);
        for (int c = 0; c < nch[d]; c++) if (hn && hn->GetBinContent(c + 1) > 0) { x[d][c] = hx->GetBinContent(c + 1) / hn->GetBinContent(c + 1); y[d][c] = hy->GetBinContent(c + 1) / hn->GetBinContent(c + 1); }
    }
    std::vector<MaskRow> rows;
    std::map<int, std::pair<int, int> > prange;
    {
        std::ifstream in(mask);
        std::string l;
        while (std::getline(in, l))
        {
            if (l.empty() || l[0] == '#') continue;
            std::istringstream s(l);
            MaskRow r;
            if (s >> r.period >> r.first >> r.last >> r.det >> r.id >> r.x >> r.y >> r.flag) { rows.push_back(r); prange[r.period] = std::make_pair(r.first, r.last); }
        }
    }
    // periods (make_mask_sample.py <prefix>_periods.txt), so periods without a
    // masked channel are drawn too; without it, those in the mask file
    int pmax = -1;
    if (periods && *periods)
    {
        std::ifstream in(periods);
        std::string l;
        while (std::getline(in, l))
        {
            if (l.empty() || l[0] == '#') continue;
            std::istringstream s(l);
            int k, a, b;
            if (s >> k >> a >> b) { prange[k] = std::make_pair(a, b); if (k > pmax) pmax = k; }
        }
    }
    for (const auto &r : rows) if (r.period > pmax) pmax = r.period;
    std::map<std::string, int> col = {{"hot", kRed + 1}, {"dead", kBlue + 3}, {"hotE", kOrange + 1}, {"coldE", kGreen + 2}};
    TCanvas c("c", "", 1500, 720);
    c.Print(Form("%s.pdf[", out));
    for (int p = 0; p <= pmax; p++)
    {
        c.Clear();
        TPad *pads[2] = {new TPad(Form("pa%d", p), "", 0.0, 0.0, 0.5, 0.93), new TPad(Form("pb%d", p), "", 0.5, 0.0, 1.0, 0.93)};
        pads[0]->Draw(); pads[1]->Draw();
        int nmask[4] = {0, 0, 0, 0};
        std::map<std::string, int> nflag;
        for (int pad = 0; pad < 2; pad++)
        {
            pads[pad]->cd();
            gPad->SetLeftMargin(0.1); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.1);
            const int d0 = pad == 0 ? 0 : 2;
            const double xr = pad == 0 ? 150 : 170, yr = pad == 0 ? 110 : 120;
            TH2D *fr = new TH2D(Form("fr%d_%d", p, pad), ";STAR x (cm);STAR y (cm)", 10, -xr, xr, 10, -yr, yr + 35);  // room above for the legend
            fr->Draw();
            for (int d = d0; d < d0 + 2; d++)
            {
                const double h = 0.5 * pitch[d];
                for (int ch = 0; ch < nch[d]; ch++)
                {
                    if (x[d][ch] > 1e8) continue;
                    TBox *b = new TBox(x[d][ch] - h, y[d][ch] - h, x[d][ch] + h, y[d][ch] + h);
                    b->SetFillStyle(0); b->SetLineColor(kGray); b->Draw("l");
                }
                for (const auto &r : rows)
                {
                    if (r.period != p || r.det != d || r.id < 0 || r.id >= nch[d]) continue;
                    // the mask row's own position (a dead channel may have none in the occ file)
                    const double cx = r.x == r.x ? r.x : x[d][r.id], cy = r.y == r.y ? r.y : y[d][r.id];
                    if (cx > 1e8) continue;
                    TBox *b = new TBox(cx - h, cy - h, cx + h, cy + h);
                    b->SetFillColor(col.count(r.flag) ? col[r.flag] : kMagenta); b->SetFillStyle(1001); b->Draw();
                    nmask[d]++; nflag[r.flag]++;
                }
            }
            TLatex t; t.SetNDC(); t.SetTextSize(0.035);
            t.DrawLatex(0.12, 0.92, Form("%s (north x < 0, south x > 0): masked north %d, south %d", pad == 0 ? "ECal" : "HCal", nmask[d0], nmask[d0 + 1]));
            if (pad == 1)
            {
                TLegend *lg = new TLegend(0.12, 0.81, 0.96, 0.89);
                lg->SetNColumns(4); lg->SetBorderSize(1); lg->SetFillColor(kWhite); lg->SetTextSize(0.033);
                for (const char *fl : {"hot", "dead", "hotE", "coldE"})
                {
                    TBox *bx = new TBox(); bx->SetFillColor(col[fl]); bx->SetFillStyle(1001);
                    lg->AddEntry(bx, Form("%s (%d)", fl, nflag.count(fl) ? nflag[fl] : 0), "f");
                }
                lg->Draw();
            }
        }
        c.cd(0);
        TLatex tt; tt.SetNDC(); tt.SetTextSize(0.03); tt.SetTextAlign(21);
        const bool known = prange.count(p);
        tt.SetTextSize(0.032);
        tt.DrawLatex(0.5, 0.95, known ? Form("Masked FCS channels, period %d (runs %d-%d)", p, prange[p].first, prange[p].second)
                                       : Form("Masked FCS channels, period %d", p));
        c.SaveAs(Form("%s_p%02d.png", out, p));
        c.Print(Form("%s.pdf", out));
    }
    c.Print(Form("%s.pdf]", out));
    printf("plot_channel_mask: %zu masked channel-periods, periods 0-%d -> %s_p*.png, %s.pdf\n", rows.size(), pmax, out, out);
}
