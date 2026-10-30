// Jet leakage, data against simulation (reco only): the point energy per live
// eta-phi area in a thin ring just outside the jet, dR = R-1.1R from the jet
// axis (all azimuth), against jet pT and jet E; with a far ring (2.0-2.5R) as
// the UE baseline. Means are ratios of sums over jets (sum E / sum live area;
// simulation weighted per ptHat). From in_event_ue v8+ outputs (lr_E, lr_A;
// "-x": lrx_*, cells within R of other reco jets left out). "-n": each ring's
// density over the same sample's far-ring density in that bin, so the UE level
// (lower in simulation) drops out and the leakage is the rise above 1; drawn
// for the 1.0-1.1R and 1.2-1.4R rings, with sim / data underneath. "-d": the
// excess density, 1.0-1.1R minus 2.0-2.5R of the same sample (errors of the
// two added in quadrature, conservative: they are positively correlated),
// one row, data and simulation. "-e": the same excess split by the highest
// FCS EM trigger level of the jet's event (EM0 only, EM1, EM2; data fired
// bits, simulation emulated bits), against jet pT and E in coarser bins.
// "-u": as -e, but both samples classified by the same re-emulated level
// (run 23007048 settings; in_event_ue v9+ "emu" branch); "-u0": the same with
// data restricted to events recorded through fcsEM0 (prescaled, so an unbiased
// sample of EM0-condition events). With "-d": "-s2 <files>" adds a second
// simulation sample to the plot, "-l1 <label>", "-l2 <label>" name the two.
//
// Build and run from this folder:
//   g++ -m64 -O2 leak_ring.cpp -o leak_ring `root-config --cflags --libs` -std=c++11
//   ./leak_ring "<data files>" "<sim files>" <out.png> [-x] [-n | -d | -e | -u | -u0]

#include <TROOT.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

static const int kNL = 6, kNear = 0, kMid = 2, kFar = 5;
static const double kEdge[kNL + 1] = {1.0, 1.1, 1.2, 1.4, 1.6, 2.0, 2.5};

static std::vector<std::string> Split(const char *s)
{
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (ss >> x) v.push_back(x);
    return v;
}
// sum E / sum A with the error from the per-jet covariance
struct Dens
{
    double w = 0, w2 = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    void Add(double x, double y, double ww) { w += ww; w2 += ww * ww; sx += ww * x; sy += ww * y; sxx += ww * x * x; syy += ww * y * y; sxy += ww * x * y; }
    double Val() const { return sy > 0 ? sx / sy : NAN; }
    double Err() const
    {
        if (!(w > 0) || !(sy > 0) || !(sx > 0)) return NAN;
        const double mx = sx / w, my = sy / w, r = mx / my;
        const double vx = sxx / w - mx * mx, vy = syy / w - my * my, c = sxy / w - mx * my;
        return std::fabs(r) * std::sqrt(std::max(0.0, vx / (mx * mx) + vy / (my * my) - 2 * c / (mx * my)) / (w * w / w2));
    }
};

// a ring's density, optionally over (norm) or minus (diff) the far ring's (errors in quadrature)
static void Value(const Dens &a, const Dens &f, bool norm, double &y, double &ey, bool diff = false)
{
    y = a.Val(); ey = a.Err();
    if (diff) { y -= f.Val(); ey = std::sqrt(ey * ey + f.Err() * f.Err()); return; }
    if (!norm) return;
    const double r = y / f.Val();
    ey = r * std::sqrt(std::pow(ey / y, 2) + std::pow(f.Err() / f.Val(), 2));
    y = r;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 4) { std::fprintf(stderr, "usage: ./leak_ring \"<data files>\" \"<sim files>\" <out.png> [-x] [-n | -d | -e | -u | -u0]\n"); return 1; }
    bool excl = false, norm = false, diff = false, bylevel = false, useemu = false, hw0only = false;
    std::string sim2, lab[3] = {"data", "simulation (reco)", "simulation 2 (reco)"};
    for (int a = 4; a < argc; a++)
    {
        if (!std::strcmp(argv[a], "-x")) excl = true;
        if (!std::strcmp(argv[a], "-n")) norm = true;
        if (!std::strcmp(argv[a], "-d")) diff = true;
        if (!std::strcmp(argv[a], "-e")) bylevel = true;
        if (!std::strcmp(argv[a], "-u")) { bylevel = true; useemu = true; }
        if (!std::strcmp(argv[a], "-u0")) { bylevel = true; useemu = true; hw0only = true; }
        if (!std::strcmp(argv[a], "-s2") && a + 1 < argc) sim2 = argv[++a];
        else if (!std::strcmp(argv[a], "-l1") && a + 1 < argc) lab[1] = argv[++a];
        else if (!std::strcmp(argv[a], "-l2") && a + 1 < argc) lab[2] = argv[++a];
    }
    const std::vector<double> ed[2] = {{2, 2.5, 3, 3.5, 4, 5, 6, 8, 12}, {10, 15, 20, 25, 30, 35, 45, 60, 80, 120}};
    // [var: pT, E][sample: data, sim][ring][bin]
    // per EM level: [var][sample][level][near, far][bin], coarser bins
    const std::vector<double> el[2] = {{2, 3, 4, 6, 12}, {10, 20, 30, 45, 60, 120}};
    const int nS = sim2.empty() ? 2 : 3;
    std::vector<Dens> dl[2][3][3][2];
    for (int v = 0; v < 2; v++) for (int s = 0; s < 3; s++) for (int l = 0; l < 3; l++) for (int r = 0; r < 2; r++) dl[v][s][l][r].assign(el[v].size() - 1, Dens());
    double nlev[3][3] = {{0}};
    std::vector<Dens> d[2][3][kNL];
    for (int v = 0; v < 2; v++) for (int s = 0; s < 3; s++) for (int q = 0; q < kNL; q++) d[v][s][q].assign(ed[v].size() - 1, Dens());

    for (int s = 0; s < nS; s++)
    {
        TChain ch("jets");
        std::vector<std::string> f = Split(s < 2 ? argv[1 + s] : sim2.c_str());
        for (size_t i = 0; i < f.size(); i++) ch.Add(f[i].c_str());
        Float_t pt, E, w, eE[kNL], eA[kNL];
        Int_t em, hw0 = 1;
        ch.SetBranchAddress("jet_pt", &pt); ch.SetBranchAddress("jet_E", &E); ch.SetBranchAddress("weight", &w);
        ch.SetBranchAddress(useemu ? "emu" : "em", &em);
        if (hw0only && s == 0) ch.SetBranchAddress("hw0", &hw0);
        ch.SetBranchAddress(excl ? "lrx_E" : "lr_E", eE); ch.SetBranchAddress(excl ? "lrx_A" : "lr_A", eA);
        for (Long64_t i = 0; i < ch.GetEntries(); i++)
        {
            ch.GetEntry(i);
            const double ww = s ? w : 1.0;
            const double x[2] = {pt, E};
            for (int v = 0; v < 2; v++)
            {
                int b = -1;
                for (size_t q = 0; q + 1 < ed[v].size(); q++) if (x[v] >= ed[v][q] && x[v] < ed[v][q + 1]) b = (int)q;
                if (b < 0) continue;
                for (int q = 0; q < kNL; q++) if (eA[q] > 0) d[v][s][q][b].Add(eE[q], eA[q], ww);
            }
            if (em < 0 || em > 2 || !hw0) continue;
            if (E >= 10 && E < 120) nlev[s][em] += ww;
            for (int v = 0; v < 2; v++)
            {
                int b = -1;
                for (size_t q = 0; q + 1 < el[v].size(); q++) if (x[v] >= el[v][q] && x[v] < el[v][q + 1]) b = (int)q;
                if (b < 0) continue;
                if (eA[kNear] > 0) dl[v][s][em][0][b].Add(eE[kNear], eA[kNear], ww);
                if (eA[kFar] > 0) dl[v][s][em][1][b].Add(eE[kFar], eA[kFar], ww);
            }
        }
        std::printf("%s: %lld jets\n", lab[s].c_str(), ch.GetEntries());
    }
    for (int v = 0; v < 2; v++)
    {
        std::printf("\nUE density (GeV per unit live eta-phi area), rings %s: data | sim, sim/data\n", excl ? "without other jets" : "all points");
        std::printf(" jet %-9s", v ? "E" : "pT");
        for (int q = 0; q < kNL; q++) std::printf("        %.1f-%.1fR          ", kEdge[q], kEdge[q + 1]);
        std::printf("\n");
        for (size_t b = 0; b + 1 < ed[v].size(); b++)
        {
            std::printf(" %5.1f-%-5.1f", ed[v][b], ed[v][b + 1]);
            for (int q = 0; q < kNL; q++)
                std::printf("  %5.2f|%5.2f %5.2f    ", d[v][0][q][b].Val(), d[v][1][q][b].Val(), d[v][1][q][b].Val() / d[v][0][q][b].Val());
            std::printf("\n");
        }
    }

    gStyle->SetOptStat(0);
    const int cd = TColor::GetColor("#222222"), cs = TColor::GetColor("#d9480f");
    if (bylevel)
    {
        const char *lname[3] = {"EM0 only", "EM1 (not EM2)", "EM2"};
        std::printf("\nshare of jets by highest EM level (10-120 GeV): data %.3f %.3f %.3f | sim %.3f %.3f %.3f\n",
                    nlev[0][0] / (nlev[0][0] + nlev[0][1] + nlev[0][2]), nlev[0][1] / (nlev[0][0] + nlev[0][1] + nlev[0][2]), nlev[0][2] / (nlev[0][0] + nlev[0][1] + nlev[0][2]),
                    nlev[1][0] / (nlev[1][0] + nlev[1][1] + nlev[1][2]), nlev[1][1] / (nlev[1][0] + nlev[1][1] + nlev[1][2]), nlev[1][2] / (nlev[1][0] + nlev[1][1] + nlev[1][2]));
        TCanvas cv("cv", "", 1700, 1000);
        cv.Divide(3, 2, 0.002, 0.002);
        std::vector<TObject *> keep;
        for (int v = 0; v < 2; v++)
            for (int l = 0; l < 3; l++)
            {
                cv.cd(v * 3 + l + 1);
                gPad->SetMargin(0.14, 0.03, 0.13, 0.09);
                if (v == 1) gPad->SetLogx();
                TH1D *fr = new TH1D(Form("fre%d%d", v, l), "", 1, el[v].front(), el[v].back()); keep.push_back(fr);
                fr->SetTitle(Form(";%s;#rho(R-1.1R) - #rho(2.0-2.5R)  (GeV per unit #eta-#phi)", v ? "jet E (GeV)" : "jet p_{T} (GeV)"));
                fr->SetMinimum(-2.0); fr->SetMaximum(2.5); fr->GetYaxis()->SetTitleOffset(1.3);
                if (v == 1) fr->GetXaxis()->SetMoreLogLabels();
                fr->Draw();
                TLine *l0 = new TLine(el[v].front(), 0, el[v].back(), 0); l0->SetLineStyle(2); l0->SetLineColor(TColor::GetColor("#868e96")); l0->Draw(); keep.push_back(l0);
                TLegend *lg = new TLegend(0.17, 0.78, 0.6, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042); keep.push_back(lg);
                std::printf("\n%s, by jet %s: excess data | sim\n", lname[l], v ? "E" : "pT");
                for (int sm = 0; sm < 2; sm++)
                {
                    TGraphErrors *g = new TGraphErrors(); keep.push_back(g);
                    for (size_t b = 0; b + 1 < el[v].size(); b++)
                    {
                        double y, ey;
                        Value(dl[v][sm][l][0][b], dl[v][sm][l][1][b], false, y, ey, true);
                        if (!std::isfinite(y) || dl[v][sm][l][0][b].w2 <= 0) continue;
                        // need a handful of jets: effective count from the weights
                        const double neff = dl[v][sm][l][0][b].w * dl[v][sm][l][0][b].w / dl[v][sm][l][0][b].w2;
                        if (neff < 5) continue;
                        const double xc = v ? std::sqrt(el[v][b] * el[v][b + 1]) : 0.5 * (el[v][b] + el[v][b + 1]);
                        const double dx = v ? xc * 0.015 * (sm ? 1 : -1) : 0.05 * (sm ? 1 : -1);
                        const int i = g->GetN();
                        g->SetPoint(i, xc + dx, y); g->SetPointError(i, 0, ey);
                        std::printf("   %s %5.1f-%-5.1f  %+.2f+-%.2f (n_eff %.0f)\n", sm ? "sim " : "data", el[v][b], el[v][b + 1], y, ey, neff);
                    }
                    g->SetMarkerStyle(sm ? 21 : 20); g->SetMarkerSize(1.3);
                    g->SetMarkerColor(sm ? cs : cd); g->SetLineColor(sm ? cs : cd);
                    g->Draw("PL SAME");
                    lg->AddEntry(g, sm ? "simulation (reco)" : "data", "pl");
                }
                lg->Draw();
                TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.05);
                t.DrawLatex(0.14, 0.93, Form("highest %s: %s%s", useemu ? "emulated trigger" : "trigger", lname[l], hw0only ? " (data: EM0-recorded)" : ""));
            }
        cv.SaveAs(argv[3]);
        return 0;
    }
    if (diff)
    {
        std::printf("\nexcess density 1.0-1.1R minus 2.0-2.5R (GeV per unit live eta-phi area): data | sim\n");
        TCanvas cv("cv", "", 1500, 620);
        cv.Divide(2, 1, 0.002, 0.002);
        std::vector<TObject *> keep;
        for (int v = 0; v < 2; v++)
        {
            cv.cd(v + 1);
            gPad->SetMargin(0.13, 0.03, 0.13, 0.09);
            if (v == 1) gPad->SetLogx();
            TGraphErrors *g[3];
            double lo = 0, hi = 0;
            for (int s = 0; s < nS; s++)
            {
                g[s] = new TGraphErrors(); keep.push_back(g[s]);
                for (size_t b = 0; b + 1 < ed[v].size(); b++)
                {
                    double y, ey;
                    Value(d[v][s][kNear][b], d[v][s][kFar][b], false, y, ey, true);
                    if (!std::isfinite(y)) continue;
                    const double xc = v ? std::sqrt(ed[v][b] * ed[v][b + 1]) : 0.5 * (ed[v][b] + ed[v][b + 1]);
                    const double dx = v ? xc * 0.012 * (s - 1) : 0.03 * (s - 1);
                    const int i = g[s]->GetN();
                    g[s]->SetPoint(i, xc + dx, y); g[s]->SetPointError(i, 0, ey);
                    lo = std::min(lo, y - ey); hi = std::max(hi, y + ey);
                    std::printf(" %-22s %s %5.1f-%-5.1f  %+.2f+-%.2f\n", lab[s].c_str(), v ? "E " : "pT", ed[v][b], ed[v][b + 1], y, ey);
                }
            }
            TH1D *fr = new TH1D(Form("frd%d", v), "", 1, ed[v].front(), ed[v].back()); keep.push_back(fr);
            fr->SetTitle(Form(";%s;#rho(R-1.1R) - #rho(2.0-2.5R)  (GeV per unit #eta-#phi)", v ? "jet E (GeV)" : "jet p_{T} (GeV)"));
            lo = std::max(lo, -2.5); hi = std::min(hi, 3.0);  // keep poorly measured bins from setting the scale
            fr->SetMinimum(lo - 0.15 * (hi - lo)); fr->SetMaximum(hi + 0.35 * (hi - lo));
            fr->GetYaxis()->SetTitleOffset(1.2);
            if (v == 1) fr->GetXaxis()->SetMoreLogLabels();
            fr->Draw();
            TLine *l = new TLine(ed[v].front(), 0, ed[v].back(), 0); l->SetLineStyle(2); l->SetLineColor(TColor::GetColor("#868e96")); l->Draw(); keep.push_back(l);
            TLegend *lg = new TLegend(0.16, nS > 2 ? 0.73 : 0.78, 0.7, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04); keep.push_back(lg);
            const int c3 = TColor::GetColor("#1c6bb0");
            for (int s = 0; s < nS; s++)
            {
                const int col = s == 0 ? cd : (s == 1 ? cs : c3);
                g[s]->SetMarkerStyle(s == 0 ? 20 : (s == 1 ? 21 : 22)); g[s]->SetMarkerSize(1.3);
                g[s]->SetMarkerColor(col); g[s]->SetLineColor(col);
                g[s]->Draw("PL SAME");
                lg->AddEntry(g[s], lab[s].c_str(), "pl");
            }
            lg->Draw();
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
            t.DrawLatex(0.13, 0.93, Form("Energy density just outside the jet minus far from it%s", excl ? ", other jets excluded" : ""));
        }
        cv.SaveAs(argv[3]);
        return 0;
    }
    TCanvas cv("cv", "", 1500, 1000);
    std::vector<TObject *> keep;
    for (int v = 0; v < 2; v++)
    {
        for (int row = 0; row < 2; row++)
        {
            cv.cd();
            const double x0 = v * 0.5, x1 = x0 + 0.5, y0 = row == 0 ? 0.36 : 0.0, y1 = row == 0 ? 1.0 : 0.36;
            TPad *pd = new TPad(Form("p%d%d", v, row), "", x0, y0, x1, y1);
            pd->SetMargin(0.13, 0.03, row == 0 ? 0.02 : 0.25, row == 0 ? 0.09 : 0.02);
            if (v == 1) pd->SetLogx();
            pd->Draw(); pd->cd(); keep.push_back(pd);
            TH1D *fr = new TH1D(Form("fr%d%d", v, row), "", 1, ed[v].front(), ed[v].back());
            fr->SetTitle(Form(";%s;%s", v ? "jet E (GeV)" : "jet p_{T} (GeV)",
                              row == 0 ? (norm ? "ring density / 2.0-2.5R density (same sample)" : "point energy per live #eta-#phi area (GeV)") : "sim / data"));
            double ymax = 0;
            for (int s = 0; s < 2; s++)
                for (size_t b = 0; b + 1 < ed[v].size(); b++)
                {
                    double y, ey; Value(d[v][s][kNear][b], d[v][s][kFar][b], norm, y, ey);
                    if (std::isfinite(y)) ymax = std::max(ymax, y);
                }
            fr->SetMinimum(row == 0 ? 0 : (norm ? 0.6 : 0.3)); fr->SetMaximum(row == 0 ? 1.35 * ymax : (norm ? 2.2 : 1.6));
            fr->GetXaxis()->SetLabelSize(row == 0 ? 0 : 0.09); fr->GetXaxis()->SetTitleSize(0.10); fr->GetXaxis()->SetTitleOffset(1.0);
            fr->GetYaxis()->SetLabelSize(row == 0 ? 0.045 : 0.08); fr->GetYaxis()->SetTitleSize(row == 0 ? 0.05 : 0.09);
            fr->GetYaxis()->SetTitleOffset(row == 0 ? 1.2 : 0.65); fr->GetYaxis()->SetNdivisions(row == 0 ? 510 : 505);
            if (v == 1) fr->GetXaxis()->SetMoreLogLabels();
            fr->Draw(); keep.push_back(fr);
            if (row == 1 || norm) { TLine *l = new TLine(ed[v].front(), 1, ed[v].back(), 1); l->SetLineStyle(2); l->Draw(); keep.push_back(l); }
            TLegend *lg = row == 0 ? new TLegend(0.16, 0.70, 0.95, 0.89) : 0;
            if (lg) { lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038); lg->SetNColumns(2); keep.push_back(lg); }
            for (int ring = 0; ring < 2; ring++)
            {
                const int q = ring == 0 ? kNear : (norm ? kMid : kFar);
                for (int s = 0; s < (row == 0 ? 2 : 1); s++)
                {
                    TGraphErrors *g = new TGraphErrors();
                    for (size_t b = 0; b + 1 < ed[v].size(); b++)
                    {
                        const double xc = v ? std::sqrt(ed[v][b] * ed[v][b + 1]) : 0.5 * (ed[v][b] + ed[v][b + 1]);
                        const double dx = v ? xc * 0.012 * (ring ? 1 : -1) : 0.03 * (ring ? 1 : -1);
                        double ya, ea, yc, ec, y, ey;
                        Value(d[v][0][q][b], d[v][0][kFar][b], norm, ya, ea);
                        Value(d[v][1][q][b], d[v][1][kFar][b], norm, yc, ec);
                        if (row == 0) { y = s ? yc : ya; ey = s ? ec : ea; }
                        else { y = yc / ya; ey = y * std::sqrt(std::pow(ec / yc, 2) + std::pow(ea / ya, 2)); }
                        if (!std::isfinite(y)) continue;
                        const int i = g->GetN();
                        g->SetPoint(i, xc + dx, y); g->SetPointError(i, 0, ey);
                    }
                    const int col = row == 1 ? (ring ? TColor::GetColor("#868e96") : cs) : (s ? cs : cd);
                    g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerSize(1.3);
                    g->SetMarkerStyle(ring == 0 ? (s ? 21 : 20) : (s ? 25 : 24));
                    if (ring) g->SetLineStyle(2);
                    g->Draw("PL SAME"); keep.push_back(g);
                    if (lg) lg->AddEntry(g, Form("%s, %.1f-%.1fR", s ? "sim" : "data", kEdge[q], kEdge[q + 1]), "pl");
                    if (row == 1 && v == 0)
                    {
                        static TLegend l2(0.40, 0.76, 0.95, 0.96);
                        l2.SetBorderSize(0); l2.SetFillStyle(0); l2.SetTextSize(0.075); l2.SetNColumns(2);
                        l2.AddEntry(g, Form("%.1f-%.1fR", kEdge[q], kEdge[q + 1]), "pl");
                        if (ring == 1) l2.Draw();
                    }
                }
            }
            if (lg) lg->Draw();
            if (row == 0)
            {
                TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
                t.DrawLatex(0.13, 0.93, Form("Energy just outside the jet (reco)%s%s", norm ? ", relative to the far ring" : "", excl ? ", other jets excluded" : ""));
            }
        }
    }
    cv.SaveAs(argv[3]);
    return 0;
}
