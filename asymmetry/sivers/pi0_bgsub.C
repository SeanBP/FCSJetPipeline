// Background-subtracted A_N of the FMS-style pi0 candidates (Sean 2026-11-16;
// AnConfig.h "FMS-style pi0 candidates"), as the FMS did it (STAR PRD 103,
// 092009, Sec. III):
//   A_sig = f A_pi0 + (1 - f) A_bg  ->  A_pi0 = (A_sig - (1 - f) A_sb) / f
// The inputs are two AnAsymmetry -I outputs of the same yields and class, one
// for the signal window (M_gg 0-0.2 GeV) and one for the sideband (0.2-0.3 GeV).
// A_bg is taken as the sideband asymmetry; the FMS assumed a constant
// background asymmetry.
// The signal fraction f comes per x_F bin from the signal output's pi0_mass_<row>
// spectrum:
//   - a Gaussian plus a quadratic background is fitted to the spectrum summed
//     over x_F in 0.03-0.3 GeV;
//   - each bin is refitted with that Gaussian mean and width fixed (free
//     amplitude and background), so sparse bins stay stable;
//   - f = Gaussian integral / counts in 0-0.2 GeV. It is the same for both
//     beams, since they share the candidates.
// Error: sqrt(dA_sig^2 + (1 - f)^2 dA_sb^2) / f (statistical; f's error and the
// small correlation of candidates sharing a jet are not included).
// Writes <out>.root (gr_AN_<beam>_<row>_pi0, gr_AN_<row>_pi0), <out>.txt and
// <out>_mass.png (the fits per x_F bin, the window boundary marked).
// Usage (compiled, ACLiC; ROOT 5's CINT cannot run it):
//   root -b -q 'pi0_bgsub.C+("an_iso_sig.root", "an_iso_sb.root", "an_iso_pi0", "all", "isolated #pi^{0} (z_{em} > 0.98)")'
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

struct Row { double lo, hi, xf, an, err, n, cosm, pol; };

// rows of an AnAsymmetry table for one level and class: key (beam, x_F bin low edge * 1000)
static std::map<std::pair<std::string, int>, Row> ReadTable(const char *root, const char *level, const char *cls)
{
    std::map<std::pair<std::string, int>, Row> out;
    TString tn(root);
    tn.ReplaceAll(".root", ".txt");
    ifstream in(tn.Data());
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        std::string l, c, b;
        Row r;
        if (!(is >> l >> c >> b >> r.lo >> r.hi >> r.xf >> r.an >> r.err >> r.n >> r.cosm >> r.pol)) continue;
        if (l != level || c != cls) continue;
        out[std::make_pair(b, (int)TMath::Nint(r.lo * 1000))] = r;
    }
    return out;
}

void pi0_bgsub(const char *sig_root, const char *sb_root, const char *out, const char *row = "all", const char *label = "#pi^{0} candidates")
{
    gStyle->SetOptStat(0);
    // AnAsymmetry -I fills the np1 rows; "all" holds the same entries
    const char *cls = "np1";
    std::map<std::pair<std::string, int>, Row> ts = ReadTable(sig_root, row, cls), tb = ReadTable(sb_root, row, cls);
    TFile fs(sig_root);
    TH2D *m2 = (TH2D *)fs.Get(Form("pi0_mass_%s", row));
    if (!m2) { printf("ERROR: no pi0_mass_%s in %s (not an AnAsymmetry -I output)\n", row, sig_root); return; }
    const double kSigHi = 0.2;
    // global shape
    TH1D *mall = m2->ProjectionY("mall");
    TF1 *fg = new TF1("fg", "gaus(0)+pol2(3)", 0.03, 0.3);
    {
        const int b0 = mall->FindBin(0.06), b1 = mall->FindBin(0.2);
        int bm = b0;
        for (int b = b0; b <= b1; b++) if (mall->GetBinContent(b) > mall->GetBinContent(bm)) bm = b;
        fg->SetParameters(mall->GetBinContent(bm), mall->GetBinCenter(bm), 0.02, 0, 0, 0);
        fg->SetParLimits(2, 0.005, 0.08);
        mall->Fit(fg, "RQN0");
    }
    const double mu = fg->GetParameter(1), sg = fabs(fg->GetParameter(2));
    printf("pi0_bgsub: %s, row %s: global fit mean %.4f GeV, sigma %.4f GeV, chi2/ndf %.1f/%d\n", label, row, mu, sg, fg->GetChisquare(), fg->GetNDF());
    // per x_F bin: f
    const int nx = m2->GetNbinsX();
    std::vector<double> fsig(nx + 1, -1);
    std::vector<TH1D *> hb(nx + 1, (TH1D *)0);
    std::vector<TF1 *> fb(nx + 1, (TF1 *)0);
    for (int i = 1; i <= nx; i++)
    {
        TH1D *h = m2->ProjectionY(Form("m_%d", i), i, i);
        const double nsig = h->Integral(1, h->FindBin(kSigHi - 1e-6));
        if (nsig < 20) { delete h; continue; }
        TF1 *f = new TF1(Form("f_%d", i), "gaus(0)+pol2(3)", 0.03, 0.3);
        f->SetParameters(h->GetMaximum(), mu, sg, 0, 0, 0);
        f->FixParameter(1, mu);
        f->FixParameter(2, sg);
        h->Fit(f, "RQN0");
        TF1 g("g", "gaus", 0, 0.3);
        g.SetParameters(f->GetParameter(0), mu, sg);
        const double s = g.Integral(0, kSigHi) / h->GetBinWidth(1);
        fsig[i] = TMath::Max(0.0, TMath::Min(1.0, s / nsig));
        hb[i] = h; fb[i] = f;
    }
    // combine
    TString tn = TString(out) + ".txt";
    FILE *tab = fopen(tn.Data(), "w");
    fprintf(tab, "# background-subtracted pi0 A_N: (A_sig - (1-f) A_sb) / f; %s; row %s; sig %s, sb %s\n", label, row, sig_root, sb_root);
    fprintf(tab, "# global mass fit: mean %.4f sigma %.4f GeV\n", mu, sg);
    fprintf(tab, "# beam xF_lo xF_hi <xF> A_pi0 stat f A_sig stat_sig A_sb stat_sb n_sig n_sb\n");
    TFile fo(TString(out) + ".root", "RECREATE");
    TGraphAsymmErrors *gboth = new TGraphAsymmErrors();
    gboth->SetName(Form("gr_AN_%s_pi0", row));
    gboth->SetTitle(Form("%s;x_{F};A_{N}", label));
    const char *beams[2] = {"blue", "yellow"};
    printf("  %-6s %13s %8s %18s %6s %18s %18s\n", "beam", "x_F bin", "<x_F>", "A_pi0 +- stat", "f", "A_sig", "A_sb");
    for (int b = 0; b < 2; b++)
    {
        TGraphAsymmErrors *g = new TGraphAsymmErrors();
        g->SetName(Form("gr_AN_%s_%s_pi0", beams[b], row));
        g->SetTitle(Form("%s beam, %s;x_{F};A_{N}", beams[b], label));
        for (std::map<std::pair<std::string, int>, Row>::iterator it = ts.begin(); it != ts.end(); ++it)
        {
            if (it->first.first != beams[b]) continue;
            const Row &s = it->second;
            std::map<std::pair<std::string, int>, Row>::iterator jt = tb.find(it->first);
            if (jt == tb.end()) continue;
            const Row &k = jt->second;
            // the stored (blue-frame) x_F bin of this row
            const double blo = b == 0 ? s.lo : -s.hi;
            const int ib = m2->GetXaxis()->FindBin(blo + 1e-6);
            if (ib < 1 || ib > nx || fsig[ib] <= 0.05) continue;
            const double f = fsig[ib];
            const double an = (s.an - (1 - f) * k.an) / f;
            const double dan = sqrt(s.err * s.err + (1 - f) * (1 - f) * k.err * k.err) / f;
            const double exl = b == 0 ? s.xf - s.lo : s.xf - s.lo, exh = s.hi - s.xf;
            g->SetPoint(g->GetN(), s.xf, an); g->SetPointError(g->GetN() - 1, exl, exh, dan, dan);
            gboth->SetPoint(gboth->GetN(), s.xf, an); gboth->SetPointError(gboth->GetN() - 1, exl, exh, dan, dan);
            fprintf(tab, "%s %.3f %.3f %.4f %+.6f %.6f %.4f %+.6f %.6f %+.6f %.6f %.0f %.0f\n", beams[b], s.lo, s.hi, s.xf, an, dan, f, s.an, s.err, k.an, k.err, s.n, k.n);
            printf("  %-6s %6.3f-%6.3f %8.4f %+9.5f +- %.5f %6.3f %+9.5f +- %.5f %+9.5f +- %.5f\n", beams[b], s.lo, s.hi, s.xf, an, dan, f, s.an, s.err, k.an, k.err);
        }
        g->Write();
    }
    gboth->Sort();
    gboth->Write();
    TNamed("pi0_bgsub", Form("%s; row %s; global fit mean %.4f sigma %.4f GeV", label, row, mu, sg)).Write();
    fo.Close();
    fclose(tab);
    // mass fits per x_F bin
    int nfit = 0;
    for (int i = 1; i <= nx; i++) if (hb[i]) nfit++;
    if (nfit == 0) { printf("no x_F bin with >= 20 candidates in the signal window\n"); return; }
    const int ncol = TMath::Min(4, nfit), nrow = (nfit + ncol - 1) / ncol;
    TCanvas *c = new TCanvas("c", "", 450 * ncol, 360 * nrow + 60);
    TPad *mp = new TPad("mp", "", 0, 0, 1, 1 - 1.15 * 60.0 / (360 * nrow + 60)); mp->Draw(); mp->Divide(ncol, nrow, 0.003, 0.003);
    int ip = 0;
    for (int i = 1; i <= nx; i++)
    {
        if (!hb[i]) continue;
        mp->cd(++ip); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13);
        TH1D *h = hb[i];
        h->SetTitle(Form("x_{F} %.3f-%.3f;M_{#gamma#gamma} (GeV);candidates", m2->GetXaxis()->GetBinLowEdge(i), m2->GetXaxis()->GetBinUpEdge(i)));
        h->SetLineColor(kBlack); h->SetMarkerStyle(20); h->SetMarkerSize(0.6);
        h->SetMinimum(0); h->SetMaximum(1.35 * h->GetMaximum());
        h->Draw("E");
        const double ymax = h->GetMaximum();
        // signal / sideband boundary (ROOT 5 batch PNGs drop transparent fills)
        TLine *lw = new TLine(0.2, 0, 0.2, ymax); lw->SetLineColor(kAzure + 2); lw->SetLineStyle(2); lw->SetLineWidth(2); lw->Draw();
        fb[i]->SetLineColor(kRed + 1); fb[i]->SetLineWidth(2); fb[i]->Draw("same");
        TF1 *bg = new TF1(Form("bg_%d", i), "pol2", 0.03, 0.3);
        for (int k = 0; k < 3; k++) bg->SetParameter(k, fb[i]->GetParameter(3 + k));
        bg->SetLineColor(kGray + 2); bg->SetLineStyle(2); bg->SetLineWidth(2); bg->Draw("same");
        TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.055);
        l.DrawLatex(0.6, 0.82, Form("f_{sig} = %.2f", fsig[i]));
    }
    c->cd();
    const double hh = 60.0 / (360 * nrow + 60);
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.34 * hh);
    t.DrawLatex(0.01, 1 - 0.42 * hh, Form("#bf{%s}: M_{#gamma#gamma} per #pi^{0} x_{F} bin; signal 0-0.2, sideband 0.2-0.3 GeV (blue line)", label));
    t.DrawLatex(0.01, 1 - 0.9 * hh, Form("red: Gaussian (#mu %.3f, #sigma %.3f GeV, fixed) + quadratic; gray: background", mu, sg));
    c->SaveAs(Form("%s_mass.png", out));
}
