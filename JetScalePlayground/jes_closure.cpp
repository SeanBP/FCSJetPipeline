// Closure test of an EM-only JES map (CalibrationMap_em_only JSON) on
// cross-section-weighted simulation: applies the map in memory exactly as
// jet_calibration/ApplyCorrections_em_only.cpp does (bilinear interpolation
// of the parameters at the jet's |x|, y; one evaluation at the jet's own
// energy; E_corr = E / scale) and reports the weighted median of
// E_corr / E_truth in bins of truth energy, for all positions and three
// radial bands. The map is built from weighted medians, so the closure has
// to be weighted too (FCSJetAnalysis/CalibrationAnalysis_em.cpp is not).
//
// With "ue" the map is taken to be built from UE-subtracted energies (the
// default of JetEnergyScaleFineGrid_em_only): reco and truth energies are both
// UE subtracted and only pairs with both estimates are used. With "noue" raw
// energies are used throughout.
//
// It also writes <output_base>.root with the 2D histograms of
// FCSJetAnalysis/CalibrationAnalysis_em.cpp, same names and binning
// (h2_response_Energy_before/after: truth vs reco energy; h2_scale_Energy_
// before/after: truth energy vs reco/truth), but cross-section weighted and
// for matched jets inside the reco fiducial region where the map is defined.
//
// With a fifth argument, the jet_calibration.root the map was built from, it
// also compares the map with the per-bin medians it was fitted to: the
// n_eff-weighted mean of scale_mu - map(x, y, response_mu) per truth-energy
// bin and radial band. A residual there is a failure of the fitted form (or of
// the smoothing); a non-closure without one comes from the inversion itself.
//
// Usage: ./jes_closure <file list> <map.json> [ue|noue] [output_base] [jet_calibration.root]
// Build: g++ -m64 -O2 jes_closure.cpp ../shared/VectorDict.cxx -o jes_closure -I../shared -I../jet_scale
//            -I../UnfoldingPlayground `root-config --cflags --libs` -std=c++11

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

#include <algorithm>
#include <numeric>
#include <sstream>
#include <limits>
#include <map>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>

#include "UnfoldingInput.h"

struct Map
{
    std::vector<double> x, y;
    std::vector<std::vector<double> > par[7];   // A B C D E0 sigmaL sigmaR
    double xmin, xmax, ymin, ymax;
};

static std::vector<double> JsonArray(const std::string &json, const std::string &key)
{
    std::vector<double> out;
    size_t pos = json.find("\"" + key + "\"");
    pos = json.find('[', pos);
    const size_t end = json.find(']', pos);
    std::stringstream ss(json.substr(pos + 1, end - pos - 1));
    double v;
    char comma;
    while (ss >> v) { out.push_back(v); ss >> comma; }
    return out;
}

static bool ReadMap(const std::string &file, Map &m)
{
    std::ifstream in(file.c_str());
    if (!in.is_open()) return false;
    std::stringstream buf;
    buf << in.rdbuf();
    const std::string json = buf.str();
    m.x = JsonArray(json, "x_grid");
    m.y = JsonArray(json, "y_grid");
    for (int k = 0; k < 7; k++)
        m.par[k].assign(m.x.size(), std::vector<double>(m.y.size(), std::numeric_limits<double>::quiet_NaN()));
    const char *names[7] = {"A", "B", "C", "D", "E0", "sigmaL", "sigmaR"};
    int imin = 1 << 30, imax = -1, jmin = 1 << 30, jmax = -1;
    size_t pos = json.find('[', json.find("\"bins\"")) + 1;
    while (true)
    {
        pos = json.find_first_of("{]", pos);
        if (pos == std::string::npos || json[pos] == ']') break;
        const size_t end = json.find('}', pos);
        const std::string obj = json.substr(pos, end - pos + 1);
        pos = end + 1;
        auto field = [&](const std::string &k) {
            size_t p = obj.find("\"" + k + "\"");
            p = obj.find(':', p) + 1;
            return std::strtod(obj.c_str() + p, nullptr);
        };
        const int i = (int)field("i"), j = (int)field("j");
        for (int k = 0; k < 7; k++) m.par[k][i][j] = field(names[k]);
        imin = std::min(imin, i); imax = std::max(imax, i); jmin = std::min(jmin, j); jmax = std::max(jmax, j);
    }
    if (imax < 0) return false;
    m.xmin = m.x[imin]; m.xmax = m.x[imax]; m.ymin = m.y[jmin]; m.ymax = m.y[jmax];
    return true;
}

static int FindIndex(const std::vector<double> &g, double v)
{
    if (v <= g.front()) return 0;
    if (v >= g.back()) return (int)g.size() - 2;
    return (int)(std::upper_bound(g.begin(), g.end(), v) - g.begin()) - 1;
}

// The map's scale at (|x|, y, E), or NaN where undefined (as ApplyCorrections_em_only)
static double Scale(const Map &m, double x, double y, double E)
{
    x = std::min(std::max(std::fabs(x), m.xmin), m.xmax);
    y = std::min(std::max(y, m.ymin), m.ymax);
    const int i = FindIndex(m.x, x), j = FindIndex(m.y, y);
    const double tx = (x - m.x[i]) / (m.x[i + 1] - m.x[i]), ty = (y - m.y[j]) / (m.y[j + 1] - m.y[j]);
    double p[7];
    for (int k = 0; k < 7; k++)
    {
        const std::vector<std::vector<double> > &Z = m.par[k];
        const double a = Z[i][j] * (1 - tx) + Z[i + 1][j] * tx, b = Z[i][j + 1] * (1 - tx) + Z[i + 1][j + 1] * tx;
        p[k] = a * (1 - ty) + b * ty;
    }
    const double e = std::max(E, 1.0);
    double g = 0.0;
    if (p[3] != 0.0 && std::isfinite(p[3]))
    {
        const double s = e < p[4] ? p[5] : p[6];
        if (s != 0.0 && std::isfinite(s)) g = p[3] * std::exp(-0.5 * (e - p[4]) * (e - p[4]) / (s * s));
    }
    return p[0] + p[1] / e + p[2] / (e * e) + g;
}

static double WeightedMedian(std::vector<std::pair<double, double> > &v)
{
    std::sort(v.begin(), v.end());
    double tot = 0;
    for (size_t i = 0; i < v.size(); i++) tot += v[i].second;
    double cum = 0;
    for (size_t i = 0; i < v.size(); i++) { cum += v[i].second; if (cum >= 0.5 * tot) return v[i].first; }
    return NAN;
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    if (argc < 3) { std::cerr << "usage: ./jes_closure <file list> <map.json> [ue|noue] [output_base]" << std::endl; return 1; }
    const bool ue = !(argc > 3 && std::string(argv[3]) == "noue");
    const std::string out_base = argc > 4 ? argv[4] : "jes_closure";

    Map m;
    if (!ReadMap(argv[2], m)) { std::cerr << "ERROR: cannot read map " << argv[2] << std::endl; return 1; }
    std::vector<JetPair> pairs;
    if (!LoadJetPairs(ReadFileList(argv[1]), kMatchPipeline, pairs)) return 1;

    // truth-energy bins: those of JetEnergyScaleFineGrid_em_only (3-100 GeV, power 1.5)
    const int nb = 20;
    std::vector<double> eb(nb + 1);
    for (int i = 0; i <= nb; i++) eb[i] = 3.0 + 97.0 * std::pow((double)i / nb, 1.5);
    const double r1 = 66.35, r2 = 87.53;
    const int nreg = 4;
    std::vector<std::vector<std::vector<std::pair<double, double> > > > before(nreg, std::vector<std::vector<std::pair<double, double> > >(nb)),
        after = before;
    // CalibrationAnalysis_em.cpp's binning: truth energy 3-100 GeV in 50 bins of
    // E0 + (Emax - E0) (i/N)^2, reco energy 0-100 GeV in 100 bins, ratio 0-3 in 50
    std::vector<double> e2(51);
    for (int i = 0; i <= 50; i++) e2[i] = 3.0 + 97.0 * std::pow(i / 50.0, 2.0);
    TFile fh((out_base + ".root").c_str(), "RECREATE");
    TH2D *h_resp[2], *h_scale[2];
    const char *when[2] = {"before", "after"};
    for (int q = 0; q < 2; q++)
    {
        h_resp[q] = new TH2D(Form("h2_response_Energy_%s", when[q]), Form("Response (Energy, %s);Truth Jet Energy [GeV];Reco Jet Energy [GeV]", when[q]),
                             50, &e2[0], 100, 0.0, 100.0);
        h_scale[q] = new TH2D(Form("h2_scale_Energy_%s", when[q]), Form("Scale (Energy, %s);Truth Jet Energy [GeV];Reco/Truth", when[q]),
                              50, &e2[0], 50, 0.0, 3.0);
        h_resp[q]->Sumw2(); h_scale[q]->Sumw2();
    }
    long n_used = 0, n_invalid = 0, n_skip = 0;
    for (size_t k = 0; k < pairs.size(); k++)
    {
        const JetPair &p = pairs[k];
        if (ue && !p.has_ue) { n_skip++; continue; }
        const double reco = ue ? p.rEsub() : p.rE, truth = ue ? p.tEsub() : p.tE;
        const double s = Scale(m, p.rx, p.ry, reco);
        if (!std::isfinite(s) || std::fabs(s) < 1e-8) { n_invalid++; continue; }
        const double rv[2] = {reco, reco / s};
        for (int q = 0; q < 2; q++)
        {
            if (truth > 0 && truth < 100 && rv[q] > 0 && rv[q] < 100) h_resp[q]->Fill(truth, rv[q], p.w);
            if (truth > 0) h_scale[q]->Fill(truth, rv[q] / truth, p.w);
        }
        if (truth < eb.front() || truth >= eb.back()) continue;
        const int b = (int)(std::upper_bound(eb.begin(), eb.end(), truth) - eb.begin()) - 1;
        const double r = std::sqrt(p.rx * p.rx + p.ry * p.ry);
        const int regs[2] = {0, r < r1 ? 1 : (r < r2 ? 2 : 3)};
        for (int q = 0; q < 2; q++)
        {
            before[regs[q]][b].push_back(std::make_pair(reco / truth, (double)p.w));
            after[regs[q]][b].push_back(std::make_pair(reco / s / truth, (double)p.w));
        }
        n_used++;
    }
    std::cout << n_used << " pairs used (" << (ue ? "UE subtracted at both levels" : "raw energies") << "), " << n_invalid
              << " where the map is undefined, " << n_skip << " without UE estimates" << std::endl;

    const char *reg_name[nreg] = {"all positions", "inner third", "middle third", "outer third"};
    TCanvas c("c", "", 1400, 1000);
    c.Divide(2, 2, 0.001, 0.001);
    std::vector<TObject *> keep;
    for (int g = 0; g < nreg; g++)
    {
        std::cout << "\n== " << reg_name[g] << " ==\n  truth E (GeV)    jets   median scale before   after" << std::endl;
        TGraph *gb = new TGraph(), *ga = new TGraph();
        double worst = 1, worst_E = 0;
        for (int b = 0; b < nb; b++)
        {
            if (after[g][b].size() < 50) continue;
            const double mb = WeightedMedian(before[g][b]), ma = WeightedMedian(after[g][b]);
            const double ec = 0.5 * (eb[b] + eb[b + 1]);
            std::printf("  %6.1f-%6.1f  %8zu   %8.4f   %8.4f\n", eb[b], eb[b + 1], after[g][b].size(), mb, ma);
            gb->SetPoint(gb->GetN(), ec, mb);
            ga->SetPoint(ga->GetN(), ec, ma);
            if (std::fabs(ma - 1) > std::fabs(worst - 1)) { worst = ma; worst_E = ec; }
        }
        std::cout << "  largest deviation after correction: " << worst << " at " << worst_E << " GeV" << std::endl;
        c.cd(g + 1);
        gPad->SetLogx(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.11);
        TH1F *fr = gPad->DrawFrame(3, 0.5, 100, 1.5);
        fr->GetXaxis()->SetTitle(ue ? "truth jet energy - UE (GeV)" : "truth jet energy (GeV)");
        fr->GetYaxis()->SetTitle("weighted median reco / truth energy");
        fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent();
        TLine *l1 = new TLine(3, 1, 100, 1); l1->SetLineStyle(2); l1->Draw();
        TLine *l2 = new TLine(3, 1.02, 100, 1.02); l2->SetLineStyle(3); l2->SetLineColor(kGray + 1); l2->Draw();
        TLine *l3 = new TLine(3, 0.98, 100, 0.98); l3->SetLineStyle(3); l3->SetLineColor(kGray + 1); l3->Draw();
        gb->SetMarkerStyle(24); gb->SetMarkerColor(kGray + 2); gb->SetLineColor(kGray + 2); gb->Draw("PL same");
        ga->SetMarkerStyle(20); ga->SetMarkerColor(kRed + 1); ga->SetLineColor(kRed + 1); ga->Draw("PL same");
        TLegend *leg = new TLegend(0.55, 0.75, 0.96, 0.9);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
        leg->AddEntry(gb, "before correction", "pl");
        leg->AddEntry(ga, "after correction", "pl");
        leg->Draw();
        TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.05); t->SetTextFont(62);
        t->DrawLatex(0.12, 0.935, reg_name[g]);
        keep.push_back(gb); keep.push_back(ga); keep.push_back(leg); keep.push_back(t);
        keep.push_back(l1); keep.push_back(l2); keep.push_back(l3);
    }
    fh.cd();
    for (int q = 0; q < 2; q++) { h_resp[q]->Write(); h_scale[q]->Write(); }
    fh.Close();

    if (argc > 5)
    {
        TFile fb(argv[5]);
        TTree *t = (TTree *)fb.Get("bins");
        double El, Eh, xc, yc, smu, rmu, neff;
        t->SetBranchAddress("E_low", &El); t->SetBranchAddress("E_high", &Eh);
        t->SetBranchAddress("x_center", &xc); t->SetBranchAddress("y_center", &yc);
        t->SetBranchAddress("scale_mu", &smu); t->SetBranchAddress("response_mu", &rmu);
        t->SetBranchAddress("n_eff", &neff);
        std::map<double, std::vector<double> > sum_w, sum_wr, sum_wabs;
        for (Long64_t k = 0; k < t->GetEntries(); k++)
        {
            t->GetEntry(k);
            if (!(smu > 0) || !(rmu > 0) || neff < 5) continue;
            const double sc = Scale(m, xc, yc, rmu);
            if (!std::isfinite(sc)) continue;
            const double r = std::sqrt(xc * xc + yc * yc);
            const int g = r < r1 ? 1 : (r < r2 ? 2 : 3);
            std::vector<double> &w = sum_w[El], &wr = sum_wr[El], &wa = sum_wabs[El];
            if (w.empty()) { w.assign(nreg, 0); wr.assign(nreg, 0); wa.assign(nreg, 0); }
            const int gs[2] = {0, g};
            for (int q = 0; q < 2; q++) { w[gs[q]] += neff; wr[gs[q]] += neff * (smu - sc); wa[gs[q]] += neff * std::fabs(smu - sc); }
        }
        std::cout << "\nMap minus fitted points, n_eff-weighted mean (and mean |difference|), per truth-energy bin:" << std::endl;
        std::cout << "  E_low      all              inner            middle           outer" << std::endl;
        for (std::map<double, std::vector<double> >::iterator it = sum_w.begin(); it != sum_w.end(); ++it)
        {
            std::printf("  %6.1f", it->first);
            for (int g = 0; g < nreg; g++)
                std::printf("   %+7.4f (%6.4f)", it->second[g] > 0 ? -sum_wr[it->first][g] / it->second[g] : 0.0,
                            it->second[g] > 0 ? sum_wabs[it->first][g] / it->second[g] : 0.0);
            std::printf("\n");
        }
    }
    c.SaveAs((out_base + ".png").c_str());
    c.SaveAs((out_base + ".pdf").c_str());
    return 0;
}
