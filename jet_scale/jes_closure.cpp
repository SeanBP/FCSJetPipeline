// Closure test of an EM-only JES map (CalibrationMap_em_only JSON) on
// cross-section-weighted simulation: applies the map in memory exactly as
// jet_calibration/ApplyCorrections_em_only.cpp does (bilinear interpolation
// of A, B at the jet's |x|, y; one evaluation at the jet's own energy
// clamped to the map's [Emin, Emax]; E_corr = E / scale) and reports the central
// value of E_corr / E_truth in bins of truth energy, for all positions and three
// radial bands, with the estimators the map is built from (arXiv:1609.05195,
// Sec. 2.1): the Gaussian-core mean (jet_scale/GaussCore.h) where the bin has
// n_eff >= core_min_neff, and below that, where the fit can lock onto a single
// heavily weighted jet (widths under 0.01, against 0.016 or more above it), the
// kernel peak (GaussCoreKernelPeak). The bootstrap spread is printed and drawn as
// the error, since the fit's own error shrinks with the same lock-on.
//
// The table also gives the core's width and the first-order non-closure
// expected from the inversion itself (their Eq. 24): C_mo ~ 1 + f''(x) / f'(x)^3
// sigma~(x)^2 / x, with f(x) = R_mo(x) x from the cores before correction (finite
// differences between bins) and sigma~ = core width x x.
//
// With "ue" the map is taken to be built from UE-subtracted energies (the
// default of JetEnergyScaleFineGrid_em_only): reco and truth energies are both
// UE subtracted and only pairs with both estimates are used. With "noue" raw
// energies are used throughout.
//
// It also prints the jet energy resolution per truth-energy bin and band, before
// and after correction, as "RES" rows at the end (plot_jes_resolution.C draws
// them): sigma_eff / central value, with sigma_eff the half-width of the shortest
// interval holding 68.27% of the weight (EffWidth) and bootstrap errors. The core's
// width is no resolution: it follows only the peak (9% at 10 GeV, against
// sigma_eff 25%) and collapses where the fit locks onto a heavily weighted jet.
//
// It also writes <output_base>.root with the 2D histograms of
// FCSJetAnalysis/CalibrationAnalysis_em.cpp, same names and binning
// (h2_response_Energy_before/after: truth vs reco energy; h2_scale_Energy_
// before/after: truth energy vs reco/truth), but cross-section weighted and
// for matched jets inside the reco fiducial region where the map is defined.
//
// With a fifth argument, the jet_calibration.root the map was built from, it
// also compares the map with the points it was fitted to: the
// n_eff-weighted mean of scale_mu - map(x, y, response_mu) per truth-energy
// bin and radial band. A residual there is a failure of the fitted form (or of
// the smoothing); a non-closure without one comes from the inversion itself.
//
// Usage: ./jes_closure <file list> <map.json> [ue|noue] [output_base] [jet_calibration.root]
// Build (from this folder): g++ -m64 -O2 jes_closure.cpp ../shared/VectorDict.cxx -o jes_closure -I../shared
//            `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
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
#include "GaussCore.h"
#include <TTree.h>
#include <TH2D.h>

#include "JetPairs.h"

struct Map
{
    std::vector<double> x, y;
    std::vector<std::vector<double> > par[2];   // A B
    double C, D, Emin, Emax;                    // the same everywhere
    double xmin, xmax, ymin, ymax;
    int imin, imax, jmin, jmax;                 // the bins' index range
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

static double JsonNumber(const std::string &json, const std::string &key)
{
    const size_t k = json.find("\"" + key + "\"");
    if (k == std::string::npos) { std::cerr << "ERROR: map has no \"" << key << "\"" << std::endl; std::exit(1); }
    const size_t pos = json.find(':', k);
    return std::strtod(json.c_str() + pos + 1, nullptr);
}

static bool ReadMap(const std::string &file, Map &m)
{
    std::ifstream in(file.c_str());
    if (!in.is_open()) return false;
    std::stringstream buf;
    buf << in.rdbuf();
    const std::string json = buf.str();
    // Log-cubic maps only, as ApplyCorrections_em_only
    {
        const size_t mk = json.find("\"model\"");
        const size_t v = mk == std::string::npos ? mk : json.find('"', json.find(':', mk));
        if (v == std::string::npos || json.compare(v, 10, "\"logcubic\"") != 0)
        {
            std::cerr << "ERROR: " << file << " is not a log-cubic map (no \"model\": \"logcubic\")" << std::endl;
            return false;
        }
    }
    m.C = JsonNumber(json, "C");
    m.D = JsonNumber(json, "D");
    m.Emin = JsonNumber(json, "Emin");
    m.Emax = JsonNumber(json, "Emax");
    m.x = JsonArray(json, "x_grid");
    m.y = JsonArray(json, "y_grid");
    for (int k = 0; k < 2; k++)
        m.par[k].assign(m.x.size(), std::vector<double>(m.y.size(), std::numeric_limits<double>::quiet_NaN()));
    const char *names[2] = {"A", "B"};
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
        for (int k = 0; k < 2; k++) m.par[k][i][j] = field(names[k]);
        imin = std::min(imin, i); imax = std::max(imax, i); jmin = std::min(jmin, j); jmax = std::max(jmax, j);
    }
    if (imax < 0) return false;
    m.xmin = m.x[imin]; m.xmax = m.x[imax]; m.ymin = m.y[jmin]; m.ymax = m.y[jmax];
    m.imin = imin; m.imax = imax; m.jmin = jmin; m.jmax = jmax;
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
    // cell kept inside the map, so a jet on the upper edge does not reach the
    // NaN cells beyond it (0 x NaN is NaN), as in ApplyCorrections_em_only
    const int i = std::min(std::max(FindIndex(m.x, x), m.imin), m.imax - 1);
    const int j = std::min(std::max(FindIndex(m.y, y), m.jmin), m.jmax - 1);
    const double tx = (x - m.x[i]) / (m.x[i + 1] - m.x[i]), ty = (y - m.y[j]) / (m.y[j + 1] - m.y[j]);
    double p[2];
    for (int k = 0; k < 2; k++)
    {
        const std::vector<std::vector<double> > &Z = m.par[k];
        const double a = Z[i][j] * (1 - tx) + Z[i + 1][j] * tx, b = Z[i][j + 1] * (1 - tx) + Z[i + 1][j + 1] * tx;
        p[k] = a * (1 - ty) + b * ty;
    }
    const double l = std::log(std::min(std::max(E, m.Emin), m.Emax));
    return p[0] + l * (p[1] + l * (m.C + l * m.D));
}

static const double core_min_neff = 250;

// The bin's central value: Gaussian core, or kernel peak below core_min_neff
static double Central(const std::vector<double> &v, const std::vector<double> &w, bool use_core)
{
    return use_core ? GaussCoreFit(v, w).mu : GaussCoreKernelPeak(v, w);
}

// Spread of the central value over bootstrap resamples of the bin (fewer for the big low-energy bins)
static double CentralBootstrap(const std::vector<double> &v, const std::vector<double> &w, bool use_core, TRandom3 &rnd)
{
    const int nboot = v.size() > 20000 ? 20 : 100;
    std::vector<double> bv(v.size()), bw(v.size());
    double s = 0, s2 = 0;
    int n = 0;
    for (int ib = 0; ib < nboot; ib++)
    {
        for (size_t i = 0; i < v.size(); i++) { const size_t j = rnd.Integer(v.size()); bv[i] = v[j]; bw[i] = w[j]; }
        const double c = Central(bv, bw, use_core);
        if (std::isfinite(c)) { s += c; s2 += c * c; n++; }
    }
    return n > 1 ? std::sqrt(std::max(0.0, s2 / n - (s / n) * (s / n))) : NAN;
}

// Half-width of the shortest interval holding 68.27% of the weight
static double EffWidth(const std::vector<double> &v, const std::vector<double> &w)
{
    std::vector<std::pair<double, double> > a(v.size());
    double tot = 0;
    for (size_t i = 0; i < v.size(); i++) { a[i] = std::make_pair(v[i], w[i]); tot += w[i]; }
    std::sort(a.begin(), a.end());
    double best = NAN, run = 0;
    size_t r = 0;
    for (size_t l = 0; l < a.size(); l++)
    {
        while (r < a.size() && run < 0.6827 * tot) run += a[r++].second;
        if (run < 0.6827 * tot) break;
        const double width = a[r - 1].first - a[l].first;
        if (!(width >= best)) best = width;
        run -= a[l].second;
    }
    return 0.5 * best;
}

// Spread of EffWidth over bootstrap resamples of the bin (fewer for the big low-energy bins)
static double EffWidthBootstrap(const std::vector<double> &v, const std::vector<double> &w, TRandom3 &rnd)
{
    const int nboot = v.size() > 20000 ? 20 : 100;
    std::vector<double> bv(v.size()), bw(v.size());
    double s = 0, s2 = 0;
    int n = 0;
    for (int ib = 0; ib < nboot; ib++)
    {
        for (size_t i = 0; i < v.size(); i++) { const size_t j = rnd.Integer(v.size()); bv[i] = v[j]; bw[i] = w[j]; }
        const double c = EffWidth(bv, bw);
        if (std::isfinite(c)) { s += c; s2 += c * c; n++; }
    }
    return n > 1 ? std::sqrt(std::max(0.0, s2 / n - (s / n) * (s / n))) : NAN;
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
    if (!LoadJetPairs(ReadFileList(argv[1]), pairs)) return 1;

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
    // central values and their errors, kept for the resolution: [region][0 before, 1 after][bin]
    std::vector<std::vector<std::vector<double> > > cen(nreg, std::vector<std::vector<double> >(2, std::vector<double>(nb, NAN))), cen_err = cen;
    for (int g = 0; g < nreg; g++)
    {
        std::cout << "\n== " << reg_name[g] << " ==\n  truth E (GeV)    jets   before    after  [core width]  Eq.24 C_mo  boot before  after  n_eff  (1 core, 2 kernel peak; core needs n_eff >= " << core_min_neff << ")" << std::endl;
        // central values before / after, and the Eq. 24 expectation
        std::vector<double> mb(nb, NAN), ma(nb, NAN), eb_b(nb, NAN), eb_a(nb, NAN), ne(nb, 0.0);
        std::vector<GaussCoreResult> cb(nb), ca(nb);
        std::vector<bool> use_core(nb, false);
        TRandom3 rnd(12345 + g);
        for (int b = 0; b < nb; b++)
        {
            if (after[g][b].size() < 50) continue;
            std::vector<double> vb, va, ww;
            for (size_t i = 0; i < before[g][b].size(); i++) { vb.push_back(before[g][b][i].first); va.push_back(after[g][b][i].first); ww.push_back(before[g][b][i].second); }
            double sw = 0, sw2 = 0;
            for (size_t i = 0; i < ww.size(); i++) { sw += ww[i]; sw2 += ww[i] * ww[i]; }
            ne[b] = sw * sw / sw2;
            use_core[b] = ne[b] >= core_min_neff;
            if (use_core[b]) { cb[b] = GaussCoreFit(vb, ww); ca[b] = GaussCoreFit(va, ww); mb[b] = cb[b].mu; ma[b] = ca[b].mu; }
            else { mb[b] = GaussCoreKernelPeak(vb, ww); ma[b] = GaussCoreKernelPeak(va, ww); }
            eb_b[b] = CentralBootstrap(vb, ww, use_core[b], rnd); eb_a[b] = CentralBootstrap(va, ww, use_core[b], rnd);
            cen[g][0][b] = mb[b]; cen[g][1][b] = ma[b]; cen_err[g][0][b] = eb_b[b]; cen_err[g][1][b] = eb_a[b];
        }
        // q = 0 core, 1 kernel peak
        TGraphErrors *gb[2], *ga[2];
        for (int q = 0; q < 2; q++) { gb[q] = new TGraphErrors(); ga[q] = new TGraphErrors(); }
        double worst = 1, worst_E = 0;
        for (int b = 0; b < nb; b++)
        {
            if (after[g][b].size() < 50) continue;
            const double ec = 0.5 * (eb[b] + eb[b + 1]);
            // f(x) = R(x) x at the neighbouring bin centres (cores before correction)
            double cmo = NAN;
            if (b > 0 && b < nb - 1 && std::isfinite(cb[b - 1].mu) && std::isfinite(cb[b].mu) && std::isfinite(cb[b + 1].mu))
            {
                const double x0 = 0.5 * (eb[b - 1] + eb[b]), x1 = ec, x2 = 0.5 * (eb[b + 1] + eb[b + 2]);
                const double f0 = cb[b - 1].mu * x0, f1 = cb[b].mu * x1, f2 = cb[b + 1].mu * x2;
                const double d1 = (f2 - f0) / (x2 - x0);
                const double d2 = 2 * ((f2 - f1) / (x2 - x1) - (f1 - f0) / (x1 - x0)) / (x2 - x0);
                const double st = cb[b].sigma * x1;
                if (d1 > 0) cmo = 1 + d2 / (d1 * d1 * d1) * st * st / x1;
            }
            std::printf("  %6.1f-%6.1f  %8zu   %8.4f %8.4f  [%6.4f]  %8.4f  %8.4f %8.4f  %9.1f  %d\n", eb[b], eb[b + 1], after[g][b].size(), mb[b], ma[b],
                        ca[b].sigma, cmo, eb_b[b], eb_a[b], ne[b], use_core[b] ? 1 : 2);
            if (!std::isfinite(mb[b]) || !std::isfinite(ma[b])) continue;
            const int q = use_core[b] ? 0 : 1;
            gb[q]->SetPoint(gb[q]->GetN(), ec, mb[b]); gb[q]->SetPointError(gb[q]->GetN() - 1, 0, eb_b[b]);
            ga[q]->SetPoint(ga[q]->GetN(), ec, ma[b]); ga[q]->SetPointError(ga[q]->GetN() - 1, 0, eb_a[b]);
            if (std::fabs(ma[b] - 1) > std::fabs(worst - 1)) { worst = ma[b]; worst_E = ec; }
        }
        std::cout << "  largest deviation after correction: " << worst << " at " << worst_E << " GeV" << std::endl;
        c.cd(g + 1);
        gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.11);
        TH1F *fr = gPad->DrawFrame(0, 0.5, 100, 1.5);
        fr->GetXaxis()->SetTitle(ue ? "truth jet energy - UE (GeV)" : "truth jet energy (GeV)");
        fr->GetYaxis()->SetTitle("reco / truth energy");
        TLine *l1 = new TLine(0, 1, 100, 1); l1->SetLineStyle(2); l1->Draw();
        TLine *l2 = new TLine(0, 1.02, 100, 1.02); l2->SetLineStyle(3); l2->SetLineColor(kGray + 1); l2->Draw();
        TLine *l3 = new TLine(0, 0.98, 100, 0.98); l3->SetLineStyle(3); l3->SetLineColor(kGray + 1); l3->Draw();
        gb[0]->SetMarkerStyle(24); gb[1]->SetMarkerStyle(26); ga[0]->SetMarkerStyle(20); ga[1]->SetMarkerStyle(22);
        TLegend *leg = new TLegend(0.5, 0.68, 0.96, 0.9);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
        for (int q = 0; q < 2; q++)
        {
            gb[q]->SetMarkerColor(kGray + 2); gb[q]->SetLineColor(kGray + 2); ga[q]->SetMarkerColor(kRed + 1); ga[q]->SetLineColor(kRed + 1);
            gb[q]->Draw("PL same"); ga[q]->Draw("PL same");
            keep.push_back(gb[q]); keep.push_back(ga[q]);
        }
        leg->AddEntry(gb[0], "core, before correction", "pl");
        leg->AddEntry(ga[0], "core, after correction", "pl");
        leg->AddEntry(gb[1], "kernel peak, before correction", "pl");
        leg->AddEntry(ga[1], "kernel peak, after correction", "pl");
        leg->Draw();
        TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.05); t->SetTextFont(62);
        t->DrawLatex(0.12, 0.935, reg_name[g]);
        keep.push_back(leg); keep.push_back(t);
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

    std::cout << "\nResolution, sigma_eff / central value of reco / truth energy (region 0 all, 1 inner, 2 middle, 3 outer):\n"
              << "RES region lo hi jets before err_before after err_after" << std::endl;
    TRandom3 rnd_res(54321);
    for (int g = 0; g < nreg; g++)
        for (int b = 0; b < nb; b++)
        {
            if (after[g][b].size() < 50) continue;
            std::vector<double> ww;
            for (size_t i = 0; i < after[g][b].size(); i++) ww.push_back(after[g][b][i].second);
            double res[2], res_err[2];
            for (int q = 0; q < 2; q++)
            {
                const std::vector<std::pair<double, double> > &src = q == 0 ? before[g][b] : after[g][b];
                std::vector<double> v;
                for (size_t i = 0; i < src.size(); i++) v.push_back(src[i].first);
                const double se = EffWidth(v, ww), se_err = EffWidthBootstrap(v, ww, rnd_res), m = cen[g][q][b];
                res[q] = se / m;
                res_err[q] = res[q] * std::sqrt(std::pow(se_err / se, 2) + std::pow(cen_err[g][q][b] / m, 2));
            }
            std::printf("RES %d %.1f %.1f %zu %.4f %.4f %.4f %.4f\n", g, eb[b], eb[b + 1], after[g][b].size(), res[0], res_err[0], res[1], res_err[1]);
        }
    c.SaveAs((out_base + ".png").c_str());
    c.SaveAs((out_base + ".pdf").c_str());
    return 0;
}
