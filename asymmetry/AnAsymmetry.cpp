// A_N from the yield histograms of AnHistograms: a left-right asymmetry per
// beam, by the cross ratio (luminosity and acceptance cancel, so no relative
// luminosity is needed):
//   eps = (sqrt(L_up R_down) - sqrt(L_down R_up)) / (sqrt(L_up R_down) + sqrt(L_down R_up))
//   A_N = eps / (P <|cos phi|>)
// with L, R the jets left and right of the polarized beam and up, down its
// spin. The other beam's spin is held fixed: eps is formed separately for its
// up and down states and the two are averaged (inverse variance). Within one
// state the other beam's own asymmetry is the same for both spins of this
// beam, so the cross ratio removes it; summing over the other beam instead
// lets it leak in wherever the spin-state luminosities differ (the toy,
// an_toy_closure.cpp, with luminosities 0.30/0.20/0.25/0.25: up to 7 sigma).
//
// Sides. STAR coordinates: +z west, +y up, +x south; reco_x is in them (FCS
// north points, det 0, all have x < 0). Left of a beam is spin x momentum
// direction with spin up (+y):
//   blue   (moving +z, towards the FCS): left = y x z = +x  (x_F > 0)
//   yellow (moving -z, away from the FCS): left = y x (-z) = -x  (x_F < 0)
// So the yellow beam's left and right are the blue beam's right and left,
// and its x_F is minus the stored (blue-frame) x_F: the yellow A_N is
// reported at negative x_F. Spin states (AnConfig.h): blue up BuYu, BuYd;
// yellow up BuYu, BdYu.
//
// A jet at azimuth phi from the beam's left carries A_N cos(phi), so a
// left-right eps / P is A_N times the jets' mean |cos phi|; <|cos phi|> per
// bin is h_abscos_* over the summed yields (AnHistograms.cpp), the same for
// both beams. Small asymmetries: N_L,up ~ 1 + P A_N <cos phi>_L.
//
// Polarization. With one yield file and P_blue, P_yellow given, P is that
// constant. With -r, the yields of many files (one run each, AnHistograms'
// "runs" entry) are summed and P per run comes from make_run_polarization.py's
// table (the fill's polarization at the middle of the run); each bin then uses
// the mean P of its jets, sum_runs N_run P_run / sum_runs N_run, as the FMS
// EM-jet analysis weighted runs by events. Runs not in the table, with no
// value for a beam, or with more than one run in a file are left out (listed).
//
// Bad runs (-r only). Runs in a bad-run list are left out. By default the list
// is D. Kapukchyan's pi0 QA list, read from his directory at every call, so
// his updates apply without changes here (kBadRunListDefault: "Runnumber |
// Flag | AdditionalNotes" lines). The list has a line for every run he
// reviewed (as of 2026-02-07, all 413 analysis runs up to 23024003; 59 with a
// flag such as LowStat, EnShift, DoublePeak), so an empty flag reads as a run
// that passed. Mode "flagged" (default) drops the flagged runs; "all" drops
// every listed run, i.e. the whole reviewed period. The
// list's path, modification time, mode and the runs it removed are printed and
// written to the output as TNamed "bad_runs". "none" uses no list.
//
// x_F. Each point sits at the bin's mean corrected x_F (h_xfsum_* over the
// summed yields; negative for yellow), its horizontal error bar spanning the
// bin.
//
// Statistical errors. Each cross ratio's error is from Poisson counts,
// d eps = sqrt(L_up R_down L_down R_up) / (sqrt(L_up R_down) + sqrt(L_down R_up))^2
//         x sqrt(1/L_up + 1/L_down + 1/R_up + 1/R_down),
// the two spin states of the other beam are combined with inverse-variance
// weights, and the result is divided by <P><|cos phi|> (dA_N = d eps /
// (<P><|cos phi|>)); the polarization and <|cos phi|> uncertainties are not in
// it (systematic, and statistically negligible). The yields must be raw
// counts: a yield histogram filled with weights is refused. The vertical
// error bars of every graph and the "stat" column of the table are these
// statistical errors only.
//
// EM levels. A jet counts in one level only (the lowest listed one it passes),
// so the levels add up: "all" is their sum and gives the final result.
//
// Outputs: <out.root> with the graphs and TNamed "bad_runs" (-r), and
// <out.root minus .root>.txt, the table of every point (level, photon class,
// beam, x_F bin, <x_F>, A_N, stat error, jets, <|cos phi|>, <P>).
// plot_an.C draws the final A_N vs x_F from <out.root>.
//
// Usage: ./AnAsymmetry <yields.root> <out.root> <P_blue> <P_yellow> [rebin = 2]
//        ./AnAsymmetry -r <yield file list> <run polarization table> <out.root> [rebin = 2]
//                         [bad-run list = kBadRunListDefault | none] [flagged | all]
//   rebin 2 gives x_F bins of 0.05, as the FMS EM-jet analysis
//   Prints, for each EM level and all levels, A_N per x_F bin of both beams
//   with all photon classes summed; writes graphs
//   gr_AN_<beam>_<EM level|all>_<photon class|all> (x_F > 0 blue, x_F < 0
//   yellow) and gr_AN_<EM level|all>_<class|all> with both beams.
// Build: g++ -m64 -O2 AnAsymmetry.cpp -o AnAsymmetry `root-config --cflags --libs` -std=c++11

#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TDatime.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "AnConfig.h"

static const char *const kBadRunListDefault = "/star/u/dkap7827/Fcs2019/FcsAsymRun22/BadRuns.txt";

enum Beam { kBeamBlue, kBeamYellow };
static const char *const kBeamName[2] = {"blue", "yellow"};

// spin states (AnConfig.h order BuYu, BuYd, BdYu, BdYd) with the beam's spin
// up / down, [beam][other beam's spin: up, down]
static const int kUp[2][2] = {{0, 1}, {0, 2}}, kDown[2][2] = {{2, 3}, {1, 3}};
// the side (AnConfig.h: 0 x <= 0, 1 x > 0) left of the beam
static const int kLeftSide[2] = {1, 0};

struct Asym { double eps, err; bool ok; };

// cross-ratio asymmetry and its Poisson error
static Asym CrossRatio(double Lu, double Ld, double Ru, double Rd)
{
    Asym a = {0, 0, false};
    if (!(Lu > 0 && Ld > 0 && Ru > 0 && Rd > 0)) return a;
    const double p = std::sqrt(Lu * Rd), m = std::sqrt(Ld * Ru);
    a.eps = (p - m) / (p + m);
    a.err = p * m / ((p + m) * (p + m)) * std::sqrt(1 / Lu + 1 / Ld + 1 / Ru + 1 / Rd);
    a.ok = true;
    return a;
}

static TString YieldName(int l, int k, int s, int p)
{
    return TString::Format("h_xF_%s_%s_%s_%s", kEmLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]);
}

// Histograms summed over the input files, and per (EM level, photon class)
// the yields summed over sides and spin states weighted by each run's P
struct Sums
{
    std::map<std::string, TH1D *> h;
    TH1D *pw[2][kNEmLevels][kNPhotonClasses] = {};
    // adds one file's histograms, its jets weighted by pol[beam]; false if any is missing
    bool Add(TFile &f, const double pol[2])
    {
        std::vector<std::string> names;
        for (int l = 0; l < kNEmLevels; l++)
            for (int k = 0; k < kNPhotonClasses; k++)
            {
                names.push_back(TString::Format("h_abscos_%s_%s", kEmLevelName[l], PhotonClassName(k)).Data());
                names.push_back(TString::Format("h_xfsum_%s_%s", kEmLevelName[l], PhotonClassName(k)).Data());
                for (int s = 0; s < kNSides; s++)
                    for (int p = 0; p < kNSpinStates; p++) names.push_back(YieldName(l, k, s, p).Data());
            }
        for (size_t i = 0; i < names.size(); i++)
        {
            TH1D *hh = (TH1D *)f.Get(names[i].c_str());
            if (!hh) { std::cerr << "ERROR: " << f.GetName() << " has no " << names[i] << std::endl; return false; }
            // yields must be raw counts, so their Poisson errors are sqrt(N)
            if (names[i].compare(0, 5, "h_xF_") == 0 && hh->GetSumw2N() > 0)
                for (int b = 1; b <= hh->GetNbinsX(); b++)
                    if (std::fabs(hh->GetSumw2()->At(b) - hh->GetBinContent(b)) > 1e-6 * (1 + hh->GetBinContent(b)))
                    {
                        std::cerr << "ERROR: " << f.GetName() << " " << names[i] << " is weighted; the statistical errors need raw counts" << std::endl;
                        return false;
                    }
        }
        for (size_t i = 0; i < names.size(); i++)
        {
            TH1D *in = (TH1D *)f.Get(names[i].c_str());
            if (h.count(names[i])) h[names[i]]->Add(in);
            else { h[names[i]] = (TH1D *)in->Clone((names[i] + "_sum").c_str()); h[names[i]]->SetDirectory(0); }
        }
        for (int l = 0; l < kNEmLevels; l++)
            for (int k = 0; k < kNPhotonClasses; k++)
            {
                TH1D *n = (TH1D *)((TH1D *)f.Get(YieldName(l, k, 0, 0)))->Clone("n_tmp");
                n->SetDirectory(0);
                n->Reset();
                for (int s = 0; s < kNSides; s++)
                    for (int p = 0; p < kNSpinStates; p++) n->Add((TH1D *)f.Get(YieldName(l, k, s, p)));
                for (int b = 0; b < 2; b++)
                {
                    if (!pw[b][l][k]) { pw[b][l][k] = (TH1D *)n->Clone(TString::Format("pw_%d_%d_%d", b, l, k)); pw[b][l][k]->SetDirectory(0); pw[b][l][k]->Reset(); }
                    pw[b][l][k]->Add(n, pol[b]);
                }
                delete n;
            }
        return true;
    }
};

// Bad-run list: run -> flag ("" when the flag column is empty); false if unreadable
static bool ReadBadRuns(const std::string &path, std::map<int, std::string> &bad)
{
    std::ifstream in(path.c_str());
    if (!in) return false;
    std::string line;
    while (std::getline(in, line))
    {
        std::istringstream is(line);
        int run = 0;
        if (!(is >> run) || run <= 0) continue;  // header, blank
        std::string rest, flag;
        std::getline(is, rest);
        const size_t a = rest.find('|');
        if (a != std::string::npos)
        {
            const size_t b = rest.find('|', a + 1);
            flag = rest.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1);
            const size_t f0 = flag.find_first_not_of(" \t"), f1 = flag.find_last_not_of(" \t");
            flag = f0 == std::string::npos ? "" : flag.substr(f0, f1 - f0 + 1);
        }
        bad[run] = flag;
    }
    return true;
}

// the run of an AnHistograms file (its "runs" entry), 0 if it holds none or several
static int RunOf(TFile &f)
{
    TNamed *r = (TNamed *)f.Get("runs");
    if (!r) return 0;
    std::istringstream is(r->GetTitle());
    int run = 0, extra = 0;
    if (!(is >> run) || (is >> extra)) return 0;
    return run;
}

int main(int argc, char **argv)
{
    const bool per_run = argc > 1 && std::string(argv[1]) == "-r";
    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0] << " <yields.root> <out.root> <P_blue> <P_yellow> [rebin = 2]\n"
                  << "       " << argv[0] << " -r <yield file list> <run polarization table> <out.root> [rebin = 2] [bad-run list | none] [flagged | all]" << std::endl;
        return 1;
    }
    const int rebin = argc > 5 ? std::atoi(argv[5]) : 2;
    if (rebin < 1 || kNXfBins % rebin) { std::cerr << "ERROR: bad rebin" << std::endl; return 1; }

    Sums sums;
    std::string pol_note, bad_note;
    if (!per_run)
    {
        const double pol[2] = {std::atof(argv[3]), std::atof(argv[4])};
        if (!(pol[0] > 0 && pol[1] > 0)) { std::cerr << "ERROR: bad polarization" << std::endl; return 1; }
        TFile in(argv[1]);
        if (in.IsZombie() || !sums.Add(in, pol)) return 1;
        pol_note = TString::Format("P_blue %.3f, P_yellow %.3f", pol[0], pol[1]).Data();
    }
    else
    {
        std::map<int, std::pair<double, double> > table;
        {
            std::ifstream t(argv[3]);
            std::string line;
            while (std::getline(t, line))
            {
                if (line.empty() || line[0] == '#') continue;
                std::istringstream is(line);
                int run, fill; double tmid, pb, dpb, py, dpy;
                if (is >> run >> fill >> tmid >> pb >> dpb >> py >> dpy) table[run] = std::make_pair(pb, py);
            }
        }
        // bad-run list, read now so the current version applies
        const std::string bad_path = argc > 6 ? argv[6] : kBadRunListDefault;
        const std::string bad_mode = argc > 7 ? argv[7] : "flagged";
        if (bad_mode != "all" && bad_mode != "flagged") { std::cerr << "ERROR: bad-run mode must be all or flagged" << std::endl; return 1; }
        std::map<int, std::string> bad;
        if (bad_path != "none")
        {
            if (!ReadBadRuns(bad_path, bad)) { std::cerr << "ERROR: cannot read bad-run list " << bad_path << " (give \"none\" to run without one)" << std::endl; return 1; }
            FileStat_t st;
            const bool have_st = gSystem->GetPathInfo(bad_path.c_str(), st) == 0;
            int n_flagged = 0;
            for (std::map<int, std::string>::const_iterator b = bad.begin(); b != bad.end(); ++b) n_flagged += !b->second.empty();
            bad_note = TString::Format("%s (modified %s), mode %s: %zu runs listed, %d with a flag", bad_path.c_str(),
                                       have_st ? TDatime((UInt_t)st.fMtime).AsSQLString() : "unknown", bad_mode.c_str(), bad.size(), n_flagged).Data();
        }
        else bad_note = "no bad-run list";
        std::printf("AnAsymmetry: bad runs: %s\n", bad_note.c_str());
        std::ifstream list(argv[2]);
        std::string f;
        int n_used = 0;
        std::vector<std::string> skipped, bad_removed;
        while (list >> f)
        {
            TFile in(f.c_str());
            const int run = in.IsZombie() ? 0 : RunOf(in);
            std::map<int, std::pair<double, double> >::const_iterator it = table.find(run);
            if (run == 0 || it == table.end() || !(it->second.first > 0 && it->second.second > 0))
            {
                skipped.push_back(f + (run == 0 ? " (not one run)" : it == table.end() ? " (run not in the polarization table)" : " (no polarization for a beam)"));
                continue;
            }
            std::map<int, std::string>::const_iterator bi = bad.find(run);
            if (bi != bad.end() && (bad_mode == "all" || !bi->second.empty()))
            {
                bad_removed.push_back(TString::Format("%d%s%s", run, bi->second.empty() ? "" : " ", bi->second.c_str()).Data());
                continue;
            }
            const double pol[2] = {it->second.first, it->second.second};
            if (!sums.Add(in, pol)) return 1;
            n_used++;
        }
        std::printf("AnAsymmetry: %d files used, %zu bad runs removed, %zu left out otherwise\n", n_used, bad_removed.size(), skipped.size());
        for (size_t i = 0; i < bad_removed.size(); i++) std::printf("  bad run: %s\n", bad_removed[i].c_str());
        for (size_t i = 0; i < skipped.size(); i++) std::printf("  left out: %s\n", skipped[i].c_str());
        bad_note += "; removed:";
        for (size_t i = 0; i < bad_removed.size(); i++) bad_note += " " + bad_removed[i];
        if (n_used == 0) return 1;
        pol_note = std::string("per-run P from ") + argv[3];
    }

    TFile out(argv[per_run ? 4 : 2], "RECREATE");
    if (per_run) TNamed("bad_runs", bad_note.c_str()).Write();
    const std::string table_name = std::string(argv[per_run ? 4 : 2]).substr(0, std::string(argv[per_run ? 4 : 2]).rfind(".root")) + ".txt";
    FILE *table = std::fopen(table_name.c_str(), "w");
    std::fprintf(table, "# A_N = eps / (<P> <|cos phi|>), %s; stat = statistical error only\n", pol_note.c_str());
    if (per_run) std::fprintf(table, "# bad runs: %s\n", bad_note.c_str());
    std::fprintf(table, "# level class beam xF_lo xF_hi <xF> A_N stat jets <|cos|> <P>\n");
    for (int L = 0; L <= kNEmLevels; L++)
    {
        // levels summed: one level, or all (L = kNEmLevels)
        const int l0 = L < kNEmLevels ? L : 0, l1 = L < kNEmLevels ? L + 1 : kNEmLevels;
        const char *lname = L < kNEmLevels ? kEmLevelName[L] : "all";
        const int l = L;
        // [photon class or kNPhotonClasses = all]: yields [side][spin state], and per bin the jets
        // summed over sides and spin states, sum of |cos phi|, sum of x_F, sum of P x jets per beam
        TH1D *y[kNPhotonClasses + 1][kNSides][kNSpinStates], *n[kNPhotonClasses + 1], *mcos[kNPhotonClasses + 1], *mxf[kNPhotonClasses + 1], *mpol[2][kNPhotonClasses + 1];
        for (int k = 0; k <= kNPhotonClasses; k++)
        {
            const int k0 = k < kNPhotonClasses ? k : 0, k1 = k < kNPhotonClasses ? k + 1 : kNPhotonClasses;
            const TString tag = TString::Format("%d_%d", l, k);
            bool first = true;
            for (int ll = l0; ll < l1; ll++)
            for (int kk = k0; kk < k1; kk++, first = false)
            {
                TH1D *c = (TH1D *)sums.h[TString::Format("h_abscos_%s_%s", kEmLevelName[ll], PhotonClassName(kk)).Data()]->Rebin(rebin, "tmp_c");
                TH1D *x = (TH1D *)sums.h[TString::Format("h_xfsum_%s_%s", kEmLevelName[ll], PhotonClassName(kk)).Data()]->Rebin(rebin, "tmp_x");
                if (first) { mcos[k] = (TH1D *)c->Clone("c_" + tag); mxf[k] = (TH1D *)x->Clone("x_" + tag); }
                else { mcos[k]->Add(c); mxf[k]->Add(x); }
                delete c; delete x;
                for (int b = 0; b < 2; b++)
                {
                    TH1D *pr = (TH1D *)sums.pw[b][ll][kk]->Rebin(rebin, "tmp_p");
                    if (first) mpol[b][k] = (TH1D *)pr->Clone(TString::Format("p%d_", b) + tag);
                    else mpol[b][k]->Add(pr);
                    delete pr;
                }
                for (int s = 0; s < kNSides; s++)
                    for (int p = 0; p < kNSpinStates; p++)
                    {
                        TH1D *r = (TH1D *)sums.h[YieldName(ll, kk, s, p).Data()]->Rebin(rebin, "tmp_y");
                        if (first) y[k][s][p] = (TH1D *)r->Clone(TString::Format("y_%d_%d_", s, p) + tag);
                        else y[k][s][p]->Add(r);
                        delete r;
                    }
            }
            n[k] = (TH1D *)y[k][0][0]->Clone("n_" + tag);
            n[k]->Reset();
            for (int s = 0; s < kNSides; s++)
                for (int p = 0; p < kNSpinStates; p++) n[k]->Add(y[k][s][p]);
            // per-bin means: <|cos phi|>, <x_F>, <P> per beam
            mcos[k]->Divide(n[k]);
            mxf[k]->Divide(n[k]);
            for (int b = 0; b < 2; b++) mpol[b][k]->Divide(n[k]);
        }
        const int nb = y[0][0][0]->GetNbinsX();

        for (int k = 0; k <= kNPhotonClasses; k++)
        {
            const char *cname = k < kNPhotonClasses ? PhotonClassName(k) : "all";
            TGraphAsymmErrors *g[2], *gboth = new TGraphAsymmErrors();
            gboth->SetName(TString::Format("gr_AN_%s_%s", lname, cname));
            gboth->SetTitle(TString::Format("%s, %s;x_{F};A_{N}", lname, cname));
            if (k == kNPhotonClasses)
                std::printf("\n%s, all photon classes: A_N = eps / (P <|cos phi|>) (%s)\n  %-6s %13s %8s %18s %10s %8s %7s\n",
                            lname, pol_note.c_str(), "beam", "x_F bin", "<x_F>", "A_N +- stat", "jets", "<|cos|>", "<P>");
            for (int b = 0; b < 2; b++)
            {
                g[b] = new TGraphAsymmErrors();
                g[b]->SetName(TString::Format("gr_AN_%s_%s_%s", kBeamName[b], lname, cname));
                g[b]->SetTitle(TString::Format("%s beam, %s, %s;x_{F};A_{N}", kBeamName[b], lname, cname));
                const int sl = kLeftSide[b], sr = 1 - sl;
                for (int i = 1; i <= nb; i++)
                {
                    // one cross ratio per spin state q of the other beam, then their weighted mean
                    double sw = 0, swe = 0;
                    for (int q = 0; q < 2; q++)
                    {
                        const double Lu = y[k][sl][kUp[b][q]]->GetBinContent(i), Ld = y[k][sl][kDown[b][q]]->GetBinContent(i);
                        const double Ru = y[k][sr][kUp[b][q]]->GetBinContent(i), Rd = y[k][sr][kDown[b][q]]->GetBinContent(i);
                        const Asym c = CrossRatio(Lu, Ld, Ru, Rd);
                        if (!c.ok || !(c.err > 0)) continue;
                        sw += 1 / (c.err * c.err);
                        swe += c.eps / (c.err * c.err);
                    }
                    const double mc = mcos[k]->GetBinContent(i), P = mpol[b][k]->GetBinContent(i);
                    if (!(sw > 0) || !(mc > 0) || !(P > 0)) continue;
                    const double eps = swe / sw, deps = 1 / std::sqrt(sw);
                    // x_F of the polarized beam's frame: the stored x_F for blue, its negative for yellow
                    const double sign = b == kBeamBlue ? 1 : -1;
                    const double lo = y[k][0][0]->GetBinLowEdge(i), hi = lo + y[k][0][0]->GetBinWidth(i);
                    const double xf = sign * mxf[k]->GetBinContent(i);
                    const double exl = b == kBeamBlue ? xf - lo : xf + hi, exh = b == kBeamBlue ? hi - xf : -lo - xf;
                    const double an = eps / (P * mc), dan = deps / (P * mc);
                    g[b]->SetPoint(g[b]->GetN(), xf, an); g[b]->SetPointError(g[b]->GetN() - 1, exl, exh, dan, dan);
                    gboth->SetPoint(gboth->GetN(), xf, an); gboth->SetPointError(gboth->GetN() - 1, exl, exh, dan, dan);
                    std::fprintf(table, "%s %s %s %.3f %.3f %.4f %+.6f %.6f %.0f %.4f %.4f\n", lname, cname, kBeamName[b],
                                 b == kBeamBlue ? lo : -hi, b == kBeamBlue ? hi : -lo, xf, an, dan, n[k]->GetBinContent(i), mc, P);
                    if (k == kNPhotonClasses)
                        std::printf("  %-6s %6.3f-%6.3f %8.4f %+9.5f +- %.5f %10.0f %8.4f %7.4f\n", kBeamName[b],
                                    b == kBeamBlue ? lo : -hi, b == kBeamBlue ? hi : -lo, xf, an, dan, n[k]->GetBinContent(i), mc, P);
                }
                g[b]->Write();
            }
            gboth->Sort();
            gboth->Write();
        }
    }
    out.Close();
    std::fclose(table);
    std::printf("\nwrote %s and %s\n", argv[per_run ? 4 : 2], table_name.c_str());
    return 0;
}
