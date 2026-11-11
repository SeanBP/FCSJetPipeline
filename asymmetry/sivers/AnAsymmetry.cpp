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
// written to the output as TNamed "bad_runs". "none" uses no list. The runs
// of AnConfig.h kSuppBadRuns (hot runs from the 20% pass QA) are left out in
// every case.
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
// Spin convention. Yield files carry TNamed "spin_convention" (AnConfig.h);
// files made before 2026-11-04 do not and used the old Spin_config mapping,
// which inverted both beams' spins: their states are relabelled p -> 3 - p
// (BuYu <-> BdYd, BuYd <-> BdYu) when summed.
//
// Trigger levels (AnConfig.h): EM0-2 and JP0-2, judged independently: a jet
// firing both is in its EM level and in its JP level, so a level's A_N has
// every jet that fired it. Within a family a jet is in one level only (the
// lowest listed), so "all" (the EM levels) and "allJP" (the JP levels) are
// sums without double counting. With -J, "all" is the EM + JP union: the EM
// levels plus the JP jets that did not fire EM at the chosen factor (the JP
// histograms' EM tags, AnConfig.h kEmTagSuffix), each jet once. Yield files keep the triggering jets per slice
// of trigger value / threshold (_t115, _t130, _t145, _t160); the slices at
// or above the factor are summed: -fe <f> for the EM levels (default
// kEmCutFactor = 1.3), -fj <f> for the JP levels (default kJpCutFactor =
// 1.3), -f <f> both (1.15, 1.3, 1.45 or 1.6). The -J union takes the JP jets
// whose EM slice is below the EM factor. -V adds the EM jets that failed only the HCal-veto margin
// (_veto). Yield files made before 2026-11-05 (no "trigger_rule") have
// neither and are refused.
//
// Pi0 mass window (-m in:<mean>:<sigma> / -m out:<mean>:<sigma>). The np2 rows
// then hold the two-photon jets whose pair mass is inside (in) or outside
// (out) mean +- 2 sigma of the Gaussian fitted to the pi0 peak (mgg_peak.C;
// Sean 2026-11-11: the window in units of the peak width). The older form
// in:<peak> / out:<peak> still takes 0.5-1.5 x <peak>. The yields come from
// AnHistograms' x_F vs mass yields (h_xFm_*, E <= sqrt(s)/2, jets at the
// default factor 1.3, -J union from _jponly): a mass bin (25 MeV)
// counts as inside if its centre is in the window, and the window actually
// used is printed. <|cos phi|>, <x_F> and the run-weighted P per bin are those
// of all np2 jets (the mass split is not stored for them; the jets in and out
// of the window share the acceptance, so the difference is small). Only the
// np2 rows are meaningful in this mode; it excludes -E, -C, -A, -V and factors
// other than 1.3.
//
// Leading-point pi0 probability range (-p <lo>:<hi>). The np1 rows (and so
// "all" and np12) then hold every jet with P_lead (reco_pi0_lead_prob) in
// [lo, hi), from AnHistograms' x_F vs P_lead yields (h_xFP_*, E <= sqrt(s)/2,
// jets at the default factor 1.3, -J union from _jponly): a 0.05 P_lead bin
// counts if its centre is in the range. Their <|cos phi|> and <x_F> come from
// h_abscosP_* / h_xfsumP_* of the same jets. The other photon-class rows are
// empty in this mode. Same restrictions as -m (which it excludes); yield files
// without h_xFP_* or with "plead_missing" (jet files before 2026-11-10) are
// refused.
//
// Pi0 candidate energy fraction (-z <lo>:<hi>, with or without -p): the same
// rows, from the x_F vs P_lead vs z yields (h_xFPz_*, h_abscosPz_*,
// h_xfsumPz_*), jets with z in [lo, hi) (bins of 0.1 by centre; the last bin
// holds z = 1, the one-point jets). z is the chosen pair's energy fraction when
// P_leadpair >= 0.5, else the leading point's (AnConfig.h kZCandPairMin). With
// -z, -p must use the P_lead class edges 0, 0.5, 0.7, 0.85, 1.
//
// Outputs: <out.root> with the graphs and TNamed "bad_runs" (-r), and
// <out.root minus .root>.txt, the table of every point (level, photon class,
// beam, x_F bin, <x_F>, A_N, stat error, jets, <|cos phi|>, <P>).
// plot_an.C draws the final A_N vs x_F from <out.root>.
//
// Jets above sqrt(s)/2. The yield files keep the jets with corrected E above
// sqrt(s)/2 apart (_hiE histograms, AnConfig.h); by default they are left
// out, as the FMS EM-jet analysis; -E adds them to every yield, |cos phi|
// and x_F sum. The choice is printed, in the table header and as TNamed
// "high_E_jets". x_F bins then reach 1.5 (kXfMax).
//
// Companion jets. Jets below every listed trigger threshold, in events where
// another jet passes one, are kept apart too (_comp histograms, AnConfig.h)
// and left out by default. -C forms the A_N of the companion jets alone (in
// place of the triggering jets), -A of both together. The choice is printed,
// in the table header and as TNamed "jet_sample". -E applies to either.
//
// Usage: ./AnAsymmetry [-E] [-C | -A] [-f factor | -fe f_EM -fj f_JP] [-J] [-V] [-m in|out:<mean>:<sigma>] <yields.root> <out.root> <P_blue> <P_yellow> [rebin = 2]
//        ./AnAsymmetry [-E] [-C | -A] [-f factor | -fe f_EM -fj f_JP] [-J] [-V] -r <yield file list> <run polarization table> <out.root> [rebin = 2]
//                         [bad-run list = kBadRunListDefault | none] [flagged | all]
//   rebin 2 gives x_F bins of 0.05, as the FMS EM-jet analysis
//   Prints, for each level, all and allJP, A_N per x_F bin of both beams
//   with all photon classes summed; writes graphs
//   gr_AN_<beam>_<level|all|allJP>_<photon class|all|np3plus|np12> (x_F > 0 blue, x_F < 0
//   yellow) and gr_AN_<level|all|allJP>_<class|all> with both beams.
// Build: g++ -m64 -O2 AnAsymmetry.cpp -o AnAsymmetry -I.. `root-config --cflags --libs` -std=c++11
//        (the report scripts copy AnConfig.h next to it, so need no -I..)

#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TDatime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
    return TString::Format("h_xF_%s_%s_%s_%s", kTrgLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]);
}

// Histograms summed over the input files, and per (EM level, photon class)
// the yields summed over sides and spin states weighted by each run's P.
// With hi_e, the jets above sqrt(s)/2 (the _hiE histograms) are added to
// the default ones; files without them (made before 2026-11-04, when those
// jets were dropped) are counted in n_no_hie. use_trig / use_comp choose the
// triggering jets (no suffix) and the companions (_comp); a file without the
// companion histograms (made before 2026-11-05) is an error when they are used.
struct Sums
{
    std::map<std::string, TH1D *> h;
    TH1D *pw[2][kNTrgLevels][kNPhotonClasses] = {};
    // JP levels, JP jets not firing EM at the factor (the -J union)
    TH1D *pwo[2][kNTrgLevels][kNPhotonClasses] = {};
    int n_old_labels = 0, n_no_hie = 0, nbins = 0;
    bool hi_e = false, use_trig = true, use_comp = false, use_veto = false;
    double factor_em = kEmCutFactor, factor_jp = kJpCutFactor;
    // -m: 0 off, 1 inside the mass window, 2 outside; np2 yields from h_xFm_*
    int mass_mode = 0;
    double mass_lo = 0, mass_hi = 0;
    // the mass bins actually in the window and their width (written as TNamed
    // "mass_window" for plot_an_mass.C)
    double mass_bin_lo = 0, mass_bin_hi = 0, mass_bin_w = 0;
    // -p: P_lead range [plead_lo, plead_hi); np1 rows = all jets in it, other classes empty
    bool plead_mode = false, z_mode = false;
    double plead_lo = 0, plead_hi = 1, z_lo = 0, z_hi = 1;
    static bool StartsWith(const std::string &a, const char *b) { return a.compare(0, std::strlen(b), b) == 0; }
    // -p: a yield (h_xF_<l>_<cls>_<side>_<spin>) or sum (h_abscos_ / h_xfsum_<l>_<cls>)
    // of the jets in the P_lead range for cls np1, zero for the other classes
    TH1D *ReadPlead(TFile &f, const std::string &name, bool jp_only)
    {
        TH1D *a = new TH1D((name + "_rdp").c_str(), "", kNXfBins, kXfMin, kXfMax);
        a->SetDirectory(0);
        const bool yield = StartsWith(name, "h_xF_");
        if (!yield) a->Sumw2();
        if (name.find("_np1") == std::string::npos) return a;
        TString n(name.c_str());
        const char *tag = z_mode ? "Pz_" : "P_";
        if (yield) { n.ReplaceAll("h_xF_", TString("h_xF") + tag); n.ReplaceAll("_np1_", "_"); }
        else { n.ReplaceAll("h_abscos_", TString("h_abscos") + tag); n.ReplaceAll("h_xfsum_", TString("h_xfsum") + tag); n.ReplaceAll("_np1", ""); }
        if (jp_only) n += "_jponly";
        if (z_mode)
        {
            TH3D *m3 = (TH3D *)f.Get(n);
            if (!m3 || m3->GetNbinsX() != a->GetNbinsX()) { delete a; return 0; }
            for (int by = 1; by <= m3->GetNbinsY(); by++)
            {
                const double c = m3->GetYaxis()->GetBinCenter(by);
                if (c < plead_lo || c >= plead_hi) continue;
                for (int bz = 1; bz <= m3->GetNbinsZ(); bz++)
                {
                    const double cz = m3->GetZaxis()->GetBinCenter(bz);
                    if (cz < z_lo || cz >= z_hi) continue;
                    for (int bx = 1; bx <= m3->GetNbinsX(); bx++) a->AddBinContent(bx, m3->GetBinContent(bx, by, bz));
                }
            }
            a->SetEntries(a->Integral());
            return a;
        }
        TH2D *m = (TH2D *)f.Get(n);
        if (!m) { delete a; return 0; }
        for (int by = 1; by <= m->GetNbinsY(); by++)
        {
            const double c = m->GetYaxis()->GetBinCenter(by);
            if (c < plead_lo || c >= plead_hi) continue;
            for (int bx = 1; bx <= m->GetNbinsX(); bx++) a->AddBinContent(bx, m->GetBinContent(bx, by));
        }
        a->SetEntries(a->Integral());
        return a;
    }
    bool InWindow(double m) const { return m >= mass_lo && m <= mass_hi; }
    // np2 yield of one (level, side, spin) from the x_F vs mass yields
    TH1D *ReadMass(TFile &f, const std::string &name, bool jp_only)
    {
        // name = h_xF_<level>_np2_<side>_<spin>
        TString n(name.c_str());
        n.ReplaceAll("h_xF_", "h_xFm_");
        n.ReplaceAll("_np2_", "_");
        if (jp_only) n += "_jponly";
        TH2D *m = (TH2D *)f.Get(n);
        if (!m) return 0;
        TH1D *a = new TH1D((name + "_rdm").c_str(), "", m->GetNbinsX(), m->GetXaxis()->GetXmin(), m->GetXaxis()->GetXmax());
        a->SetDirectory(0);
        a->Sumw2(false);
        for (int by = 1; by <= m->GetNbinsY(); by++)
        {
            const bool in = InWindow(m->GetYaxis()->GetBinCenter(by));
            if ((mass_mode == 1) != in) continue;
            for (int bx = 1; bx <= m->GetNbinsX(); bx++) a->AddBinContent(bx, m->GetBinContent(bx, by));
        }
        a->SetEntries(a->Integral());
        return a;
    }
    // -m: per x_F bin, the fraction of a file's two-photon jets of level l in
    // (out of) the mass window: window yield / all np2 yield, 0 where empty
    TH1D *WindowFraction(TFile &f, int l, bool jp_only, bool &have_hie)
    {
        int k2 = -1;
        for (int k = 0; k < kNPhotonClasses; k++)
            if (std::string(PhotonClassName(k)) == "np2") k2 = k;
        TH1D *win = 0, *all = 0;
        const int mm = mass_mode;
        for (int s = 0; s < kNSides; s++)
            for (int p = 0; p < kNSpinStates; p++)
            {
                mass_mode = mm;
                TH1D *a = Read(f, YieldName(l, k2, s, p).Data(), have_hie, jp_only);
                mass_mode = 0;
                TH1D *b = Read(f, YieldName(l, k2, s, p).Data(), have_hie, jp_only);
                mass_mode = mm;
                if (!win) { win = (TH1D *)a->Clone("wf_win"); win->SetDirectory(0); all = (TH1D *)b->Clone("wf_all"); all->SetDirectory(0); }
                else { win->Add(a); all->Add(b); }
                delete a; delete b;
            }
        for (int bx = 0; bx <= win->GetNbinsX() + 1; bx++)
        {
            const double na = all->GetBinContent(bx);
            win->SetBinContent(bx, na > 0 ? win->GetBinContent(bx) / na : 0);
            win->SetBinError(bx, 0);
        }
        delete all;
        return win;
    }
    // the suffixes summed for one histogram name: the triggering slices at or
    // above the factor (JP levels: with every EM tag, or with jp_only the EM
    // tags below the factor), companions, veto failures (EM levels only)
    std::vector<std::string> Suffixes(bool jp, bool jp_only) const
    {
        std::vector<std::string> v;
        if (use_trig)
            for (int g = 0; g < kNTrgSlices; g++)
            {
                if (kTrgSliceEdge[g] < (jp ? factor_jp : factor_em) - 1e-9) continue;
                if (!jp) { v.push_back(kTrgSliceSuffix[g]); continue; }
                for (int t = 0; t < kNEmTags; t++)
                    if (!jp_only || t == 0 || kTrgSliceEdge[t - 1] < factor_em - 1e-9)
                        v.push_back(std::string(kTrgSliceSuffix[g]) + kEmTagSuffix[t]);
            }
        if (use_comp) v.push_back(kCompSuffix);
        if (use_veto && !jp) v.push_back(kVetoSuffix);
        return v;
    }
    static bool IsJpName(const std::string &name) { return name.find("_JP") != std::string::npos; }
    TH1D *Read(TFile &f, const std::string &name, bool &have_hie, bool jp_only = false)
    {
        if (mass_mode && name.compare(0, 5, "h_xF_") == 0 && name.find("_np2_") != std::string::npos)
            return ReadMass(f, name, jp_only);
        if (plead_mode && (StartsWith(name, "h_xF_") || StartsWith(name, "h_abscos_") || StartsWith(name, "h_xfsum_")))
            return ReadPlead(f, name, jp_only);
        TH1D *a = 0;
        const std::vector<std::string> sfx = Suffixes(IsJpName(name), jp_only);
        for (size_t c = 0; c < sfx.size(); c++)
        {
            const std::string base = name + sfx[c];
            TH1D *b = (TH1D *)f.Get(base.c_str());
            if (!b) { delete a; return 0; }
            if (!a) { a = (TH1D *)b->Clone((name + "_rd").c_str()); a->SetDirectory(0); }
            else a->Add(b);
            if (hi_e)
            {
                TH1D *e = (TH1D *)f.Get((base + kHiESuffix).c_str());
                if (e) a->Add(e);
                else have_hie = false;
            }
        }
        return a;
    }
    // adds one file's histograms, its jets weighted by pol[beam]; false if any is missing
    bool Add(TFile &f, const double pol[2])
    {
        // files without "spin_convention" were filled with the old Spin_config
        // mapping (AnConfig.h), whose state p is the physical state 3 - p
        if (!f.Get("trigger_rule"))
        {
            std::cerr << "ERROR: " << f.GetName() << " has no trigger_rule: a yield file made before 2026-11-05 (pT-based EM levels, no JP); remake it" << std::endl;
            return false;
        }
        const bool old_labels = !f.Get("spin_convention");
        if (mass_mode && (old_labels || !f.Get("h_xFm_EM0_xpos_BuYu")))
        {
            std::cerr << "ERROR: " << f.GetName() << " has no h_xFm_* (or old spin labels); -m needs yield files from 2026-11-05 on" << std::endl;
            return false;
        }
        if (z_mode && (!f.Get("h_xFPz_EM0_xpos_BuYu") || f.Get("pz_missing")))
        {
            std::cerr << "ERROR: " << f.GetName() << " has no filled h_xFPz_* (jet files without reco_pi0_lead_z1, or an older AnHistograms); -z needs yield files from 2026-11-10 on" << std::endl;
            return false;
        }
        if (plead_mode && (old_labels || !f.Get("h_xFP_EM0_xpos_BuYu") || f.Get("plead_missing")))
        {
            std::cerr << "ERROR: " << f.GetName() << " has no filled h_xFP_* (jet files without reco_pi0_lead_prob, or an older AnHistograms); -p needs yield files from 2026-11-10 on" << std::endl;
            return false;
        }
        if (old_labels) n_old_labels++;
        bool have_hie = true;
        std::vector<std::string> names;
        for (int l = 0; l < kNTrgLevels; l++)
            for (int k = 0; k < kNPhotonClasses; k++)
            {
                names.push_back(TString::Format("h_abscos_%s_%s", kTrgLevelName[l], PhotonClassName(k)).Data());
                names.push_back(TString::Format("h_xfsum_%s_%s", kTrgLevelName[l], PhotonClassName(k)).Data());
                for (int s = 0; s < kNSides; s++)
                    for (int p = 0; p < kNSpinStates; p++) names.push_back(YieldName(l, k, s, p).Data());
            }
        std::vector<std::string> check;
        for (size_t i = 0; i < names.size(); i++)
        {
            const std::vector<std::string> sfx = Suffixes(IsJpName(names[i]), false);
            for (size_t c = 0; c < sfx.size(); c++) check.push_back(names[i] + sfx[c]);
        }
        for (size_t i = 0; i < check.size(); i++)
        {
            TH1D *hh = (TH1D *)f.Get(check[i].c_str());
            if (!hh)
            {
                std::cerr << "ERROR: " << f.GetName() << " has no " << check[i]
                          << std::endl;
                return false;
            }
            // every file's x_F binning must agree (files before 2026-11-04: 40 bins to 1; since: 60 to 1.5)
            if (nbins == 0) nbins = hh->GetNbinsX();
            if (hh->GetNbinsX() != nbins)
            {
                std::cerr << "ERROR: " << f.GetName() << " " << check[i] << " has " << hh->GetNbinsX() << " x_F bins, earlier files "
                          << nbins << " (yield files from before and after 2026-11-04 cannot be summed)" << std::endl;
                return false;
            }
            // yields must be raw counts, so their Poisson errors are sqrt(N)
            if (check[i].compare(0, 5, "h_xF_") == 0 && hh->GetSumw2N() > 0)
                for (int b = 1; b <= hh->GetNbinsX(); b++)
                    if (std::fabs(hh->GetSumw2()->At(b) - hh->GetBinContent(b)) > 1e-6 * (1 + hh->GetBinContent(b)))
                    {
                        std::cerr << "ERROR: " << f.GetName() << " " << check[i] << " is weighted; the statistical errors need raw counts" << std::endl;
                        return false;
                    }
        }
        for (size_t i = 0; i < names.size(); i++)
        {
            // the summed histogram this one adds to: for old files, the yield
            // of state p goes to physical state 3 - p
            std::string target = names[i];
            if (old_labels && target.compare(0, 5, "h_xF_") == 0)
                for (int p = 0; p < kNSpinStates; p++)
                {
                    const std::string from = std::string("_") + kSpinStateName[p];
                    if (target.size() > from.size() && target.compare(target.size() - from.size(), from.size(), from) == 0)
                    {
                        target = target.substr(0, target.size() - from.size()) + "_" + kSpinStateName[kNSpinStates - 1 - p];
                        break;
                    }
                }
            for (int o = 0; o < (IsJpName(names[i]) ? 2 : 1); o++)
            {
                TH1D *in = Read(f, names[i], have_hie, o == 1);
                // -m: h_abscos / h_xfsum of np2 hold every two-photon jet, while the
                // yields hold only those in (out of) the window; scale the sums by
                // this file's window fraction per x_F bin, so that dividing by the
                // window yield gives <|cos phi|> and <x_F> of all np2 jets
                if (mass_mode)
                    for (int l = 0; l < kNTrgLevels; l++)
                    {
                        const TString c = TString::Format("h_abscos_%s_np2", kTrgLevelName[l]), x = TString::Format("h_xfsum_%s_np2", kTrgLevelName[l]);
                        if (names[i] != c.Data() && names[i] != x.Data()) continue;
                        TH1D *wf = WindowFraction(f, l, o == 1, have_hie);
                        in->Multiply(wf);
                        delete wf;
                    }
                const std::string key = target + (o ? "#jponly" : "");
                if (h.count(key)) h[key]->Add(in);
                else { h[key] = (TH1D *)in->Clone((key + "_sum").c_str()); h[key]->SetDirectory(0); }
                delete in;
            }
        }
        for (int l = 0; l < kNTrgLevels; l++)
            for (int k = 0; k < kNPhotonClasses; k++)
            {
                for (int o = 0; o < (IsJpLevel(l) ? 2 : 1); o++)
                {
                    TH1D *n = Read(f, YieldName(l, k, 0, 0).Data(), have_hie, o == 1);
                    n->SetName("n_tmp");
                    n->Reset();
                    for (int s = 0; s < kNSides; s++)
                        for (int p = 0; p < kNSpinStates; p++)
                        {
                            TH1D *r = Read(f, YieldName(l, k, s, p).Data(), have_hie, o == 1);
                            n->Add(r);
                            delete r;
                        }
                    TH1D *(&w)[2][kNTrgLevels][kNPhotonClasses] = o ? pwo : pw;
                    for (int b = 0; b < 2; b++)
                    {
                        if (!w[b][l][k]) { w[b][l][k] = (TH1D *)n->Clone(TString::Format("pw%s_%d_%d_%d", o ? "o" : "", b, l, k)); w[b][l][k]->SetDirectory(0); w[b][l][k]->Reset(); }
                        w[b][l][k]->Add(n, pol[b]);
                    }
                    delete n;
                }
            }
        if (hi_e && !have_hie) n_no_hie++;
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
    // -E (anywhere among the arguments): add the jets above sqrt(s)/2;
    // -C: companion jets alone; -A: triggering and companion jets
    // -f <factor>: trigger-value cut factor, both families (a slice edge);
    // -fe / -fj: EM / JP alone; -J: JP levels in
    // "all"; -V: add the EM jets failing only the veto margin
    bool hi_e = false, comp_only = false, comp_add = false, jp_in_all = false, add_veto = false;
    double factor_em = kEmCutFactor, factor_jp = kJpCutFactor;
    std::string mass_opt, plead_opt, z_opt;
    {
        int m = 1;
        for (int i = 1; i < argc; i++)
        {
            const std::string a = argv[i];
            if (a == "-E") { hi_e = true; continue; }
            if (a == "-C") { comp_only = true; continue; }
            if (a == "-A") { comp_add = true; continue; }
            if (a == "-J") { jp_in_all = true; continue; }
            if (a == "-V") { add_veto = true; continue; }
            if (a == "-f" && i + 1 < argc) { factor_em = factor_jp = std::atof(argv[++i]); continue; }
            if (a == "-fe" && i + 1 < argc) { factor_em = std::atof(argv[++i]); continue; }
            if (a == "-fj" && i + 1 < argc) { factor_jp = std::atof(argv[++i]); continue; }
            if (a == "-m" && i + 1 < argc) { mass_opt = argv[++i]; continue; }
            if (a == "-p" && i + 1 < argc) { plead_opt = argv[++i]; continue; }
            if (a == "-z" && i + 1 < argc) { z_opt = argv[++i]; continue; }
            argv[m++] = argv[i];
        }
        argc = m;
    }
    {
        bool edge_em = false, edge_jp = false;
        for (int g = 0; g < kNTrgSlices; g++)
        {
            edge_em |= std::fabs(kTrgSliceEdge[g] - factor_em) < 1e-6;
            edge_jp |= std::fabs(kTrgSliceEdge[g] - factor_jp) < 1e-6;
        }
        if (!edge_em || !edge_jp) { std::cerr << "ERROR: -f / -fe / -fj must be a slice edge (1.15, 1.3, 1.45, 1.6)" << std::endl; return 1; }
    }
    if (comp_only && comp_add) { std::cerr << "ERROR: -C and -A exclude each other" << std::endl; return 1; }
    const bool per_run = argc > 1 && std::string(argv[1]) == "-r";
    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0] << " [-E] [-C | -A] [-f factor | -fe f_EM -fj f_JP] [-J] [-V] <yields.root> <out.root> <P_blue> <P_yellow> [rebin = 2]\n"
                  << "       " << argv[0] << " [-E] [-C | -A] [-f factor | -fe f_EM -fj f_JP] [-J] [-V] -r <yield file list> <run polarization table> <out.root> [rebin = 2] [bad-run list | none] [flagged | all]\n"
                  << "  -E: include the jets with corrected E above sqrt(s)/2 (default: left out, as the FMS EM-jet analysis)\n"
                  << "  -C: companion jets alone (below the trigger threshold, in events with a jet above it)\n"
                  << "  -A: triggering and companion jets together (default: triggering jets alone)\n"
                  << "  -f: trigger value cut, factor x nominal threshold, both families: 1.15, 1.3, 1.45 or 1.6\n"
                  << "  -fe / -fj: the factor for EM levels / JP levels alone (default 1.3 for both)\n"
                  << "  -J: \"all\" = EM + JP union, each jet once (default: EM levels only; allJP is always written)\n"
                  << "  -V: add EM jets that failed only the HCal-veto margin\n"
                  << "  -m in:<mean>:<sigma> | out:<mean>:<sigma>: np2 rows = two-photon jets with pair mass inside / outside mean +- 2 sigma (GeV;\n"
                  << "     in:<peak> / out:<peak>: 0.5-1.5 x <peak>)\n"
                  << "  -p <lo>:<hi>: np1 rows (and all, np12) = every jet with leading-point pi0 probability P_lead in [lo, hi)\n"
                  << "  -z <lo>:<hi>: the same rows, jets whose pi0 candidate carries a jet energy fraction z in [lo, hi) (with -p: class edges 0, 0.5, 0.7, 0.85, 1)" << std::endl;
        return 1;
    }
    const int rebin = argc > 5 ? std::atoi(argv[5]) : 2;
    if (rebin < 1) { std::cerr << "ERROR: bad rebin" << std::endl; return 1; }

    Sums sums;
    sums.hi_e = hi_e;
    sums.use_trig = !comp_only;
    sums.use_comp = comp_only || comp_add;
    sums.use_veto = add_veto;
    sums.factor_em = factor_em;
    sums.factor_jp = factor_jp;
    std::string mass_note;
    if (!mass_opt.empty())
    {
        const size_t c = mass_opt.find(':');
        const std::string w = mass_opt.substr(0, c);
        const size_t c2 = c == std::string::npos ? std::string::npos : mass_opt.find(':', c + 1);
        const double peak = c == std::string::npos ? 0 : std::atof(mass_opt.substr(c + 1).c_str());
        const double sig = c2 == std::string::npos ? 0 : std::atof(mass_opt.substr(c2 + 1).c_str());
        if ((w != "in" && w != "out") || !(peak > 0) || (c2 != std::string::npos && !(sig > 0)))
        {
            std::cerr << "ERROR: -m takes in:<mean>:<sigma> or out:<mean>:<sigma> (GeV; or in:<peak> / out:<peak>)" << std::endl;
            return 1;
        }
        // the window: mean +- kMassWindowSigma sigma, or 0.5-1.5 x peak
        const double kMassWindowSigma = 2.0;
        const double wlo = sig > 0 ? peak - kMassWindowSigma * sig : 0.5 * peak;
        const double whi = sig > 0 ? peak + kMassWindowSigma * sig : 1.5 * peak;
        const std::string wtext = sig > 0 ? TString::Format("%.4f +- 2 x %.4f GeV", peak, sig).Data()
                                          : TString::Format("0.5-1.5 x %.4f GeV", peak).Data();
        if (hi_e || comp_only || comp_add || add_veto || std::fabs(factor_em - 1.3) > 1e-6 || std::fabs(factor_jp - 1.3) > 1e-6)
        {
            std::cerr << "ERROR: -m works with the default jets only (no -E, -C, -A, -V; factors 1.3)" << std::endl;
            return 1;
        }
        sums.mass_mode = w == "in" ? 1 : 2;
        // the window in whole mass bins (a bin is in if its centre is), on the
        // mass axis of the yield files (25 MeV before 2026-11-11, AnConfig.h
        // kNMassBins after), read from the first file
        int nmb = kNMassBins;
        double mmin = kMassMin, mmax = kMassMax;
        {
            std::string first = per_run ? "" : argv[1];
            if (per_run)
            {
                std::ifstream lst(argv[2]);
                std::string line;
                while (first.empty() && std::getline(lst, line))
                    if (!line.empty() && line[0] != '#') first = line;
            }
            TFile *f0 = first.empty() ? 0 : TFile::Open(first.c_str());
            TH2 *h0 = f0 && !f0->IsZombie() ? (TH2 *)f0->Get("h_xFm_EM0_xpos_BuYu") : 0;
            if (h0) { nmb = h0->GetNbinsY(); mmin = h0->GetYaxis()->GetXmin(); mmax = h0->GetYaxis()->GetXmax(); }
            delete f0;
        }
        const double bw = (mmax - mmin) / nmb;
        double lo = 1e9, hi = -1e9;
        for (int b = 0; b < nmb; b++)
        {
            const double cen = mmin + (b + 0.5) * bw;
            if (cen >= wlo && cen <= whi) { lo = std::min(lo, cen - 0.5 * bw); hi = std::max(hi, cen + 0.5 * bw); }
        }
        sums.mass_bin_lo = lo; sums.mass_bin_hi = hi; sums.mass_bin_w = bw;
        sums.mass_lo = wlo; sums.mass_hi = whi;
        mass_note = TString::Format("np2 = two-photon jets with M_gg %s %s (bins %.3f-%.3f GeV); <|cos|>, <x_F>, P of all np2", w == "in" ? "inside" : "outside", wtext.c_str(), lo, hi).Data();
        std::printf("AnAsymmetry: %s\n", mass_note.c_str());
    }
    if (!z_opt.empty() && plead_opt.empty()) plead_opt = "0:1";
    if (!plead_opt.empty())
    {
        const size_t c = plead_opt.find(':');
        const double lo = std::atof(plead_opt.substr(0, c).c_str()), hi = c == std::string::npos ? -1 : std::atof(plead_opt.substr(c + 1).c_str());
        if (c == std::string::npos || !(lo >= 0 && hi > lo && hi <= 1.0 + 1e-9)) { std::cerr << "ERROR: -p takes <lo>:<hi> with 0 <= lo < hi <= 1" << std::endl; return 1; }
        if (!mass_opt.empty()) { std::cerr << "ERROR: -p and -m exclude each other" << std::endl; return 1; }
        if (hi_e || comp_only || comp_add || add_veto || std::fabs(factor_em - 1.3) > 1e-6 || std::fabs(factor_jp - 1.3) > 1e-6)
        {
            std::cerr << "ERROR: -p works with the default jets only (no -E, -C, -A, -V; factors 1.3)" << std::endl;
            return 1;
        }
        sums.plead_mode = true;
        sums.plead_lo = lo; sums.plead_hi = hi + (hi >= 1.0 - 1e-9 ? 1e-6 : 0);
        mass_note = TString::Format("np1 rows (all, np12) = every jet with P_lead in [%.2f, %.2f%s (0.05 bins by centre); other photon classes empty", lo, hi, hi >= 1.0 - 1e-9 ? "]" : ")").Data();
        if (!z_opt.empty())
        {
            const size_t cz = z_opt.find(':');
            const double zl = std::atof(z_opt.substr(0, cz).c_str()), zh = cz == std::string::npos ? -1 : std::atof(z_opt.substr(cz + 1).c_str());
            if (cz == std::string::npos || !(zl >= 0 && zh > zl && zh <= 1.0 + 1e-9)) { std::cerr << "ERROR: -z takes <lo>:<hi> with 0 <= lo < hi <= 1" << std::endl; return 1; }
            bool elo = false, ehi = false;
            for (int k = 0; k <= kNPleadClsZ; k++) { elo |= std::fabs(kPleadEdgesZ[k] - lo) < 1e-6; ehi |= std::fabs(std::min(kPleadEdgesZ[k], 1.0) - hi) < 1e-6; }
            if (!elo || !ehi) { std::cerr << "ERROR: with -z, -p must use the P_lead class edges 0, 0.5, 0.7, 0.85, 1" << std::endl; return 1; }
            sums.z_mode = true;
            sums.z_lo = zl; sums.z_hi = zh + (zh >= 1.0 - 1e-9 ? 1e-3 : 0);
            mass_note += TString::Format("; pi0 candidate z in [%.2f, %.2f%s", zl, zh, zh >= 1.0 - 1e-9 ? "] (z = 1: one-point jets)" : ")").Data();
            std::printf("AnAsymmetry: pi0 candidate z in [%.2f, %.2f%s\n", zl, zh, zh >= 1.0 - 1e-9 ? "]" : ")");
        }
        std::printf("AnAsymmetry: %s\n", mass_note.c_str());
    }
    const std::string trg_note = TString::Format("trigger value >= %.2f (EM) / %.2f (JP) x threshold; levels judged per family (EM, JP); all = %s%s", factor_em, factor_jp,
                                                 jp_in_all ? "EM levels + JP jets not firing EM (-J union)" : "EM levels",
                                                 add_veto ? "; EM veto failures added (-V)" : "").Data() + (mass_note.empty() ? std::string() : "; " + mass_note);
    std::printf("AnAsymmetry: %s\n", trg_note.c_str());
    const std::string hie_note = hi_e ? "jets above sqrt(s)/2 included (-E)" : "jets above sqrt(s)/2 left out (default, as the FMS)";
    const std::string sample_note = comp_only ? "companion jets alone (-C)" : comp_add ? "triggering and companion jets (-A)" : "triggering jets (default)";
    std::printf("AnAsymmetry: %s; %s\n", hie_note.c_str(), sample_note.c_str());
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
        bad_note += TString::Format(" + %d runs of AnConfig.h kSuppBadRuns", kNSuppBadRuns).Data();
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
            if (const char *why = SuppBadRunReason(run))
            {
                bad_removed.push_back(TString::Format("%d %s (AnConfig.h)", run, why).Data());
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

    if (sums.nbins % rebin) { std::cerr << "ERROR: rebin " << rebin << " does not divide the " << sums.nbins << " x_F bins" << std::endl; return 1; }
    if (sums.n_no_hie > 0)
        std::printf("AnAsymmetry: WARNING: %d yield file(s) have no _hiE histograms (made before 2026-11-04, jets above sqrt(s)/2 dropped there)\n", sums.n_no_hie);
    if (sums.n_old_labels > 0)
        std::printf("AnAsymmetry: %d yield file(s) without spin_convention (old Spin_config mapping): spin states relabelled p -> 3 - p\n", sums.n_old_labels);
    TFile out(argv[per_run ? 4 : 2], "RECREATE");
    if (per_run) TNamed("bad_runs", bad_note.c_str()).Write();
    TNamed("high_E_jets", hie_note.c_str()).Write();
    TNamed("jet_sample", sample_note.c_str()).Write();
    TNamed("trigger_selection", trg_note.c_str()).Write();
    if (sums.mass_mode) TNamed("mass_window", TString::Format("%.4f %.4f %.4f", sums.mass_bin_lo, sums.mass_bin_hi, sums.mass_bin_w).Data()).Write();
    const std::string table_name = std::string(argv[per_run ? 4 : 2]).substr(0, std::string(argv[per_run ? 4 : 2]).rfind(".root")) + ".txt";
    FILE *table = std::fopen(table_name.c_str(), "w");
    std::fprintf(table, "# A_N = eps / (<P> <|cos phi|>), %s; stat = statistical error only\n", pol_note.c_str());
    std::fprintf(table, "# %s; %s; %s\n", hie_note.c_str(), sample_note.c_str(), trg_note.c_str());
    if (per_run) std::fprintf(table, "# bad runs: %s\n", bad_note.c_str());
    std::fprintf(table, "# level class beam xF_lo xF_hi <xF> A_N stat jets <|cos|> <P>\n");
    for (int L = 0; L <= kNTrgLevels + 1; L++)
    {
        // levels summed: one level, all (L = kNTrgLevels: EM, with -J also JP),
        // allJP (L = kNTrgLevels + 1)
        const int l0 = L < kNTrgLevels ? L : L == kNTrgLevels ? 0 : kNEmLevels;
        const int l1 = L < kNTrgLevels ? L + 1 : L == kNTrgLevels ? (jp_in_all ? kNTrgLevels : kNEmLevels) : kNTrgLevels;
        const char *lname = L < kNTrgLevels ? kTrgLevelName[L] : L == kNTrgLevels ? "all" : "allJP";
        const int l = L;
        // [photon class, kNPhotonClasses = all, kNPhotonClasses + 1 = np3plus (3 or
        // more points: np3 + np4 + np5plus), kNPhotonClasses + 2 = np12 (1 or 2
        // points: np1 + np2)]: yields [side][spin state], and per bin the jets
        // summed over sides and spin states, sum of |cos phi|, sum of x_F, sum of P x jets per beam
        const int nK = kNPhotonClasses + 3;
        TH1D *y[nK][kNSides][kNSpinStates], *n[nK], *mcos[nK], *mxf[nK], *mpol[2][nK];
        for (int k = 0; k < nK; k++)
        {
            const int k0 = k < kNPhotonClasses ? k : (k == kNPhotonClasses ? 0 : k == kNPhotonClasses + 1 ? 3 : 1);
            const int k1 = k < kNPhotonClasses ? k + 1 : (k == kNPhotonClasses + 2 ? 3 : kNPhotonClasses);
            const TString tag = TString::Format("%d_%d", l, k);
            bool first = true;
            for (int ll = l0; ll < l1; ll++)
            for (int kk = k0; kk < k1; kk++, first = false)
            {
                // the -J union ("all" with JP levels): JP jets that did not fire EM at the factor
                const bool jpo = L == kNTrgLevels && IsJpLevel(ll);
                const std::string ok_ = jpo ? "#jponly" : "";
                TH1D *c = (TH1D *)sums.h[TString::Format("h_abscos_%s_%s", kTrgLevelName[ll], PhotonClassName(kk)).Data() + ok_]->Rebin(rebin, "tmp_c");
                TH1D *x = (TH1D *)sums.h[TString::Format("h_xfsum_%s_%s", kTrgLevelName[ll], PhotonClassName(kk)).Data() + ok_]->Rebin(rebin, "tmp_x");
                if (first) { mcos[k] = (TH1D *)c->Clone("c_" + tag); mxf[k] = (TH1D *)x->Clone("x_" + tag); }
                else { mcos[k]->Add(c); mxf[k]->Add(x); }
                delete c; delete x;
                for (int b = 0; b < 2; b++)
                {
                    TH1D *pr = (TH1D *)(jpo ? sums.pwo : sums.pw)[b][ll][kk]->Rebin(rebin, "tmp_p");
                    if (first) mpol[b][k] = (TH1D *)pr->Clone(TString::Format("p%d_", b) + tag);
                    else mpol[b][k]->Add(pr);
                    delete pr;
                }
                for (int s = 0; s < kNSides; s++)
                    for (int p = 0; p < kNSpinStates; p++)
                    {
                        TH1D *r = (TH1D *)sums.h[YieldName(ll, kk, s, p).Data() + ok_]->Rebin(rebin, "tmp_y");
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

        for (int k = 0; k < nK; k++)
        {
            const char *cname = k < kNPhotonClasses ? PhotonClassName(k) : (k == kNPhotonClasses ? "all" : k == kNPhotonClasses + 1 ? "np3plus" : "np12");
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
