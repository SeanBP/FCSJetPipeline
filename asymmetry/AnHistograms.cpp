// A_N yield histograms from calibrated jet files (after jet finding, the
// underlying-event step and ApplyCorrections_em_only), so the jet and
// SimpleTree files need not be kept.
//
// Per event: no FCS LED (fcs_led), at least one EM or JP trigger listed
// (AnConfig.h) and a spin state (Spin_config); events failing any are
// dropped. In data the listed triggers are the ones that recorded the event
// (MuDst nominal trigger ids: condition fired and the prescale accepted it).
// Per jet: a valid calibration (reco_E_corr > 0), x_F > 0, corrected pT =
// reco_E_corr / cosh(eta) above kJetPtCorrMin, and a trigger level the jet
// fired (AnConfig.h "Trigger levels"): its rebuilt trigger value (reco_trg_*)
// above a factor times the nominal threshold of a listed level, with the EM
// HCal-veto margin and, for JP0/JP1, the jet in a listed region. EM and JP
// are judged independently: a jet firing both fills its EM level and its JP
// level (the lowest listed level within each family), so each level's A_N
// has every jet that fired it; the JP histograms carry the same jet's EM
// result (_noem, _emt115 ... _emt160 after the slice suffix) for an EM + JP
// union without double counting (AnAsymmetry -J). The yields are kept per slice of
// trigger value / threshold, the factor's choice left to AnAsymmetry -f
// (default kEmCutFactor = kJpCutFactor = 1.3; slices from 1.15). Each jet fills x_F
// (reco_x_F_corr) in
//   h_xF_<level>_<photon class>_<side>_<spin state><sample>
// with <level> EM0, EM1, EM2, JP0, JP1, JP2, e.g. h_xF_EM1_np2_xpos_BuYd_t130,
// and with weight |cos phi| (phi = atan2(reco_y, reco_x), the jet's azimuth;
// the same for both beams, whose left is +x or -x)
//   h_abscos_<level>_<photon class><sample>
// summed over sides and spin states: divided by the summed yields it is the
// bin's mean |cos phi|, which a left-right asymmetry is diluted by
// (AnAsymmetry.cpp); likewise, with weight x_F,
//   h_xfsum_<level>_<photon class><sample>
// for the bin's mean corrected x_F, where its A_N point is placed. <sample>:
//   _t115, _t130, _t145, _t160  triggering jets with trigger value / threshold
//                               in [1.15, 1.3), [1.3, 1.45), [1.45, 1.6), >= 1.6
//   _veto                       EM jets above 1.15 x an EM threshold failing only
//                               the HCal-veto margin (EM levels only)
//   _comp                       companion jets: every jet cut passed but no
//                               trigger level at 1.15 in either family (and no
//                               veto failure), in an event with a jet at >=
//                               the default factor x threshold; at the highest EM
//                               level among those jets, else the highest JP
// each followed by _hiE for jets with corrected E above sqrt(s)/2 (the FMS
// EM-jet cut; AnAsymmetry -E adds them). Two-photon jets (two point
// constituents), E <= sqrt(s)/2 only, also fill
//   h_xFm_<level>_<side>_<spin state>[_t115 | _comp | _jponly]
// x_F vs the pair's invariant mass (no suffix: triggering jets at >= the
// default factor x threshold (1.3); _t115: jets in 1.15-1.3;
// _jponly, JP levels: jets at >= 1.3 in JP but not in EM, for an
// EM + JP union), for a pi0 mass window chosen later. Every jet with
// E <= sqrt(s)/2 and a leading-point pi0 probability (reco_pi0_lead_prob >= 0;
// jet files from 2026-11-10 on) likewise fills
//   h_xFP_<level>_<side>_<spin state>[_t115 | _comp | _jponly]
// x_F vs P_lead, and h_abscosP_<level>[...] / h_xfsumP_<level>[...] (x_F vs
// P_lead weighted with |cos phi| / x_F, summed over sides and spin states),
// for P_lead ranges chosen later (AnAsymmetry -p). Without the branch they stay
// empty and TNamed "plead_missing" is written. The Collins yields (h_xFg_*,
// h_absg_*, h_xfsumg_*; TNamed "collins_missing" without their branches) are
// filled in the same loop by collins/CollinsHistograms.h; see collins/README.md.
// A_N vs the pi0 candidate's energy fraction: every jet with P_lead and z
// (AnConfig.h kZCandPairMin), E <= sqrt(s)/2, default (no suffix) or _jponly, fills
//   h_xFPz_<level>_<side>_<spin state>[_jponly]
// (TH3D: x_F as h_xF_*, P_lead in 0-0.5-0.7-0.85-1, z in 10 bins to 1), and
// h_abscosPz_<level>[_jponly] / h_xfsumPz_<level>[_jponly] (the same axes,
// weighted with |cos phi| / x_F, summed over sides and spin states), for
// AnAsymmetry -p ... -z.
// Also written:
//   h_cutflow                    jets (and events) at each cut
//   h_events_<level>             events per spin state, by the highest listed
//                                EM level, or JP level if no EM is listed
//   h_ptcorr_<level>, h_ptraw_<level>  corrected / uncorrected pT before the
//                                trigger cut, by the same event level
//   h_trgratio_<level>           trigger value / threshold of the jets given
//                                that level (triggering and _veto), before
//                                the slice cut
//   runs                         TNamed, the run numbers of the input
//                                (event_run, else from the file names,
//                                jet_output_<run>_<file>.root)
//   spin_convention              TNamed, the Spin_config -> spin state
//                                mapping used (AnConfig.h)
// and sanity distributions, qa_<variable>_<stage>, at two stages: "evt",
// every jet of the events passing the event cuts (no LED, an EM or JP
// trigger, a spin state), and "sel", the jets entering the default yields
// (>= the default factor x threshold), those above sqrt(s)/2
// included: eta, phi, energy (uncorrected reco_E, UE reco_ue_E, corrected),
// pT (uncorrected, corrected), corrected x_F, points in the jet, jet position
// (reco_x, reco_y), and per event the jets and selected jets, the vertex z
// (event_vertex_z, -9999 = none, and its source), the listed EM and JP triggers and
// sqrt(s) (qa_sqrt_s). qa_xy_sel_<side> is split by side for the left/right
// acceptance. The same jet distributions per level,
// qa_<variable>_<stage>_<level>: at "evt" by the event's highest listed
// level, at "sel" by the jet's level (the one its yields go to); jets per
// event by the event's highest listed level at both stages. A third stage,
// "comp", has the companion jets (by their level). Each stage also has, for
// its two-photon jets, the pair's invariant mass (qa_mgg_*), its orientation
// on the ECal face (qa_pairang_*: angle of the line joining the two points,
// [0, 180) degrees, 0 = side by side in x, 90 = one above the other), mass vs
// corrected jet E (qa_mggE_*) and mass vs orientation (qa_mggang_*). The pair
// uses the constituents' momenta (points, energy along the direction from the
// nominal vertex at the origin), so the mass carries the point energy scale
// and no vertex correction.
// Output files of one run, or of a fill, add with hadd.
//
// Usage: ./AnHistograms <jet file list> <output.root>
// Build: g++ -m64 -O2 AnHistograms.cpp VectorDict.cxx -o AnHistograms -Icollins
//            `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//        (the jobs copy collins/*.h next to it, so need no -Icollins)

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "VectorDict.h"
#include "AnConfig.h"
#include "CollinsHistograms.h"

// The run number in a jet_output_<run>_<file>.root name, 0 if none
static int RunFromFileName(const std::string &path)
{
    const std::string base = path.substr(path.find_last_of('/') + 1);
    int run = 0, file_number = 0;
    return std::sscanf(base.c_str(), "jet_output_%d_%d.root", &run, &file_number) == 2 ? run : 0;
}

// kTrg115: jets firing a trigger level (EM or JP) at >= 1.15 x threshold;
// kTrgDefault: at >= the default factor x threshold (the default yields);
// kBelowHalfS: those with E <= sqrt(s)/2. Not steps of the flow: kEmDefault,
// kJpDefault, kBothDefault (jets at the default factor in the EM family, the
// JP family, both), kVeto (_veto jets) and kCompanion (_comp jets).
enum Cut { kAllJets, kNoLed, kTrgListed, kSpin, kCalibrated, kXfPositive, kPtMin, kTrg115, kTrgDefault, kBelowHalfS,
           kEmDefault, kJpDefault, kBothDefault, kVeto, kCompanion, kNCuts };
static const char *const kCutName[kNCuts] = {"all jets", "no LED", "EM/JP trigger", "spin state", "calibrated", "x_F > 0",
                                             "pT > 2", "trg >= 1.15 thr", "trg >= default", "E <= sqrt(s)/2",
                                             "EM >= default", "JP >= default", "EM and JP dflt", "EM veto failed", "companion"};

// histogram samples: the triggering slices, companions, EM veto failures
static const int kNSamp = kNTrgSlices + 2, kSampComp = kNTrgSlices, kSampVeto = kNTrgSlices + 1;
static const char *SampSuffix(int g) { return g < kNTrgSlices ? kTrgSliceSuffix[g] : g == kSampComp ? kCompSuffix : kVetoSuffix; }
// the histograms of level l, sample g, EM tag t exist for: EM levels t = 0
// (no tag), every sample; JP levels, triggering slices with every EM tag
// (kEmTagSuffix after the slice suffix), companions t = 0 (no tag), no _veto
static bool Booked(int l, int g, int t)
{
    if (!IsJpLevel(l)) return t == 0;
    if (g < kNTrgSlices) return true;
    return g == kSampComp && t == 0;
}
static TString Suffix(int l, int g, int t)
{
    return TString(SampSuffix(g)) + (IsJpLevel(l) && g < kNTrgSlices ? kEmTagSuffix[t] : "");
}

// JP region masks (reco_trg_jp_mask bit ns*5 + patch): A = patch 0, BC = 1, 2, DE = 3, 4
static const int kJpMaskA = (1 << 0) | (1 << 5), kJpMaskBC = (1 << 1) | (1 << 2) | (1 << 6) | (1 << 7),
                 kJpMaskDE = (1 << 3) | (1 << 4) | (1 << 8) | (1 << 9);

// Invariant mass (GeV) and orientation on the ECal face ([0, 180) degrees)
// of a jet's two constituents a and b; mass -1 if not computable
static void PairOf(double ax, double ay, double az, double ae, double bx, double by, double bz, double be, double &mass, double &angle)
{
    mass = -1; angle = -1;
    const double ra = std::sqrt(ax * ax + ay * ay + az * az), rb = std::sqrt(bx * bx + by * by + bz * bz);
    if (!(ra > 0 && rb > 0 && az > 0 && bz > 0)) return;
    double c = (ax * bx + ay * by + az * bz) / (ra * rb);
    if (c > 1) c = 1;
    const double m2 = 2.0 * ae * be * (1.0 - c);
    mass = m2 > 0 ? std::sqrt(m2) : 0;
    // positions where the directions cross the plane z = z_proj
    const double dx = z_proj * (bx / bz - ax / az), dy = z_proj * (by / bz - ay / az);
    angle = std::atan2(dy, dx) * 180.0 / M_PI;
    if (angle < 0) angle += 180.0;
    if (angle >= 180.0) angle -= 180.0;
}

// The jet sanity distributions of one stage (and EM level)
struct QaSet
{
    TH1D *eta, *phi, *e, *ue, *ecorr, *pt, *ptcorr, *xf, *np, *njet, *mgg, *pang;
    TH2D *xy, *mggE, *mggang;
    void Book(const TString &tag)
    {
        eta = new TH1D("qa_eta_" + tag, ";jet #eta;jets", 60, 2.0, 5.0);
        phi = new TH1D("qa_phi_" + tag, ";jet #phi (reco_phi, 0 to 2#pi);jets", 72, 0, 2 * M_PI);
        e = new TH1D("qa_E_" + tag, ";uncorrected jet E (GeV);jets", 160, 0, 400);
        ue = new TH1D("qa_ueE_" + tag, ";UE energy in the jet (GeV);jets", 100, 0, 10);
        ecorr = new TH1D("qa_Ecorr_" + tag, ";corrected jet E (GeV);jets", 160, 0, 400);
        pt = new TH1D("qa_pt_" + tag, ";uncorrected jet p_{T} (GeV);jets", 100, 0, 20);
        ptcorr = new TH1D("qa_ptcorr_" + tag, ";corrected jet p_{T} (GeV);jets", 100, 0, 20);
        xf = new TH1D("qa_xF_" + tag, ";corrected x_{F};jets", kNXfBins, kXfMin, kXfMax);
        np = new TH1D("qa_np_" + tag, ";points in the jet;jets", 20, 0, 20);
        njet = new TH1D("qa_njet_" + tag, ";jets per event;events", 10, 0, 10);
        xy = new TH2D("qa_xy_" + tag, ";jet x (cm);jet y (cm)", 60, -150, 150, 44, -110, 110);
        mgg = new TH1D("qa_mgg_" + tag, ";two-photon jets: M_{#gamma#gamma} (GeV);jets", kNMassBins, kMassMin, kMassMax);
        pang = new TH1D("qa_pairang_" + tag, ";two-photon jets: pair orientation (deg, 0 = along x, 90 = along y);jets", 36, 0, 180);
        mggE = new TH2D("qa_mggE_" + tag, ";corrected jet E (GeV);M_{#gamma#gamma} (GeV)", 50, 0, 250, kNMassBins, kMassMin, kMassMax);
        mggang = new TH2D("qa_mggang_" + tag, ";pair orientation (deg);M_{#gamma#gamma} (GeV)", 18, 0, 180, kNMassBins, kMassMin, kMassMax);
    }
    void FillPair(double mass, double angle, double ecorr_)
    {
        if (mass < 0) return;
        mgg->Fill(mass); pang->Fill(angle); mggE->Fill(ecorr_, mass); mggang->Fill(angle, mass);
    }
    void Fill(double eta_, double phi_, double e_, double ue_, double ecorr_, double pt_, double ptcorr_, double xf_, int np_, double x_, double y_)
    {
        eta->Fill(eta_); phi->Fill(phi_); e->Fill(e_); ue->Fill(ue_); ecorr->Fill(ecorr_);
        pt->Fill(pt_); ptcorr->Fill(ptcorr_); xf->Fill(xf_); np->Fill(np_); xy->Fill(x_, y_);
    }
};

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <jet file list> <output.root>" << std::endl;
        return 1;
    }
    TChain t("jetTree");
    {
        std::ifstream in(argv[1]);
        std::string f;
        while (in >> f) t.Add(f.c_str());
    }
    if (t.GetNtrees() == 0) { std::cerr << "ERROR: no files in " << argv[1] << std::endl; return 1; }
    if (!t.GetBranch("reco_E_corr")) { std::cerr << "ERROR: no reco_E_corr; run ApplyCorrections_em_only first" << std::endl; return 1; }
    if (!t.GetBranch("reco_trg_em_et")) { std::cerr << "ERROR: no reco_trg_em_et; jet files made before 2026-11-05 lack the trigger values (remake them)" << std::endl; return 1; }

    Int_t spin_config = -1, led = 0, em[kNEmLevels] = {0}, em_tpc[kNEmLevels] = {0};
    // JP flags: JP2; JPA1, JPBC1, JPDE1; JPA0, JPBC0, JPDE0
    Int_t jp2 = 0, jp1[3] = {0}, jp0[3] = {0};
    Double_t sqrt_s = 0;
    Float_t vz = -9999;
    Int_t vz_source = 0;
    Int_t event_run = 0;  // jet files from 2026-11-05 on (merged files need it)
    std::vector<float> *eta = 0, *phi = 0, *e_raw = 0, *e_ue = 0, *x = 0, *y = 0, *pt_raw = 0, *e_corr = 0, *xf_corr = 0;
    std::vector<float> *trg_em = 0, *trg_hmax = 0, *trg_jp = 0;
    std::vector<int> *trg_jpmask = 0;
    // Photons in a jet = its point constituents (reco_nconst), the FMS EM-jet
    // definition (BrightSTAR: jet->GetNumberOfTowers(), FMS points). Needs
    // point jets (reco_from_points = 1); reco_nphotons is not read, since in
    // files made before 2026-11-05 it is the deltaR < R point count (differs
    // on ~2.5% of jets).
    std::vector<int> *nconst = 0, *c_jet = 0;
    Int_t from_points = -1;
    std::vector<float> *c_px = 0, *c_py = 0, *c_pz = 0, *c_e = 0;
    std::vector<float> *plead = 0, *gup = 0, *zpair = 0, *pleadpair = 0, *z1 = 0;
    t.SetBranchAddress("Spin_config", &spin_config);
    t.SetBranchAddress("fcs_led", &led);
    t.SetBranchAddress("event_sqrt_s_gev", &sqrt_s);
    t.SetBranchAddress("fcsEM0", &em[0]);
    t.SetBranchAddress("fcsEM1", &em[1]);
    t.SetBranchAddress("fcsEM2", &em[2]);
    // the _tpc copies: the same EM conditions under their own trigger ids
    for (int l = 0; l < kNEmLevels; l++)
    {
        const TString b = TString::Format("fcsEM%d_tpc", l);
        if (t.GetBranch(b)) t.SetBranchAddress(b, &em_tpc[l]);
    }
    t.SetBranchAddress("fcsJP2", &jp2);
    t.SetBranchAddress("fcsJPA1", &jp1[0]);
    t.SetBranchAddress("fcsJPBC1", &jp1[1]);
    t.SetBranchAddress("fcsJPDE1", &jp1[2]);
    t.SetBranchAddress("fcsJPA0", &jp0[0]);
    t.SetBranchAddress("fcsJPBC0", &jp0[1]);
    t.SetBranchAddress("fcsJPDE0", &jp0[2]);
    const bool have_run = t.GetBranch("event_run") != 0;
    if (have_run) t.SetBranchAddress("event_run", &event_run);
    if (t.GetBranch("event_vertex_z")) t.SetBranchAddress("event_vertex_z", &vz);
    if (t.GetBranch("event_vertex_source")) t.SetBranchAddress("event_vertex_source", &vz_source);
    t.SetBranchAddress("reco_eta", &eta);
    t.SetBranchAddress("reco_phi", &phi);
    t.SetBranchAddress("reco_E", &e_raw);
    t.SetBranchAddress("reco_ue_E", &e_ue);
    t.SetBranchAddress("reco_x", &x);
    t.SetBranchAddress("reco_y", &y);
    t.SetBranchAddress("reco_pt", &pt_raw);
    t.SetBranchAddress("reco_E_corr", &e_corr);
    t.SetBranchAddress("reco_x_F_corr", &xf_corr);
    t.SetBranchAddress("reco_trg_em_et", &trg_em);
    t.SetBranchAddress("reco_trg_em_hmax", &trg_hmax);
    t.SetBranchAddress("reco_trg_jp_et", &trg_jp);
    t.SetBranchAddress("reco_trg_jp_mask", &trg_jpmask);
    if (!t.GetBranch("reco_from_points")) { std::cerr << "ERROR: no reco_from_points; photon counting needs point jets" << std::endl; return 1; }
    t.SetBranchAddress("reco_from_points", &from_points);
    t.SetBranchAddress("reco_nconst", &nconst);
    t.SetBranchAddress("reco_const_jetIndex", &c_jet);
    t.SetBranchAddress("reco_const_px", &c_px);
    t.SetBranchAddress("reco_const_py", &c_py);
    t.SetBranchAddress("reco_const_pz", &c_pz);
    t.SetBranchAddress("reco_const_E", &c_e);
    const bool have_plead = t.GetBranch("reco_pi0_lead_prob") != 0;
    if (have_plead) t.SetBranchAddress("reco_pi0_lead_prob", &plead);
    else std::cerr << "Warning: no reco_pi0_lead_prob in the jet files; h_xFP_* stay empty" << std::endl;
    const bool have_gup = t.GetBranch("reco_pi0_lead_gup") != 0 && t.GetBranch("reco_pi0_lead_zpair") != 0 && t.GetBranch("reco_pi0_lead_pair") != 0;
    if (have_gup)
    {
        t.SetBranchAddress("reco_pi0_lead_gup", &gup); t.SetBranchAddress("reco_pi0_lead_zpair", &zpair);
        t.SetBranchAddress("reco_pi0_lead_pair", &pleadpair);
    }
    else std::cerr << "Warning: no reco_pi0_lead_gup in the jet files; h_xFg_* stay empty" << std::endl;
    const bool have_z = have_gup && t.GetBranch("reco_pi0_lead_z1") != 0;
    if (have_z) t.SetBranchAddress("reco_pi0_lead_z1", &z1);
    else std::cerr << "Warning: no reco_pi0_lead_z1 in the jet files; h_xFPz_* stay empty" << std::endl;

    TFile out(argv[2], "RECREATE");
    // [sample][EM tag][e]: samples SampSuffix(), EM tags for JP levels
    // (Booked()); e 0 E <= sqrt(s)/2, 1 above it (_hiE)
    static TH1D *h[kNSamp][kNEmTags][2][kNTrgLevels][kNPhotonClasses][kNSides][kNSpinStates];
    static TH1D *h_cos[kNSamp][kNEmTags][2][kNTrgLevels][kNPhotonClasses], *h_xfs[kNSamp][kNEmTags][2][kNTrgLevels][kNPhotonClasses];
    // pair yields, E <= sqrt(s)/2: [0] triggering >= the default factor, [1]
    // triggering 1.15 to the default factor, [2] companions, [3] JP levels: jets
    // >= the default factor in JP that are not >= the default factor in EM (for an EM +
    // JP union of [0] without double counting)
    static TH2D *h_m[4][kNTrgLevels][kNSides][kNSpinStates];
    // x_F vs P_lead, the same four samples as h_m, every jet with E <= sqrt(s)/2
    static TH2D *h_p[4][kNTrgLevels][kNSides][kNSpinStates], *h_pcos[4][kNTrgLevels], *h_pxfs[4][kNTrgLevels];
    // Collins (collins/CollinsHistograms.h): sample 0 default, 1 _jponly
    static CollinsHistograms collins;
    // A_N vs P_lead class vs z: [0] default, [1] _jponly
    static TH3D *h_pz[2][kNTrgLevels][kNSides][kNSpinStates], *h_pzcos[2][kNTrgLevels], *h_pzxfs[2][kNTrgLevels];
    const char *const msfx[4] = {"", "_t115", kCompSuffix, "_jponly"};
    TH1D *h_events[kNTrgLevels], *h_pt[kNTrgLevels], *h_ptraw[kNTrgLevels], *h_ratio[kNTrgLevels];
    for (int l = 0; l < kNTrgLevels; l++)
    {
        for (int g = 0; g < kNSamp; g++)
            for (int tg = 0; tg < kNEmTags; tg++)
            {
                if (!Booked(l, g, tg)) continue;
                for (int e = 0; e < 2; e++)
                {
                    const TString sfx = Suffix(l, g, tg) + (e ? kHiESuffix : "");
                    for (int k = 0; k < kNPhotonClasses; k++)
                        for (int s = 0; s < kNSides; s++)
                            for (int p = 0; p < kNSpinStates; p++)
                            {
                                TString name = TString::Format("h_xF_%s_%s_%s_%s", kTrgLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]) + sfx;
                                h[g][tg][e][l][k][s][p] = new TH1D(name, name + ";x_{F};jets", kNXfBins, kXfMin, kXfMax);
                            }
                    for (int k = 0; k < kNPhotonClasses; k++)
                    {
                        TString name = TString::Format("h_abscos_%s_%s", kTrgLevelName[l], PhotonClassName(k)) + sfx;
                        h_cos[g][tg][e][l][k] = new TH1D(name, name + ";x_{F};sum of |cos #phi|", kNXfBins, kXfMin, kXfMax);
                        name = TString::Format("h_xfsum_%s_%s", kTrgLevelName[l], PhotonClassName(k)) + sfx;
                        h_xfs[g][tg][e][l][k] = new TH1D(name, name + ";x_{F};sum of x_{F}", kNXfBins, kXfMin, kXfMax);
                    }
                }
            }
        for (int m = 0; m < 4; m++)
        {
            if (m == 3 && !IsJpLevel(l)) continue;
            for (int s = 0; s < kNSides; s++)
                for (int p = 0; p < kNSpinStates; p++)
                {
                    TString name = TString::Format("h_xFm_%s_%s_%s", kTrgLevelName[l], kSideName[s], kSpinStateName[p]) + msfx[m];
                    h_m[m][l][s][p] = new TH2D(name, name + ";x_{F};M_{#gamma#gamma} (GeV)", kNXfBins, kXfMin, kXfMax, kNMassBins, kMassMin, kMassMax);
                    name = TString::Format("h_xFP_%s_%s_%s", kTrgLevelName[l], kSideName[s], kSpinStateName[p]) + msfx[m];
                    h_p[m][l][s][p] = new TH2D(name, name + ";x_{F};P_{lead}", kNXfBins, kXfMin, kXfMax, kNPleadBins, kPleadMin, kPleadMax);
                }
            TString name = TString::Format("h_abscosP_%s", kTrgLevelName[l]) + msfx[m];
            h_pcos[m][l] = new TH2D(name, name + ";x_{F};P_{lead}", kNXfBins, kXfMin, kXfMax, kNPleadBins, kPleadMin, kPleadMax);
            name = TString::Format("h_xfsumP_%s", kTrgLevelName[l]) + msfx[m];
            h_pxfs[m][l] = new TH2D(name, name + ";x_{F};P_{lead}", kNXfBins, kXfMin, kXfMax, kNPleadBins, kPleadMin, kPleadMax);
            if (m != 0 && m != 3) continue;
            const int mg = m == 0 ? 0 : 1;
            {
                double xa[kNXfBins + 1], za[kNZBinsA + 1];
                for (int i = 0; i <= kNXfBins; i++) xa[i] = kXfMin + i * (kXfMax - kXfMin) / kNXfBins;
                for (int i = 0; i <= kNZBinsA; i++) za[i] = i * kZMaxA / kNZBinsA;
                for (int s = 0; s < kNSides; s++)
                    for (int p = 0; p < kNSpinStates; p++)
                    {
                        TString zn = TString::Format("h_xFPz_%s_%s_%s", kTrgLevelName[l], kSideName[s], kSpinStateName[p]) + msfx[m];
                        h_pz[mg][l][s][p] = new TH3D(zn, zn + ";x_{F};P_{lead};z", kNXfBins, xa, kNPleadClsZ, kPleadEdgesZ, kNZBinsA, za);
                    }
                TString zn = TString::Format("h_abscosPz_%s", kTrgLevelName[l]) + msfx[m];
                h_pzcos[mg][l] = new TH3D(zn, zn + ";x_{F};P_{lead};z", kNXfBins, xa, kNPleadClsZ, kPleadEdgesZ, kNZBinsA, za);
                zn = TString::Format("h_xfsumPz_%s", kTrgLevelName[l]) + msfx[m];
                h_pzxfs[mg][l] = new TH3D(zn, zn + ";x_{F};P_{lead};z", kNXfBins, xa, kNPleadClsZ, kPleadEdgesZ, kNZBinsA, za);
            }
            collins.Book(mg, l, msfx[m]);
        }
        h_events[l] = new TH1D(TString::Format("h_events_%s", kTrgLevelName[l]), ";spin state;events", kNSpinStates, 0, kNSpinStates);
        for (int p = 0; p < kNSpinStates; p++) h_events[l]->GetXaxis()->SetBinLabel(p + 1, kSpinStateName[p]);
        h_pt[l] = new TH1D(TString::Format("h_ptcorr_%s", kTrgLevelName[l]), ";corrected jet p_{T} (GeV);jets", 100, 0, 20);
        h_ptraw[l] = new TH1D(TString::Format("h_ptraw_%s", kTrgLevelName[l]), ";uncorrected jet p_{T} (GeV);jets", 100, 0, 20);
        h_ratio[l] = new TH1D(TString::Format("h_trgratio_%s", kTrgLevelName[l]), TString::Format(";%s trigger value / %.1f GeV;jets", kTrgLevelName[l], kTrgThresholdGeV[l]), 60, 0, 6);
    }
    // sanity distributions: [0] jets of events passing the event cuts, [1]
    // selected jets, [2] companion jets; all levels together, and per level
    const int kNStages = 3;
    const char *const stage[kNStages] = {"evt", "sel", "comp"};
    static QaSet qa[kNStages], qa_lvl[kNStages][kNTrgLevels];
    for (int g = 0; g < kNStages; g++)
    {
        qa[g].Book(stage[g]);
        for (int l = 0; l < kNTrgLevels; l++) qa_lvl[g][l].Book(TString::Format("%s_%s", stage[g], kTrgLevelName[l]));
    }
    TH1D *q_sqrts = new TH1D("qa_sqrt_s", ";event_sqrt_s_gev (GeV);events", 200, 400, 600);
    TH1D *q_vz = new TH1D("qa_vz", ";event vertex z (cm), events with a vertex;events", 80, -80, 80);
    TH1D *q_vzsrc = new TH1D("qa_vz_source", ";vertex source;events", 4, 0, 4);
    const char *const vz_name[4] = {"none", "TPC", "VPD", "BBC"};
    for (int i = 0; i < 4; i++) q_vzsrc->GetXaxis()->SetBinLabel(i + 1, vz_name[i]);
    TH1D *q_trig = new TH1D("qa_em_listed", ";listed EM triggers;events", 8, 0, 8);
    for (int i = 0; i < 8; i++)
        q_trig->GetXaxis()->SetBinLabel(i + 1, i == 0 ? "none" : TString::Format("%s%s%s", i & 1 ? "0" : "", i & 2 ? "1" : "", i & 4 ? "2" : ""));
    TH1D *q_jptrig = new TH1D("qa_jp_listed", ";listed JP levels;events", 8, 0, 8);
    for (int i = 0; i < 8; i++)
        q_jptrig->GetXaxis()->SetBinLabel(i + 1, i == 0 ? "none" : TString::Format("%s%s%s", i & 1 ? "0" : "", i & 2 ? "1" : "", i & 4 ? "2" : ""));
    TH1D *q_spin = new TH1D("qa_spin_config", ";Spin_config;events", 17, -1, 16);
    TH2D *q_xy_side[kNSides];
    for (int s = 0; s < kNSides; s++)
        q_xy_side[s] = new TH2D(TString::Format("qa_xy_sel_%s", kSideName[s]), ";jet x (cm);jet y (cm)", 60, -150, 150, 44, -110, 110);
    TH1D *h_cut = new TH1D("h_cutflow", ";;jets", kNCuts + 1, 0, kNCuts + 1);
    for (int c = 0; c < kNCuts; c++) h_cut->GetXaxis()->SetBinLabel(c + 1, kCutName[c]);
    h_cut->GetXaxis()->SetBinLabel(kNCuts + 1, "events read");

    std::set<int> runs;
    int tree_number = -1, file_run = 0;
    const Long64_t n = t.GetEntries();
    for (Long64_t i = 0; i < n; i++)
    {
        t.GetEntry(i);
        // the run: event_run, else from the file name (jet_output_<run>_<file>.root)
        if (t.GetTreeNumber() != tree_number) file_run = RunFromFileName(t.GetFile()->GetName());
        tree_number = t.GetTreeNumber();
        const int cur_run = have_run && event_run > 0 ? event_run : file_run;
        runs.insert(cur_run);
        h_cut->Fill(kNCuts);
        // listed levels; for JP0/JP1 the listed regions' patch masks
        bool listed[kNTrgLevels];
        for (int l = 0; l < kNEmLevels; l++) listed[l] = em[l] || em_tpc[l];
        int jp_region[2] = {0, 0};
        for (int r = 0; r < 3; r++)
        {
            const int m = r == 0 ? kJpMaskA : r == 1 ? kJpMaskBC : kJpMaskDE;
            if (jp0[r]) jp_region[0] |= m;
            if (jp1[r]) jp_region[1] |= m;
        }
        // runs whose JP triggers are not used (AnConfig.h kJpBadRuns)
        int jp2_ = jp2;
        if (IsJpBadRun(cur_run)) { jp_region[0] = jp_region[1] = 0; jp2_ = 0; }
        listed[3] = jp_region[0] != 0;
        listed[4] = jp_region[1] != 0;
        listed[5] = jp2_ != 0;
        // the event's level: highest listed EM, else highest listed JP
        int level = -1;
        for (int l = 0; l < kNEmLevels; l++) if (listed[l]) level = l;
        if (level < 0) for (int l = kNEmLevels; l < kNTrgLevels; l++) if (listed[l]) level = l;
        const int spin = SpinState(spin_config);
        const size_t nj = eta->size();
        if (from_points != 1) { std::cerr << "ERROR: " << t.GetFile()->GetName() << " has hit jets (reco_from_points = " << from_points << "); photon counting needs point jets" << std::endl; return 1; }
        if (nconst->size() != nj || trg_em->size() != nj || trg_jp->size() != nj) { std::cerr << "ERROR: jet branch sizes differ at entry " << i << std::endl; return 1; }
        h_cut->Fill(kAllJets, nj);
        if (led) continue;
        h_cut->Fill(kNoLed, nj);
        q_trig->Fill((listed[0] ? 1 : 0) + (listed[1] ? 2 : 0) + (listed[2] ? 4 : 0));
        q_jptrig->Fill((listed[3] ? 1 : 0) + (listed[4] ? 2 : 0) + (listed[5] ? 4 : 0));
        if (level < 0) continue;
        h_cut->Fill(kTrgListed, nj);
        q_spin->Fill(spin_config);
        if (spin < 0) continue;
        h_cut->Fill(kSpin, nj);
        h_events[level]->Fill(spin);
        q_vzsrc->Fill(vz_source >= 0 && vz_source < 4 ? vz_source : 0);
        if (vz > -9000) q_vz->Fill(vz);
        q_sqrts->Fill(sqrt_s);
        qa[0].njet->Fill(nj);
        qa_lvl[0][level].njet->Fill(nj);
        // the two constituents of each two-photon jet: pair mass and orientation (-1 if none)
        std::vector<double> mass(nj, -1), pang(nj, -1);
        {
            std::vector<int> first(nj, -1), second(nj, -1);
            for (size_t c = 0; c < c_jet->size(); c++)
            {
                const int jj = (*c_jet)[c];
                if (jj < 0 || jj >= (int)nj) continue;
                if (first[jj] < 0) first[jj] = c;
                else if (second[jj] < 0) second[jj] = c;
            }
            for (size_t j = 0; j < nj; j++)
                if ((*nconst)[j] == 2 && second[j] >= 0)
                {
                    const int a = first[j], b = second[j];
                    PairOf((*c_px)[a], (*c_py)[a], (*c_pz)[a], (*c_e)[a], (*c_px)[b], (*c_py)[b], (*c_pz)[b], (*c_e)[b], mass[j], pang[j]);
                }
        }
        // per jet, every cut but the trigger passed (base); then, independently,
        // its EM result (level, slice; or veto failure) and JP result (level,
        // slice), -1 if none
        std::vector<int> base(nj, 0), emL(nj, -1), emS(nj, -1), emVeto(nj, 0), jpL(nj, -1), jpS(nj, -1);
        int event_em = -1, event_jp = -1;  // highest EM / JP level among jets at >= the default factor x threshold
        for (size_t j = 0; j < nj; j++)
        {
            const double E = (*e_corr)[j];
            // sanity distributions of every jet here, by the event's level; the selected jets' below
            const double ptc0 = E > 0 ? E / std::cosh((*eta)[j]) : -1;
            qa[0].Fill((*eta)[j], (*phi)[j], (*e_raw)[j], (*e_ue)[j], E, (*pt_raw)[j], ptc0, (*xf_corr)[j], (*nconst)[j], (*x)[j], (*y)[j]);
            qa_lvl[0][level].Fill((*eta)[j], (*phi)[j], (*e_raw)[j], (*e_ue)[j], E, (*pt_raw)[j], ptc0, (*xf_corr)[j], (*nconst)[j], (*x)[j], (*y)[j]);
            qa[0].FillPair(mass[j], pang[j], E);
            qa_lvl[0][level].FillPair(mass[j], pang[j], E);
            if (!(E > 0)) continue;
            h_cut->Fill(kCalibrated);
            if (!((*xf_corr)[j] > 0)) continue;
            h_cut->Fill(kXfPositive);
            const double pt = E / std::cosh((*eta)[j]);
            if (pt <= kJetPtCorrMin) continue;
            h_cut->Fill(kPtMin);
            h_pt[level]->Fill(pt);
            h_ptraw[level]->Fill((*pt_raw)[j]);
            base[j] = 1;
            // EM: the lowest listed EM level at >= 1.15 x threshold, and the veto margin
            const double em_et = (*trg_em)[j];
            for (int l = 0; l < kNEmLevels && emL[j] < 0; l++)
                if (listed[l] && em_et >= kTrgSliceEdge[0] * kTrgThresholdGeV[l]) emL[j] = l;
            if (emL[j] >= 0)
            {
                const double ratio = em_et / kTrgThresholdGeV[emL[j]];
                h_ratio[emL[j]]->Fill(ratio);
                if ((*trg_hmax)[j] >= 0 && (*trg_hmax)[j] * 4 * kEmVetoMargin <= em_et)
                {
                    emS[j] = 0;
                    for (int sl = 1; sl < kNTrgSlices; sl++) if (ratio >= kTrgSliceEdge[sl]) emS[j] = sl;
                    if (ratio >= kEmCutFactor - 1e-9 && emL[j] > event_em) event_em = emL[j];
                }
                else emVeto[j] = 1;
            }
            // JP: the lowest listed JP level at >= 1.15 x threshold (JP0/JP1 in a listed region)
            const double jp_et = (*trg_jp)[j];
            for (int l = kNEmLevels; l < kNTrgLevels && jpL[j] < 0; l++)
            {
                if (!listed[l]) continue;
                if (l < 5 && !((*trg_jpmask)[j] & jp_region[l - kNEmLevels])) continue;
                if (jp_et >= kTrgSliceEdge[0] * kTrgThresholdGeV[l]) jpL[j] = l;
            }
            if (jpL[j] >= 0)
            {
                const double ratio = jp_et / kTrgThresholdGeV[jpL[j]];
                h_ratio[jpL[j]]->Fill(ratio);
                jpS[j] = 0;
                for (int sl = 1; sl < kNTrgSlices; sl++) if (ratio >= kTrgSliceEdge[sl]) jpS[j] = sl;
                if (ratio >= kJpCutFactor - 1e-9 && jpL[j] > event_jp) event_jp = jpL[j];
            }
        }
        // companions go to the highest EM level of the event's default jets, else the highest JP level
        const int event_jl = event_em >= 0 ? event_em : event_jp;
        int nsel = 0, ncomp = 0;
        for (size_t j = 0; j < nj; j++)
        {
            if (!base[j]) continue;
            const double E = (*e_corr)[j];
            const double pt = E / std::cosh((*eta)[j]);
            // above sqrt(s)/2: kept apart (the FMS cut), in the _hiE histograms
            const int hi = E > 0.5 * sqrt_s ? 1 : 0;
            const int side = computeSide((*x)[j]) > 0 ? 1 : 0;
            const int cls = PhotonClass((*nconst)[j]);
            const double xf = (*xf_corr)[j];
            const double acos_ = std::fabs((*x)[j]) / std::hypot((*x)[j], (*y)[j]);
            const bool em_on = emS[j] >= 0, jp_on = jpS[j] >= 0;
            const bool em_dflt = em_on && kTrgSliceEdge[emS[j]] >= kEmCutFactor - 1e-9;
            const bool jp_dflt = jp_on && kTrgSliceEdge[jpS[j]] >= kJpCutFactor - 1e-9;
            // fill one (sample, tag, level)
            auto fill = [&](int g, int tg, int l)
            {
                h[g][tg][hi][l][cls][side][spin]->Fill(xf);
                h_xfs[g][tg][hi][l][cls]->Fill(xf, xf);
                h_cos[g][tg][hi][l][cls]->Fill(xf, acos_);
            };
            if (em_on) fill(emS[j], 0, emL[j]);
            else if (emVeto[j]) fill(kSampVeto, 0, emL[j]);
            if (jp_on) fill(jpS[j], em_on ? 1 + emS[j] : 0, jpL[j]);
            const bool any_on = em_on || jp_on, dflt = em_dflt || jp_dflt;
            const bool comp = !any_on && !emVeto[j] && event_jl >= 0;
            if (comp) { fill(kSampComp, 0, event_jl); h_cut->Fill(kCompanion); ncomp++; }
            if (emVeto[j]) h_cut->Fill(kVeto);
            if (any_on) h_cut->Fill(kTrg115);
            if (dflt)
            {
                h_cut->Fill(kTrgDefault);
                if (!hi) h_cut->Fill(kBelowHalfS);
                if (em_dflt) h_cut->Fill(kEmDefault);
                if (jp_dflt) h_cut->Fill(kJpDefault);
                if (em_dflt && jp_dflt) h_cut->Fill(kBothDefault);
                nsel++;
            }
            // pair yields (E <= sqrt(s)/2): per family level, as the yields; [3] JP not EM at the default
            if (mass[j] >= 0 && !hi)
            {
                if (em_on) h_m[em_dflt ? 0 : 1][emL[j]][side][spin]->Fill(xf, mass[j]);
                if (jp_on) h_m[jp_dflt ? 0 : 1][jpL[j]][side][spin]->Fill(xf, mass[j]);
                if (jp_dflt && !em_dflt) h_m[3][jpL[j]][side][spin]->Fill(xf, mass[j]);
                if (comp) h_m[2][event_jl][side][spin]->Fill(xf, mass[j]);
            }
            // x_F vs P_lead yields (E <= sqrt(s)/2), the same samples
            const double pl = have_plead && j < plead->size() ? (*plead)[j] : -1;
            if (pl >= 0 && !hi)
            {
                auto pfill = [&](int m, int l)
                {
                    h_p[m][l][side][spin]->Fill(xf, pl);
                    h_pcos[m][l]->Fill(xf, pl, acos_);
                    h_pxfs[m][l]->Fill(xf, pl, xf);
                };
                if (em_on) pfill(em_dflt ? 0 : 1, emL[j]);
                if (jp_on) pfill(jp_dflt ? 0 : 1, jpL[j]);
                if (jp_dflt && !em_dflt) pfill(3, jpL[j]);
                if (comp) pfill(2, event_jl);
            }
            // A_N vs P_lead vs the pi0 candidate's z: default jets and JP-not-EM (_jponly)
            if (have_z && pl >= 0 && !hi)
            {
                const double lp = (*pleadpair)[j], zp = (*zpair)[j], zl = (*z1)[j];
                const double zc = (lp >= kZCandPairMin && zp >= 0) ? zp : zl;
                if (zc >= 0)
                {
                    auto zfill = [&](int mg, int l)
                    {
                        h_pz[mg][l][side][spin]->Fill(xf, pl, zc);
                        h_pzcos[mg][l]->Fill(xf, pl, zc, acos_);
                        h_pzxfs[mg][l]->Fill(xf, pl, zc, xf);
                    };
                    if (em_dflt) zfill(0, emL[j]);
                    if (jp_dflt) zfill(0, jpL[j]);
                    if (jp_dflt && !em_dflt) zfill(1, jpL[j]);
                }
            }
            // Collins yields: default jets (EM and JP levels) and JP-not-EM (_jponly)
            const double gg = have_gup && j < gup->size() ? (*gup)[j] : -9;
            const double plp = have_gup && j < pleadpair->size() ? (*pleadpair)[j] : -1;
            if (CollinsHistograms::Has(gg, plp) && !hi)
            {
                const double zz = (*zpair)[j];
                if (em_dflt) collins.Fill(0, emL[j], side, spin, xf, gg, zz, plp);
                if (jp_dflt) collins.Fill(0, jpL[j], side, spin, xf, gg, zz, plp);
                if (jp_dflt && !em_dflt) collins.Fill(1, jpL[j], side, spin, xf, gg, zz, plp);
            }
            // QA: default triggering jets ("sel", by the EM level, else the JP level) and companions ("comp")
            const int st = comp ? 2 : dflt ? 1 : -1;
            if (st < 0) continue;
            const int ql = comp ? event_jl : em_dflt ? emL[j] : jpL[j];
            qa[st].Fill((*eta)[j], (*phi)[j], (*e_raw)[j], (*e_ue)[j], E, (*pt_raw)[j], pt, xf, (*nconst)[j], (*x)[j], (*y)[j]);
            qa_lvl[st][ql].Fill((*eta)[j], (*phi)[j], (*e_raw)[j], (*e_ue)[j], E, (*pt_raw)[j], pt, xf, (*nconst)[j], (*x)[j], (*y)[j]);
            qa[st].FillPair(mass[j], pang[j], E);
            qa_lvl[st][ql].FillPair(mass[j], pang[j], E);
            if (st == 1) q_xy_side[side]->Fill((*x)[j], (*y)[j]);
        }
        qa[1].njet->Fill(nsel);
        qa_lvl[1][level].njet->Fill(nsel);
        // companions per event, in events with a triggering jet
        if (event_jl >= 0) { qa[2].njet->Fill(ncomp); qa_lvl[2][event_jl].njet->Fill(ncomp); }
    }

    std::ostringstream rs;
    for (std::set<int>::const_iterator r = runs.begin(); r != runs.end(); ++r) rs << (r == runs.begin() ? "" : " ") << *r;
    TNamed("runs", rs.str().c_str()).Write();
    TNamed("spin_convention", kSpinConventionTag).Write();
    if (!have_plead) TNamed("plead_missing", "jet files without reco_pi0_lead_prob: h_xFP_* empty").Write();
    if (!have_gup) TNamed("collins_missing", "jet files without reco_pi0_lead_gup: h_xFg_* empty").Write();
    if (!have_z) TNamed("pz_missing", "jet files without reco_pi0_lead_z1: h_xFPz_* empty").Write();
    TNamed("trigger_rule", TString::Format("trigger values (reco_trg_*), slices from %.2f, default factor EM %.2f JP %.2f, EM veto margin %.2f", kTrgSliceEdge[0], kEmCutFactor, kJpCutFactor, kEmVetoMargin).Data()).Write();
    out.Write();
    printf("AnHistograms: %lld events, runs %s\n", n, rs.str().c_str());
    for (int c = 0; c < kNCuts; c++) printf("  %-16s %.0f jets\n", kCutName[c], h_cut->GetBinContent(c + 1));
    out.Close();
    return 0;
}
