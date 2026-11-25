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
// (default kEmCutFactor = kJpCutFactor = 1.15 since 2026-11-17; slices from 1.15). Each jet fills x_F
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
// FMS-style pi0 candidates (AnConfig.h "FMS-style pi0 candidates"): every
// photon pair of a default jet (the A_N jet cuts, E <= sqrt(s)/2), with
// Z_gg < 0.7, both photons half a tower inside the ECal, pT > 2 GeV (scaled) and
// pT >= 1.3 x the threshold of the level it fills (the host jet's level),
// M_gg < 0.3 GeV and pi0 x_F > 0, fills
//   h_xFpi_<level>_<side>_<spin state>[_jponly]
// (TH3D: the pi0's x_F (h_xF_* binning), M_gg in 5 MeV bins to 0.3, both on
// the pi0 energy scale (pi0_scale.txt table, else AnConfig.h kPi0EnergyScale;
// TNamed "pi0_scale" records the factor per run), z_em =
// E_pair / raw jet E in the classes 0-0.9-0.98-1). _jponly: pi0s passing a JP
// level but not an EM level (the EM + JP union per pi0). It also fills
// h_abscospi_<level>[_jponly] and h_xfsumpi_<level>[_jponly], with the same
// axes, weighted with the pair's |cos phi| and x_F. Side and phi are the
// pair's own. These feed AnAsymmetry -I and sivers/pi0_bgsub.C.
// Also written:
//   h_cutflow                    jets (and events) at each cut
//   qa_cut_vz, qa_cut_vz_source, qa_cut_bbc_e_ntiles, qa_cut_tof_mult
//                                each event-level cut variable, N-1 (every
//                                other event cut passed), with TNamed
//                                event_cut_lines (the cut and variation
//                                positions) for sivers/plot_event_cuts.C
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
// (event_vertex_z and its source: EPD, BBC TAC or nominal), the listed EM and JP triggers and
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
// Collision cut (2026-11-17, AnConfig.h PassCollisionCut: east BBC >= 1 tile,
// TOF multiplicity > 2, as the FMS): events failing it fill none of the above;
// they go only to the _vbbcveto copies below. Trees without the branches are
// not cut.
// Systematic variations (AnConfig.h kVarSuffix): copies of the default
// photon-class yields (h_xF_ / h_abscos_ / h_xfsum_, JP levels split by
// kVarEmTag) and of the pi0-candidate yields (h_xFpi_ / h_abscospi_ /
// h_xfsumpi_, default and _jponly) with the variation's cut, suffix at the end:
// fiducial 25 / 30 cm (distance of the jet axis from the ECal edge at z_proj,
// from reco_x / reco_y), |vz| < 60, measured vertex only, collision-cut
// failures; and _vjesup / _vjesdn (photon classes only): every jet's E, pT and
// x_F scaled by 1 +- delta(E), the JES map's error (--jes-error=<map json>) at
// the map's lookup energy: the calibration repeated with the map shifted to
// F / (1 +- delta) (identical to ApplyCorrections_em_only --jes-shift=+1 / -1
// and a nominal pass), same events and raw trigger values. AnAsymmetry -v reads
// them; the systematic is max(|A_up - A|, |A_dn - A|) per bin, symmetric.
// Point-threshold copies _vthr150 / _vthr200 (Sean 2026-11-24): a second call
// on the job's jets rebuilt at FCS_POINT_EMIN = 1.5 / 2.0 GeV and calibrated
// with that threshold's map, with --thr-copy=150 / 200: it opens the nominal
// call's output file (UPDATE) and adds only that copy's histograms (photon
// classes and pi0 candidates, filled from its default jets as the other copies
// are), plus a TNamed thr_copy_<thr>. The nominal call leaves them out.
// Trigger reference tallies (AnConfig.h "Trigger reference tallies"): with a
// third argument, a list of the job's trigger reference files, every
// reference event of those files fills h_trgref_* (no other cut), for
// TriggerPlayground/trg_nominal.cpp.
// Output files of one run, or of a fill, add with hadd.
//
// Usage: ./AnHistograms <jet file list> <output.root> [trigger reference file list]
// Build: g++ -m64 -O2 AnHistograms.cpp VectorDict.cxx -o AnHistograms -Icollins
//            `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//        (the jobs copy collins/*.h next to it, so need no -Icollins)

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TList.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "VectorDict.h"
#include "AnConfig.h"
#include "CollinsHistograms.h"

// The run number in a jet_output_<run>_<file>.root name, 0 if none
// pi0 energy scale table (Sean 2026-11-16: global, per run or per period,
// chosen after a pass and applied by a rehist): $PI0_SCALE_TABLE, else
// ./pi0_scale.txt (shipped by an_job / rehist_job). Lines "default <factor>"
// and "<run> <factor>", '#' comments. Runs not listed take the table's default,
// else AnConfig.h kPi0EnergyScale; with no table, everything takes
// kPi0EnergyScale.
struct Pi0Scale
{
    std::map<int, double> run;
    double dflt = kPi0EnergyScale;
    std::string source = "AnConfig.h kPi0EnergyScale (no table)";
    void Load()
    {
        const char *env = std::getenv("PI0_SCALE_TABLE");
        const std::string path = env && *env ? env : "pi0_scale.txt";
        std::ifstream in(path.c_str());
        if (!in) return;
        std::string line;
        int n = 0;
        while (std::getline(in, line))
        {
            const size_t h = line.find('#');
            if (h != std::string::npos) line = line.substr(0, h);
            std::istringstream is(line);
            std::string key;
            double f = 0;
            if (!(is >> key >> f) || !(f > 0)) continue;
            if (key == "default") dflt = f;
            else { run[std::atoi(key.c_str())] = f; n++; }
        }
        source = path + TString::Format(" (%d runs, default %.5f)", n, dflt).Data();
    }
    double Of(int r) const
    {
        std::map<int, double>::const_iterator it = run.find(r);
        return it == run.end() ? dflt : it->second;
    }
};

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
enum Cut { kAllJets, kNoLed, kTrgListed, kSpin, kVertex, kCollision, kCalibrated, kXfPositive, kPtMin, kTrg115, kTrgDefault, kBelowHalfS,
           kEmDefault, kJpDefault, kBothDefault, kVeto, kCompanion, kNCuts };
static const char *const kCutName[kNCuts] = {"all jets", "no LED", "EM/JP trigger", "spin state", "|vz| < 100", "BBC/TOF", "calibrated", "x_F > 0",
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

// The JES map's error vs energy (Sean 2026-11-18): the "jes_error_E" and
// "jes_error_total" arrays of the map JSON (relative, statistical (+)
// systematic; jet_scale/jes_map_compare.py error). delta(E) is linear in E
// between the points and held at the end values outside them.
struct JesError
{
    std::vector<double> e, d;
    bool ok() const { return e.size() >= 2 && e.size() == d.size(); }
    double at(double E) const
    {
        if (E <= e.front()) return d.front();
        if (E >= e.back()) return d.back();
        size_t i = 1;
        while (e[i] < E) i++;
        const double t = (E - e[i - 1]) / (e[i] - e[i - 1]);
        return d[i - 1] + t * (d[i] - d[i - 1]);
    }
};
static std::vector<double> JsonArray(const std::string &js, const std::string &key)
{
    std::vector<double> v;
    size_t p = js.find("\"" + key + "\"");
    if (p == std::string::npos) return v;
    p = js.find('[', p);
    const size_t q = js.find(']', p);
    if (p == std::string::npos || q == std::string::npos) return v;
    std::string body = js.substr(p + 1, q - p - 1);
    for (size_t i = 0; i < body.size(); i++) if (body[i] == ',') body[i] = ' ';
    std::istringstream in(body);
    double x;
    while (in >> x) v.push_back(x);
    return v;
}
static JesError ReadJesError(const std::string &path)
{
    JesError je;
    std::ifstream in(path.c_str());
    if (!in) return je;
    std::stringstream ss;
    ss << in.rdbuf();
    je.e = JsonArray(ss.str(), "jes_error_E");
    je.d = JsonArray(ss.str(), "jes_error_total");
    for (size_t i = 1; i < je.e.size(); i++) if (!(je.e[i] > je.e[i - 1])) { je.e.clear(); break; }
    return je;
}

int main(int argc, char **argv)
{
    // --jes-error=<map json>: the map's error vs energy for the _vjesup /
    // _vjesdn copies (an_job passes the JES map it applied)
    std::string jes_err_path;
    // --thr-copy=150 / 200: fill only the _vthr150 / _vthr200 copy into an existing output file
    int thr_v = -1;
    {
        int n = 1;
        for (int i = 1; i < argc; i++)
        {
            const std::string a = argv[i];
            if (a.compare(0, 12, "--jes-error=") == 0) jes_err_path = a.substr(12);
            else if (a.compare(0, 11, "--thr-copy=") == 0)
            {
                const std::string t = a.substr(11);
                if (t == "150") thr_v = kVarThr150;
                else if (t == "200") thr_v = kVarThr200;
                else { std::cerr << "ERROR: --thr-copy must be 150 or 200" << std::endl; return 1; }
            }
            else argv[n++] = argv[i];
        }
        argc = n;
    }
    const JesError jes_err = jes_err_path.empty() ? JesError() : ReadJesError(jes_err_path);
    const bool have_jes_err = jes_err.ok();
    if (!have_jes_err && thr_v < 0)
        std::cerr << "Warning: no JES error (" << (jes_err_path.empty() ? "no --jes-error" : "no jes_error_E / jes_error_total arrays in " + jes_err_path)
                  << "): _vjesup / _vjesdn copies empty" << std::endl;
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <jet file list> <output.root> [trigger reference file list] [--jes-error=<JES map json>] [--thr-copy=150|200]" << std::endl;
        return 1;
    }
    TChain t("jetTree");
    {
        std::ifstream in(argv[1]);
        std::string f;
        while (in >> f) t.Add(f.c_str());
    }
    if (t.GetNtrees() == 0) { std::cerr << "ERROR: no files in " << argv[1] << std::endl; return 1; }
    // trigger reference files (base names): their reference events fill h_trgref_*
    std::set<std::string> ref_files;
    if (argc > 3)
    {
        std::ifstream in(argv[3]);
        std::string f;
        while (in >> f) ref_files.insert(f.substr(f.find_last_of('/') + 1));
    }
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
    // the vertex cut (JetParameters.h pass_vertex_cut) needs the 2026-11-16 vertex
    // (EPD / BBC TAC calibrated / nominal); older trees, without epd_source, are not cut
    const bool vertex_cut = t.GetBranch("epd_source") != 0;
    // collision cut variables (jet trees from 2026-11-17 on; -1 = not cut)
    Int_t bbc_e_ntiles = -1, tof_mult = -1;
    if (t.GetBranch("bbc_e_ntiles")) t.SetBranchAddress("bbc_e_ntiles", &bbc_e_ntiles);
    if (t.GetBranch("tof_mult")) t.SetBranchAddress("tof_mult", &tof_mult);
    else std::cerr << "Warning: no bbc_e_ntiles / tof_mult (jet trees before 2026-11-17): no collision cut, _vbbcveto empty" << std::endl;
    // trigger reference flags and the hardware bits (as trg_nominal.cpp)
    UInt_t tcu = 0;
    if (t.GetBranch("fcs_tcu_bits")) t.SetBranchAddress("fcs_tcu_bits", &tcu);
    const char *ref_em_names[] = {"fcsJP2", "fcsJPA1", "fcsJPBC1", "fcsJPDE1", "fcsJPA0", "fcsJPBC0", "fcsJPDE0", "fcsDiJP", "fcsDiJPAsy",
                                  "fcsHad0", "fcsHad1", "fcsHad2", "fcsHHT_N", "fcsHHT_S", "fcsHTOT_N", "fcsHTOT_S", "fcs_HcalTot"};
    const char *ref_jp_names[] = {"fcsEM0", "fcsEM1", "fcsEM0_tpc", "fcsEM1_tpc", "fcsDY", "fcsDYAsy", "fcsJPsi", "fcsELE2"};
    const int n_ref_em = sizeof(ref_em_names) / sizeof(ref_em_names[0]), n_ref_jp = sizeof(ref_jp_names) / sizeof(ref_jp_names[0]);
    std::vector<Int_t> ref_em(n_ref_em, 0), ref_jp(n_ref_jp, 0);
    if (!ref_files.empty())
    {
        // the JP flags are already bound above to the A_N variables: read them from those
        for (int k = 7; k < n_ref_em; k++) if (t.GetBranch(ref_em_names[k])) t.SetBranchAddress(ref_em_names[k], &ref_em[k]);
        // fcsEM0/1 and their _tpc copies are bound to em[] / em_tpc[] above
        for (int k = 4; k < n_ref_jp; k++) if (t.GetBranch(ref_jp_names[k])) t.SetBranchAddress(ref_jp_names[k], &ref_jp[k]);
    }
    if (!vertex_cut) std::cerr << "WARNING: no epd_source branch (jet trees before the 2026-11-16 vertex code): no vertex cut" << std::endl;
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
    // the vertex the jets were made from (pi0 photon positions at z_proj); 0 in
    // trees before 2026-11-16 (jets from the origin)
    Float_t vx_used = 0, vy_used = 0, vz_used = 0;
    if (t.GetBranch("reco_vz_used")) t.SetBranchAddress("reco_vz_used", &vz_used);
    if (t.GetBranch("reco_vx_used")) t.SetBranchAddress("reco_vx_used", &vx_used);
    if (t.GetBranch("reco_vy_used")) t.SetBranchAddress("reco_vy_used", &vy_used);
    const bool have_z = have_gup && t.GetBranch("reco_pi0_lead_z1") != 0;
    if (have_z) t.SetBranchAddress("reco_pi0_lead_z1", &z1);
    else std::cerr << "Warning: no reco_pi0_lead_z1 in the jet files; h_xFPz_* stay empty" << std::endl;

    if (thr_v >= 0 && gSystem->AccessPathName(argv[2])) { std::cerr << "ERROR: --thr-copy adds to an existing output file; " << argv[2] << " not found" << std::endl; return 1; }
    TFile out(argv[2], thr_v >= 0 ? "UPDATE" : "RECREATE");
    if (out.IsZombie()) { std::cerr << "ERROR: cannot open " << argv[2] << std::endl; return 1; }
    // [sample][EM tag][e]: samples SampSuffix(), EM tags for JP levels
    // (Booked()); e 0 E <= sqrt(s)/2, 1 above it (_hiE)
    static TH1D *h[kNSamp][kNEmTags][2][kNTrgLevels][kNPhotonClasses][kNSides][kNSpinStates];
    static TH1D *h_cos[kNSamp][kNEmTags][2][kNTrgLevels][kNPhotonClasses], *h_xfs[kNSamp][kNEmTags][2][kNTrgLevels][kNPhotonClasses];
    // pair yields, E <= sqrt(s)/2: [0] triggering >= the default factor, [1]
    // triggering 1.15 to the default factor (empty since the default is 1.15),
    // [2] companions, [3] JP levels: jets >= the default factor in JP that are
    // not >= the default factor in EM (for an EM + JP union of [0] without
    // double counting), [4] / [5] the same as [0] / [3] at the trigger-
    // systematic factor (AnConfig.h kTrgSysFactor, suffix kSysSuffix)
    static TH2D *h_m[6][kNTrgLevels][kNSides][kNSpinStates];
    // x_F vs P_lead, the same six samples as h_m, every jet with E <= sqrt(s)/2
    static TH2D *h_p[6][kNTrgLevels][kNSides][kNSpinStates], *h_pcos[6][kNTrgLevels], *h_pxfs[6][kNTrgLevels];
    // Collins (collins/CollinsHistograms.h): sample 0 default, 1 _jponly, 2 / 3 the same at kTrgSysFactor
    static CollinsHistograms collins;
    // A_N vs P_lead class vs z: [0] default, [1] _jponly, [2] / [3] the same at kTrgSysFactor
    static TH3D *h_pz[4][kNTrgLevels][kNSides][kNSpinStates], *h_pzcos[4][kNTrgLevels], *h_pzxfs[4][kNTrgLevels];
    // FMS-style pi0 candidates, x_F x M_gg x z_em class: [0] default, [1] _jponly, [2] / [3] at kTrgSysFactor
    static TH3D *h_pi[4][kNTrgLevels][kNSides][kNSpinStates], *h_picos[4][kNTrgLevels], *h_pixfs[4][kNTrgLevels];
    const TString msfx[6] = {"", "_t115", kCompSuffix, "_jponly", kSysSuffix, TString("_jponly") + kSysSuffix};
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
        for (int m = 0; m < 6; m++)
        {
            if ((m == 3 || m == 5) && !IsJpLevel(l)) continue;
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
            if (m == 1 || m == 2) continue;
            const int mg = m == 0 ? 0 : m == 3 ? 1 : m == 4 ? 2 : 3;
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
                // FMS-style pi0 candidates
                double ma[kNPi0MassBins + 1];
                for (int i = 0; i <= kNPi0MassBins; i++) ma[i] = i * kPi0MassMax / kNPi0MassBins;
                const char *ax = ";#pi^{0} x_{F};M_{#gamma#gamma} (GeV);z_{em}";
                for (int s = 0; s < kNSides; s++)
                    for (int p = 0; p < kNSpinStates; p++)
                    {
                        TString pn = TString::Format("h_xFpi_%s_%s_%s", kTrgLevelName[l], kSideName[s], kSpinStateName[p]) + msfx[m];
                        h_pi[mg][l][s][p] = new TH3D(pn, pn + ax, kNXfBins, xa, kNPi0MassBins, ma, kNPi0ZCls, kPi0ZEdges);
                    }
                TString pn = TString::Format("h_abscospi_%s", kTrgLevelName[l]) + msfx[m];
                h_picos[mg][l] = new TH3D(pn, pn + ax, kNXfBins, xa, kNPi0MassBins, ma, kNPi0ZCls, kPi0ZEdges);
                pn = TString::Format("h_xfsumpi_%s", kTrgLevelName[l]) + msfx[m];
                h_pixfs[mg][l] = new TH3D(pn, pn + ax, kNXfBins, xa, kNPi0MassBins, ma, kNPi0ZCls, kPi0ZEdges);
            }
            collins.Book(mg, l, msfx[m]);
        }
        h_events[l] = new TH1D(TString::Format("h_events_%s", kTrgLevelName[l]), ";spin state;events", kNSpinStates, 0, kNSpinStates);
        for (int p = 0; p < kNSpinStates; p++) h_events[l]->GetXaxis()->SetBinLabel(p + 1, kSpinStateName[p]);
        h_pt[l] = new TH1D(TString::Format("h_ptcorr_%s", kTrgLevelName[l]), ";corrected jet p_{T} (GeV);jets", 100, 0, 20);
        h_ptraw[l] = new TH1D(TString::Format("h_ptraw_%s", kTrgLevelName[l]), ";uncorrected jet p_{T} (GeV);jets", 100, 0, 20);
        h_ratio[l] = new TH1D(TString::Format("h_trgratio_%s", kTrgLevelName[l]), TString::Format(";%s trigger value / the period's measured threshold (nominal %.1f GeV);jets", kTrgLevelName[l], kTrgThresholdGeV[l]), 60, 0, 6);
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
    TH1D *q_vz = new TH1D("qa_vz", ";event vertex z (cm);events", 120, -300, 300);
    TH1D *q_vzsrc = new TH1D("qa_vz_source", ";vertex source;events", 4, 0, 4);
    const char *const vz_name[4] = {"not set", "EPD", "BBC TAC", "nominal"};
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
    // event-level cut variables (Sean 2026-11-17), N-1: events with no LED, a
    // listed trigger and a spin state, passing every other event cut; drawn by
    // sivers/plot_event_cuts.C with the cut (solid) and its variations (dotted)
    // from the TNamed "event_cut_lines"
    TH1D *q_cut_vz = new TH1D("qa_cut_vz", "measured vertex (EPD or BBC TAC), collision cut passed;event vertex z (cm);events", 120, -300, 300);
    TH1D *q_cut_vzsrc = new TH1D("qa_cut_vz_source", "collision cut passed;vertex source;events", 4, 0, 4);
    for (int i = 0; i < 4; i++) q_cut_vzsrc->GetXaxis()->SetBinLabel(i + 1, vz_name[i]);
    TH1D *q_cut_bbc = new TH1D("qa_cut_bbc_e_ntiles", "vertex and TOF cuts passed;east BBC small tiles fired (100 < TDC < 2400, ADC > 50);events", 18, -1.5, 16.5);
    TH1D *q_cut_tof = new TH1D("qa_cut_tof_mult", "vertex and BBC cuts passed;TOF trigger multiplicity;events", 101, -1.5, 99.5);
    // systematic-variation copies (AnConfig.h kVarSuffix): photon classes [var][JP EM tag][level]...,
    // pi0 candidates [var][0 default, 1 _jponly][level]...
    static TH1D *hv[kNVar][2][kNTrgLevels][kNPhotonClasses][kNSides][kNSpinStates];
    static TH1D *hv_cos[kNVar][2][kNTrgLevels][kNPhotonClasses], *hv_xfs[kNVar][2][kNTrgLevels][kNPhotonClasses];
    static TH3D *hv_pi[kNVar][2][kNTrgLevels][kNSides][kNSpinStates], *hv_picos[kNVar][2][kNTrgLevels], *hv_pixfs[kNVar][2][kNTrgLevels];
    {
        double xa[kNXfBins + 1], ma[kNPi0MassBins + 1];
        for (int i = 0; i <= kNXfBins; i++) xa[i] = kXfMin + i * (kXfMax - kXfMin) / kNXfBins;
        for (int i = 0; i <= kNPi0MassBins; i++) ma[i] = i * kPi0MassMax / kNPi0MassBins;
        const char *ax = ";#pi^{0} x_{F};M_{#gamma#gamma} (GeV);z_{em}";
        for (int v = 0; v < kNVar; v++)
            for (int l = 0; l < kNTrgLevels; l++)
                for (int tg = 0; tg < (IsJpLevel(l) ? 2 : 1); tg++)
                {
                    const TString sfx = TString(IsJpLevel(l) ? kVarEmTag[tg] : "") + kVarSuffix[v];
                    for (int k = 0; k < kNPhotonClasses; k++)
                    {
                        for (int s = 0; s < kNSides; s++)
                            for (int p = 0; p < kNSpinStates; p++)
                            {
                                TString name = TString::Format("h_xF_%s_%s_%s_%s", kTrgLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p]) + sfx;
                                hv[v][tg][l][k][s][p] = new TH1D(name, name + ";x_{F};jets", kNXfBins, kXfMin, kXfMax);
                            }
                        TString name = TString::Format("h_abscos_%s_%s", kTrgLevelName[l], PhotonClassName(k)) + sfx;
                        hv_cos[v][tg][l][k] = new TH1D(name, name + ";x_{F};sum of |cos #phi|", kNXfBins, kXfMin, kXfMax);
                        name = TString::Format("h_xfsum_%s_%s", kTrgLevelName[l], PhotonClassName(k)) + sfx;
                        hv_xfs[v][tg][l][k] = new TH1D(name, name + ";x_{F};sum of x_{F}", kNXfBins, kXfMin, kXfMax);
                    }
                    // pi0: tg 0 default, 1 _jponly (JP levels); none for the JES copy
                    if (!VarHasPi0(v)) continue;
                    const TString psfx = TString(tg ? "_jponly" : "") + kVarSuffix[v];
                    for (int s = 0; s < kNSides; s++)
                        for (int p = 0; p < kNSpinStates; p++)
                        {
                            TString pn = TString::Format("h_xFpi_%s_%s_%s", kTrgLevelName[l], kSideName[s], kSpinStateName[p]) + psfx;
                            hv_pi[v][tg][l][s][p] = new TH3D(pn, pn + ax, kNXfBins, xa, kNPi0MassBins, ma, kNPi0ZCls, kPi0ZEdges);
                        }
                    TString pn = TString::Format("h_abscospi_%s", kTrgLevelName[l]) + psfx;
                    hv_picos[v][tg][l] = new TH3D(pn, pn + ax, kNXfBins, xa, kNPi0MassBins, ma, kNPi0ZCls, kPi0ZEdges);
                    pn = TString::Format("h_xfsumpi_%s", kTrgLevelName[l]) + psfx;
                    hv_pixfs[v][tg][l] = new TH3D(pn, pn + ax, kNXfBins, xa, kNPi0MassBins, ma, kNPi0ZCls, kPi0ZEdges);
                }
    }
    // trigger reference tallies (only with a reference file list)
    TH2D *tr_emn[kNEmLevels] = {0}, *tr_emk[kNEmLevels] = {0};
    TH1D *tr_jpn[3] = {0}, *tr_jpk[3] = {0};
    if (!ref_files.empty())
        for (int l = 0; l < 3; l++)
        {
            tr_emn[l] = new TH2D(TString::Format("h_trgref_em%d_n", l), TString::Format("EM%d reference events;E_{4x4} / %.1f GeV;r = 4 H_{max} / E_{4x4}", l, kTrgThresholdGeV[l]),
                                 kTrgRefNE, 0, kTrgRefNE * kTrgRefDE, kTrgRefNR, 0, kTrgRefNR * kTrgRefDR);
            tr_emk[l] = (TH2D *)tr_emn[l]->Clone(TString::Format("h_trgref_em%d_k", l));
            tr_emk[l]->SetTitle(TString::Format("EM%d reference events, bit set", l));
            tr_jpn[l] = new TH1D(TString::Format("h_trgref_jp%d_n", l), TString::Format("JP%d reference events;patch / %.1f GeV;events", l, kTrgThresholdGeV[3 + l]),
                                 kTrgRefNE, 0, kTrgRefNE * kTrgRefDE);
            tr_jpk[l] = (TH1D *)tr_jpn[l]->Clone(TString::Format("h_trgref_jp%d_k", l));
        }
    long n_ref_files_read = 0;

    std::set<int> runs;
    Pi0Scale pi0sc;
    pi0sc.Load();
    std::map<int, double> pi0_used;
    std::cout << "AnHistograms: pi0 energy scale from " << pi0sc.source << std::endl;
    int tree_number = -1, file_run = 0;
    bool is_ref = false;
    const Long64_t n = t.GetEntries();
    for (Long64_t i = 0; i < n; i++)
    {
        t.GetEntry(i);
        // the run: event_run, else from the file name (jet_output_<run>_<file>.root)
        if (t.GetTreeNumber() != tree_number)
        {
            file_run = RunFromFileName(t.GetFile()->GetName());
            const std::string fname = t.GetFile()->GetName();
            is_ref = ref_files.count(fname.substr(fname.find_last_of('/') + 1)) > 0;
            if (is_ref) n_ref_files_read++;
        }
        tree_number = t.GetTreeNumber();
        const int cur_run = have_run && event_run > 0 ? event_run : file_run;
        const double pi0_scale = pi0sc.Of(cur_run);
        pi0_used[cur_run] = pi0_scale;
        runs.insert(cur_run);
        h_cut->Fill(kNCuts);
        // trigger reference tallies: reference events of the reference file, no other cut
        // (trg_nominal.cpp: leading E_4x4 jet with an HCal sum and its r; largest patch)
        if (is_ref)
        {
            bool eref = jp2 || jp1[0] || jp1[1] || jp1[2] || jp0[0] || jp0[1] || jp0[2];
            for (int k = 7; k < n_ref_em; k++) eref |= ref_em[k] != 0;
            bool jref = em[0] || em[1] || em_tpc[0] || em_tpc[1];
            for (int k = 4; k < n_ref_jp; k++) jref |= ref_jp[k] != 0;
            double lE = -1, lr = -1, ljp = -1;
            for (size_t j = 0; j < trg_em->size(); j++)
            {
                if ((*trg_em)[j] > lE && (*trg_hmax)[j] >= 0) { lE = (*trg_em)[j]; lr = 4 * (*trg_hmax)[j] / (*trg_em)[j]; }
                ljp = std::max(ljp, (double)(*trg_jp)[j]);
            }
            const double emax = kTrgRefNE * kTrgRefDE - 1e-6, rmax = kTrgRefNR * kTrgRefDR - 1e-6;
            if (eref && lE > 0)
                for (int l = 0; l < 3; l++)
                {
                    const double xe = std::min(lE / kTrgThresholdGeV[l], emax), xr = std::min(lr, rmax);
                    tr_emn[l]->Fill(xe, xr);
                    if ((tcu >> (3 + l)) & 1) tr_emk[l]->Fill(xe, xr);
                }
            if (jref && cur_run != 23070019 && ljp > 0)  // 23070019: hot JP run
                for (int l = 0; l < 3; l++)
                {
                    const int tr = 3 + l;
                    const bool hw = tr == 5 ? ((tcu >> 6) & 1) : (((tcu >> (tr == 3 ? 10 : 7)) & 7) != 0);
                    const double xj = std::min(ljp / kTrgThresholdGeV[tr], emax);
                    tr_jpn[l]->Fill(xj);
                    if (hw) tr_jpk[l]->Fill(xj);
                }
        }
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
        {
            // N-1 distributions of the event-level cut variables
            const bool vz_ok = !vertex_cut || pass_vertex_cut(vz, vz_source);
            const bool bbc_ok = bbc_e_ntiles < 0 || bbc_e_ntiles >= kBbcEastMinTiles;
            const bool tof_ok = tof_mult < 0 || tof_mult >= kTofMultMin;
            if (bbc_ok && tof_ok)
            {
                q_cut_vzsrc->Fill(vz_source >= 0 && vz_source < 4 ? vz_source : 0);
                if ((vz_source == 1 || vz_source == 2) && vz > -9000) q_cut_vz->Fill(vz);
            }
            if (vz_ok && tof_ok && bbc_e_ntiles >= 0) q_cut_bbc->Fill(bbc_e_ntiles);
            if (vz_ok && bbc_ok && tof_mult >= 0) q_cut_tof->Fill(tof_mult);
        }
        if (vertex_cut && !pass_vertex_cut(vz, vz_source)) continue;
        h_cut->Fill(kVertex, nj);
        // collision cut: failing events fill only the _vbbcveto copies
        const bool evt_nominal = PassCollisionCut(bbc_e_ntiles, tof_mult);
        bool var_evt[kNVar];
        var_evt[kVarFid25] = var_evt[kVarFid30] = evt_nominal;
        var_evt[kVarVz60] = evt_nominal && (vz_source <= 0 || std::fabs(vz) < kVarVzCut);
        var_evt[kVarVtxMeas] = evt_nominal && (vz_source == 1 || vz_source == 2);
        var_evt[kVarBbcVeto] = !evt_nominal;
        var_evt[kVarJesUp] = var_evt[kVarJesDn] = evt_nominal && have_jes_err;
        // point-threshold copies: filled only in their own call (--thr-copy), which fills nothing else
        var_evt[kVarThr150] = var_evt[kVarThr200] = false;
        if (thr_v >= 0)
        {
            for (int v = 0; v < kNVar; v++) var_evt[v] = false;
            var_evt[thr_v] = evt_nominal;
        }
        if (evt_nominal)
        {
            h_cut->Fill(kCollision, nj);
            h_events[level]->Fill(spin);
            q_vzsrc->Fill(vz_source >= 0 && vz_source < 4 ? vz_source : 0);
            if (vz > -9000) q_vz->Fill(vz);
            q_sqrts->Fill(sqrt_s);
            qa[0].njet->Fill(nj);
            qa_lvl[0][level].njet->Fill(nj);
        }
        // the two constituents of each two-photon jet: pair mass and orientation (-1 if none)
        std::vector<double> mass(nj, -1), pang(nj, -1);
        // every jet's constituents (for the FMS-style pi0 pairs)
        std::vector<std::vector<int> > cons(nj);
        for (size_t c = 0; c < c_jet->size(); c++)
        {
            const int jj = (*c_jet)[c];
            if (jj >= 0 && jj < (int)nj) cons[jj].push_back(c);
        }
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
        // baseJ[0 / 1]: the same with the energy scaled by 1 + / - delta(E) (the
        // _vjesup / _vjesdn copies; E, pT and x_F scale alike at fixed eta);
        // a jet passing only there gets its trigger results too, nothing nominal
        std::vector<int> base(nj, 0), emL(nj, -1), emS(nj, -1), emVeto(nj, 0), jpL(nj, -1), jpS(nj, -1);
        std::vector<int> baseJ[2] = {std::vector<int>(nj, 0), std::vector<int>(nj, 0)};
        std::vector<double> fJ[2] = {std::vector<double>(nj, 1.0), std::vector<double>(nj, 1.0)};
        int event_em = -1, event_jp = -1;  // highest EM / JP level among jets at >= the default factor x threshold
        for (size_t j = 0; j < nj; j++)
        {
            const double E = (*e_corr)[j];
            // sanity distributions of every jet here, by the event's level; the selected jets' below
            const double ptc0 = E > 0 ? E / std::cosh((*eta)[j]) : -1;
            if (evt_nominal)
            {
                qa[0].Fill((*eta)[j], (*phi)[j], (*e_raw)[j], (*e_ue)[j], E, (*pt_raw)[j], ptc0, (*xf_corr)[j], (*nconst)[j], (*x)[j], (*y)[j]);
                qa_lvl[0][level].Fill((*eta)[j], (*phi)[j], (*e_raw)[j], (*e_ue)[j], E, (*pt_raw)[j], ptc0, (*xf_corr)[j], (*nconst)[j], (*x)[j], (*y)[j]);
                qa[0].FillPair(mass[j], pang[j], E);
                qa_lvl[0][level].FillPair(mass[j], pang[j], E);
            }
            if (var_evt[kVarJesUp] && E > 0)
            {
                // delta at the map's lookup energy (UE-subtracted reco E, held as in
                // ApplyCorrections_em_only), so the copy equals the calibration rerun
                // with the shifted map (ApplyCorrections_em_only --jes-shift=+1 / -1)
                const double e_map = std::min(std::max((double)(*e_raw)[j] - std::max((double)(*e_ue)[j], 0.0), kJesEvalMinGeV), 0.5 * sqrt_s);
                const double dl = jes_err.at(e_map);
                for (int u = 0; u < 2; u++)
                {
                    fJ[u][j] = u == 0 ? 1 + dl : 1 - dl;
                    baseJ[u][j] = fJ[u][j] > 0 && (*xf_corr)[j] > 0 && fJ[u][j] * E / std::cosh((*eta)[j]) > kJetPtCorrMin;
                }
            }
            // the nominal chain (cutflow); a jet failing it goes on only for baseJ
            bool b1 = E > 0;
            if (b1 && evt_nominal) h_cut->Fill(kCalibrated);
            b1 = b1 && (*xf_corr)[j] > 0;
            if (b1 && evt_nominal) h_cut->Fill(kXfPositive);
            const double pt = E > 0 ? E / std::cosh((*eta)[j]) : -1;
            b1 = b1 && pt > kJetPtCorrMin;
            if (b1 && evt_nominal)
            {
                h_cut->Fill(kPtMin);
                h_pt[level]->Fill(pt);
                h_ptraw[level]->Fill((*pt_raw)[j]);
            }
            base[j] = b1;
            if (!b1 && !baseJ[0][j] && !baseJ[1][j]) continue;
            // EM: the lowest listed EM level at >= 1.15 x threshold; its slice from the
            // combined EM factor f = min(E_4x4 / threshold, veto edge / r) (AnConfig.h
            // EmTrgFactor, 2026-11-17: threshold and veto move together as one cut)
            const double em_et = (*trg_em)[j];
            for (int l = 0; l < kNEmLevels && emL[j] < 0; l++)
                if (listed[l] && em_et >= kTrgSliceEdge[0] * TrgThresholdGeV(l, cur_run)) emL[j] = l;
            if (emL[j] >= 0)
            {
                const double ratio = em_et / TrgThresholdGeV(emL[j], cur_run);
                if (evt_nominal && b1) h_ratio[emL[j]]->Fill(ratio);
                const double fem = EmTrgFactor(em_et, (*trg_hmax)[j], emL[j], cur_run);
                if (fem >= kTrgSliceEdge[0] - 1e-9)
                {
                    emS[j] = 0;
                    for (int sl = 1; sl < kNTrgSlices; sl++) if (fem >= kTrgSliceEdge[sl] - 1e-9) emS[j] = sl;
                    if (b1 && fem >= kEmCutFactor - 1e-9 && emL[j] > event_em) event_em = emL[j];
                }
                else emVeto[j] = 1;
            }
            // JP: the lowest listed JP level at >= 1.15 x threshold (JP0/JP1 in a listed region)
            const double jp_et = (*trg_jp)[j];
            for (int l = kNEmLevels; l < kNTrgLevels && jpL[j] < 0; l++)
            {
                if (!listed[l]) continue;
                if (l < 5 && !((*trg_jpmask)[j] & jp_region[l - kNEmLevels])) continue;
                if (jp_et >= kTrgSliceEdge[0] * TrgThresholdGeV(l, cur_run)) jpL[j] = l;
            }
            if (jpL[j] >= 0)
            {
                const double ratio = jp_et / TrgThresholdGeV(jpL[j], cur_run);
                if (evt_nominal && b1) h_ratio[jpL[j]]->Fill(ratio);
                jpS[j] = 0;
                for (int sl = 1; sl < kNTrgSlices; sl++) if (ratio >= kTrgSliceEdge[sl]) jpS[j] = sl;
                if (b1 && ratio >= kJpCutFactor - 1e-9 && jpL[j] > event_jp) event_jp = jpL[j];
            }
        }
        // companions go to the highest EM level of the event's default jets, else the highest JP level
        const int event_jl = event_em >= 0 ? event_em : event_jp;
        int nsel = 0, ncomp = 0;
        for (size_t j = 0; j < nj; j++)
        {
            if (!base[j] && !baseJ[0][j] && !baseJ[1][j]) continue;
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
            // at the trigger-systematic factor (the "_f130" copies of the default-only histograms)
            const bool em_sys = em_on && kTrgSliceEdge[emS[j]] >= kTrgSysFactor - 1e-9;
            const bool jp_sys = jp_on && kTrgSliceEdge[jpS[j]] >= kTrgSysFactor - 1e-9;
            // the JES copies: this jet with E, pT, x_F x (1 +- delta(E)) (default
            // trigger factor, E <= sqrt(s)/2 with the scaled energy)
            for (int u = 0; u < 2; u++)
            {
                if (!baseJ[u][j] || !(em_dflt || jp_dflt) || fJ[u][j] * E > 0.5 * sqrt_s) continue;
                const int vj = u == 0 ? kVarJesUp : kVarJesDn;
                const double xfj = fJ[u][j] * xf;
                auto jfill = [&](int tg, int l)
                {
                    hv[vj][tg][l][cls][side][spin]->Fill(xfj);
                    hv_xfs[vj][tg][l][cls]->Fill(xfj, xfj);
                    hv_cos[vj][tg][l][cls]->Fill(xfj, acos_);
                };
                if (em_dflt) jfill(0, emL[j]);
                if (jp_dflt) jfill(em_dflt ? 1 : 0, jpL[j]);
            }
            if (!base[j]) continue;  // the rest is nominal
            // FMS-style pi0 candidates of a default jet (E <= sqrt(s)/2): every photon pair, built once
            // for the nominal and the variation yields
            struct PiPair { double pxf, mgg, zem, pcos, ppt_raw; int pside; };
            std::vector<PiPair> pairs;
            if ((em_dflt || jp_dflt) && !hi && (*e_raw)[j] > 0)
            {
                const std::vector<int> &cj = cons[j];
                // photon at the ECal face (z_proj), from the vertex the jets were made from;
                // inside the ECal rectangle by half a tower
                std::vector<char> fid(cj.size(), 0);
                for (size_t a = 0; a < cj.size(); a++)
                {
                    const int c = cj[a];
                    if (!((*c_pz)[c] > 0)) continue;
                    const double sc = (z_proj - vz_used) / (*c_pz)[c];
                    const double px_ = vx_used + sc * (*c_px)[c], py_ = vy_used + sc * (*c_py)[c];
                    const FiducialRect &rr = kFiducialRectEcal;
                    fid[a] = std::fabs(px_) >= rr.x_inner + kPi0EdgeMargin && std::fabs(px_) <= rr.x_outer - kPi0EdgeMargin &&
                             py_ >= rr.y_min + kPi0EdgeMargin && py_ <= rr.y_max - kPi0EdgeMargin;
                }
                for (size_t a = 0; a < cj.size(); a++)
                    for (size_t b = a + 1; b < cj.size(); b++)
                    {
                        if (!fid[a] || !fid[b]) continue;
                        const int ca = cj[a], cb = cj[b];
                        const double ea = (*c_e)[ca], eb = (*c_e)[cb];
                        if (!(ea > 0 && eb > 0) || std::fabs(ea - eb) / (ea + eb) >= kPi0ZggMax) continue;
                        const double ppx = (*c_px)[ca] + (*c_px)[cb], ppy = (*c_py)[ca] + (*c_py)[cb], ppz = (*c_pz)[ca] + (*c_pz)[cb], pe = ea + eb;
                        const double m2 = pe * pe - ppx * ppx - ppy * ppy - ppz * ppz;
                        // mass, pT and x_F on the pi0 energy scale (AnConfig.h
                        // kPi0EnergyScale or the pi0_scale.txt table: the data peak moved to the pi0 mass); raw
                        // pT for the trigger cut (the trigger value is uncorrected)
                        const double mgg = pi0_scale * (m2 > 0 ? std::sqrt(m2) : 0);
                        if (mgg >= kPi0MassMax) continue;
                        const double ppt_raw = std::hypot(ppx, ppy), ppt = pi0_scale * ppt_raw;
                        if (ppt <= kPi0PtMin) continue;
                        const double pxf = 2 * pi0_scale * ppz / sqrt_s;
                        if (!(pxf > 0)) continue;
                        PiPair pp;
                        pp.pxf = pxf; pp.mgg = mgg; pp.ppt_raw = ppt_raw; pp.pcos = std::fabs(ppx) / ppt;
                        // raw over raw: one scale for all points cancels
                        pp.zem = pe / (*e_raw)[j];
                        pp.pside = computeSide(vx_used + (z_proj - vz_used) * ppx / ppz) > 0 ? 1 : 0;
                        pairs.push_back(pp);
                    }
            }
            // the pi0 itself above the level's cut (FMS: pi0 pT above the trigger
            // threshold; Sean 2026-11-16: factor x threshold, as the host jet's
            // trigger value); the union (_jponly) is per pi0: JP pi0s not passing in EM
            auto pi0Ok = [&](const PiPair &pp, bool &em_ok, bool &jp_ok)
            {
                em_ok = em_dflt && pp.ppt_raw >= kEmCutFactor * TrgThresholdGeV(emL[j], cur_run);
                jp_ok = jp_dflt && pp.ppt_raw >= kJpCutFactor * TrgThresholdGeV(jpL[j], cur_run);
            };
            // systematic-variation copies: default jets, E <= sqrt(s)/2 (AnConfig.h kVarSuffix)
            if ((em_dflt || jp_dflt) && !hi)
            {
                // distance of the jet axis from the ECal edge at z_proj (nominal >= 20 cm, RecoJets)
                const FiducialRect &fr = kFiducialRectEcal;
                const double ax_ = std::fabs((*x)[j]), ay_ = (*y)[j];
                const double edge = std::min(std::min(ax_ - fr.x_inner, fr.x_outer - ax_), std::min(ay_ - fr.y_min, fr.y_max - ay_));
                for (int v = 0; v < kNVar; v++)
                {
                    if (!var_evt[v] || IsJesVar(v)) continue;  // _vjesup / _vjesdn: filled above
                    if (v == kVarFid25 && edge < kVarFidCm[0]) continue;
                    if (v == kVarFid30 && edge < kVarFidCm[1]) continue;
                    auto vfill = [&](int tg, int l)
                    {
                        hv[v][tg][l][cls][side][spin]->Fill(xf);
                        hv_xfs[v][tg][l][cls]->Fill(xf, xf);
                        hv_cos[v][tg][l][cls]->Fill(xf, acos_);
                    };
                    if (em_dflt) vfill(0, emL[j]);
                    if (jp_dflt) vfill(em_dflt ? 1 : 0, jpL[j]);
                    for (size_t q = 0; q < pairs.size(); q++)
                    {
                        const PiPair &pp = pairs[q];
                        bool em_ok, jp_ok;
                        pi0Ok(pp, em_ok, jp_ok);
                        auto vpi = [&](int tg, int l)
                        {
                            hv_pi[v][tg][l][pp.pside][spin]->Fill(pp.pxf, pp.mgg, pp.zem);
                            hv_picos[v][tg][l]->Fill(pp.pxf, pp.mgg, pp.zem, pp.pcos);
                            hv_pixfs[v][tg][l]->Fill(pp.pxf, pp.mgg, pp.zem, pp.pxf);
                        };
                        if (em_ok) vpi(0, emL[j]);
                        if (jp_ok) vpi(0, jpL[j]);
                        if (jp_ok && !em_ok) vpi(1, jpL[j]);
                    }
                }
            }
            if (!evt_nominal) continue;  // collision-cut failures: only the _vbbcveto copies
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
                if (em_sys) h_m[4][emL[j]][side][spin]->Fill(xf, mass[j]);
                if (jp_sys) h_m[4][jpL[j]][side][spin]->Fill(xf, mass[j]);
                if (jp_sys && !em_sys) h_m[5][jpL[j]][side][spin]->Fill(xf, mass[j]);
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
                if (em_sys) pfill(4, emL[j]);
                if (jp_sys) pfill(4, jpL[j]);
                if (jp_sys && !em_sys) pfill(5, jpL[j]);
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
                    if (em_sys) zfill(2, emL[j]);
                    if (jp_sys) zfill(2, jpL[j]);
                    if (jp_sys && !em_sys) zfill(3, jpL[j]);
                }
            }
            // FMS-style pi0 candidates: every photon pair of a default jet (and _jponly)
            for (size_t q = 0; q < pairs.size(); q++)
            {
                const PiPair &pp = pairs[q];
                auto pifill = [&](int mg, int l)
                {
                    h_pi[mg][l][pp.pside][spin]->Fill(pp.pxf, pp.mgg, pp.zem);
                    h_picos[mg][l]->Fill(pp.pxf, pp.mgg, pp.zem, pp.pcos);
                    h_pixfs[mg][l]->Fill(pp.pxf, pp.mgg, pp.zem, pp.pxf);
                };
                bool em_ok, jp_ok;
                pi0Ok(pp, em_ok, jp_ok);
                if (em_ok) pifill(0, emL[j]);
                if (jp_ok) pifill(0, jpL[j]);
                if (jp_ok && !em_ok) pifill(1, jpL[j]);
                // the same at the trigger-systematic factor (host jet and pi0 pT)
                const bool em_oks = em_sys && pp.ppt_raw >= kTrgSysFactor * TrgThresholdGeV(emL[j], cur_run);
                const bool jp_oks = jp_sys && pp.ppt_raw >= kTrgSysFactor * TrgThresholdGeV(jpL[j], cur_run);
                if (em_oks) pifill(2, emL[j]);
                if (jp_oks) pifill(2, jpL[j]);
                if (jp_oks && !em_oks) pifill(3, jpL[j]);
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
                if (em_sys) collins.Fill(2, emL[j], side, spin, xf, gg, zz, plp);
                if (jp_sys) collins.Fill(2, jpL[j], side, spin, xf, gg, zz, plp);
                if (jp_sys && !em_sys) collins.Fill(3, jpL[j], side, spin, xf, gg, zz, plp);
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
    // only the histograms of this call's copies are written: the nominal call
    // leaves out the _vthr copies (filled by their own calls), a --thr-copy call
    // writes only its copy (and no metadata but thr_copy_<thr>)
    {
        std::vector<TObject *> drop;
        TIter it(out.GetList());
        while (TObject *o = it())
        {
            const TString nm = o->GetName();
            bool keep;
            if (thr_v >= 0) keep = nm.EndsWith(kVarSuffix[thr_v]);
            else keep = !nm.EndsWith(kVarSuffix[kVarThr150]) && !nm.EndsWith(kVarSuffix[kVarThr200]);
            if (!keep) drop.push_back(o);
        }
        // detached, not deleted (some are owned by QaSet / CollinsHistograms)
        for (size_t i = 0; i < drop.size(); i++)
        {
            if (drop[i]->InheritsFrom(TH1::Class())) ((TH1 *)drop[i])->SetDirectory(0);
            else out.GetList()->Remove(drop[i]);
        }
    }
    if (thr_v >= 0)
    {
        const char *const thr_name = thr_v == kVarThr150 ? "150" : "200";
        const char *env = std::getenv("FCS_POINT_EMIN");
        TNamed(TString::Format("thr_copy_%s", thr_name).Data(), TString::Format("%s copy: %lld events, runs %s; jets made with FCS_POINT_EMIN=%s", kVarSuffix[thr_v], n, rs.str().c_str(), env ? env : "(unset in this call)").Data()).Write();
        out.Write();
        printf("AnHistograms --thr-copy=%s: %lld events, runs %s\n", thr_name, n, rs.str().c_str());
        for (int c = 0; c < kNCuts; c++) printf("  %-16s %.0f jets\n", kCutName[c], h_cut->GetBinContent(c + 1));
        out.Close();
        return 0;
    }
    TNamed("runs", rs.str().c_str()).Write();
    TNamed("spin_convention", kSpinConventionTag).Write();
    {
        // the pi0 energy scale used: "<source>; <run>=<factor> ..." (pi0_peak_runs.C reads it)
        std::ostringstream ps;
        ps << pi0sc.source << ";";
        for (std::map<int, double>::const_iterator r = pi0_used.begin(); r != pi0_used.end(); ++r) ps << " " << r->first << "=" << TString::Format("%.6f", r->second).Data();
        TNamed("pi0_scale", ps.str().c_str()).Write();
    }
    if (!ref_files.empty()) TNamed("trgref_files", TString::Format("%ld of %zu listed trigger reference files read", n_ref_files_read, ref_files.size()).Data()).Write();
    TNamed("variations", "systematic-variation copies (AnConfig.h kVarSuffix): _vfid25 _vfid30 _vvz60 _vvtxmeas _vbbcveto _vjesup _vjesdn (JES: photon classes only) _vthr150 _vthr200 (added by --thr-copy calls: see thr_copy_*); collision cut east BBC >= 1, TOF > 2").Write();
    if (have_jes_err)
    {
        std::ostringstream je;
        je << jes_err_path << ": delta(E)";
        for (size_t i = 0; i < jes_err.e.size(); i++) je << " " << jes_err.e[i] << "=" << jes_err.d[i];
        TNamed("jes_error", je.str().c_str()).Write();
    }
    else TNamed("jes_error_missing", "no JES map error (--jes-error): _vjesup / _vjesdn copies empty").Write();
    // cut positions for plot_event_cuts.C: "<hist>:<nominal edges>:<variation edges>[:<variation label>];..."
    TNamed("event_cut_lines", TString::Format("qa_cut_vz:%g,%g:%g,%g:|v_{z}| < %g cm;qa_cut_bbc_e_ntiles:%g::;qa_cut_tof_mult:%g::",
                                              -kVertexZCut, kVertexZCut, -kVarVzCut, kVarVzCut, kVarVzCut, kBbcEastMinTiles - 0.5, kTofMultMin - 0.5).Data()).Write();
    if (!have_plead) TNamed("plead_missing", "jet files without reco_pi0_lead_prob: h_xFP_* empty").Write();
    if (!have_gup) TNamed("collins_missing", "jet files without reco_pi0_lead_gup: h_xFg_* empty").Write();
    if (!have_z) TNamed("pz_missing", "jet files without reco_pi0_lead_z1: h_xFPz_* empty").Write();
    TNamed("trigger_rule", TString::Format("trigger values (reco_trg_*) / measured per-period thresholds (x50 x nominal, AnConfig.h kTrgX50), slices from %.2f, default factor EM %.2f JP %.2f, EM factor = min(E/threshold, veto edge/r), veto edges per period %.3f / %.3f / %.3f", kTrgSliceEdge[0], kEmCutFactor, kJpCutFactor, kEmVetoEdge[0], kEmVetoEdge[1], kEmVetoEdge[2]).Data()).Write();
    out.Write();
    printf("AnHistograms: %lld events, runs %s\n", n, rs.str().c_str());
    for (int c = 0; c < kNCuts; c++) printf("  %-16s %.0f jets\n", kCutName[c], h_cut->GetBinContent(c + 1));
    out.Close();
    return 0;
}
