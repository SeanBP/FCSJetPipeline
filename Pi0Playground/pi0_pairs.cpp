// Two-point (pi0) invariant mass from FCS ECal points in SimpleTrees, data or
// simulation, to compare the point energy scale and two-photon separation of
// data and simulation.
//
// Events: no FCS LED and (unless "all") an EM trigger fired (fcsEM0/1/2;
// simulation: emulated flags). Points: ECal (detid 0, 1), E > point_energy_min,
// inside the ECal rectangle (kFiducialRectEcal) shrunk by half a tower
// (kPointEdgeBuffer), as the FMS photon fiducial cut. Pairs: both points in the
// same detector half; the mass uses directions from the origin (no vertex).
// The FMS pi0 selection keeps Z = |E1 - E2| / (E1 + E2) < 0.7.
//
// Output (pair energy E = E1 + E2 bins kEEdges, point distance d at the ECal
// bins kDEdges, asymmetry bins kZEdges):
//   h3_m_E_d          mass, Z < 0.7, against E bin and d bin
//   h3_m_E_Z          mass against E bin and Z bin (all Z)
//   h3_m_E_c          mass, Z < 0.7, against E bin and cluster class: 0 = both
//                     points from one cluster (split by the shower-shape fit),
//                     1 = from two clusters (Cal_pt_parentClusterId)
//   h3_m_E_tail       as h3_m_E_d (E bin and d bin of the plain energies),
//                     with each point's energy plus the shower tail the
//                     reclustering threshold dropped: stored ECal hits with
//                     kTailMin < E < kTowerThr, each added to the nearest
//                     selected point within kTailCm (same detector)
//   h2_tail_E         point E against that added fraction (tail / point E)
//   h2_liang_m_E      the ECal gain calibration's own selection
//                     (StFcsPi0FinderForEcal, X. Liang): per event the
//                     highest-energy cluster pair of one detector half,
//                     both clusters above 1 GeV, Z < 0.7, events with at
//                     most 40 ECal hits above 1 GeV; cluster positions;
//                     mass against pair energy
//   h3_liang_m_E_ring_<class>
//                     the same against pair energy and ring of the pair's
//                     energy-weighted radius (kRingEdges, 0 = inner), per
//                     trigger class (kClsName)
//   h3_liangf_m_E_ring_<class>
//                     as h3_liang_m_E_ring_<class> with a common 45 MeV tower
//                     floor: each cluster loses its 10-45 MeV hits within
//                     kTailCm (simulation keeps ECal hits from 2 ADC counts,
//                     ~13 MeV; production data from ~45 MeV)
// Overlap study (calibration selection, classes "all" and "em", 45 MeV floor
// applied as in h3_liangf; tower = 5.542 cm):
//   h3_lsep_m_E_sep_<all|em>    mass, pair E, cluster separation in towers
//   h3_lph_m_Eg_Z_<all|em>      mass, single-photon (cluster) energy -- filled
//                               once per photon of the pair -- and Z
//   h2_leta_m_Eg_<all|em>       mass up to 0.8 GeV (eta) against the lower
//                               photon energy of the pair
//   h2_clwidth_E                energy-weighted RMS distance of a cluster's
//                               hits (within kTailCm) from its centre (cm),
//                               against cluster energy; both pair clusters
//   h_isohit_E                  energy of isolated single ECal towers (no
//                               ECal hit above 10 MeV within 1.5 towers),
//                               events with no ECal cluster above 1 GeV
//                               within 3 towers: the MIP peak
//   p_liang_sumf2_E   mean sum of squared tower energy fractions of the
//                     pair's clusters, against pair energy
//   h_maxclet_<all|em0|em1|em2>  highest ECal cluster E_T per event, all
//                     events and those with each EM bit (trigger turn-on)
//   h3_mclu_E_d       as h3_m_E_d with clusters in place of points (every
//                     pair, same cuts; cluster positions Cal_clus_x/y/z)
// The mass range reaches the eta (0.548 GeV), whose photons are about four
// times further apart than a pi0's of the same energy.
//   h_ptE             point energy spectrum (selected points)
//   h_npts            selected points per event
// data only, with an event vertex (event_vertex_source > 0):
//   h3_m_E_d_vtx      as h3_m_E_d, directions from the vertex
// simulation only (truth from mcpart_*): each point is matched to the nearest
// GEANT photon (geid 1, E > 0.2 GeV) at the point's z, within kMatchCm.
// Pairs whose points match the two photons of one pi0 (geid 7, its end
// vertex = their start vertex) are "tagged":
//   h3_m_E_d_tag      mass of tagged pairs
//   h3_mt_E_d_tag     the same with the true photon energies at the reco
//                     positions (position effect only)
//   h3_mp_E_d_tag     reco energies at the true positions (energy effect only)
//   h2_resp_E         point E / photon E against photon E (matched points)
//   h2_dpos_E         reco - true position at the ECal (cm) against photon E
//   h2_dsep_E_d       tagged pairs: reco d - true d (cm) against E and d bin
//                     (TH3: E bin, true d bin, difference)
// Simulation (files ..._pt<ptHat>_<job>.root) is weighted by sigma / events of
// its ptHat bin, absolute (mb per event): a job must hold every file of each
// ptHat bin it reads, so that outputs add with hadd.
//
// "simuw" is simulation with every event at weight 1 (more statistics for
// detector-level peak positions in narrow bins; the truth histograms too).
//
// With a 5th argument <dir>, the EM bits fcsEM0-2 are read from the friend
// tree "trgemu" in <dir>/<input basename> (trigemu/reemulate.C) instead of
// the SimpleTree's own ("-" for none). A 6th argument k scales every ECal
// (det 0/1) hit, cluster and point energy at read-in, standing in for an
// ECal energy scale applied before clustering (exact enough with the default
// clustering's 0.01 / 0.05 GeV thresholds; the selection cuts see the scaled
// energies).
//
// Usage: ./pi0_pairs <SimpleTree list> <out.root> <data|sim|simuw> [trig|all] [<trgemu dir> | -] [k = 1]
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TVector3.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "JetParameters.h"
#include "JetXSecWeights.h"

static const int kNE = 7;
static const double kEEdges[kNE + 1] = {5, 10, 15, 20, 30, 45, 70, 100};
static const int kND = 7;
static const double kDEdges[kND + 1] = {0, 6, 9, 12, 16, 22, 30, 60};
static const int kNZ = 5;
static const double kZEdges[kNZ + 1] = {0, 0.2, 0.4, 0.55, 0.7, 1.0};
static const double kZMax = 0.7;
static const double kPointEdgeBuffer = 0.5 * 5.542; // half an ECal tower (cm)
static const double kMatchCm = 3.0;
static const double kTowerThr = 0.4, kTailMin = 0.05, kTailCm = 8.0; // GeV, GeV, cm
static const int kMaxHits = 4000;
// Trigger classes of the calibration-selection histograms (run with "all" to
// keep untriggered events): every event, any EM trigger, no EM trigger, fcsEM2,
// EM-triggered with the pair inside / outside the event's highest-E_T 4x4
// tower window (an offline stand-in for the firing patch: E_T summed over
// hits within kPatchHalf of each hit, the best window's centre; inside = the
// pair's energy-weighted position within kPatchHalf of that centre)
// "emet<N>": any EM bit and the event's highest ECal cluster E_T above N GeV
// (an offline trigger condition applied alike to data and simulation; removes
// data's EM-triggered events with no hard ECal cluster)
static const int kNCls = 9;
static const char *kClsName[kNCls] = {"all", "em", "noem", "em2", "inpatch", "outpatch", "emet2", "emet3", "emet4"};
static const double kPatchHalf = 2 * 5.542; // cm, half of a 4x4 tower window
static const int kNRing = 3;
static const double kRingEdges[kNRing + 1] = {0, 66.35, 87.53, 1e9}; // cm, thirds of the ECal radius (as jes_closure)
static const int kNM = 160;
static const double kMMax = 0.8;
static const int kMaxPts = 200, kMaxMc = 10000;

static int Bin(double v, const double *edges, int n)
{
    for (int k = 0; k < n; k++) if (v >= edges[k] && v < edges[k + 1]) return k;
    return -1;
}

static int GroupOf(const std::string &f)
{
    const size_t p = f.rfind("_pt");
    return p == std::string::npos ? -1 : std::atoi(f.c_str() + p + 3);
}

static double Mass(double e1, const TVector3 &p1, double e2, const TVector3 &p2)
{
    const double c = p1.Unit().Dot(p2.Unit());
    const double m2 = 2 * e1 * e2 * (1 - c);
    return m2 > 0 ? std::sqrt(m2) : 0;
}

int main(int argc, char **argv)
{
    if (argc < 4) { std::cerr << "usage: ./pi0_pairs <SimpleTree list> <out.root> <data|sim> [trig|all]" << std::endl; return 1; }
    const bool unweighted = std::string(argv[3]) == "simuw";
    const bool sim = std::string(argv[3]) == "sim" || unweighted;
    const bool need_trig = !(argc > 4 && std::string(argv[4]) == "all");
    const std::string emu_dir = argc > 5 && std::string(argv[5]) != "-" ? argv[5] : "";
    const double ecal_k = argc > 6 ? std::atof(argv[6]) : 1.0;
    if (ecal_k != 1.0) std::printf("pi0_pairs: ECal energies scaled by k = %.3f\n", ecal_k);
    std::vector<std::string> files;
    {
        std::ifstream in(argv[1]);
        std::string f;
        while (in >> f) files.push_back(f);
    }

    std::map<std::string, int> file_group;
    std::map<int, double> group_N, group_sigmaN;
    if (sim)
    {
        for (size_t i = 0; i < files.size(); i++)
        {
            TFile *tf = TFile::Open(files[i].c_str());
            TTree *t = tf && !tf->IsZombie() ? (TTree *)tf->Get("data") : 0;
            if (t && t->GetEntries() > 0)
            {
                Double_t sg = -1; Int_t n = -1;
                t->SetBranchStatus("*", 0); t->SetBranchStatus("sigma_job", 1); t->SetBranchStatus("N_job", 1);
                t->SetBranchAddress("sigma_job", &sg); t->SetBranchAddress("N_job", &n);
                t->GetEntry(0);
                const double s = DecodePythia6Sigma(sg);
                const int g = GroupOf(files[i]);
                if (s > 0 && n > 0 && g >= 0) { file_group[files[i]] = g; group_N[g] += n; group_sigmaN[g] += n * s; }
            }
            delete tf;
        }
        for (std::map<int, double>::iterator it = group_N.begin(); it != group_N.end(); ++it)
            std::printf("ptHat %3d: %.0f events, sigma %.4g mb\n", it->first, it->second, group_sigmaN[it->first] / it->second);
    }

    TFile out(argv[2], "RECREATE");
    TH1::SetDefaultSumw2(kTRUE);
    TH3D *h3 = new TH3D("h3_m_E_d", ";m_{#gamma#gamma} (GeV);E bin;d bin", kNM, 0, kMMax, kNE, 0, kNE, kND, 0, kND);
    TH3D *h3z = new TH3D("h3_m_E_Z", ";m_{#gamma#gamma} (GeV);E bin;Z bin", kNM, 0, kMMax, kNE, 0, kNE, kNZ, 0, kNZ);
    TH3D *h3c = new TH3D("h3_m_E_c", ";m_{#gamma#gamma} (GeV);E bin;cluster class", kNM, 0, kMMax, kNE, 0, kNE, 2, 0, 2);
    TH3D *h3tail = new TH3D("h3_m_E_tail", ";m_{#gamma#gamma} (GeV);E bin;d bin", kNM, 0, kMMax, kNE, 0, kNE, kND, 0, kND);
    TH2D *h_tail = new TH2D("h2_tail_E", ";point E (GeV);dropped tail / point E", 40, 0, 40, 100, 0, 0.5);
    TH3D *h3clu = new TH3D("h3_mclu_E_d", ";m_{#gamma#gamma} (GeV);E bin;d bin", kNM, 0, kMMax, kNE, 0, kNE, kND, 0, kND);
    TH2D *h_liang = new TH2D("h2_liang_m_E", ";m_{#gamma#gamma} (GeV);cluster pair E (GeV)", 80, 0, 0.4, 20, 0, 40);
    TH3D *h_lr[kNCls];
    for (int k = 0; k < kNCls; k++)
        h_lr[k] = new TH3D(Form("h3_liang_m_E_ring_%s", kClsName[k]), ";m_{#gamma#gamma} (GeV);cluster pair E (GeV);ring",
                           80, 0, 0.4, 20, 0, 40, kNRing, 0, kNRing);
    TH3D *h_lf[kNCls];
    for (int k = 0; k < kNCls; k++)
        h_lf[k] = new TH3D(Form("h3_liangf_m_E_ring_%s", kClsName[k]), ";m_{#gamma#gamma} (GeV);cluster pair E (GeV);ring",
                           80, 0, 0.4, 20, 0, 40, kNRing, 0, kNRing);
    TH3D *h_lsep[2], *h_lph[2];
    TH2D *h_leta[2];
    for (int k = 0; k < 2; k++)
    {
        const char *cn = k == 0 ? "all" : "em";
        h_lsep[k] = new TH3D(Form("h3_lsep_m_E_sep_%s", cn), ";m_{#gamma#gamma} (GeV);pair E (GeV);separation (towers)", 80, 0, 0.4, 20, 0, 40, 20, 0, 10);
        h_lph[k] = new TH3D(Form("h3_lph_m_Eg_Z_%s", cn), ";m_{#gamma#gamma} (GeV);photon E (GeV);Z", 80, 0, 0.4, 40, 0, 20, 7, 0, 0.7);
        h_leta[k] = new TH2D(Form("h2_leta_m_Eg_%s", cn), ";m_{#gamma#gamma} (GeV);lower photon E (GeV)", 160, 0, 0.8, 40, 0, 20);
    }
    TH2D *h_clw = new TH2D("h2_clwidth_E", ";cluster E (GeV);energy-weighted RMS radius (cm)", 40, 0, 20, 60, 0, 6);
    TH1D *h_iso = new TH1D("h_isohit_E", ";isolated ECal tower E (GeV);towers", 100, 0, 1);
    TProfile *p_sumf2 = new TProfile("p_liang_sumf2_E", ";cluster pair E (GeV);#LT#Sigma f_{i}^{2}#GT of the pair's clusters", 20, 0, 40);
    TH1D *h_turn[4];
    for (int k = 0; k < 4; k++)
        h_turn[k] = new TH1D(Form("h_maxclet_%s", k == 0 ? "all" : Form("em%d", k - 1)), ";highest ECal cluster E_{T} (GeV);events", 80, 0, 8);
    TH1D *h_ptE = new TH1D("h_ptE", ";point E (GeV);points", 100, 0, 100);
    TH1D *h_npts = new TH1D("h_npts", ";selected points;events", 30, -0.5, 29.5);
    TH3D *h3v = sim ? 0 : new TH3D("h3_m_E_d_vtx", ";m_{#gamma#gamma} (GeV);E bin;d bin", kNM, 0, kMMax, kNE, 0, kNE, kND, 0, kND);
    TH3D *h3t = 0, *h3mt = 0, *h3mp = 0, *h3sep = 0;
    TH2D *h_resp = 0, *h_dpos = 0;
    if (sim)
    {
        h3t = new TH3D("h3_m_E_d_tag", ";m_{#gamma#gamma} (GeV);E bin;d bin", kNM, 0, kMMax, kNE, 0, kNE, kND, 0, kND);
        h3mt = new TH3D("h3_mt_E_d_tag", ";m (true E, reco positions) (GeV);E bin;d bin", kNM, 0, kMMax, kNE, 0, kNE, kND, 0, kND);
        h3mp = new TH3D("h3_mp_E_d_tag", ";m (reco E, true positions) (GeV);E bin;d bin", kNM, 0, kMMax, kNE, 0, kNE, kND, 0, kND);
        h3sep = new TH3D("h3_dsep_E_d", ";reco d - true d (cm);E bin;true d bin", 80, -10, 10, kNE, 0, kNE, kND, 0, kND);
        h_resp = new TH2D("h2_resp_E", ";photon E (GeV);point E / photon E", 50, 0, 100, 150, 0, 1.5);
        h_dpos = new TH2D("h2_dpos_E", ";photon E (GeV);reco - true position (cm)", 50, 0, 100, 60, 0, 6);
    }
    TH1D *h_ev = new TH1D("h_events", ";;events (weighted)", 3, 0, 3);
    h_ev->GetXaxis()->SetBinLabel(1, "read"); h_ev->GetXaxis()->SetBinLabel(2, "no LED, triggered"); h_ev->GetXaxis()->SetBinLabel(3, "with vertex");

    const FiducialRect &rect = kFiducialRectEcal;
    long n_files = 0;
    for (size_t i = 0; i < files.size(); i++)
    {
        double w = 1.0;
        if (sim)
        {
            std::map<std::string, int>::iterator it = file_group.find(files[i]);
            if (it == file_group.end()) continue;
            w = unweighted ? 1.0 : group_sigmaN[it->second] / group_N[it->second] / group_N[it->second];
        }
        TFile *tf = TFile::Open(files[i].c_str());
        TTree *t = tf && !tf->IsZombie() ? (TTree *)tf->Get("data") : 0;
        if (!t) { delete tf; continue; }
        n_files++;
        static Int_t nhit, hdet[kMaxHits];
        static Float_t hE[kMaxHits], hx[kMaxHits], hy[kMaxHits];
        static Int_t nclu, cdet[kMaxPts];
        static Float_t cE[kMaxPts], cx[kMaxPts], cy[kMaxPts], cz[kMaxPts];
        static Int_t npts, det[kMaxPts], clu[kMaxPts], led, em[3], vsrc, nmc, geid[kMaxMc], vtx[kMaxMc], vtxend[kMaxMc];
        static Float_t pE[kMaxPts], px[kMaxPts], py[kMaxPts], pz[kMaxPts], vz;
        static Float_t mpx[kMaxMc], mpy[kMaxMc], mpz[kMaxMc], mE[kMaxMc], mvx[kMaxMc], mvy[kMaxMc], mvz[kMaxMc];
        t->SetBranchStatus("*", 0);
        const char *br[] = {"Cal_nclus", "Cal_clus_detid", "Cal_clus_energy", "Cal_clus_x", "Cal_clus_y", "Cal_clus_z", "Cal_nhits", "Cal_detid", "Cal_hit_energy", "Cal_hit_posx", "Cal_hit_posy", "Cal_npts", "Cal_pt_detid", "Cal_pt_parentClusterId", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "fcs_led", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (size_t b = 0; b < sizeof(br) / sizeof(br[0]); b++) t->SetBranchStatus(br[b], 1);
        t->SetBranchAddress("Cal_npts", &npts); t->SetBranchAddress("Cal_pt_detid", det); t->SetBranchAddress("Cal_pt_parentClusterId", clu); t->SetBranchAddress("Cal_pt_energy", pE);
        t->SetBranchAddress("Cal_pt_x", px); t->SetBranchAddress("Cal_pt_y", py); t->SetBranchAddress("Cal_pt_z", pz);
        t->SetBranchAddress("fcs_led", &led);
        t->SetBranchAddress("Cal_nclus", &nclu); t->SetBranchAddress("Cal_clus_detid", cdet); t->SetBranchAddress("Cal_clus_energy", cE);
        t->SetBranchAddress("Cal_clus_x", cx); t->SetBranchAddress("Cal_clus_y", cy); t->SetBranchAddress("Cal_clus_z", cz);
        t->SetBranchAddress("Cal_nhits", &nhit); t->SetBranchAddress("Cal_detid", hdet); t->SetBranchAddress("Cal_hit_energy", hE);
        t->SetBranchAddress("Cal_hit_posx", hx); t->SetBranchAddress("Cal_hit_posy", hy);
        if (emu_dir.empty())
        {
            t->SetBranchAddress("fcsEM0", &em[0]); t->SetBranchAddress("fcsEM1", &em[1]); t->SetBranchAddress("fcsEM2", &em[2]);
        }
        else
        {
            const std::string base = files[i].substr(files[i].rfind('/') + 1);
            t->AddFriend("trgemu", (emu_dir + "/" + base).c_str());
            t->SetBranchAddress("trgemu.fcsEM0", &em[0]); t->SetBranchAddress("trgemu.fcsEM1", &em[1]); t->SetBranchAddress("trgemu.fcsEM2", &em[2]);
        }
        if (sim)
        {
            const char *mb[] = {"mcpart_num", "mcpart_geid", "mcpart_idVtx", "mcpart_idVtxEnd", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
            for (size_t b = 0; b < sizeof(mb) / sizeof(mb[0]); b++) t->SetBranchStatus(mb[b], 1);
            t->SetBranchAddress("mcpart_num", &nmc); t->SetBranchAddress("mcpart_geid", geid);
            t->SetBranchAddress("mcpart_idVtx", vtx); t->SetBranchAddress("mcpart_idVtxEnd", vtxend);
            t->SetBranchAddress("mcpart_px", mpx); t->SetBranchAddress("mcpart_py", mpy); t->SetBranchAddress("mcpart_pz", mpz);
            t->SetBranchAddress("mcpart_E", mE);
            t->SetBranchAddress("mcpart_Vtx_x", mvx); t->SetBranchAddress("mcpart_Vtx_y", mvy); t->SetBranchAddress("mcpart_Vtx_z", mvz);
        }
        else
        {
            t->SetBranchStatus("event_vertex_z", 1); t->SetBranchStatus("event_vertex_source", 1);
            t->SetBranchAddress("event_vertex_z", &vz); t->SetBranchAddress("event_vertex_source", &vsrc);
        }
        for (Long64_t ev = 0; ev < t->GetEntries(); ev++)
        {
            t->GetEntry(ev);
            if (ecal_k != 1.0)
            {
                for (int k = 0; k < std::min((int)nhit, kMaxHits); k++) if (hdet[k] <= 1) hE[k] *= ecal_k;
                for (int k = 0; k < std::min((int)nclu, kMaxPts); k++) if (cdet[k] <= 1) cE[k] *= ecal_k;
                for (int k = 0; k < std::min((int)npts, kMaxPts); k++) if (det[k] <= 1) pE[k] *= ecal_k;
            }
            h_ev->Fill(0.0, w);
            if (led) continue;
            if (need_trig && !(em[0] || em[1] || em[2])) continue;
            h_ev->Fill(1, w);
            const bool has_vtx = !sim && vsrc > 0;
            if (has_vtx) h_ev->Fill(2, w);

            // MIP peak: isolated single ECal towers away from any cluster above 1 GeV
            for (int k = 0; k < nhit && k < kMaxHits; k++)
            {
                if (hdet[k] > 1 || hE[k] < 0.045) continue;
                bool iso = true;
                for (int l = 0; l < nhit && l < kMaxHits && iso; l++)
                    if (l != k && hdet[l] == hdet[k] && hE[l] > 0.010 && std::hypot(hx[l] - hx[k], hy[l] - hy[k]) < 1.5 * 5.542) iso = false;
                for (int a = 0; a < nclu && a < kMaxPts && iso; a++)
                    if (cdet[a] == hdet[k] && cE[a] > 1.0 && std::hypot(cx[a] - hx[k], cy[a] - hy[k]) < 3 * 5.542) iso = false;
                if (iso) h_iso->Fill(std::min((double)hE[k], 0.999), w);
            }

            // Trigger turn-on: highest ECal cluster E_T, all events and per EM bit
            double maxet = 0;
            {
                for (int a = 0; a < nclu && a < kMaxPts; a++)
                    if (cdet[a] <= 1) maxet = std::max(maxet, (double)cE[a] * std::hypot(cx[a], cy[a]) / std::sqrt(cx[a] * cx[a] + cy[a] * cy[a] + cz[a] * cz[a]));
                h_turn[0]->Fill(std::min(maxet, 7.99), w);
                for (int k = 0; k < 3; k++) if (em[k]) h_turn[k + 1]->Fill(std::min(maxet, 7.99), w);
            }

            // Gain-calibration selection (StFcsPi0FinderForEcal)
            {
                int nhot = 0;
                for (int k = 0; k < nhit && k < kMaxHits; k++) if (hdet[k] <= 1 && hE[k] > 1.0) nhot++;
                if (nhot <= 40)
                {
                    double bestE = -1, bestM = -1, bestR = -1, bestX = 0, bestY = 0;
                    int bestA = -1, bestB = -1;
                    for (int a = 0; a < nclu && a < kMaxPts; a++)
                    {
                        if (cdet[a] > 1 || cE[a] < 1.0) continue;
                        for (int b = a + 1; b < nclu && b < kMaxPts; b++)
                        {
                            if (cdet[b] != cdet[a] || cE[b] < 1.0) continue;
                            if (std::fabs(cE[a] - cE[b]) / (cE[a] + cE[b]) > 0.7) continue;
                            if (cE[a] + cE[b] > bestE)
                            {
                                bestE = cE[a] + cE[b];
                                bestM = Mass(cE[a], TVector3(cx[a], cy[a], cz[a]), cE[b], TVector3(cx[b], cy[b], cz[b]));
                                bestX = (cE[a] * cx[a] + cE[b] * cx[b]) / bestE; bestY = (cE[a] * cy[a] + cE[b] * cy[b]) / bestE;
                                bestR = std::hypot(bestX, bestY); bestA = a; bestB = b;
                            }
                        }
                    }
                    if (bestE > 0)
                    {
                        h_liang->Fill(std::min(bestM, 0.3999), std::min(bestE, 39.99), w);
                        const int ir = Bin(bestR, kRingEdges, kNRing);
                        const bool anyem = em[0] || em[1] || em[2];
                        // offline firing-patch stand-in: highest-E_T window of hits
                        double bestEt = -1, px0 = 0, py0 = 0;
                        for (int k = 0; k < nhit && k < kMaxHits; k++)
                        {
                            if (hdet[k] > 1 || hE[k] < 0.1) continue;
                            double et = 0;
                            for (int l = 0; l < nhit && l < kMaxHits; l++)
                                if (hdet[l] <= 1 && std::fabs(hx[l] - hx[k]) < kPatchHalf && std::fabs(hy[l] - hy[k]) < kPatchHalf)
                                    et += hE[l] * std::hypot(hx[l], hy[l]) / 724.0;
                            if (et > bestEt) { bestEt = et; px0 = hx[k]; py0 = hy[k]; }
                        }
                        const bool inpatch = std::fabs(bestX - px0) < kPatchHalf && std::fabs(bestY - py0) < kPatchHalf;
                        const bool cls[kNCls] = {true, anyem, !anyem, em[2] != 0, anyem && inpatch, anyem && !inpatch,
                                                 anyem && maxet > 2, anyem && maxet > 3, anyem && maxet > 4};
                        for (int k = 0; k < kNCls; k++)
                            if (cls[k]) h_lr[k]->Fill(std::min(bestM, 0.3999), std::min(bestE, 39.99), ir + 0.5, w);
                        // common 45 MeV tower floor
                        double ef[2];
                        const int cc[2] = {bestA, bestB};
                        for (int j = 0; j < 2; j++)
                        {
                            const int c = cc[j];
                            double low = 0;
                            for (int k = 0; k < nhit && k < kMaxHits; k++)
                                if (hdet[k] == cdet[c] && hE[k] > 0.010 && hE[k] < 0.045 && std::hypot(hx[k] - cx[c], hy[k] - cy[c]) < kTailCm) low += hE[k];
                            ef[j] = cE[c] - low;
                        }
                        const double mf = bestM * std::sqrt(ef[0] * ef[1] / ((double)cE[bestA] * cE[bestB])), Ef = ef[0] + ef[1];
                        for (int k = 0; k < kNCls; k++)
                            if (cls[k]) h_lf[k]->Fill(std::min(mf, 0.3999), std::min(Ef, 39.99), ir + 0.5, w);
                        // overlap study, on the floored energies
                        const double sep = std::hypot(cx[bestA] - cx[bestB], cy[bestA] - cy[bestB]) / 5.542;
                        const double zf = std::fabs(ef[0] - ef[1]) / Ef;
                        for (int k = 0; k < 2; k++)
                        {
                            if (k == 1 && !anyem) continue;
                            h_lsep[k]->Fill(std::min(mf, 0.3999), std::min(Ef, 39.99), std::min(sep, 9.99), w);
                            for (int j = 0; j < 2; j++) h_lph[k]->Fill(std::min(mf, 0.3999), std::min(ef[j], 19.99), std::min(zf, 0.699), w);
                            h_leta[k]->Fill(std::min(mf, 0.7999), std::min(std::min(ef[0], ef[1]), 19.99), w);
                        }
                        for (int j = 0; j < 2; j++)
                        {
                            const int c = cc[j];
                            double s0 = 0, s2 = 0;
                            for (int k = 0; k < nhit && k < kMaxHits; k++)
                            {
                                if (hdet[k] != cdet[c] || hE[k] < 0.045) continue;
                                const double d = std::hypot(hx[k] - cx[c], hy[k] - cy[c]);
                                if (d < kTailCm) { s0 += hE[k]; s2 += hE[k] * d * d; }
                            }
                            if (s0 > 0) h_clw->Fill(std::min(ef[j], 19.99), std::min(std::sqrt(s2 / s0), 5.99), w);
                        }
                        // tower energy sharing of the two clusters: sum of (hit E / cluster E)^2,
                        // hits within kTailCm of the cluster position
                        for (int c : {bestA, bestB})
                        {
                            double s1 = 0, s2 = 0;
                            for (int k = 0; k < nhit && k < kMaxHits; k++)
                                if (hdet[k] == cdet[c] && std::hypot(hx[k] - cx[c], hy[k] - cy[c]) < kTailCm) { s1 += hE[k]; s2 += hE[k] * hE[k]; }
                            if (s1 > 0) p_sumf2->Fill(std::min(bestE, 39.99), s2 / (s1 * s1), w);
                        }
                    }
                }
            }

            // Every cluster pair, with the point pairs' cuts
            {
                std::vector<int> cs;
                for (int a = 0; a < nclu && a < kMaxPts; a++)
                {
                    if (cdet[a] > 1 || !(cE[a] > point_energy_min)) continue;
                    const double ax = std::fabs(cx[a]);
                    if (ax < rect.x_inner + kPointEdgeBuffer || ax > rect.x_outer - kPointEdgeBuffer) continue;
                    if (cy[a] < rect.y_min + kPointEdgeBuffer || cy[a] > rect.y_max - kPointEdgeBuffer) continue;
                    cs.push_back(a);
                }
                for (size_t i = 0; i < cs.size(); i++)
                    for (size_t j = i + 1; j < cs.size(); j++)
                    {
                        const int a = cs[i], b = cs[j];
                        if (cdet[a] != cdet[b]) continue;
                        const double E = cE[a] + cE[b];
                        if (std::fabs(cE[a] - cE[b]) / E >= kZMax) continue;
                        const int ie = Bin(E, kEEdges, kNE);
                        const TVector3 u1(cx[a], cy[a], cz[a]), u2(cx[b], cy[b], cz[b]);
                        const int id = Bin((u1 - u2).Mag(), kDEdges, kND);
                        if (ie < 0 || id < 0) continue;
                        h3clu->Fill(std::min(Mass(cE[a], u1, cE[b], u2), kMMax - 1e-6), ie + 0.5, id + 0.5, w);
                    }
            }

            std::vector<int> sel;
            for (int p = 0; p < npts && p < kMaxPts; p++)
            {
                if (det[p] != 0 && det[p] != 1) continue;
                if (!(pE[p] > point_energy_min)) continue;
                const double ax = std::fabs(px[p]);
                if (ax < rect.x_inner + kPointEdgeBuffer || ax > rect.x_outer - kPointEdgeBuffer) continue;
                if (py[p] < rect.y_min + kPointEdgeBuffer || py[p] > rect.y_max - kPointEdgeBuffer) continue;
                sel.push_back(p);
                h_ptE->Fill(pE[p], w);
            }
            h_npts->Fill(sel.size(), w);
            if (sel.size() < 2) continue;

            // Dropped shower tails: sub-threshold hits to the nearest selected point
            std::vector<double> tail(sel.size(), 0.0);
            for (int k = 0; k < nhit && k < kMaxHits; k++)
            {
                if (hdet[k] > 1 || hE[k] <= kTailMin || hE[k] >= kTowerThr) continue;
                int best = -1; double bd = kTailCm;
                for (size_t a = 0; a < sel.size(); a++)
                {
                    if (det[sel[a]] != hdet[k]) continue;
                    const double dd = std::hypot(hx[k] - px[sel[a]], hy[k] - py[sel[a]]);
                    if (dd < bd) { bd = dd; best = (int)a; }
                }
                if (best >= 0) tail[best] += hE[k];
            }
            for (size_t a = 0; a < sel.size(); a++) h_tail->Fill(pE[sel[a]], tail[a] / pE[sel[a]], w);

            // Simulation: nearest GEANT photon at each point's z
            std::vector<int> match(sel.size(), -1);
            std::vector<TVector3> tpos(sel.size());
            if (sim)
            {
                for (size_t a = 0; a < sel.size(); a++)
                {
                    const int p = sel[a];
                    double best = kMatchCm;
                    for (int k = 0; k < nmc && k < kMaxMc; k++)
                    {
                        if (geid[k] != 1 || mE[k] < 0.2 || mpz[k] <= 0) continue;
                        const double s = (pz[p] - mvz[k]) / mpz[k];
                        if (s <= 0) continue;
                        const double x = mvx[k] + s * mpx[k], y = mvy[k] + s * mpy[k];
                        const double dd = std::hypot(x - px[p], y - py[p]);
                        if (dd < best) { best = dd; match[a] = k; tpos[a].SetXYZ(x, y, pz[p]); }
                    }
                    if (match[a] >= 0)
                    {
                        h_resp->Fill(mE[match[a]], pE[p] / mE[match[a]], w);
                        h_dpos->Fill(mE[match[a]], best, w);
                    }
                }
            }

            for (size_t a = 0; a < sel.size(); a++)
                for (size_t b = a + 1; b < sel.size(); b++)
                {
                    const int p = sel[a], q = sel[b];
                    if (det[p] != det[q]) continue;
                    const double E = pE[p] + pE[q], Z = std::fabs(pE[p] - pE[q]) / E;
                    const int ie = Bin(E, kEEdges, kNE);
                    if (ie < 0) continue;
                    const TVector3 v1(px[p], py[p], pz[p]), v2(px[q], py[q], pz[q]);
                    const double d = (v1 - v2).Mag(), m = Mass(pE[p], v1, pE[q], v2);
                    const int iz = Bin(Z, kZEdges, kNZ);
                    h3z->Fill(std::min(m, kMMax - 1e-6), ie + 0.5, iz + 0.5, w);
                    if (Z >= kZMax) continue;
                    h3c->Fill(std::min(m, kMMax - 1e-6), ie + 0.5, (clu[p] >= 0 && clu[p] == clu[q]) ? 0.5 : 1.5, w);
                    const int id = Bin(d, kDEdges, kND);
                    if (id < 0) continue;
                    h3->Fill(std::min(m, kMMax - 1e-6), ie + 0.5, id + 0.5, w);
                    h3tail->Fill(std::min(Mass(pE[p] + tail[a], v1, pE[q] + tail[b], v2), kMMax - 1e-6), ie + 0.5, id + 0.5, w);
                    if (has_vtx)
                    {
                        const TVector3 o(0, 0, vz);
                        h3v->Fill(std::min(Mass(pE[p], v1 - o, pE[q], v2 - o), kMMax - 1e-6), ie + 0.5, id + 0.5, w);
                    }
                    if (sim && match[a] >= 0 && match[b] >= 0 && match[a] != match[b] && vtx[match[a]] == vtx[match[b]])
                    {
                        // both photons from one pi0: a pi0 ending where they start
                        bool pi0 = false;
                        for (int k = 0; k < nmc && k < kMaxMc && !pi0; k++) pi0 = geid[k] == 7 && vtxend[k] == vtx[match[a]];
                        if (!pi0) continue;
                        const int ka = match[a], kb = match[b];
                        h3t->Fill(std::min(m, kMMax - 1e-6), ie + 0.5, id + 0.5, w);
                        h3mt->Fill(std::min(Mass(mE[ka], v1, mE[kb], v2), kMMax - 1e-6), ie + 0.5, id + 0.5, w);
                        h3mp->Fill(std::min(Mass(pE[p], tpos[a], pE[q], tpos[b]), kMMax - 1e-6), ie + 0.5, id + 0.5, w);
                        const double dt = (tpos[a] - tpos[b]).Mag();
                        const int idt = Bin(dt, kDEdges, kND);
                        if (idt >= 0) h3sep->Fill(std::max(-9.999, std::min(9.999, d - dt)), ie + 0.5, idt + 0.5, w);
                    }
                }
        }
        delete tf;
    }
    std::printf("%ld files read; events read %.4g, triggered %.4g, with vertex %.4g\n", n_files,
                h_ev->GetBinContent(1), h_ev->GetBinContent(2), h_ev->GetBinContent(3));
    out.cd();
    out.Write();
    out.Close();
    return 0;
}
