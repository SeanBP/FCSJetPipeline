// In-event UE estimate for FCS EM jets: cones at the jet's eta, rotated in
// phi by 180, +-135 and +-90 deg, in the jet's own event. Partly dead cones
// (beam hole, gap between the halves, outer edge) are scaled to the jet cone
// two ways:
//   area      E_live * A(jet cone, live) / A(rotated cone, live)
//   template  E_live * T(jet cone) / T(rotated cone, live), with T the sum
//             over live cells of the mean UE density map rho(cell) * area:
//             the map, from this job's events with cells >= 2R from every
//             reco jet, carries the eta dependence of the UE and the edge
//             response, so a cone missing its high-eta side is not
//             under-corrected
// Recoil test: the same cones from truth photons (simulation; no acceptance,
// full cones) and their dependence on jet energy, 180 against 90 deg.
//
// Input: pairs "<jet_output> <SimpleTree>" (one data run or one ptHat bin per
// job, as proxy_dr_profile). Jets: EM-triggered events, reco pT >= 2 GeV,
// inside the reco fiducial region, reco E 10-120 GeV.
// Output tree "jets", one row per jet:
//   jet_pt, jet_eta, jet_phi, jet_E, jet_x, jet_y, em (highest EM level, stored
//   bits: data hardware, simulation its own emulation; -1 none), emu (highest
//   re-emulated level, -1 none or no emulation dir), hw0 (stored fcsEM0 bit),
//   group (run or -1-ptHatMin), weight (simulation sigma / N per ptHat bin, data 1)
//   a_jet, t_jet        live area (fraction of pi R^2) and template sum of the jet cone
//   e[k], a[k], t[k]    rotated cone k: live point energy, live area fraction,
//                       template sum over its live cells
//   te[k], te_jet       simulation: truth-photon energy (> 0.2 GeV, full cone)
//   s[k]                eta-strip estimate of the jet cone's UE from cone k: each
//                       point weighted by A(jet cone, live; its eta strip) /
//                       A(cone k, live; same strip), strips 0.1 wide; strips
//                       of cone k under 20% live are uncovered, and the
//                       covered strips' template share of the jet cone scales
//                       the sum up to the whole live jet cone
//   closure (simulation, truth photons projected straight onto the ECal plane):
//   tf[k]               truth energy in the jet cone's live footprint moved to
//                       cone k (the target: same eta distribution, no jet)
//   tla[k], tlt[k], tls[k]  its estimates from the truth photons on live cells
//                       of cone k, by area, template and eta-strip scaling
//   k = 0: 180 deg, 1: +135, 2: -135, 3: +90, 4: -90
//   annulus (in_event_geo.h), v = 0: whole ring, 1: transverse 60-120 deg on
//   both sides, 2: away side > 120 deg:
//   an_e[v]             eta-weighted annulus estimate of the jet cone's UE
//   an_area[v], an_fill[v]  live annulus area (eta-phi) and uncovered-strip scale-up
//   an_tl[v]            the same from the truth photons on live cells
//   an_tf[v]            its truth target: the jet cone's live footprint moved to
//                       the annulus's rotations (10 deg steps, >= 2R), averaged
//   v + 3               the same three without the other-jet exclusion
//   separation scan, whole ring with other jets excluded (from truth photons on
//   live cells: *_tl):
//   sh_e[q], sh_fill[q], sh_tl[q]  shells 0.7-0.8-...-2.0-2.5-3.3 in dR from the jet
//                       (strips covered at >= 3% of a cone slice)
//   dm_e[q], dm_fill[q], dm_tl[q]  dR >= 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.75
//   pixel ring (in_event_geo.h: ring cells, each point weighted by the ratio
//   of the jet's to the ring's live area in its own eta bin, matching the
//   ring's area-per-eta to the jet's, weighted-average density times the jet's live
//   area; no jet dropped), p = 0 other jets excluded, 1 not:
//   px_e[p], px_tl[p]   UE estimate (reco points; truth photons on live cells)
//   px_rho[p]           the weighted-average UE density (reco)
//   px_area[p], px_cov[p], px_wmax[p]  live ring area, share of the jet area in
//                       eta bins the ring reaches, largest cell weight
//   px_fb[p]            1: the ring reaches none of the jet's eta bins; the
//                       template map's expectation t_jet is used
//   leakage rings, dR = 1.0-1.1-1.2-1.4-1.6-2.0-2.5 R from the jet axis, all azimuth:
//   lr_E[q], lr_A[q]    point energy and live eta-phi area (sub-pixels)
//   lrx_E[q], lrx_A[q]  the same without what lies within R of other reco jets
// and h_map (mean point energy per unit eta-phi area per event, per cell)
// and h_rho_eta (the same against eta).
//
// Build and run from this folder:
//   g++ -m64 -O2 in_event_ue.cpp ../../shared/VectorDict.cxx -o in_event_ue -I../../shared -I../../underlying_event -I../../jet_scale `root-config --cflags --libs` -std=c++11
//   ./in_event_ue <pairs.list> <out.root> [max jets = 2000000] [emulation dir]
// The emulation dir holds Pi0Playground/trigemu/reemulate.C friend files
// (tree "trgemu", named as the SimpleTrees): each event's re-emulated EM
// level is then kept (emu), and events are kept if either their stored or
// their re-emulated EM level fired.

#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "JetXSecWeights.h"
#include "in_event_geo.h"

static const int kNC = 5;
static const double kRot[kNC] = {180, 135, -135, 90, -90};
static const double kMapExcl = 2 * R; // map cells this far from every reco jet

struct Pt { float eta, phi, E; int cell; };
struct Jet { float pt, eta, phi, x, y; };
struct Event { int group, em, emu, hw0; int pt_first, pt_n, jet_first, jet_n, t_first, t_n; };

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 3) { std::cerr << "usage: ./in_event_ue <pairs.list> <out.root> [max jets] [emulation dir]" << std::endl; return 1; }
    const long max_jets = argc > 3 ? std::atol(argv[3]) : 2000000;
    const std::string emu_dir = argc > 4 && std::string(argv[4]) != "-" ? argv[4] : "";

    std::vector<std::pair<std::string, std::string> > pairs;
    {
        std::ifstream in(argv[1]);
        std::string a, b;
        while (in >> a >> b) pairs.push_back(std::make_pair(a, b));
    }
    if (pairs.empty()) { std::cerr << "ERROR: no pairs" << std::endl; return 1; }

    Geometry geo; geo.Build();
    std::printf("in_event_ue: %zu ECal cells of %.0f cm\n", geo.cells.size(), kCell);

    std::vector<Event> events;
    std::vector<Pt> pts;
    std::vector<Jet> jets;
    std::vector<float> T_eta, T_phi, T_E;
    std::map<int, double> group_N, group_sigmaN;
    const int MAX_PTS = 200, MAX_MC = 10000;
    for (size_t f = 0; f < pairs.size(); f++)
    {
        const std::string &jf = pairs[f].first;
        const std::string base = jf.substr(jf.find_last_of('/') + 1);
        int group = 0, a = 0, b = 0;
        if (std::sscanf(base.c_str(), "jet_output_%d_%d.root", &a, &b) == 2) group = a;
        else { size_t p = base.rfind("_pt"); if (p != std::string::npos && std::sscanf(base.c_str() + p, "_pt%d_%d.root", &a, &b) == 2) group = -1 - a; }

        TFile *fj = TFile::Open(jf.c_str());
        TFile *fs = TFile::Open(pairs[f].second.c_str());
        TFile *fe = 0;
        TTree *te = 0;
        Int_t m0 = 0, m1 = 0, m2 = 0;
        if (!emu_dir.empty())
        {
            const std::string sb = pairs[f].second.substr(pairs[f].second.find_last_of('/') + 1);
            fe = TFile::Open((emu_dir + "/" + sb).c_str());
            te = (fe && !fe->IsZombie()) ? (TTree *)fe->Get("trgemu") : 0;
            if (!te) std::cerr << "no emulation file for " << sb << std::endl;
            else { te->SetBranchAddress("fcsEM0", &m0); te->SetBranchAddress("fcsEM1", &m1); te->SetBranchAddress("fcsEM2", &m2); }
        }
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        const bool has_entry = jt && jt->GetBranch("simpletree_entry");
        if (!jt || !st || (!has_entry && jt->GetEntries() != st->GetEntries())) { std::cerr << "skip " << jf << std::endl; delete fj; delete fs; continue; }

        std::vector<float> *reco_pt = 0, *reco_eta = 0, *reco_phi = 0, *reco_x = 0, *reco_y = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        Long64_t simpletree_entry = -1;
        Double_t sqrt_s = kSqrtS_MainRun22_GeV, sigma_job = -1;
        Int_t N_job = -1;
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int i = 0; i < 8; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &reco_pt); jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi); jt->SetBranchAddress("reco_x", &reco_x); jt->SetBranchAddress("reco_y", &reco_y);
        jt->SetBranchAddress("fcsEM0", &em0); jt->SetBranchAddress("fcsEM1", &em1); jt->SetBranchAddress("fcsEM2", &em2);
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &simpletree_entry); }
        if (jt->GetBranch("event_sqrt_s_gev")) { jt->SetBranchStatus("event_sqrt_s_gev", 1); jt->SetBranchAddress("event_sqrt_s_gev", &sqrt_s); }
        if (group < 0 && jt->GetBranch("sigma_job"))
        {
            jt->SetBranchStatus("sigma_job", 1); jt->SetBranchStatus("N_job", 1);
            jt->SetBranchAddress("sigma_job", &sigma_job); jt->SetBranchAddress("N_job", &N_job);
        }
        Int_t Cal_npts = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        Int_t Cal_pt_detid[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pb[i], 1);
        st->SetBranchAddress("Cal_npts", &Cal_npts); st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x); st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z); st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);
        static Int_t mc_n = 0, mc_geid[MAX_MC], mc_idVtx[MAX_MC];
        static Float_t mc_px[MAX_MC], mc_py[MAX_MC], mc_pz[MAX_MC], mc_E[MAX_MC], mc_vx[MAX_MC], mc_vy[MAX_MC], mc_vz[MAX_MC];
        const bool file_truth = group < 0 && st->GetBranch("mcpart_num");
        if (file_truth)
        {
            const char *mb[] = {"mcpart_num", "mcpart_geid", "mcpart_idVtx", "mcpart_px", "mcpart_py", "mcpart_pz", "mcpart_E", "mcpart_Vtx_x", "mcpart_Vtx_y", "mcpart_Vtx_z"};
            for (int i = 0; i < 10; i++) st->SetBranchStatus(mb[i], 1);
            st->SetBranchAddress("mcpart_num", &mc_n); st->SetBranchAddress("mcpart_geid", mc_geid); st->SetBranchAddress("mcpart_idVtx", mc_idVtx);
            st->SetBranchAddress("mcpart_px", mc_px); st->SetBranchAddress("mcpart_py", mc_py); st->SetBranchAddress("mcpart_pz", mc_pz);
            st->SetBranchAddress("mcpart_E", mc_E); st->SetBranchAddress("mcpart_Vtx_x", mc_vx); st->SetBranchAddress("mcpart_Vtx_y", mc_vy);
            st->SetBranchAddress("mcpart_Vtx_z", mc_vz);
        }

        bool counted = false;
        for (Long64_t ev = 0; ev < jt->GetEntries(); ev++)
        {
            if (jt->GetEntry(ev) <= 0 || reco_pt->empty()) continue;
            if (!counted && group < 0)
            {
                const double sg = DecodePythia6Sigma(sigma_job);
                if (sg > 0 && N_job > 0) { group_N[group] += N_job; group_sigmaN[group] += N_job * sg; }
                counted = true;
            }
            const int em_level = ProxyEmLevel(em0, em1, em2);
            int emu_level = -1;
            if (te && te->GetEntry(has_entry ? simpletree_entry : ev) > 0) emu_level = ProxyEmLevel(m0, m1, m2);
            if (em_level < 0 && emu_level < 0) continue;
            if (st->GetEntry(has_entry ? simpletree_entry : ev) <= 0) continue;
            Event e;
            e.group = group; e.em = em_level; e.emu = emu_level; e.hw0 = em0 ? 1 : 0;
            e.jet_first = (int)jets.size(); e.pt_first = (int)pts.size(); e.t_first = (int)T_E.size();
            for (size_t j = 0; j < reco_pt->size(); j++)
            {
                Jet jet = {reco_pt->at(j), reco_eta->at(j), reco_phi->at(j), reco_x->at(j), reco_y->at(j)};
                jets.push_back(jet);
            }
            for (int i = 0; i < std::min((int)Cal_npts, MAX_PTS); i++)
            {
                if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                if (!pass_point_energy_cut(Cal_pt_energy[i], sqrt_s)) continue;
                const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                Pt p = {(float)std::asinh(Cal_pt_z[i] / rho), (float)std::atan2(Cal_pt_y[i], Cal_pt_x[i]), Cal_pt_energy[i], geo.CellOf(Cal_pt_x[i], Cal_pt_y[i])};
                pts.push_back(p);
            }
            // truth photons for the full cones: no acceptance, no eta range beyond 1.5-5.5
            float pvx, pvy, pvz;
            if (file_truth && FindEmTruthPrimaryVertex(std::min((int)mc_n, MAX_MC), mc_idVtx, mc_vx, mc_vy, mc_vz, pvx, pvy, pvz))
                for (int i = 0; i < std::min((int)mc_n, MAX_MC); i++)
                {
                    if (mc_E[i] < 0.2 || !IsEmTruthParticle(mc_geid[i], mc_vx[i], mc_vy[i], mc_vz[i], pvx, pvy, pvz)) continue;
                    const double ptp = std::sqrt(mc_px[i] * mc_px[i] + mc_py[i] * mc_py[i]);
                    if (!(ptp > 0)) continue;
                    const double teta = std::asinh(mc_pz[i] / ptp);
                    if (teta < 1.5 || teta > 5.5) continue;
                    T_eta.push_back(teta); T_phi.push_back(std::atan2(mc_py[i], mc_px[i])); T_E.push_back(mc_E[i]);
                }
            e.jet_n = (int)jets.size() - e.jet_first; e.pt_n = (int)pts.size() - e.pt_first; e.t_n = (int)T_E.size() - e.t_first;
            events.push_back(e);
        }
        fj->Close(); delete fj; fs->Close(); delete fs;
        if (fe) { fe->Close(); delete fe; }
    }
    std::printf("in_event_ue: %zu pairs, %zu EM-triggered events, %zu jets, %zu points\n", pairs.size(), events.size(), jets.size(), pts.size());

    // The UE density map: point energy and exposure per cell over cells >= 2R from every reco jet
    const size_t ncell = geo.cells.size();
    std::vector<double> m_E(ncell, 0.0), m_n(ncell, 0.0);
    std::vector<char> near(ncell, 0);
    std::vector<int> cl;
    for (size_t ie = 0; ie < events.size(); ie++)
    {
        const Event &e = events[ie];
        std::vector<int> marked;
        for (int j = e.jet_first; j < e.jet_first + e.jet_n; j++)
        {
            geo.Near(jets[j].eta, jets[j].phi, kMapExcl, cl);
            for (size_t q = 0; q < cl.size(); q++) if (!near[cl[q]]) { near[cl[q]] = 1; marked.push_back(cl[q]); }
        }
        for (size_t c = 0; c < ncell; c++) if (!near[c]) m_n[c] += 1.0;
        for (int i = e.pt_first; i < e.pt_first + e.pt_n; i++)
            if (pts[i].cell >= 0 && !near[pts[i].cell]) m_E[pts[i].cell] += pts[i].E;
        for (size_t q = 0; q < marked.size(); q++) near[marked[q]] = 0;
    }
    std::vector<double> rho(ncell, 0.0);
    // cells never exposed (always near a jet) take the mean of exposed cells at similar eta
    TH1D h_rho_eta("h_rho_eta", "mean UE density away from jets;#eta;GeV per unit #eta-#phi area per event", 40, 2.0, 4.6);
    TH1D h_rho_eta_n("h_rho_eta_n", "", 40, 2.0, 4.6);
    for (size_t c = 0; c < ncell; c++)
        if (m_n[c] > 0) { rho[c] = m_E[c] / m_n[c] / geo.cells[c].area; h_rho_eta.Fill(geo.cells[c].eta, rho[c] * m_n[c]); h_rho_eta_n.Fill(geo.cells[c].eta, m_n[c]); }
    h_rho_eta.Divide(&h_rho_eta_n);
    int n_fill = 0;
    for (size_t c = 0; c < ncell; c++)
        if (!(m_n[c] > 0)) { rho[c] = h_rho_eta.GetBinContent(h_rho_eta.FindBin(geo.cells[c].eta)); n_fill++; }
    std::printf("in_event_ue: map filled from %zu events; %d cells without exposure take the eta mean\n", events.size(), n_fill);

    // Per-group weights (simulation sigma / N, absolute; data 1)
    std::map<int, double> w_group;
    for (std::map<int, double>::iterator it = group_N.begin(); it != group_N.end(); ++it)
        w_group[it->first] = group_sigmaN[it->first] / it->second / it->second;

    TFile out(argv[2], "RECREATE");
    TTree tj("jets", "in-event UE cones per jet");
    Float_t t_pt, t_eta, t_phi, t_E, t_x, t_y, t_w, a_jet, t_jet, te_jet, e_[kNC], a_[kNC], tt_[kNC], te_[kNC];
    Float_t s_[kNC], tf_[kNC], tla_[kNC], tlt_[kNC], tls_[kNC];
    Int_t t_em, t_group, t_emu, t_hw0;
    tj.Branch("jet_pt", &t_pt, "jet_pt/F"); tj.Branch("jet_eta", &t_eta, "jet_eta/F"); tj.Branch("jet_phi", &t_phi, "jet_phi/F");
    tj.Branch("jet_E", &t_E, "jet_E/F"); tj.Branch("jet_x", &t_x, "jet_x/F"); tj.Branch("jet_y", &t_y, "jet_y/F");
    tj.Branch("em", &t_em, "em/I"); tj.Branch("emu", &t_emu, "emu/I"); tj.Branch("hw0", &t_hw0, "hw0/I"); tj.Branch("group", &t_group, "group/I"); tj.Branch("weight", &t_w, "weight/F");
    tj.Branch("a_jet", &a_jet, "a_jet/F"); tj.Branch("t_jet", &t_jet, "t_jet/F"); tj.Branch("te_jet", &te_jet, "te_jet/F");
    tj.Branch("e", e_, Form("e[%d]/F", kNC)); tj.Branch("a", a_, Form("a[%d]/F", kNC));
    tj.Branch("t", tt_, Form("t[%d]/F", kNC)); tj.Branch("te", te_, Form("te[%d]/F", kNC));
    tj.Branch("s", s_, Form("s[%d]/F", kNC)); tj.Branch("tf", tf_, Form("tf[%d]/F", kNC));
    tj.Branch("tla", tla_, Form("tla[%d]/F", kNC)); tj.Branch("tlt", tlt_, Form("tlt[%d]/F", kNC)); tj.Branch("tls", tls_, Form("tls[%d]/F", kNC));
    // v + kNA: the same ranges without the other-jet exclusion
    const int kNAV = 2 * kNA;
    Float_t an_e[kNAV], an_area[kNAV], an_fill[kNAV], an_tl[kNAV], an_tf[kNAV];
    tj.Branch("an_e", an_e, Form("an_e[%d]/F", kNAV)); tj.Branch("an_area", an_area, Form("an_area[%d]/F", kNAV));
    tj.Branch("an_fill", an_fill, Form("an_fill[%d]/F", kNAV)); tj.Branch("an_tl", an_tl, Form("an_tl[%d]/F", kNAV));
    tj.Branch("an_tf", an_tf, Form("an_tf[%d]/F", kNAV));
    // truth targets of the annulus: the jet cone's live footprint moved to these
    // rotations (10 deg steps, at least 2R from the jet), averaged
    std::vector<double> an_rot[kNA];
    for (int v = 0; v < kNA; v++)
        for (int d = -180; d < 180; d += 10)
        {
            const double r = std::fabs(d) * M_PI / 180.0;
            if (r >= kAnnJetExcl && r >= kAnnPhiLo[v] - 1e-9 && r <= kAnnPhiHi[v] + 1e-9) an_rot[v].push_back(d * M_PI / 180.0);
        }
    // separation scan (whole ring, other jets excluded): eta-weighted estimates
    // from shells d_lo <= dR < d_hi, and from the whole ring with dR >= d_min
    const int kNSh = 15, kNDm = 8;
    const double sh_edge[kNSh + 1] = {0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.3};
    const double dmin[kNDm] = {0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.75};
    const double kShellMinLive = 0.03;
    Float_t sh_e[kNSh], sh_fill[kNSh], sh_tl[kNSh], dm_e[kNDm], dm_fill[kNDm], dm_tl[kNDm];
    tj.Branch("sh_e", sh_e, Form("sh_e[%d]/F", kNSh)); tj.Branch("sh_fill", sh_fill, Form("sh_fill[%d]/F", kNSh)); tj.Branch("sh_tl", sh_tl, Form("sh_tl[%d]/F", kNSh));
    tj.Branch("dm_e", dm_e, Form("dm_e[%d]/F", kNDm)); tj.Branch("dm_fill", dm_fill, Form("dm_fill[%d]/F", kNDm)); tj.Branch("dm_tl", dm_tl, Form("dm_tl[%d]/F", kNDm));
    // pixel ring (in_event_geo.h), >= 2R from the jet; p = 0 other jets excluded, 1 not
    Float_t px_e[2], px_tl[2], px_rho[2], px_area[2], px_cov[2], px_wmax[2];
    Int_t px_fb[2];
    tj.Branch("px_e", px_e, "px_e[2]/F"); tj.Branch("px_tl", px_tl, "px_tl[2]/F"); tj.Branch("px_rho", px_rho, "px_rho[2]/F");
    tj.Branch("px_area", px_area, "px_area[2]/F"); tj.Branch("px_cov", px_cov, "px_cov[2]/F"); tj.Branch("px_wmax", px_wmax, "px_wmax[2]/F");
    tj.Branch("px_fb", px_fb, "px_fb[2]/I");
    // leakage rings: point energy and live area (sub-pixels) at dR = [a, b) R from the jet axis,
    // all azimuth; lr_* all points, lrx_* without the cells within R of other jets
    const int kNL = 6;
    const double lr_edge[kNL + 1] = {1.0, 1.1, 1.2, 1.4, 1.6, 2.0, 2.5};
    Float_t lr_E[kNL], lr_A[kNL], lrx_E[kNL], lrx_A[kNL];
    tj.Branch("lr_E", lr_E, Form("lr_E[%d]/F", kNL)); tj.Branch("lr_A", lr_A, Form("lr_A[%d]/F", kNL));
    tj.Branch("lrx_E", lrx_E, Form("lrx_E[%d]/F", kNL)); tj.Branch("lrx_A", lrx_A, Form("lrx_A[%d]/F", kNL));
    PixelRing pr;
    Annulus an;
    std::vector<std::pair<float, float> > others, none;
    double strip_full[kNS];
    for (int q = 0; q < kNS; q++) strip_full[q] = DiskSlice(-R + q * kStrip, -R + (q + 1) * kStrip);

    const double cone_area = M_PI * R * R;
    long n_out = 0;
    for (size_t ie = 0; ie < events.size() && n_out < max_jets; ie++)
    {
        const Event &e = events[ie];
        const double w = e.group < 0 ? (w_group.count(e.group) ? w_group[e.group] : 0.0) : 1.0;
        if (e.group < 0 && !(w > 0)) continue;
        for (int j = e.jet_first; j < e.jet_first + e.jet_n; j++)
        {
            const Jet &jet = jets[j];
            const double E = jet.pt * std::cosh(jet.eta);
            if (jet.pt < reco_jet_pt_min_em || E < 10 || E >= 120 || !ProxyInFiducial(jet.x, jet.y)) continue;
            t_pt = jet.pt; t_eta = jet.eta; t_phi = jet.phi; t_E = E; t_x = jet.x; t_y = jet.y;
            t_em = e.em; t_emu = e.emu; t_hw0 = e.hw0; t_group = e.group; t_w = w;
            geo.Near(jet.eta, jet.phi, R, cl);
            double A = 0, T = 0, Ajs[kNS] = {0}, Tjs[kNS] = {0};
            for (size_t q = 0; q < cl.size(); q++)
            {
                const Cell &c = geo.cells[cl[q]];
                A += c.area; T += rho[cl[q]] * c.area;
                const int sq = StripOf(c.eta - jet.eta);
                if (sq >= 0) { Ajs[sq] += c.area; Tjs[sq] += rho[cl[q]] * c.area; }
            }
            a_jet = A / cone_area; t_jet = T;
            double tej = 0;
            for (int i = e.t_first; i < e.t_first + e.t_n; i++) if (ProxyDeltaR(T_eta[i], T_phi[i], jet.eta, jet.phi) < R) tej += T_E[i];
            te_jet = tej;
            for (int k = 0; k < kNC; k++)
            {
                const double cphi = std::remainder(jet.phi + kRot[k] * M_PI / 180.0, 2 * M_PI);
                geo.Near(jet.eta, cphi, R, cl);
                double Ak = 0, Tk = 0, Ek = 0, TEk = 0, Aks[kNS] = {0};
                for (size_t q = 0; q < cl.size(); q++)
                {
                    const Cell &c = geo.cells[cl[q]];
                    Ak += c.area; Tk += rho[cl[q]] * c.area;
                    const int sq = StripOf(c.eta - jet.eta);
                    if (sq >= 0) Aks[sq] += c.area;
                }
                // strip weights and the template share of the jet cone the covered strips carry
                double ws[kNS], tcov = 0;
                for (int q = 0; q < kNS; q++)
                {
                    const bool cov = strip_full[q] > 0 && Aks[q] >= kStripMinLive * strip_full[q];
                    ws[q] = cov ? Ajs[q] / Aks[q] : 0;
                    if (cov) tcov += Tjs[q];
                }
                const double sfill = tcov > 0 ? T / tcov : 0;
                double Sk = 0;
                for (int i = e.pt_first; i < e.pt_first + e.pt_n; i++)
                    if (ProxyDeltaR(pts[i].eta, pts[i].phi, jet.eta, cphi) < R)
                    {
                        Ek += pts[i].E;
                        const int sq = StripOf(pts[i].eta - jet.eta);
                        if (sq >= 0) Sk += pts[i].E * ws[sq];
                    }
                // truth: full cone, live part of the cone, and the jet cone's live footprint moved here
                double TF = 0, TL = 0, TLS = 0;
                const double rot = kRot[k] * M_PI / 180.0;
                for (int i = e.t_first; i < e.t_first + e.t_n; i++)
                {
                    if (ProxyDeltaR(T_eta[i], T_phi[i], jet.eta, cphi) >= R) continue;
                    TEk += T_E[i];
                    if (CellOfDirection(geo, T_eta[i], T_phi[i]) >= 0)
                    {
                        TL += T_E[i];
                        const int sq = StripOf(T_eta[i] - jet.eta);
                        if (sq >= 0) TLS += T_E[i] * ws[sq];
                    }
                    if (CellOfDirection(geo, T_eta[i], T_phi[i] - rot) >= 0) TF += T_E[i];
                }
                e_[k] = Ek; a_[k] = Ak / cone_area; tt_[k] = Tk; te_[k] = TEk;
                s_[k] = Sk * sfill;
                tf_[k] = TF;
                tla_[k] = Ak > 0 ? TL * A / Ak : 0;
                tlt_[k] = Tk > 0 ? TL * T / Tk : 0;
                tls_[k] = TLS * sfill;
            }
            others.clear();
            for (int o = e.jet_first; o < e.jet_first + e.jet_n; o++)
                if (o != j && jets[o].pt >= reco_jet_pt_min_em) others.push_back(std::make_pair(jets[o].eta, jets[o].phi));
            for (int vv = 0; vv < kNAV; vv++)
            {
                const int v = vv % kNA;
                BuildAnnulus(geo, rho, jet.eta, jet.phi, vv < kNA ? others : none, Ajs, T, Tjs, v, an);
                double Ea = 0, TLa = 0;
                for (int i = e.pt_first; i < e.pt_first + e.pt_n; i++) if (pts[i].cell >= 0) Ea += pts[i].E * an.w[pts[i].cell];
                for (int i = e.t_first; i < e.t_first + e.t_n; i++)
                {
                    const int c = CellOfDirection(geo, T_eta[i], T_phi[i]);
                    if (c >= 0) TLa += T_E[i] * an.w[c];
                }
                double TFa = 0;
                for (size_t q = 0; q < an_rot[v].size(); q++)
                {
                    const double rot = an_rot[v][q];
                    for (int i = e.t_first; i < e.t_first + e.t_n; i++)
                        if (ProxyDeltaR(T_eta[i], T_phi[i], jet.eta, jet.phi + rot) < R && CellOfDirection(geo, T_eta[i], T_phi[i] - rot) >= 0) TFa += T_E[i];
                }
                an_e[vv] = Ea * an.fill; an_tl[vv] = TLa * an.fill; an_fill[vv] = an.fill;
                an_area[vv] = an.area_live; an_tf[vv] = an_rot[v].empty() ? 0 : TFa / an_rot[v].size();
            }
            for (int q = 0; q < kNSh + kNDm; q++)
            {
                const bool shell = q < kNSh;
                if (shell) BuildRegion(geo, rho, jet.eta, jet.phi, others, Ajs, T, Tjs, 0, M_PI, sh_edge[q], sh_edge[q + 1], kShellMinLive, an);
                else BuildRegion(geo, rho, jet.eta, jet.phi, others, Ajs, T, Tjs, 0, M_PI, dmin[q - kNSh], 1e9, kStripMinLive, an);
                double Ea = 0, TLa = 0;
                for (int i = e.pt_first; i < e.pt_first + e.pt_n; i++) if (pts[i].cell >= 0) Ea += pts[i].E * an.w[pts[i].cell];
                for (int i = e.t_first; i < e.t_first + e.t_n; i++)
                {
                    const int c = CellOfDirection(geo, T_eta[i], T_phi[i]);
                    if (c >= 0) TLa += T_E[i] * an.w[c];
                }
                if (shell) { sh_e[q] = Ea * an.fill; sh_tl[q] = TLa * an.fill; sh_fill[q] = an.fill; }
                else { dm_e[q - kNSh] = Ea * an.fill; dm_tl[q - kNSh] = TLa * an.fill; dm_fill[q - kNSh] = an.fill; }
            }
            for (int q = 0; q < kNL; q++) lr_E[q] = lr_A[q] = lrx_E[q] = lrx_A[q] = 0;
            {
                const double dmax = lr_edge[kNL] * R;
                for (size_t q = 0; q < geo.pix.size(); q++)
                {
                    const Pix &px = geo.pix[q];
                    if (std::fabs(px.eta - jet.eta) >= dmax) continue;
                    const double d = ProxyDeltaR(px.eta, px.phi, jet.eta, jet.phi) / R;
                    if (d < lr_edge[0] || d >= lr_edge[kNL]) continue;
                    int b = 0;
                    while (d >= lr_edge[b + 1]) b++;
                    lr_A[b] += px.area;
                    bool near = false;
                    for (size_t o = 0; o < others.size() && !near; o++) near = ProxyDeltaR(px.eta, px.phi, others[o].first, others[o].second) < R;
                    if (!near) lrx_A[b] += px.area;
                }
                for (int i = e.pt_first; i < e.pt_first + e.pt_n; i++)
                {
                    if (pts[i].cell < 0) continue;
                    const double d = ProxyDeltaR(pts[i].eta, pts[i].phi, jet.eta, jet.phi) / R;
                    if (d < lr_edge[0] || d >= lr_edge[kNL]) continue;
                    int b = 0;
                    while (d >= lr_edge[b + 1]) b++;
                    lr_E[b] += pts[i].E;
                    bool near = false;
                    for (size_t o = 0; o < others.size() && !near; o++) near = ProxyDeltaR(pts[i].eta, pts[i].phi, others[o].first, others[o].second) < R;
                    if (!near) lrx_E[b] += pts[i].E;
                }
            }
            for (int p = 0; p < 2; p++)
            {
                const std::vector<std::pair<float, float> > &oth = p == 0 ? others : none;
                BuildPixelRing(geo, jet.eta, jet.phi, oth, kAnnJetExcl, pr);
                px_area[p] = pr.area_ring; px_wmax[p] = pr.w_max;
                px_cov[p] = pr.area_jet > 0 ? pr.den / pr.area_jet : 0;
                px_fb[p] = pr.den > 0 ? 0 : 1;
                if (px_fb[p]) { px_e[p] = T; px_tl[p] = T; px_rho[p] = 0; continue; }  // no ring at the jet's eta: the template map's expectation
                double Ep = 0, TLp = 0;
                for (int i = e.pt_first; i < e.pt_first + e.pt_n; i++) Ep += pts[i].E * pr.Weight(pts[i].cell, pts[i].eta);
                for (int i = e.t_first; i < e.t_first + e.t_n; i++)
                {
                    const int c = CellOfDirection(geo, T_eta[i], T_phi[i]);
                    TLp += T_E[i] * pr.Weight(c, T_eta[i]);
                }
                px_rho[p] = Ep / pr.den;
                px_e[p] = px_rho[p] * pr.area_jet; px_tl[p] = TLp / pr.den * pr.area_jet;
            }
            tj.Fill();
            n_out++;
        }
    }
    // the map as a histogram over the detector plane
    TH2D h_map("h_map", "mean UE density away from jets;x (cm);y (cm)", 2 * geo.nx, -kFiducialRectEcal.x_outer, kFiducialRectEcal.x_outer, geo.ny, kFiducialRectEcal.y_min, kFiducialRectEcal.y_max);
    for (size_t c = 0; c < ncell; c++) h_map.Fill(geo.cells[c].x, geo.cells[c].y, rho[c]);
    out.cd(); tj.Write(); h_map.Write(); h_rho_eta.Write();
    out.Close();
    std::printf("in_event_ue: %ld jets written to %s\n", n_out, argv[2]);
    return 0;
}
