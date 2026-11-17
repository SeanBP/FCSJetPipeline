// Photon-class response from simulation, for the photon-class
// misidentification systematic (Sean 2026-11-17; as the FMS EM-jet analysis).
// Reads EM-only point-jet files of a sim_to_jet production (jet_output_*:
// JetMatcher, truth jets matched to reco jets) and fills, per x_F bin, the
// cross-section-weighted counts of reco photon class (reco_nconst, AnConfig.h
// PhotonClass: np0 ... np5plus) x true photon class (photons among the
// matched truth jet's constituents with E >= --truth-emin, default 1 GeV as
// the reco points; JetMatcher's truth_nphotons has no energy cut).
//
// Selection (as close to the A_N jets as the simulation allows):
//   - reco point jet (reco_from_points = 1) matched to a truth jet (the
//     pairs match_reco_idx[k] / match_truth_idx[k]); the files hold only reco jets inside the 20 cm
//     ECal fiducial (RecoJets / JetMatcher);
//   - truth pT > kJetPtCorrMin and truth E <= sqrt(s)/2: the JES map makes
//     the corrected reco energy the truth energy on average, so the truth
//     cuts stand in for the corrected-pT and E cuts (the raw sim files carry no
//     reco_E_corr); x_F bins are truth x_F;
//   - trigger (default on): the offline part of the A_N rule on the jet's
//     trigger values (reco_trg_*), at the lowest levels: EmTrgFactor(EM0) >=
//     kEmCutFactor or reco_trg_jp_et >= kJpCutFactor x the JP0 threshold (EM +
//     JP union). The hardware flags are not required: the JES production's
//     emulated flags and fcs_tcu_bits are all 0 (2026-11-17); in data the bit
//     is set for >= ~98% of the jets above 1.15 x x50. Which level recorded a
//     data event (prescales) is not emulated. Thresholds of --trg-run
//     (default 23066056, the period the simulation's trigger test was pinned
//     to); --no-trigger drops the cut (the inclusive class mix).
// Weights: JetXSecWeights.h (sigma_bin / sum N_job per ptHat bin, as the JES).
//
// Output: <out>.root with TH3D h_resp_w (x_F, reco class, true class;
// weighted, kNXfBins of AnConfig.h) and h_resp_n (unweighted), and h_unmatched_w
// (x_F of the reco jets with no truth match: their share per reco class is
// in the table); <out>.txt per x_F bin (merged to --xf-step, default 0.05)
// and reco class: the effective entries, the matched fraction, and the
// fraction of each true class (f(true m | reco n)), the purity f(n | n) last;
// <out>_merged.txt the same for the merged classes np12 (1-2) and np3plus (>= 3).
//
// Usage: ./misid_response <jet file list> <out prefix> [--truth-emin 1.0]
//                         [--trg-run 23066056 | --no-trigger] [--xf-step 0.05]
// Build: g++ -m64 -O2 -I../../shared -I.. -I../../jet_scale misid_response.cpp ../../shared/VectorDict.cxx
//            -o misid_response `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "JetParameters.h"
#include "VectorDict.h"
#include "AnConfig.h"
#include "JetXSecWeights.h"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <jet file list> <out prefix> [--truth-emin 1.0] [--trg-run 23066056 | --no-trigger] [--xf-step 0.05]" << std::endl;
        return 1;
    }
    double truth_emin = 1.0, xf_step = 0.05;
    int trg_run = 23066056;
    bool use_trg = true;
    for (int i = 3; i < argc; i++)
    {
        const std::string a = argv[i];
        if (a == "--truth-emin" && i + 1 < argc) truth_emin = std::atof(argv[++i]);
        else if (a == "--trg-run" && i + 1 < argc) trg_run = std::atoi(argv[++i]);
        else if (a == "--no-trigger") use_trg = false;
        else if (a == "--xf-step" && i + 1 < argc) xf_step = std::atof(argv[++i]);
        else { std::cerr << "unknown option " << a << std::endl; return 1; }
    }
    std::vector<std::string> files;
    {
        std::ifstream in(argv[1]);
        std::string f;
        while (in >> f) files.push_back(f);
    }
    if (files.empty()) { std::cerr << "ERROR: no files in " << argv[1] << std::endl; return 1; }
    // truncated sim jobs (sigma_job = -1: the production's post-processing drops them) are left out
    {
        std::vector<std::string> ok;
        for (size_t i = 0; i < files.size(); i++)
        {
            TFile f(files[i].c_str());
            TTree *t = f.IsZombie() ? 0 : (TTree *)f.Get("jetTree");
            if (!t || t->GetEntries() == 0 || !t->GetBranch("sigma_job")) continue;
            Double_t sj = -1;
            t->SetBranchStatus("*", 0); t->SetBranchStatus("sigma_job", 1); t->SetBranchAddress("sigma_job", &sj);
            t->GetEntry(0);
            if (sj == -1) { std::cerr << "skip (truncated job, sigma_job = -1) " << files[i] << std::endl; continue; }
            ok.push_back(files[i]);
        }
        files.swap(ok);
    }
    std::vector<float> wfile;
    if (!ComputeXSecWeights(files, "jetTree", wfile)) return 1;

    const int NC = kNPhotonClasses;
    TFile out((std::string(argv[2]) + ".root").c_str(), "RECREATE");
    double xa[kNXfBins + 1], ca[NC + 1];
    for (int i = 0; i <= kNXfBins; i++) xa[i] = kXfMin + i * (kXfMax - kXfMin) / kNXfBins;
    for (int i = 0; i <= NC; i++) ca[i] = i;
    TH3D *hw = new TH3D("h_resp_w", "photon-class response (weighted);truth x_{F};reco class;true class", kNXfBins, xa, NC, ca, NC, ca);
    TH3D *hn = new TH3D("h_resp_n", "photon-class response (entries);truth x_{F};reco class;true class", kNXfBins, xa, NC, ca, NC, ca);
    TH3D *hw2 = new TH3D("h_resp_w2", "sum of weights^2;truth x_{F};reco class;true class", kNXfBins, xa, NC, ca, NC, ca);
    TH1D *hum[NC];
    for (int k = 0; k < NC; k++)
        hum[k] = new TH1D(Form("h_unmatched_w_%s", PhotonClassName(k)), "reco jets without a truth match (weighted);reco x_{F} (raw);jets", kNXfBins, kXfMin, kXfMax);
    for (int k = 1; k <= NC; k++) { hw->GetYaxis()->SetBinLabel(k, PhotonClassName(k - 1)); hw->GetZaxis()->SetBinLabel(k, PhotonClassName(k - 1)); }

    long n_jets = 0, n_sel = 0, n_unm = 0;
    for (size_t fi = 0; fi < files.size(); fi++)
    {
        if (wfile[fi] < 0) continue;
        TFile f(files[fi].c_str());
        TTree *t = f.IsZombie() ? 0 : (TTree *)f.Get("jetTree");
        if (!t) continue;
        std::vector<float> *r_E = 0, *r_xf = 0, *t_pt = 0, *t_E = 0, *t_xf = 0, *trg_em = 0, *trg_hmax = 0, *trg_jp = 0;
        std::vector<float> *tc_E = 0;
        std::vector<int> *nconst = 0, *mt = 0, *mr = 0, *tc_jet = 0, *tc_pid = 0;
        Int_t from_points = -1;
        t->SetBranchAddress("reco_E", &r_E);
        t->SetBranchAddress("reco_x_F", &r_xf);
        t->SetBranchAddress("truth_pt", &t_pt);
        t->SetBranchAddress("truth_E", &t_E);
        t->SetBranchAddress("truth_x_F", &t_xf);
        t->SetBranchAddress("reco_nconst", &nconst);
        t->SetBranchAddress("match_truth_idx", &mt);
        t->SetBranchAddress("match_reco_idx", &mr);
        t->SetBranchAddress("truth_const_E", &tc_E);
        t->SetBranchAddress("truth_const_jetIndex", &tc_jet);
        t->SetBranchAddress("truth_const_pid", &tc_pid);
        t->SetBranchAddress("reco_from_points", &from_points);
        t->SetBranchAddress("reco_trg_em_et", &trg_em);
        t->SetBranchAddress("reco_trg_em_hmax", &trg_hmax);
        t->SetBranchAddress("reco_trg_jp_et", &trg_jp);
        const double w = wfile[fi];
        for (Long64_t i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);
            if (from_points != 1) { std::cerr << "ERROR: " << files[fi] << " has hit jets; photon counting needs point jets" << std::endl; return 1; }
            // true photons per truth jet, E >= truth_emin
            std::vector<int> ntrue(t_E->size(), 0);
            for (size_t c = 0; c < tc_jet->size(); c++)
            {
                const int j = (*tc_jet)[c];
                if (j >= 0 && j < (int)ntrue.size() && (*tc_pid)[c] == 22 && (*tc_E)[c] >= truth_emin) ntrue[j]++;
            }
            // reco jet -> matched truth jet
            std::vector<int> match(r_E->size(), -1);
            for (size_t k = 0; k < mr->size() && k < mt->size(); k++)
                if ((*mr)[k] >= 0 && (*mr)[k] < (int)match.size()) match[(*mr)[k]] = (*mt)[k];
            for (size_t j = 0; j < r_E->size(); j++)
            {
                n_jets++;
                // trigger: the A_N rule (EM + JP union) on the emulated values and flags
                if (use_trg)
                {
                    const bool em_on = EmTrgFactor((*trg_em)[j], (*trg_hmax)[j], 0, trg_run) >= kEmCutFactor - 1e-9;
                    const bool jp_on = (*trg_jp)[j] >= kJpCutFactor * TrgThresholdGeV(3, trg_run) - 1e-9;
                    if (!em_on && !jp_on) continue;
                }
                const int cls = PhotonClass((*nconst)[j]);
                const int ti = match[j];
                if (ti < 0 || ti >= (int)t_E->size())
                {
                    hum[cls]->Fill((*r_xf)[j], w);
                    n_unm++;
                    continue;
                }
                const double tpt = (*t_pt)[ti], tE = (*t_E)[ti], txf = (*t_xf)[ti];
                if (!(tpt > kJetPtCorrMin) || !(txf > 0)) continue;
                // sqrt(s) from the truth jet itself: x_F = 2 pz / sqrt(s)
                const double pz = std::sqrt(std::max(0.0, tE * tE - tpt * tpt));
                const double sqrts = pz > 0 ? 2 * pz / txf : 0;
                if (sqrts > 0 && tE > 0.5 * sqrts) continue;
                const int tcls = PhotonClass(ntrue[ti]);
                hw->Fill(txf, cls + 0.5, tcls + 0.5, w);
                hw2->Fill(txf, cls + 0.5, tcls + 0.5, w * w);
                hn->Fill(txf, cls + 0.5, tcls + 0.5);
                n_sel++;
            }
        }
    }

    // table: per merged x_F bin and reco class, f(true m | reco n)
    FILE *tx = std::fopen((std::string(argv[2]) + ".txt").c_str(), "w");
    std::fprintf(tx, "# photon-class response (misid_response.cpp): %zu files, truth photons E >= %.2f GeV, trigger %s, x_F = truth x_F\n",
                 files.size(), truth_emin, use_trg ? Form("offline A_N rule, EM0 factor >= %.2f or JP0 patch >= %.2f x threshold (run %d periods)", kEmCutFactor, kJpCutFactor, trg_run) : "none");
    std::fprintf(tx, "# f_<m> = weighted fraction of the reco class's matched jets with true class m (purity = f of its own class); n_eff = (sum w)^2 / sum w^2;\n");
    std::fprintf(tx, "# unmatched = weighted fraction of the reco class's jets (that bin's raw reco x_F) without a truth match, not in the f_<m>\n");
    std::fprintf(tx, "# xF_lo xF_hi reco_class n_entries n_eff unmatched");
    for (int m = 0; m < NC; m++) std::fprintf(tx, " f_%s", PhotonClassName(m));
    std::fprintf(tx, " purity\n");
    const int step = std::max(1, (int)std::lround(xf_step / ((kXfMax - kXfMin) / kNXfBins)));
    for (int b0 = 1; b0 <= kNXfBins; b0 += step)
    {
        const int b1 = std::min(kNXfBins, b0 + step - 1);
        for (int n = 0; n < NC; n++)
        {
            double sw = 0, sw2 = 0, cnt = 0, fm[NC] = {0}, um = 0;
            for (int b = b0; b <= b1; b++)
            {
                for (int m = 0; m < NC; m++)
                {
                    const double v = hw->GetBinContent(b, n + 1, m + 1);
                    fm[m] += v; sw += v; sw2 += hw2->GetBinContent(b, n + 1, m + 1); cnt += hn->GetBinContent(b, n + 1, m + 1);
                }
                um += hum[n]->GetBinContent(b);
            }
            if (cnt < 1) continue;
            std::fprintf(tx, "%.3f %.3f %s %.0f %.1f %.4f", hw->GetXaxis()->GetBinLowEdge(b0), hw->GetXaxis()->GetBinUpEdge(b1), PhotonClassName(n), cnt,
                         sw2 > 0 ? sw * sw / sw2 : 0, sw + um > 0 ? um / (sw + um) : 0);
            for (int m = 0; m < NC; m++) std::fprintf(tx, " %.4f", sw > 0 ? fm[m] / sw : 0);
            std::fprintf(tx, " %.4f\n", sw > 0 ? fm[n] / sw : 0);
        }
    }
    std::fclose(tx);
    // merged classes (Sean 2026-11-17: the systematic is reported for the merged
    // classes np12 = 1-2 photons and np3plus = 3 or more; np0 left out)
    {
        FILE *tm = std::fopen((std::string(argv[2]) + "_merged.txt").c_str(), "w");
        std::fprintf(tm, "# merged photon classes (np12 = 1-2, np3plus = >= 3; np0 left out), same selection as %s.txt\n", argv[2]);
        std::fprintf(tm, "# xF_lo xF_hi reco_class n_entries n_eff f_true_np12 f_true_np3plus purity (f of the true np0: the rest)\n");
        const char *gname[2] = {"np12", "np3plus"};
        auto grp = [](int c) { return c <= 0 ? -1 : c <= 2 ? 0 : 1; };
        for (int b0 = 1; b0 <= kNXfBins; b0 += step)
        {
            const int b1 = std::min(kNXfBins, b0 + step - 1);
            for (int g = 0; g < 2; g++)
            {
                double sw = 0, sw2 = 0, cnt = 0, fg[2] = {0, 0};
                for (int b = b0; b <= b1; b++)
                    for (int n = 0; n < NC; n++)
                    {
                        if (grp(n) != g) continue;
                        for (int m = 0; m < NC; m++)
                        {
                            const double v = hw->GetBinContent(b, n + 1, m + 1);
                            sw += v; sw2 += hw2->GetBinContent(b, n + 1, m + 1); cnt += hn->GetBinContent(b, n + 1, m + 1);
                            if (grp(m) >= 0) fg[grp(m)] += v;
                        }
                    }
                if (cnt < 1) continue;
                std::fprintf(tm, "%.3f %.3f %s %.0f %.1f %.4f %.4f %.4f\n", hw->GetXaxis()->GetBinLowEdge(b0), hw->GetXaxis()->GetBinUpEdge(b1), gname[g], cnt,
                             sw2 > 0 ? sw * sw / sw2 : 0, sw > 0 ? fg[0] / sw : 0, sw > 0 ? fg[1] / sw : 0, sw > 0 ? fg[g] / sw : 0);
            }
        }
        std::fclose(tm);
    }
    out.Write();
    out.Close();
    printf("misid_response: %zu files, %ld reco jets, %ld selected and matched, %ld unmatched (after the trigger cut) -> %s.{root,txt}\n",
           files.size(), n_jets, n_sel, n_unm, argv[2]);
    return 0;
}
