// Runs the near off-axis proxy-event UE estimate (../underlying_event/
// ProxySearch.h, as the pipeline step) jet by jet and reports its cost and
// behaviour.
//
// Input: a folder of sim_to_jet -E 1 -s 1 output (SimpleTree_<s>.root next to
// jet_output_<s>.root), treated as one group. With "data", data_to_jet output
// (jet_output_<run>_<file>.root next to SimpleTree_mudst_<run>_<file>.root,
// EM-triggered): weight 1, grouped by run, files in time order. Without
// simpletree_entry, jetTree and SimpleTree entries are matched by repeating
// data_to_jet's event tests; files that then disagree are skipped.
//
// Build and run from this folder:
//   g++ -m64 -O2 proxy_ue.cpp ../shared/VectorDict.cxx -o proxy_ue -I. -I../shared -I../jet_scale -I../underlying_event \
//       `root-config --cflags --libs` -std=c++11
//   ./proxy_ue <input_dir/> [n_proxy] [data] [files_per_group]
// files_per_group > 0 restricts each jet's proxy events to its own group of
// that many consecutive files, to see what a job limited to its own files
// could do; it replaces the same-run requirement, so use it on one run.

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"
#include "ProxySearch.h"

static double Percentile(std::vector<double> &v, double p)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[(size_t)(p * (v.size() - 1))];
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2)
    {
        std::cerr << "usage: ./proxy_ue <input_dir/> [n_proxy]" << std::endl;
        return 1;
    }
    std::string input_dir = argv[1];
    if (input_dir[input_dir.size() - 1] != '/') input_dir += '/';
    const int n_proxy_arg = argc > 2 ? std::atoi(argv[2]) : proxy_n_max;
    const int n_proxy = n_proxy_arg > 0 ? n_proxy_arg : 2000000000; // 0 = no limit
    const bool is_data = argc > 3 && std::string(argv[3]) == "data";
    // Optional: proxy events only from the same group of this many consecutive
    // files (in time order), as if each group were one job that sees nothing
    // else. 0 = the whole run.
    const int files_per_group = argc > 4 ? std::atoi(argv[4]) : 0;

    std::vector<std::string> jet_files, simple_files;
    std::vector<int> file_run;
    {
        TSystemDirectory dir("input_dir", input_dir.c_str());
        TList *files = dir.GetListOfFiles();
        if (!files) { std::cerr << "ERROR: cannot list " << input_dir << std::endl; return 1; }
        TIter next(files);
        TSystemFile *file;
        std::vector<std::string> names;
        while ((file = (TSystemFile *)next()))
        {
            const std::string n = file->GetName();
            if (n.compare(0, 10, "jet_output") == 0 && n.size() > 5 && n.compare(n.size() - 5, 5, ".root") == 0)
                names.push_back(n);
        }
        std::sort(names.begin(), names.end());
        if (is_data)
        {
            // jet_output_<run>_<file>.root: order by run, then by the file's
            // sequence number (the last five digits of <file>), which follows time
            std::vector<std::pair<std::pair<int, int>, std::string> > keyed;
            for (size_t i = 0; i < names.size(); i++)
            {
                int run = 0, file_number = 0;
                // runMudst.C names st_fwd files <last four digits of the run>_raw_<file>
                if (std::sscanf(names[i].c_str(), "jet_output_%d_raw_%d.root", &run, &file_number) != 2 &&
                    std::sscanf(names[i].c_str(), "jet_output_%d_%d.root", &run, &file_number) != 2) continue;
                keyed.push_back(std::make_pair(std::make_pair(run, file_number % 100000), names[i]));
            }
            std::sort(keyed.begin(), keyed.end());
            names.clear();
            for (size_t i = 0; i < keyed.size(); i++) { names.push_back(keyed[i].second); file_run.push_back(keyed[i].first.first); }
        }
        for (size_t i = 0; i < names.size(); i++)
        {
            jet_files.push_back(input_dir + names[i]);
            simple_files.push_back(input_dir + (is_data ? "SimpleTree_mudst" : "SimpleTree") + names[i].substr(10));
        }
    }
    std::cout << "Found " << jet_files.size() << " file pairs in " << input_dir << (is_data ? " (data)" : " (simulation)") << std::endl;

    std::vector<float> file_weight;
    if (is_data)
        file_weight.assign(jet_files.size(), 1.0f);
    else if (!ComputeXSecWeights(jet_files, "jetTree", file_weight))
        return 1;
    if (!is_data) file_run.assign(jet_files.size(), 0);

    // ---------------- pass 1: read and index ----------------
    const clock_t t_start = clock();
    ProxyPool pool;
    std::vector<float> event_w;                  // weight of each pool event
    const int MAX_HITS = 10000; // as in data_to_jet/RecoJets.cpp
    int n_files_used = 0, n_files_skipped = 0, n_files_with_entry = 0;

    const int MAX_PTS = 200; // matches StSimpleReaderMaker.h's Cal_pt_* array size
    Long64_t n_read = 0, n_simple_read = 0;
    for (size_t f = 0; f < jet_files.size(); f++)
    {
        if (file_weight[f] < 0) continue;
        TFile *fj = TFile::Open(jet_files[f].c_str());
        TFile *fs = TFile::Open(simple_files[f].c_str());
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        // Newer jetTrees say which SimpleTree entry each of theirs came from
        const bool has_entry = jt && jt->GetBranch("simpletree_entry");
        bool ok = jt && st && st->GetBranch("Cal_npts") && (is_data || has_entry || jt->GetEntries() == st->GetEntries());
        if (ok && has_entry) { n_files_with_entry++; n_simple_read += st->GetEntries(); }

        // Data: the SimpleTree entries that have a jetTree entry, in order
        std::vector<Long64_t> st_entry;
        if (ok && is_data && !has_entry)
        {
            Int_t nhits = 0, s0 = 0, s1 = 0, s2 = 0;
            st->SetBranchStatus("*", 0);
            const char *fb[] = {"Cal_nhits", "fcsEM0", "fcsEM1", "fcsEM2"};
            for (int i = 0; i < 4; i++) st->SetBranchStatus(fb[i], 1);
            st->SetBranchAddress("Cal_nhits", &nhits);
            st->SetBranchAddress("fcsEM0", &s0);
            st->SetBranchAddress("fcsEM1", &s1);
            st->SetBranchAddress("fcsEM2", &s2);
            const Long64_t ns = st->GetEntries();
            for (Long64_t i = 0; i < ns; i++)
            {
                st->GetEntry(i);
                if (nhits <= MAX_HITS && (s0 || s1 || s2)) st_entry.push_back(i);
            }
            n_simple_read += ns;
            if ((Long64_t)st_entry.size() != jt->GetEntries())
            {
                std::cerr << "WARNING: " << jet_files[f] << " has " << jt->GetEntries() << " jetTree entries but "
                          << st_entry.size() << " SimpleTree entries pass the same filter" << std::endl;
                ok = false;
            }
            st->ResetBranchAddresses();
        }
        if (!ok)
        {
            std::cerr << "WARNING: skipping " << jet_files[f] << std::endl;
            n_files_skipped++;
            if (fj) { fj->Close(); delete fj; }
            if (fs) { fs->Close(); delete fs; }
            continue;
        }
        n_files_used++;

        std::vector<float> *reco_pt = 0, *reco_eta = 0, *reco_phi = 0, *reco_x = 0, *reco_y = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_x", "reco_y", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int i = 0; i < 8; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &reco_pt);
        jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi);
        jt->SetBranchAddress("reco_x", &reco_x);
        jt->SetBranchAddress("reco_y", &reco_y);
        jt->SetBranchAddress("fcsEM0", &em0);
        jt->SetBranchAddress("fcsEM1", &em1);
        jt->SetBranchAddress("fcsEM2", &em2);
        Long64_t simpletree_entry = -1;
        if (has_entry) { jt->SetBranchStatus("simpletree_entry", 1); jt->SetBranchAddress("simpletree_entry", &simpletree_entry); }

        Int_t Cal_npts = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        Int_t Cal_pt_detid[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pb[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pb[i], 1);
        st->SetBranchAddress("Cal_npts", &Cal_npts);
        st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x);
        st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z);
        st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);

        const double sqrt_s = is_data ? computeSqrtSForRun(file_run[f]) : kSqrtS_MainRun22_GeV;
        const Long64_t n = jt->GetEntries();
        for (Long64_t ev = 0; ev < n; ev++)
        {
            if (jt->GetEntry(ev) <= 0) continue;
            n_read++;
            if (reco_pt->empty()) continue;
            if (st->GetEntry(has_entry ? simpletree_entry : is_data ? st_entry[ev] : ev) <= 0) continue;

            const int em_level = ProxyEmLevel(em0, em1, em2);
            std::vector<ProxyJet> ev_jets;
            for (size_t j = 0; j < reco_pt->size(); j++)
            {
                ProxyJet jet = {reco_pt->at(j), reco_eta->at(j), reco_phi->at(j), reco_x->at(j), reco_y->at(j)};
                ev_jets.push_back(jet);
            }
            // Points: same selection and direction convention as the reco jets (use_points)
            std::vector<float> ev_eta, ev_phi, ev_E;
            for (int i = 0; i < Cal_npts; i++)
            {
                if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                if (!pass_point_energy_cut(Cal_pt_energy[i], sqrt_s)) continue;
                const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                ev_eta.push_back(std::asinh(Cal_pt_z[i] / rho)); ev_phi.push_back(std::atan2(Cal_pt_y[i], Cal_pt_x[i]));
                ev_E.push_back(Cal_pt_energy[i]);
            }
            pool.AddEvent(files_per_group > 0 ? (int)f / files_per_group : file_run[f], em_level, ev_jets, ev_eta, ev_phi, ev_E);
            event_w.push_back(file_weight[f]);
        }
        fj->Close(); delete fj;
        fs->Close(); delete fs;
    }
    const double t_read = double(clock() - t_start) / CLOCKS_PER_SEC;
    if (is_data)
        std::cout << n_files_used << " file pairs used (" << n_files_with_entry << " with simpletree_entry), " << n_files_skipped << " skipped; " << n_simple_read
                  << " recorded events, of which " << n_read << " fired an EM trigger" << std::endl;
    std::cout << "Read " << n_read << " events; kept " << pool.NEvents() << " with a reco jet (" << pool.NJets() << " jets, "
              << pool.NPoints() << " points) in " << pool.MemoryMB() << " MB; " << pool.NIndexBins() << " index bins; " << t_read << " s" << std::endl;

    // ---------------- pass 2: jet by jet ----------------
    const int kNPt = 5;
    const double kPt[kNPt + 1] = {2.0, 3.0, 4.0, 6.0, 8.0, 12.0};
    const int kNSub = 3; const int kSub[kNSub] = {2, 5, 10}; // UE also from the first few proxies only
    double sw[kNPt] = {0}, sw_full[kNPt] = {0}, sw_none[kNPt] = {0}, sw_found[kNPt] = {0};
    double oc_w[kNPt][4] = {{0}}, oc_window[kNPt] = {0}, fin_w[kNPt] = {0}, fin_wx[kNPt] = {0}; // outcomes of the search rule
    double lv_w[kNPt][3] = {{0}}, lv_full[kNPt][3] = {{0}}, lv_found[kNPt][3] = {{0}}; // the same by EM level
    double ue_w[kNPt][kNSub] = {{0}}, ue_wx[kNPt][kNSub] = {{0}}, ue_wxx[kNPt][kNSub] = {{0}};
    long n_jets_pt[kNPt] = {0};
    std::vector<double> scanned[kNPt], reach[kNPt];
    long n_done = 0; double total_scanned = 0;

    const clock_t t_loop = clock();
    for (int ie = 0; ie < pool.NEvents(); ie++)
    {
        const ProxyPool::Event &pev = pool.GetEvent(ie);
        if (pev.em_level < 0) continue;
        struct { double w; int em_level; } ev = {event_w[ie], pev.em_level};
        for (int k_jet = 0; k_jet < pev.jet_n; k_jet++)
        {
            const ProxyJet &jet = pool.GetJet(pev.jet_first + k_jet);
            int pb = -1;
            for (int k = 0; k < kNPt; k++) if (jet.pt >= kPt[k] && jet.pt < kPt[k + 1]) pb = k;
            if (pb < 0 || !ProxyInFiducial(jet.x, jet.y)) continue;

            // The search itself is ProxySearch.h's, the same one the pipeline runs
            const ProxyResult res = pool.Estimate(ie, k_jet, n_proxy);
            const std::vector<double> &ue_list = res.E_list;
            const long n_scanned = res.n_scanned;
            const int last_distance = res.reach, n_nominal = res.n_nominal;
            const double window = res.window;
            const bool widened = n_nominal < std::min(proxy_n_min, n_proxy);

            // Outcome: 0 = target reached, 1 = fewer but at least the minimum,
            // 2 = minimum reached by widening, 3 = still below the minimum
            const int outcome = n_nominal >= n_proxy ? 0 : !widened ? 1 : (int)ue_list.size() >= std::min(proxy_n_min, n_proxy) ? 2 : 3;
            oc_w[pb][outcome] += ev.w;
            if (outcome == 2) { oc_window[pb] += ev.w * window; }
            {
                double mean = 0.0;
                for (size_t i = 0; i < ue_list.size(); i++) mean += ue_list[i];
                if (!ue_list.empty()) { mean /= ue_list.size(); fin_w[pb] += ev.w; fin_wx[pb] += ev.w * mean; }
            }
            n_done++; total_scanned += n_scanned;
            n_jets_pt[pb]++;
            sw[pb] += ev.w;
            if ((int)ue_list.size() == n_proxy) { sw_full[pb] += ev.w; scanned[pb].push_back((double)n_scanned); reach[pb].push_back((double)last_distance); }
            if (ue_list.empty()) sw_none[pb] += ev.w;
            sw_found[pb] += ev.w * ue_list.size();
            lv_w[pb][ev.em_level] += ev.w; lv_found[pb][ev.em_level] += ev.w * ue_list.size();
            if ((int)ue_list.size() == n_proxy) lv_full[pb][ev.em_level] += ev.w;
            for (int s = 0; s < kNSub; s++)
            {
                if ((int)ue_list.size() < kSub[s]) continue;
                double mean = 0.0;
                for (int i = 0; i < kSub[s]; i++) mean += ue_list[i];
                mean /= kSub[s];
                ue_w[pb][s] += ev.w; ue_wx[pb][s] += ev.w * mean; ue_wxx[pb][s] += ev.w * mean * mean;
            }
        }
    }
    const double t_jets = double(clock() - t_loop) / CLOCKS_PER_SEC;

    if (files_per_group > 0)
        std::cout << "Proxy events restricted to groups of " << files_per_group << " consecutive files" << std::endl;
    std::cout << "\nProxy search for " << n_done << " jets (EM-triggered, fiducial, pT 2-12 GeV), up to " << n_proxy << " proxy events each: "
              << t_jets << " s, " << 1e6 * t_jets / std::max(1L, n_done) << " microseconds per jet, "
              << total_scanned / std::max(1L, n_done) << " candidates examined per jet on average" << std::endl;
    printf("\n  %-10s %8s %10s %10s %10s | %-22s | %-24s | %s\n", "jet pT", "jets", "all found", "none", "mean n",
           "candidates examined", "reach in events", "UE per cone (GeV): mean and jet-to-jet spread");
    printf("  %-10s %8s %10s %10s %10s | %-22s | %-24s |", "(GeV)", "", "", "", "", "median / 90%", "median / 90%");
    for (int s = 0; s < kNSub; s++) printf(" first %-2d          ", kSub[s]);
    printf("\n");
    for (int k = 0; k < kNPt; k++)
    {
        if (sw[k] <= 0) continue;
        const double sc50 = Percentile(scanned[k], 0.5), sc90 = Percentile(scanned[k], 0.9);
        const double re50 = Percentile(reach[k], 0.5), re90 = Percentile(reach[k], 0.9);
        printf("  %4.0f-%-5.0f %8ld %9.1f%% %9.1f%% %10.1f | %8.0f / %-11.0f | %9.0f / %-12.0f |", kPt[k], kPt[k + 1], n_jets_pt[k],
               100.0 * sw_full[k] / sw[k], 100.0 * sw_none[k] / sw[k], sw_found[k] / sw[k], sc50, sc90, re50, re90);
        for (int s = 0; s < kNSub; s++)
        {
            if (ue_w[k][s] <= 0) { printf(" %-18s", "-"); continue; }
            const double mean = ue_wx[k][s] / ue_w[k][s];
            printf(" %5.2f +- %-9.2f", mean, std::sqrt(std::max(0.0, ue_wxx[k][s] / ue_w[k][s] - mean * mean)));
        }
        printf("\n");
    }
    printf("\n  Outcome of the rule (target %d, minimum %d, pT window +-%.0f%% widened in %.0f%% steps to at most +-%.0f%%):\n",
           n_proxy, proxy_n_min, 100.0 * proxy_pt_window_fraction, 100.0 * proxy_pt_window_step, 100.0 * proxy_pt_window_max);
    printf("  %-10s %12s %14s %18s %16s %20s %18s\n", "jet pT", "target", "min to target", "min by widening", "below min", "mean widened window", "UE, all jets (GeV)");
    for (int k = 0; k < kNPt; k++)
    {
        if (sw[k] <= 0) continue;
        printf("  %4.0f-%-5.0f %11.2f%% %13.2f%% %17.2f%% %15.2f%%", kPt[k], kPt[k + 1], 100.0 * oc_w[k][0] / sw[k], 100.0 * oc_w[k][1] / sw[k],
               100.0 * oc_w[k][2] / sw[k], 100.0 * oc_w[k][3] / sw[k]);
        if (oc_w[k][2] > 0) printf(" %18.0f%%", 100.0 * oc_window[k] / oc_w[k][2]); else printf(" %19s", "-");
        printf(" %18.2f\n", fin_w[k] > 0 ? fin_wx[k] / fin_w[k] : 0.0);
    }
    printf("\n  By highest EM level: share of the jets, share with all %d found, mean number found\n", n_proxy);
    printf("  %-10s | %-28s | %-28s | %-28s\n", "jet pT", "level 0 (fcsEM0)", "level 1 (fcsEM1)", "level 2 (fcsEM2)");
    for (int k = 0; k < kNPt; k++)
    {
        if (sw[k] <= 0) continue;
        printf("  %4.0f-%-5.0f", kPt[k], kPt[k + 1]);
        for (int l = 0; l < 3; l++)
        {
            if (lv_w[k][l] <= 0) { printf(" | %-28s", "-"); continue; }
            printf(" | %5.1f%%  %5.1f%%  %4.1f         ", 100.0 * lv_w[k][l] / sw[k], 100.0 * lv_full[k][l] / lv_w[k][l], lv_found[k][l] / lv_w[k][l]);
        }
        printf("\n");
    }
    std::cout << "\n(\"all found\", \"none\" and the UE columns are cross-section weighted; the jet counts and the search statistics are raw."
                 " Reach is the distance in kept events from the jet's own event to its last proxy.)" << std::endl;
    return 0;
}
