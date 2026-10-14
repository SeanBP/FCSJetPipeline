// Maps of the mean UE energy estimate against the jet's position on the
// detector, in ranges of jet pT.
//
// Input: as proxy_ue.cpp (a folder of jet_output and SimpleTree files; third
// argument "data" for data_to_jet output). Every EM-triggered reco jet gets
// its UE estimate from ../underlying_event/ProxySearch.h, with proxy events
// drawn from the whole run, and the estimate is averaged over the jets whose
// axis falls in each position bin. Bins with fewer than min_jets jets are
// left empty.
//
// Build and run from this folder:
//   g++ -m64 -O2 ue_position_map.cpp ../shared/VectorDict.cxx -o ue_position_map -I. -I../shared -I../jet_scale -I../underlying_event \
//       `root-config --cflags --libs` -std=c++11
//   ./ue_position_map <input_dir/> [output_base] [data] [bin_cm] [min_jets]
//
// Writes <output_base>.pdf/.png/.root (default ue_position_map; bins of 6 cm, min_jets 5).

#include <TROOT.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TLatex.h>
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
        std::cerr << "usage: ./ue_position_map <input_dir/> [output_base] [data] [bin_cm] [min_jets]" << std::endl;
        return 1;
    }
    std::string input_dir = argv[1];
    if (input_dir[input_dir.size() - 1] != '/') input_dir += '/';
    const std::string out_base = argc > 2 ? argv[2] : "ue_position_map";
    const double bin_cm = argc > 4 ? std::atof(argv[4]) : 6.0;
    const int min_jets = argc > 5 ? std::atoi(argv[5]) : 5;
    const bool is_data = argc > 3 && std::string(argv[3]) == "data";
    const int files_per_group = 0; // proxy events from the whole run

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

    // ---------------- pass 2: estimate and map ----------------
    const int kNPt = 6;
    const double kLo[kNPt] = {2.0, 3.0, 4.0, 5.0, 6.0, 2.0}, kHi[kNPt] = {3.0, 4.0, 5.0, 6.0, 12.0, 12.0}; // the last is all of them
    const FiducialRect &rect = kFiducialRectEcal;
    const double buf = reco_fiducial_buffer_ecal;
    const double x_hi = rect.x_outer - buf, y_lo = rect.y_min + buf, y_hi = rect.y_max - buf;
    // Bins of about bin_cm that line up with the edges of the fiducial region
    // of each half; the gap between the halves is one (always empty) bin
    const double x_in = rect.x_inner + buf;
    const int nx_half = (int)((x_hi - x_in) / bin_cm + 0.5), ny = (int)((y_hi - y_lo) / bin_cm + 0.5);
    const int nx = 2 * nx_half + 1;
    std::vector<double> x_edges;
    for (int i = 0; i <= nx_half; i++) x_edges.push_back(-x_hi + i * (x_hi - x_in) / nx_half);
    for (int i = 0; i <= nx_half; i++) x_edges.push_back(x_in + i * (x_hi - x_in) / nx_half);
    std::vector<double> y_edges;
    for (int j = 0; j <= ny; j++) y_edges.push_back(y_lo + j * (y_hi - y_lo) / ny);

    TFile fout((out_base + ".root").c_str(), "RECREATE");
    TH2D *h_sum[kNPt], *h_n[kNPt], *h_mean[kNPt];
    for (int k = 0; k < kNPt; k++)
    {
        h_sum[k] = new TH2D(Form("ue_sum_%d", k), "", nx, &x_edges[0], ny, &y_edges[0]);
        h_n[k] = new TH2D(Form("n_jets_%d", k), ";jet x (cm);jet y (cm);jets", nx, &x_edges[0], ny, &y_edges[0]);
        h_n[k]->Sumw2();
        h_mean[k] = new TH2D(Form("ue_mean_%d", k), ";Jet x (cm);Jet y (cm);Mean UE energy per cone (GeV)", nx, &x_edges[0], ny, &y_edges[0]);
    }
    double tot_w[kNPt] = {0}, tot_wx[kNPt] = {0}; long tot_n[kNPt] = {0};

    for (int ie = 0; ie < pool.NEvents(); ie++)
    {
        const ProxyPool::Event &pev = pool.GetEvent(ie);
        if (pev.em_level < 0) continue;
        for (int k_jet = 0; k_jet < pev.jet_n; k_jet++)
        {
            const ProxyJet &jet = pool.GetJet(pev.jet_first + k_jet);
            if (jet.pt < kLo[kNPt - 1] || jet.pt >= kHi[kNPt - 1]) continue;
            const ProxyResult res = pool.Estimate(ie, k_jet);
            if (res.n < proxy_n_min) continue; // no usable estimate
            for (int k = 0; k < kNPt; k++)
            {
                if (jet.pt < kLo[k] || jet.pt >= kHi[k]) continue;
                h_sum[k]->Fill(jet.x, jet.y, event_w[ie] * res.ue_E);
                h_n[k]->Fill(jet.x, jet.y, event_w[ie]);
                tot_w[k] += event_w[ie]; tot_wx[k] += event_w[ie] * res.ue_E; tot_n[k]++;
            }
        }
    }

    // In data every jet has weight 1, so the bin content of h_n is the number of jets
    double z_max = 0.0;
    for (int k = 0; k < kNPt; k++)
    {
        int n_filled = 0;
        for (int bx = 1; bx <= nx; bx++)
            for (int by = 1; by <= ny; by++)
            {
                const double n = h_n[k]->GetBinContent(bx, by);
                if (n <= 0 || h_n[k]->GetBinError(bx, by) <= 0) continue;
                // effective number of jets in the bin
                const double n_eff = n * n / (h_n[k]->GetBinError(bx, by) * h_n[k]->GetBinError(bx, by));
                if (n_eff < min_jets) continue;
                const double mean = h_sum[k]->GetBinContent(bx, by) / n;
                h_mean[k]->SetBinContent(bx, by, std::max(mean, 1e-4)); // a true zero still gets drawn
                if (k < kNPt - 1) z_max = std::max(z_max, mean);
                n_filled++;
            }
        printf("jet pT %4.0f-%-4.0f GeV: %7ld jets with an estimate, mean UE %.2f GeV, %d of %d position bins filled\n",
               kLo[k], kHi[k], tot_n[k], tot_w[k] > 0 ? tot_wx[k] / tot_w[k] : 0.0, n_filled, 2 * nx_half * ny);
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    gStyle->SetNumberContours(99);
    TCanvas c("c", "", 2700, 1700);
    c.Divide(3, 2);
    const double z_top = 6.0; // one colour scale for every panel
    for (int k = 0; k < kNPt; k++)
    {
        c.cd(k + 1);
        gPad->SetRightMargin(0.17); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.09);
        TH2D *h = h_mean[k];
        h->SetMinimum(0.0); h->SetMaximum(z_top);
        h->GetXaxis()->SetTitleSize(0.048); h->GetYaxis()->SetTitleSize(0.048); h->GetZaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.042); h->GetYaxis()->SetLabelSize(0.042); h->GetZaxis()->SetLabelSize(0.042);
        h->GetZaxis()->SetTitleOffset(1.15); h->GetYaxis()->SetTitleOffset(1.15);
        h->Draw("COLZ");
        for (int side = -1; side <= 1; side += 2)
        {
            TBox *fid = new TBox(side * x_in, y_lo, side * x_hi, y_hi);
            fid->SetFillStyle(0); fid->SetLineWidth(2); fid->Draw("l");
        }
        TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.046);
        t.DrawLatex(0.12, 0.93, Form("Jet p_{T} %.0f#minus%.0f GeV", kLo[k], kHi[k]));
        t.SetTextFont(42); t.SetTextSize(0.036); t.SetTextAlign(31);
        t.DrawLatex(0.83, 0.93, Form("%ld jets, mean %.2f GeV", tot_n[k], tot_w[k] > 0 ? tot_wx[k] / tot_w[k] : 0.0));
    }
    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());
    std::cout << "largest bin mean outside the all-pT panel: " << z_max << " GeV (colour scale tops at " << z_top << ")" << std::endl;

    fout.cd();
    for (int k = 0; k < kNPt; k++) { h_mean[k]->Write(); h_n[k]->Write(); }
    fout.Close();
    return 0;
}
