// What does the trigger requirement on a proxy event cost, and does the UE care?
//
// Input: as proxy_pt_window.cpp. Uses the simulated FCS trigger bits in the
// jetTree (no prescales) and single-reco-jet events.
//
// For jets in bins of pT it prints, per trigger: how often it fired, the UE
// energy away from the jet (measured as in proxy_pt_window.cpp) when it fired
// and when it did not, and where the jets that fired it sit. It also counts
// the distinct trigger-bit patterns, which is how finely an "all triggers
// must match" requirement would divide the proxy events.
//
// Build and run from this folder:
//   g++ -m64 -O2 proxy_trigger.cpp ../shared/VectorDict.cxx -o proxy_trigger -I. -I../shared -I../jet_scale -I../underlying_event \
//       `root-config --cflags --libs` -std=c++11
//   ./proxy_trigger <input_dir/>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <cstdio>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"
#include "ProxySelection.h"

static const double kBandHalfWidthDeg = 45.0;  // azimuthal band about 90 deg where the UE is sampled
static const double kMinArea = 0.3;            // least eta-phi area the detector must offer there
static const double kCellCm = 4.0;             // cells used to measure that area


static const int kNPt = 4;                     // jet pT bins
static const double kPt[kNPt + 1] = {2.0, 3.0, 4.0, 6.0, 12.0};

static const int kNTrig = 11;                  // the first kNJP are the jet-patch triggers
static const int kNJP = 7;
static const char *kTrig[kNTrig] = {"fcsJP2", "fcsJPA1", "fcsJPBC1", "fcsJPDE1", "fcsJPA0", "fcsJPBC0", "fcsJPDE0",
                                    "fcsEM0", "fcsEM1", "fcsEM2", "fcsEM3"};

struct Cell { double eta, phi, area; };

// Weighted mean with its error
struct Mean
{
    double w, wx, wxx, w2; long n;
    Mean() : w(0), wx(0), wxx(0), w2(0), n(0) {}
    void Fill(double x, double wt) { w += wt; wx += wt * x; wxx += wt * x * x; w2 += wt * wt; n++; }
    double Value() const { return w > 0 ? wx / w : 0.0; }
    double Error() const
    {
        if (w <= 0 || n < 2) return 0.0;
        const double m = wx / w, var = std::max(0.0, wxx / w - m * m);
        return std::sqrt(var * w2) / w;
    }
};

// Summary of how the jets of one pT bin spread over trigger-bit patterns
static void PatternSummary(const std::map<int, double> &pat, const char *label)
{
    double tot = 0.0;
    std::vector<double> v;
    for (std::map<int, double>::const_iterator it = pat.begin(); it != pat.end(); ++it) { tot += it->second; v.push_back(it->second); }
    std::sort(v.begin(), v.end());
    std::reverse(v.begin(), v.end());
    double sum_p2 = 0.0, run = 0.0; int n90 = 0;
    for (size_t i = 0; i < v.size(); i++)
    {
        sum_p2 += (v[i] / tot) * (v[i] / tot);
        if (run < 0.9 * tot) n90++;
        run += v[i];
    }
    printf("    %-22s %3d patterns seen, %2d hold 90%% of jets, most common %.0f%%; a random jet shares another jet's exact pattern %.0f%% of the time\n",
           label, (int)v.size(), n90, 100.0 * v[0] / tot, 100.0 * sum_p2);
}

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2)
    {
        std::cerr << "usage: ./proxy_trigger <input_dir/>" << std::endl;
        return 1;
    }
    std::string input_dir = argv[1];
    if (input_dir[input_dir.size() - 1] != '/') input_dir += '/';

    std::vector<std::string> jet_files, simple_files;
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
        for (size_t i = 0; i < names.size(); i++)
        {
            jet_files.push_back(input_dir + names[i]);
            simple_files.push_back(input_dir + "SimpleTree" + names[i].substr(10));
        }
    }
    std::cout << "Found " << jet_files.size() << " file pairs in " << input_dir << std::endl;

    std::vector<float> file_weight;
    if (!ComputeXSecWeights(jet_files, "jetTree", file_weight))
        return 1;

    // The whole ECal face as cells, to measure the area available away from a jet
    std::vector<Cell> cells;
    {
        const FiducialRect &rect = kFiducialRectEcal;
        for (int side = -1; side <= 1; side += 2)
            for (double x = rect.x_inner; x < rect.x_outer; x += kCellCm)
                for (double y = rect.y_min; y < rect.y_max; y += kCellCm)
                {
                    const double dx = std::min(kCellCm, rect.x_outer - x), dy = std::min(kCellCm, rect.y_max - y);
                    const double cx = side * (x + dx / 2), cy = y + dy / 2, r2 = cx * cx + cy * cy;
                    Cell c;
                    computeEtaPhi(cx, cy, c.eta, c.phi);
                    c.area = dx * dy * z_proj / (r2 * std::sqrt(r2 + z_proj * z_proj));
                    cells.push_back(c);
                }
    }

    double jets_w[kNPt] = {0};
    Mean ue_all[kNPt], ue_none[kNPt], ue_fired[kNPt][kNTrig], ue_not[kNPt][kNTrig], jet_y[kNPt][kNTrig], jet_ax[kNPt][kNTrig];
    double fired_w[kNPt][kNTrig] = {{0}};
    std::map<int, double> pat_jp[kNPt], pat_all[kNPt];

    const int MAX_PTS = 200; // matches StSimpleReaderMaker.h's Cal_pt_* array size
    Long64_t n_events = 0;

    for (size_t f = 0; f < jet_files.size(); f++)
    {
        if (file_weight[f] < 0) continue;
        const double w = file_weight[f];

        TFile *fj = TFile::Open(jet_files[f].c_str());
        TFile *fs = TFile::Open(simple_files[f].c_str());
        TTree *jt = (fj && !fj->IsZombie()) ? (TTree *)fj->Get("jetTree") : 0;
        TTree *st = (fs && !fs->IsZombie()) ? (TTree *)fs->Get("data") : 0;
        if (!jt || !st || jt->GetEntries() != st->GetEntries() || !st->GetBranch("Cal_npts"))
        {
            std::cerr << "WARNING: skipping " << jet_files[f] << std::endl;
            if (fj) { fj->Close(); delete fj; }
            if (fs) { fs->Close(); delete fs; }
            continue;
        }

        std::vector<float> *reco_pt = 0, *reco_eta = 0, *reco_phi = 0, *reco_x = 0, *reco_y = 0;
        Int_t trig[kNTrig] = {0};
        jt->SetBranchStatus("*", 0);
        const char *jb[] = {"reco_pt", "reco_eta", "reco_phi", "reco_x", "reco_y"};
        for (int i = 0; i < 5; i++) jt->SetBranchStatus(jb[i], 1);
        jt->SetBranchAddress("reco_pt", &reco_pt);
        jt->SetBranchAddress("reco_eta", &reco_eta);
        jt->SetBranchAddress("reco_phi", &reco_phi);
        jt->SetBranchAddress("reco_x", &reco_x);
        jt->SetBranchAddress("reco_y", &reco_y);
        for (int i = 0; i < kNTrig; i++) { jt->SetBranchStatus(kTrig[i], 1); jt->SetBranchAddress(kTrig[i], &trig[i]); }

        Int_t Cal_npts = 0;
        Float_t Cal_pt_energy[MAX_PTS], Cal_pt_x[MAX_PTS], Cal_pt_y[MAX_PTS], Cal_pt_z[MAX_PTS];
        Int_t Cal_pt_detid[MAX_PTS];
        st->SetBranchStatus("*", 0);
        const char *pt_branches[] = {"Cal_npts", "Cal_pt_energy", "Cal_pt_x", "Cal_pt_y", "Cal_pt_z", "Cal_pt_detid"};
        for (int i = 0; i < 6; i++) st->SetBranchStatus(pt_branches[i], 1);
        st->SetBranchAddress("Cal_npts", &Cal_npts);
        st->SetBranchAddress("Cal_pt_energy", Cal_pt_energy);
        st->SetBranchAddress("Cal_pt_x", Cal_pt_x);
        st->SetBranchAddress("Cal_pt_y", Cal_pt_y);
        st->SetBranchAddress("Cal_pt_z", Cal_pt_z);
        st->SetBranchAddress("Cal_pt_detid", Cal_pt_detid);

        const Long64_t n = jt->GetEntries();
        for (Long64_t ev = 0; ev < n; ev++)
        {
            if (jt->GetEntry(ev) <= 0 || st->GetEntry(ev) <= 0) continue;
            n_events++;
            if (reco_pt->size() != 1) continue; // single-jet events: the jet is the only candidate for what fired

            const double pt = reco_pt->at(0), eta = reco_eta->at(0), phi = reco_phi->at(0);
            int pb = -1;
            for (int k = 0; k < kNPt; k++) if (pt >= kPt[k] && pt < kPt[k + 1]) pb = k;
            if (pb < 0) continue;

            int bits_jp = 0, bits_all = 0;
            for (int i = 0; i < kNTrig; i++)
                if (trig[i]) { bits_all |= 1 << i; if (i < kNJP) bits_jp |= 1 << i; }
            jets_w[pb] += w;
            pat_jp[pb][bits_jp] += w;
            pat_all[pb][bits_all] += w;
            for (int i = 0; i < kNTrig; i++)
                if (trig[i]) { fired_w[pb][i] += w; jet_y[pb][i].Fill(reco_y->at(0), w); jet_ax[pb][i].Fill(std::fabs(reco_x->at(0)), w); }

            double area = 0.0;
            for (size_t c = 0; c < cells.size(); c++)
                if (ProxyJetAllowed(eta, phi, kBandHalfWidthDeg, cells[c].eta, cells[c].phi))
                    area += cells[c].area;
            if (area < kMinArea) continue;

            double e_sum = 0.0;
            for (int i = 0; i < Cal_npts; i++)
            {
                if (Cal_pt_detid[i] != 0 && Cal_pt_detid[i] != 1) continue;
                if (!pass_point_energy_cut(Cal_pt_energy[i], kSqrtS_MainRun22_GeV)) continue;
                const double rho = std::sqrt((double)Cal_pt_x[i] * Cal_pt_x[i] + (double)Cal_pt_y[i] * Cal_pt_y[i]);
                const double p_eta = std::asinh(Cal_pt_z[i] / rho), p_phi = std::atan2(Cal_pt_y[i], Cal_pt_x[i]);
                if (ProxyJetAllowed(eta, phi, kBandHalfWidthDeg, p_eta, p_phi))
                    e_sum += Cal_pt_energy[i];
            }
            const double ue = e_sum * M_PI * R * R / area; // scaled to one cone
            ue_all[pb].Fill(ue, w);
            if (!bits_all) ue_none[pb].Fill(ue, w);
            for (int i = 0; i < kNTrig; i++) (trig[i] ? ue_fired[pb][i] : ue_not[pb][i]).Fill(ue, w);
        }
        fj->Close(); delete fj;
        fs->Close(); delete fs;
    }
    std::cout << "Read " << n_events << " events; single-reco-jet events only; simulated trigger bits, no prescales" << std::endl;

    for (int k = 0; k < kNPt; k++)
    {
        printf("\nJet pT %.0f-%.0f GeV: UE per cone, all jets %.2f +- %.2f GeV; no trigger fired %.2f +- %.2f GeV\n",
               kPt[k], kPt[k + 1], ue_all[k].Value(), ue_all[k].Error(), ue_none[k].Value(), ue_none[k].Error());
        printf("    %-9s %8s  %-22s %-22s %10s %10s\n", "trigger", "fired", "UE if fired (GeV)", "UE if not fired (GeV)", "mean y", "mean |x|");
        for (int i = 0; i < kNTrig; i++)
            printf("    %-9s %7.1f%%  %6.2f +- %-12.2f %6.2f +- %-12.2f %8.1f cm %7.1f cm\n", kTrig[i], 100.0 * fired_w[k][i] / jets_w[k],
                   ue_fired[k][i].Value(), ue_fired[k][i].Error(), ue_not[k][i].Value(), ue_not[k][i].Error(),
                   jet_y[k][i].Value(), jet_ax[k][i].Value());
        PatternSummary(pat_jp[k], "7 jet-patch bits:");
        PatternSummary(pat_all[k], "all 11 bits:");
    }
    return 0;
}
