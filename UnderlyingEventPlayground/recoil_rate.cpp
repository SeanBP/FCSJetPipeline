// How often does a truth jet appear back-to-back in azimuth with a reco jet
// inside the FCS, and how much more often than it would appear at that same
// place by chance?
//
// Input: a folder of sim_to_jet -E 1 JetTrees (jet_output_*.root). Every
// event is weighted by its cross-section weight (jet_scale/JetXSecWeights.h).
//
// Same event. For every reco jet J, every truth jet T in the same event
// other than the one matched to J is entered at |dphi(T, J)|. Dividing by the
// (weighted) number of reco jets gives the mean number of truth jets per reco
// jet in each |dphi| bin; the last bin, 175-180 deg, is the "recoil window",
// within 5 deg of back-to-back.
//
// Baseline. The chance of a truth jet being at that same place when there is
// no corresponding reco jet: for the reco jet's azimuth phi_J, the mean
// number of truth jets per event, at each azimuth, over all events that have
// NO reco jet within 5 deg of phi_J. Reading that off at phi_J +- dphi gives
// the baseline for the same |dphi| bins. It carries the detector acceptance
// (truth jets can only be where the FCS is) but no correlation with J.
//
// The ratio same / baseline in the recoil window is the enhancement: how much
// more likely a truth jet is there because a reco jet is opposite.
//
// That baseline ("generic") averages over every kind of event, and in most of
// them the truth jet found is the event's own main jet, whereas the same-event
// count leaves out the jet matched to J and so only sees additional jets. A
// second, like-for-like baseline ("reco jet elsewhere") is therefore also
// made: the mean number of truth jets, other than the matched one, at that
// place in events whose reco jet is NOT within 5 deg of phi_J. The ratio to
// it isolates the effect of the reco jet being opposite.
//
// Truth jets are the ones JetMatcher stores: pT >= truth_jet_pt_min_em, axis
// inside the reco fiducial region dilated by R in eta-phi. Results are also
// given for truth jets with pT >= 2 GeV and for bins of reco jet pT.
//
// Build and run from this folder:
//   g++ -m64 recoil_rate.cpp ../shared/VectorDict.cxx -o recoil_rate -I. -I../shared -I../jet_scale \
//       `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//   ./recoil_rate <input_dir/> [output_base]
//
// Writes <output_base>.pdf/.png/.root (default recoil_rate).

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "JetParameters.h"
#include "JetXSecWeights.h"

struct Jet { float phi, pt; int match; }; // match: index of the matched jet of the other kind, or -1
struct Event { float w; std::vector<Jet> reco, truth; };

static const int kNPhi = 360;      // 1 deg azimuth bins for the baseline
static const int kNDphi = 36;      // 5 deg |dphi| bins
static const double kNear = 5.0;   // "no reco jet within kNear deg" for the baseline, and the recoil window half-width

static int PhiBin(double phi_rad)
{
    double d = std::fmod(phi_rad * 180.0 / M_PI, 360.0);
    if (d < 0) d += 360.0;
    return std::min(kNPhi - 1, (int)d);
}
static double AbsDphiDeg(double a, double b)
{
    double d = std::fabs(a - b) * 180.0 / M_PI;
    d = std::fmod(d, 360.0);
    return d > 180.0 ? 360.0 - d : d;
}

// One selection: reco jet pT range and truth jet pT threshold
struct Selection
{
    const char *name;
    double reco_lo, reco_hi, truth_min;
    double sw, sw2;                       // weighted number of reco jets
    long n_reco;
    std::vector<double> same, same_w2, base, base2;
    std::vector<long> same_n;
    Selection(const char *n, double lo, double hi, double tm)
        : name(n), reco_lo(lo), reco_hi(hi), truth_min(tm), sw(0), sw2(0), n_reco(0),
          same(kNDphi, 0.0), same_w2(kNDphi, 0.0), base(kNDphi, 0.0), base2(kNDphi, 0.0), same_n(kNDphi, 0) {}
};

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    if (argc < 2)
    {
        std::cerr << "usage: ./recoil_rate <input_dir/> [output_base]" << std::endl;
        return 1;
    }
    std::string input_dir = argv[1];
    if (input_dir[input_dir.size() - 1] != '/') input_dir += '/';
    const std::string out_base = argc > 2 ? argv[2] : "recoil_rate";

    std::vector<std::string> jet_files;
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
        for (size_t i = 0; i < names.size(); i++) jet_files.push_back(input_dir + names[i]);
    }
    std::cout << "Found " << jet_files.size() << " JetTree files in " << input_dir << std::endl;

    std::vector<float> file_weight;
    if (!ComputeXSecWeights(jet_files, "jetTree", file_weight))
        return 1;

    // ---------------- read everything into memory ----------------
    std::vector<Event> events;
    for (size_t f = 0; f < jet_files.size(); f++)
    {
        if (file_weight[f] < 0) continue;
        TFile *fj = TFile::Open(jet_files[f].c_str());
        if (!fj || fj->IsZombie()) continue;
        TTree *jt = (TTree *)fj->Get("jetTree");
        if (!jt) { fj->Close(); delete fj; continue; }

        std::vector<float> *reco_pt = 0, *reco_phi = 0, *truth_pt = 0, *truth_phi = 0;
        std::vector<int> *match_truth_idx = 0, *match_reco_idx = 0;
        jt->SetBranchStatus("*", 0);
        const char *br[] = {"reco_pt", "reco_phi", "truth_pt", "truth_phi", "match_truth_idx", "match_reco_idx"};
        for (int i = 0; i < 6; i++) jt->SetBranchStatus(br[i], 1);
        jt->SetBranchAddress("reco_pt", &reco_pt);
        jt->SetBranchAddress("reco_phi", &reco_phi);
        jt->SetBranchAddress("truth_pt", &truth_pt);
        jt->SetBranchAddress("truth_phi", &truth_phi);
        jt->SetBranchAddress("match_truth_idx", &match_truth_idx);
        jt->SetBranchAddress("match_reco_idx", &match_reco_idx);

        const Long64_t n = jt->GetEntries();
        for (Long64_t ev = 0; ev < n; ev++)
        {
            if (jt->GetEntry(ev) <= 0) continue;
            Event e; e.w = file_weight[f];
            for (size_t j = 0; j < reco_phi->size(); j++) { Jet a = {reco_phi->at(j), reco_pt->at(j), -1}; e.reco.push_back(a); }
            for (size_t j = 0; j < truth_phi->size(); j++) { Jet a = {truth_phi->at(j), truth_pt->at(j), -1}; e.truth.push_back(a); }
            for (size_t m = 0; m < match_truth_idx->size(); m++)
            {
                const int t = match_truth_idx->at(m), r = match_reco_idx->at(m);
                if (r >= 0 && r < (int)e.reco.size() && t >= 0 && t < (int)e.truth.size()) { e.reco[r].match = t; e.truth[t].match = r; }
            }
            events.push_back(e);
        }
        fj->Close();
        delete fj;
    }
    std::cout << "Read " << events.size() << " events" << std::endl;

    std::vector<Selection> sels;
    sels.push_back(Selection("all reco jets, truth pT >= 1 GeV", 0, 1e9, 0.0));
    sels.push_back(Selection("all reco jets, truth pT >= 2 GeV", 0, 1e9, 2.0));
    sels.push_back(Selection("reco pT 2-3 GeV, truth pT >= 1 GeV", 2, 3, 0.0));
    sels.push_back(Selection("reco pT 3-5 GeV, truth pT >= 1 GeV", 3, 5, 0.0));
    sels.push_back(Selection("reco pT > 5 GeV, truth pT >= 1 GeV", 5, 1e9, 0.0));

    for (size_t s = 0; s < sels.size(); s++)
    {
        Selection &S = sels[s];

        // ---- baseline ingredients, for this truth-jet threshold ----
        // T[phi_t]: weighted truth jets per azimuth bin over all events
        // Wr[phi_r]: weighted events with a reco jet in azimuth bin phi_r (any pT)
        // H[phi_r][phi_t]: weighted truth jets at phi_t in events with a reco jet at phi_r
        std::vector<double> T(kNPhi, 0.0), Wr(kNPhi, 0.0);
        std::vector<std::vector<double> > H(kNPhi, std::vector<double>(kNPhi, 0.0));
        // Like-for-like baseline: Wr2[phi_r] and H2[phi_r][phi_t] use the reco jets of
        // this selection and only the truth jets not matched to that reco jet
        std::vector<double> Wr2(kNPhi, 0.0);
        std::vector<std::vector<double> > H2(kNPhi, std::vector<double>(kNPhi, 0.0));
        double Wtot = 0.0;
        for (size_t i = 0; i < events.size(); i++)
        {
            const Event &e = events[i];
            Wtot += e.w;
            for (size_t r = 0; r < e.reco.size(); r++)
            {
                if (e.reco[r].pt < S.reco_lo || e.reco[r].pt >= S.reco_hi) continue;
                const int br = PhiBin(e.reco[r].phi);
                Wr2[br] += e.w;
                for (size_t t = 0; t < e.truth.size(); t++)
                    if ((int)t != e.reco[r].match && e.truth[t].pt >= S.truth_min) H2[br][PhiBin(e.truth[t].phi)] += e.w;
            }
            for (size_t t = 0; t < e.truth.size(); t++)
                if (e.truth[t].pt >= S.truth_min) T[PhiBin(e.truth[t].phi)] += e.w;
            for (size_t r = 0; r < e.reco.size(); r++)
            {
                const int br = PhiBin(e.reco[r].phi);
                Wr[br] += e.w;
                for (size_t t = 0; t < e.truth.size(); t++)
                    if (e.truth[t].pt >= S.truth_min) H[br][PhiBin(e.truth[t].phi)] += e.w;
            }
        }

        // base_j[k]: for a reco jet in azimuth bin j, the mean number of truth
        // jets per event in |dphi| bin k, over events with no reco jet within kNear deg of j
        std::vector<std::vector<double> > base_j(kNPhi, std::vector<double>(kNDphi, 0.0));
        std::vector<std::vector<double> > base2_j(kNPhi, std::vector<double>(kNDphi, 0.0));
        double W2tot = 0.0;
        std::vector<double> T2(kNPhi, 0.0);
        for (int jr = 0; jr < kNPhi; jr++)
        {
            W2tot += Wr2[jr];
            for (int pt = 0; pt < kNPhi; pt++) T2[pt] += H2[jr][pt];
        }
        for (int j = 0; j < kNPhi; j++)
        {
            double w_excl = 0.0, w2_excl = 0.0;
            std::vector<double> t_excl(kNPhi, 0.0), t2_excl(kNPhi, 0.0);
            for (int d = -(int)kNear; d <= (int)kNear; d++)
            {
                const int jr = ((j + d) % kNPhi + kNPhi) % kNPhi;
                w_excl += Wr[jr]; w2_excl += Wr2[jr];
                for (int pt = 0; pt < kNPhi; pt++) { t_excl[pt] += H[jr][pt]; t2_excl[pt] += H2[jr][pt]; }
            }
            const double denom = Wtot - w_excl, denom2 = W2tot - w2_excl;
            for (int pt = 0; pt < kNPhi; pt++)
            {
                int d = std::abs(pt - j);
                if (d > 180) d = 360 - d;
                const int k = std::min(kNDphi - 1, d / 5);
                base_j[j][k] += (T[pt] - t_excl[pt]) / denom;
                base2_j[j][k] += (T2[pt] - t2_excl[pt]) / denom2;
            }
        }

        // ---- same event, and the baseline averaged over the same reco jets ----
        for (size_t i = 0; i < events.size(); i++)
        {
            const Event &e = events[i];
            for (size_t r = 0; r < e.reco.size(); r++)
            {
                if (e.reco[r].pt < S.reco_lo || e.reco[r].pt >= S.reco_hi) continue;
                S.sw += e.w; S.sw2 += (double)e.w * e.w; S.n_reco++;
                const int j = PhiBin(e.reco[r].phi);
                for (int k = 0; k < kNDphi; k++) { S.base[k] += e.w * base_j[j][k]; S.base2[k] += e.w * base2_j[j][k]; }
                for (size_t t = 0; t < e.truth.size(); t++)
                {
                    if ((int)t == e.reco[r].match || e.truth[t].pt < S.truth_min) continue;
                    const int k = std::min(kNDphi - 1, (int)(AbsDphiDeg(e.truth[t].phi, e.reco[r].phi) / 5.0));
                    S.same[k] += e.w; S.same_w2[k] += (double)e.w * e.w; S.same_n[k]++;
                }
            }
        }

        const int kb = kNDphi - 1; // 175-180 deg
        const double p_same = S.same[kb] / S.sw, e_same = std::sqrt(S.same_w2[kb]) / S.sw, p_base = S.base[kb] / S.sw;
        std::cout << "\n" << S.name << ": " << S.n_reco << " reco jets (" << S.sw * S.sw / S.sw2 << " effective)" << std::endl;
        std::cout << "  recoil window 175-180 deg: truth jets per reco jet, same event = " << 100.0 * p_same << " +- " << 100.0 * e_same
                  << " %  (" << S.same_n[kb] << " raw pairs)" << std::endl;
        const double p_base2 = S.base2[kb] / S.sw;
        std::cout << "  generic baseline (any event, no reco jet opposite)          = " << 100.0 * p_base << " %   ratio " << p_same / p_base << " +- " << e_same / p_base << std::endl;
        std::cout << "  like-for-like baseline (reco jet elsewhere, extra jets only) = " << 100.0 * p_base2 << " %   ratio " << p_same / p_base2 << " +- " << e_same / p_base2 << std::endl;
    }

    // ---------------- plot the first two selections ----------------
    gStyle->SetOptStat(0);
    TFile fout((out_base + ".root").c_str(), "RECREATE");
    TCanvas c("c", "", 1200, 1100);
    TPad top("top", "", 0, 0.36, 1, 1), bot("bot", "", 0, 0, 1, 0.36);
    top.SetBottomMargin(0.02); top.SetLeftMargin(0.13); top.SetRightMargin(0.04); top.SetLogy();
    bot.SetTopMargin(0.03); bot.SetBottomMargin(0.28); bot.SetLeftMargin(0.13); bot.SetRightMargin(0.04);
    top.Draw(); bot.Draw();

    const Selection &P = sels[0];
    TH1D *hs = new TH1D("same", "", kNDphi, 0, 180);
    TH1D *hb1 = new TH1D("base_generic", "", kNDphi, 0, 180), *hb2 = new TH1D("base_elsewhere", "", kNDphi, 0, 180);
    TH1D *hr1 = new TH1D("ratio_generic", "", kNDphi, 0, 180), *hr2 = new TH1D("ratio_elsewhere", "", kNDphi, 0, 180);
    double ymin = 1e9, ymax = 0;
    for (int k = 0; k < kNDphi; k++)
    {
        const double ps = P.same[k] / P.sw, es = std::sqrt(P.same_w2[k]) / P.sw, p1 = P.base[k] / P.sw, p2 = P.base2[k] / P.sw;
        hs->SetBinContent(k + 1, ps); hs->SetBinError(k + 1, es);
        hb1->SetBinContent(k + 1, p1); hb2->SetBinContent(k + 1, p2);
        if (p1 > 0) { hr1->SetBinContent(k + 1, ps / p1); hr1->SetBinError(k + 1, es / p1); }
        if (p2 > 0) { hr2->SetBinContent(k + 1, ps / p2); hr2->SetBinError(k + 1, es / p2); }
        if (ps > 0) ymin = std::min(ymin, ps);
        if (p2 > 0) ymin = std::min(ymin, p2);
        ymax = std::max(ymax, std::max(ps, p1));
    }
    const int c_same = kBlack, c_gen = kAzure + 1, c_else = kOrange + 7;
    hs->SetLineColor(c_same); hs->SetMarkerColor(c_same); hs->SetMarkerStyle(20); hs->SetLineWidth(2);
    hb1->SetLineColor(c_gen); hb1->SetLineStyle(2); hb1->SetLineWidth(3);
    hb2->SetLineColor(c_else); hb2->SetLineWidth(3);
    hr1->SetLineColor(c_gen); hr1->SetMarkerColor(c_gen); hr1->SetMarkerStyle(24); hr1->SetLineWidth(2);
    hr2->SetLineColor(c_else); hr2->SetMarkerColor(c_else); hr2->SetMarkerStyle(20); hr2->SetLineWidth(2);
    hs->Write(); hb1->Write(); hb2->Write(); hr1->Write(); hr2->Write();

    top.cd();
    hs->SetTitle(";;Truth jets per reco jet, per 5#circ bin");
    hs->GetYaxis()->SetRangeUser(0.5 * ymin, 12.0 * ymax);
    hs->GetYaxis()->SetTitleSize(0.050); hs->GetYaxis()->SetLabelSize(0.045); hs->GetYaxis()->SetTitleOffset(1.2);
    hs->GetXaxis()->SetLabelSize(0);
    hs->Draw("E1");
    hb1->Draw("HIST SAME"); hb2->Draw("HIST SAME"); hs->Draw("E1 SAME");
    TLegend leg(0.16, 0.68, 0.94, 0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.037);
    leg.AddEntry(hs, "Same event as the reco jet (jet matched to it excluded)", "lp");
    leg.AddEntry(hb2, "Baseline: events whose reco jet is elsewhere, extra jets only", "l");
    leg.AddEntry(hb1, "Baseline: any event with no reco jet there, all truth jets", "l");
    leg.Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.043); t.SetTextFont(62);
    t.DrawLatex(0.16, 0.92, "Truth jets (p_{T} #geq 1 GeV) around a reco jet in the FCS, simulation");

    bot.cd();
    hr2->SetTitle(";|#Delta#phi| between truth jet and reco jet (deg);Same / baseline");
    hr2->GetYaxis()->SetRangeUser(0, 2.0);
    hr2->GetYaxis()->SetTitleSize(0.085); hr2->GetYaxis()->SetLabelSize(0.075); hr2->GetYaxis()->SetTitleOffset(0.65);
    hr2->GetXaxis()->SetTitleSize(0.095); hr2->GetXaxis()->SetLabelSize(0.080);
    hr2->GetYaxis()->SetNdivisions(506);
    hr2->Draw("E1");
    TLine one(0, 1, 180, 1); one.SetLineStyle(2); one.SetLineColor(kGray + 2); one.Draw();
    hr1->Draw("E1 SAME"); hr2->Draw("E1 SAME");

    c.SaveAs((out_base + ".pdf").c_str());
    c.SaveAs((out_base + ".png").c_str());
    fout.Close();
    return 0;
}
