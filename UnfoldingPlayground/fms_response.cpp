// The response objects of the STAR FMS EM-jet analysis, made from our FCS
// simulation: BrightSTAR starSim/src/FmsSimMakeResponseMatrix.cxx and the
// particle-level correction of emJetAnalysis/src/EjCreateBinnedHistExtended.cxx
// (see UnderlyingEventPlayground/fms_emjet_reference.txt, sections 7 and 8).
//
// FMS made, from matched detector/particle jets of triggered simulated events:
//   - an energy response matrix (detector E along x, particle E along y; raw
//     energies, no UE subtraction) and its profile, the mean particle energy
//     per detector-energy bin, fitted with a 6th-order polynomial up to 60 GeV
//     and a straight line above: the 1D energy correction E -> E'
//   - the same for pT, with UE subtracted at both levels (pT - dpT), fitted
//     with a cubic up to 5.8 GeV and a line above: pT -> pT'
//   - photon-multiplicity migration matrices (detector vs particle photon
//     count), used to correct A_N per multiplicity (fms_multiplicity_unfold.cpp)
// The energy and pT corrections are one factor per jet, position-independent:
// no matrix unfolding of E, pT or x_F was done.
//
// Here, with our definitions: UE subtracted from energy AND pT at both levels
// (proxy-event method), cross-section weights (FMS filled unweighted), and our
// reco fiducial region. Analysis selection by default, as for data: an EM
// trigger fired and UE-subtracted reco pT >= 2 GeV (FMS applied 2 GeV after
// its correction, and a trigger-dependent threshold).
//
// Usage: ./fms_response <file list> [output_base] [pipeline|fms matching] [analysis|all]
//        [E_break_GeV=60]
// Build: g++ -m64 -O2 fms_response.cpp ../shared/VectorDict.cxx -o fms_response -I../shared -I../jet_scale
//            `root-config --cflags --libs` -std=c++11

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>

#include <algorithm>

#include "UnfoldingInput.h"

// photon multiplicity bins 1, 2, 3, 4, >= 5 (FMS dropped jets above 5)
const int n_nph = 5;
inline int NphBin(int n) { return std::min(std::max(n, 1), n_nph); }

const int n_ebin = 6;
const double e_edges[n_ebin + 1] = {0, 20, 40, 60, 80, 100, 150};   // FMS: 20 GeV wide, 0-100 GeV
const int n_ptbin = 9;
const double pt_edges[n_ptbin + 1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0};  // FMS pT bins

int main(int argc, char **argv)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(55);
    if (argc < 2)
    {
        std::cerr << "usage: ./fms_response <file list> [output_base] [pipeline|fms] [analysis|all] [E_break_GeV]" << std::endl;
        return 1;
    }
    const std::string out_base = argc > 2 ? argv[2] : "fms_response";
    const MatchMode mode = (argc > 3 && std::string(argv[3]) == "fms") ? kMatchFms : kMatchPipeline;
    const bool analysis = !(argc > 4 && std::string(argv[4]) == "all");
    const double e_break = argc > 5 ? std::atof(argv[5]) : 60.0;

    std::vector<JetPair> pairs;
    if (!LoadJetPairs(ReadFileList(argv[1]), mode, pairs)) return 1;

    TFile fout((out_base + ".root").c_str(), "RECREATE");
    TH2D hE("h2EngResMat", ";E^{det}_{jet} - UE (GeV);E^{part}_{jet} - UE (GeV)", 60, 0, 150, 60, 0, 150);
    TH2D hPt("h2PtResMat", ";p_{T}^{det} - UE (GeV);p_{T}^{part} - UE (GeV)", 40, 2.0, 6.0, 40, 2.0, 6.0);
    TProfile pE("h2EngResMatProf", ";E^{det}_{jet} - UE (GeV);mean E^{part}_{jet} - UE (GeV)", 60, 0, 150);
    TProfile pPt("h2PtResMatProf", ";p_{T}^{det} - UE (GeV);mean p_{T}^{part} - UE (GeV)", 32, 2.0, 10.0);
    TH2D hNph("h2nPhResMat", ";detector photon multiplicity;particle photon multiplicity", n_nph, 0.5, n_nph + 0.5, n_nph, 0.5, n_nph + 0.5);
    std::vector<TH2D *> hNphBin;
    for (int ie = 0; ie < n_ebin; ie++)
        for (int ip = 0; ip < n_ptbin; ip++)
            hNphBin.push_back(new TH2D(Form("hNph2D_E%d_pt%d", ie, ip),
                                       Form("E %g-%g GeV, p_{T} %g-%g GeV;detector photon multiplicity;particle photon multiplicity",
                                            e_edges[ie], e_edges[ie + 1], pt_edges[ip], pt_edges[ip + 1]),
                                       n_nph, 0.5, n_nph + 0.5, n_nph, 0.5, n_nph + 0.5));
    hE.Sumw2(); hPt.Sumw2(); hNph.Sumw2();

    long n_used = 0, n_no_ue = 0;
    for (size_t i = 0; i < pairs.size(); i++)
    {
        const JetPair &p = pairs[i];
        if (!p.has_ue) { n_no_ue++; continue; }
        if (analysis && (p.em_level < 0 || p.rptsub() < 2.0)) continue;
        n_used++;
        hE.Fill(p.rEsub(), p.tEsub(), p.w);
        pE.Fill(p.rEsub(), p.tEsub(), p.w);
        hPt.Fill(p.rptsub(), p.tptsub(), p.w);
        pPt.Fill(p.rptsub(), p.tptsub(), p.w);
        hNph.Fill(NphBin(p.rnph), NphBin(p.tnph), p.w);
        const int ie = (int)(std::upper_bound(e_edges, e_edges + n_ebin + 1, p.rEsub()) - e_edges) - 1;
        const int ip = (int)(std::upper_bound(pt_edges, pt_edges + n_ptbin + 1, p.rptsub()) - pt_edges) - 1;
        if (ie >= 0 && ie < n_ebin && ip >= 0 && ip < n_ptbin) hNphBin[ie * n_ptbin + ip]->Fill(NphBin(p.rnph), NphBin(p.tnph), p.w);
    }
    std::cout << pairs.size() << " matched fiducial pairs (" << (mode == kMatchFms ? "FMS" : "pipeline") << " matching), " << n_no_ue
              << " without UE estimates; " << n_used << " used (" << (analysis ? "EM-triggered, UE-subtracted reco pT >= 2 GeV" : "all")
              << ")" << std::endl;

    // FMS-form 1D corrections, fitted to the profiles
    const double e_lo = 5.0, e_hi = 150.0;
    TF1 fE_low("fE_low", "pol6", e_lo, e_break), fE_high("fE_high", "pol1", e_break, e_hi);
    pE.Fit(&fE_low, "QR0"); pE.Fit(&fE_high, "QR0+");
    TF1 fPt_low("fPt_low", "pol3", 2.0, 5.8), fPt_high("fPt_high", "pol1", 5.8, 10.0);
    pPt.Fit(&fPt_low, "QR0"); pPt.Fit(&fPt_high, "QR0+");
    std::cout << "\nFMS-form energy correction (UE-subtracted E, cross-section weighted):\n  E <= " << e_break << " GeV: E' =";
    for (int k = 0; k <= 6; k++) std::cout << " " << (fE_low.GetParameter(k) >= 0 && k ? "+" : "") << fE_low.GetParameter(k) << (k ? Form(" E^%d", k) : "");
    std::cout << "   (chi2/ndf " << fE_low.GetChisquare() << "/" << fE_low.GetNDF() << ")" << std::endl;
    std::cout << "  E >  " << e_break << " GeV: E' = " << fE_high.GetParameter(0) << " + " << fE_high.GetParameter(1) << " E   (chi2/ndf "
              << fE_high.GetChisquare() << "/" << fE_high.GetNDF() << ")" << std::endl;
    std::cout << "  continuity at the break: " << fE_low.Eval(e_break) << " vs " << fE_high.Eval(e_break) << " GeV" << std::endl;
    std::cout << "FMS-form pT correction:\n  pT <= 5.8 GeV: pT' =";
    for (int k = 0; k <= 3; k++) std::cout << " " << fPt_low.GetParameter(k) << (k ? Form(" pT^%d", k) : "");
    std::cout << "   (chi2/ndf " << fPt_low.GetChisquare() << "/" << fPt_low.GetNDF() << ")" << std::endl;
    std::cout << "  pT >  5.8 GeV: pT' = " << fPt_high.GetParameter(0) << " + " << fPt_high.GetParameter(1) << " pT" << std::endl;
    std::cout << "  (FMS, 200 GeV: E' = -8.482 + 4.069E - 0.3063E^2 + ... up to 60 GeV, 8.248 + 0.9365E above;"
                 " pT' = 1.486 + 0.07337pT + 0.215pT^2 - 0.01525pT^3 up to 5.8 GeV, 0.4478 + 0.9685pT above)" << std::endl;
    std::cout << "\nCorrection factor E'/E at E = 10, 20, 40, 60, 100 GeV:";
    const double es[5] = {10, 20, 40, 60, 100};
    for (int k = 0; k < 5; k++) std::cout << " " << (es[k] <= e_break ? fE_low.Eval(es[k]) : fE_high.Eval(es[k])) / es[k];
    std::cout << "\nCorrection factor pT'/pT at pT = 2, 3, 4, 6, 8 GeV:";
    const double ps[5] = {2, 3, 4, 6, 8};
    for (int k = 0; k < 5; k++) std::cout << " " << (ps[k] <= 5.8 ? fPt_low.Eval(ps[k]) : fPt_high.Eval(ps[k])) / ps[k];
    std::cout << std::endl;

    // Photon multiplicity: purity (fraction of a detector bin that is the
    // same particle bin) and stability (fraction of a particle bin that stays)
    std::cout << "\nPhoton multiplicity, all bins together (1, 2, 3, 4, >=5):\n  purity   ";
    for (int b = 1; b <= n_nph; b++)
    {
        double col = 0;
        for (int t = 1; t <= n_nph; t++) col += hNph.GetBinContent(b, t);
        std::cout << " " << (col > 0 ? hNph.GetBinContent(b, b) / col : 0);
    }
    std::cout << "\n  stability";
    for (int t = 1; t <= n_nph; t++)
    {
        double row = 0;
        for (int b = 1; b <= n_nph; b++) row += hNph.GetBinContent(b, t);
        std::cout << " " << (row > 0 ? hNph.GetBinContent(t, t) / row : 0);
    }
    std::cout << std::endl;

    // ---- plots ----
    TCanvas c("c", "", 1500, 1000);
    c.Divide(3, 2, 0.001, 0.001);
    std::vector<TObject *> keep;
    TLatex tl;
    tl.SetNDC(); tl.SetTextSize(0.05); tl.SetTextFont(62);
    c.cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.13); hE.Draw("colz"); tl.DrawLatex(0.1, 0.93, "energy response");
    c.cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.13); hPt.Draw("colz"); tl.DrawLatex(0.1, 0.93, "p_{T} response");
    c.cd(3); gPad->SetLogz(); gPad->SetRightMargin(0.13); hNph.Draw("colz text"); tl.DrawLatex(0.1, 0.93, "photon multiplicity, all bins");
    c.cd(4);
    pE.SetMarkerStyle(20); pE.SetMarkerSize(0.6); pE.Draw();
    fE_low.SetLineColor(kRed + 1); fE_low.Draw("same");
    fE_high.SetLineColor(kAzure + 2); fE_high.Draw("same");
    TLine *d1 = new TLine(0, 0, 150, 150); d1->SetLineStyle(2); d1->Draw(); keep.push_back(d1);
    tl.DrawLatex(0.1, 0.93, "mean particle vs detector energy");
    c.cd(5);
    pPt.SetMarkerStyle(20); pPt.SetMarkerSize(0.6); pPt.Draw();
    fPt_low.SetLineColor(kRed + 1); fPt_low.Draw("same");
    fPt_high.SetLineColor(kAzure + 2); fPt_high.Draw("same");
    TLine *d2 = new TLine(2, 2, 10, 10); d2->SetLineStyle(2); d2->Draw(); keep.push_back(d2);
    tl.DrawLatex(0.1, 0.93, "mean particle vs detector p_{T}");
    c.cd(6); gPad->SetLogz(); gPad->SetRightMargin(0.13);
    TH2D *rep = hNphBin[3 * n_ptbin + 2];
    rep->Draw("colz text");
    tl.DrawLatex(0.1, 0.93, Form("multiplicity, E %g-%g, p_{T} %g-%g GeV", e_edges[3], e_edges[4], pt_edges[2], pt_edges[3]));
    c.SaveAs((out_base + ".png").c_str());
    c.SaveAs((out_base + ".pdf").c_str());

    fout.cd();
    hE.Write(); hPt.Write(); pE.Write(); pPt.Write(); hNph.Write();
    fE_low.Write(); fE_high.Write(); fPt_low.Write(); fPt_high.Write();
    for (size_t k = 0; k < hNphBin.size(); k++) hNphBin[k]->Write();
    fout.Close();
    return 0;
}
