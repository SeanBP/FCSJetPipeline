// Generator-level recoil test for the in-event UE cones: Pythia6 with the
// simulation's "FMS tune" (PyTune 370, PARP(90) = 0.213, PARP(91) = 1,
// MSEL = 1, sqrt(s) = 508.4 GeV), optionally without multi-parton
// interactions (MSTP(81) = 20), no detector.
//
// Truth jets: anti-kT R on stable EM particles (photons, electrons; E > 0.2
// GeV, 1.5 < eta < 5.5, pi0/eta decayed, ctau > 10 mm stable), as the
// in-event truth photons. Selected jets: pT >= 2 GeV, 2.5 < eta < 4.0,
// E bins 10-20-35-60-120 GeV. Per jet:
//   h_band_b       EM energy at |eta - eta_jet| < R against |dphi| (sum over jets)
//   p_cone_b       mean EM energy in cones at the jet's eta rotated by
//                  bin 1: 180, 2: +135, 3: -135, 4: +90, 5: -90
//   h_rec_deta_b   eta(recoil parton) - eta(jet), for jets within R of one of
//                  the two outgoing hard-process partons (MSTI(7), MSTI(8));
//                  the recoil parton is the other one
//   h_rec_dphi_b   |dphi| of the recoil parton
//   h_rec_where_b  1: matched jets, 2: recoil within the eta band |deta| < R,
//                  3: recoil within R of the 180 cone, 4: of a +-90 cone
//   h_njet         jets per bin; info: sigma (mb), events
//
// Build and run inside the SL7 container (stardev):
//   source build.csh   (Pythia 6.4.28 from the STAR source without its dummy
//   PDFSET/STRUCTM, so LHAPDF supplies them, as in root4star)
//   ./recoil_truth <ckin3> <ckin4> <events> <seed> <mpi 0|1> <out.root>

#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "fastjet/ClusterSequence.hh"

extern "C" {
struct PyJets { int n, npad, k[5][4000]; double p[5][4000], v[5][4000]; }; extern PyJets pyjets_;
struct PyPars { int mstp[200]; double parp[200]; int msti[200]; double pari[200]; }; extern PyPars pypars_;
struct PySubs { int msel, mselpd, msub[500], kfin[81][2]; double ckin[200]; }; extern PySubs pysubs_;
struct PyDat1 { int mstu[200]; double paru[200]; int mstj[200]; double parj[200]; }; extern PyDat1 pydat1_;
struct PyDatr { int mrpy[6]; double rrpy[100]; }; extern PyDatr pydatr_;
void pyinit_(const char *, const char *, const char *, double *, long, long, long);
void pyevnt_();
void pytune_(int *);
void pystat_(int *);
}

static const double R = 0.7;
static const int kNB = 4;
static const double kEdge[kNB + 1] = {10, 20, 35, 60, 120};
static const int kNC = 5;
static const double kRot[kNC] = {180, 135, -135, 90, -90};

static double DPhi(double a, double b) { return std::fabs(std::remainder(a - b, 2 * M_PI)); }
static double DR(double e1, double p1, double e2, double p2) { const double dp = DPhi(p1, p2), de = e1 - e2; return std::sqrt(de * de + dp * dp); }

int main(int argc, char **argv)
{
    if (argc < 7) { std::fprintf(stderr, "usage: ./recoil_truth <ckin3> <ckin4> <events> <seed> <mpi 0|1> <out.root>\n"); return 1; }
    const double ckin3 = std::atof(argv[1]), ckin4 = std::atof(argv[2]);
    const long nev = std::atol(argv[3]);
    const int seed = std::atoi(argv[4]), mpi = std::atoi(argv[5]);

    pydatr_.mrpy[0] = seed;
    int tune = 370; pytune_(&tune);
    pypars_.parp[89] = 0.213;
    pypars_.parp[90] = 1.0;
    pysubs_.msel = 1;
    pysubs_.ckin[2] = ckin3; pysubs_.ckin[3] = ckin4;
    pydat1_.mstj[21] = 2; pydat1_.parj[70] = 10.0;  // decay only ctau < 10 mm
    if (!mpi) pypars_.mstp[80] = 20;
    pydat1_.mstu[11] = 12345; pydat1_.mstu[121] = 0;  // no banner / init listing
    double ecm = 508.4;
    pyinit_("CMS", "p", "p", &ecm, 3, 1, 1);

    TFile out(argv[6], "RECREATE");
    TH1D h_njet("h_njet", "jets;E bin", kNB, 0, kNB);
    std::vector<TH1D *> h_band, h_rdeta, h_rdphi, h_rwhere;
    std::vector<TProfile *> p_cone;
    for (int b = 0; b < kNB; b++)
    {
        h_band.push_back(new TH1D(Form("h_band_%d", b), ";|#Delta#phi|;EM energy (GeV)", 36, 0, M_PI));
        p_cone.push_back(new TProfile(Form("p_cone_%d", b), ";cone;EM energy (GeV)", kNC, 0, kNC));
        h_rdeta.push_back(new TH1D(Form("h_rec_deta_%d", b), ";#eta_{recoil} - #eta_{jet}", 120, -10, 2));
        h_rdphi.push_back(new TH1D(Form("h_rec_dphi_%d", b), ";|#Delta#phi_{recoil}|", 36, 0, M_PI));
        h_rwhere.push_back(new TH1D(Form("h_rec_where_%d", b), "", 4, 0.5, 4.5));
    }

    std::vector<double> pe, pp, pE;
    std::vector<fastjet::PseudoJet> in;
    const fastjet::JetDefinition jdef(fastjet::antikt_algorithm, R);
    for (long iev = 0; iev < nev; iev++)
    {
        pyevnt_();
        pe.clear(); pp.clear(); pE.clear(); in.clear();
        for (int i = 0; i < pyjets_.n; i++)
        {
            if (pyjets_.k[0][i] != 1) continue;
            const int kf = std::abs(pyjets_.k[1][i]);
            if (kf != 22 && kf != 11) continue;
            const double px = pyjets_.p[0][i], py = pyjets_.p[1][i], pz = pyjets_.p[2][i], E = pyjets_.p[3][i];
            const double pt = std::sqrt(px * px + py * py);
            if (E < 0.2 || !(pt > 0)) continue;
            const double eta = std::asinh(pz / pt);
            if (eta < 1.5 || eta > 5.5) continue;
            pe.push_back(eta); pp.push_back(std::atan2(py, px)); pE.push_back(E);
            in.push_back(fastjet::PseudoJet(px, py, pz, E));
        }
        if (in.empty()) continue;
        fastjet::ClusterSequence cs(in, jdef);
        const std::vector<fastjet::PseudoJet> jets = cs.inclusive_jets(2.0);
        if (jets.empty()) continue;
        // outgoing hard-process partons
        const int i7 = pypars_.msti[6] - 1, i8 = pypars_.msti[7] - 1;
        double reta[2] = {0, 0}, rphi[2] = {0, 0};
        bool rok[2] = {false, false};
        const int ih[2] = {i7, i8};
        for (int q = 0; q < 2; q++)
        {
            if (ih[q] < 0 || ih[q] >= pyjets_.n) continue;
            const double px = pyjets_.p[0][ih[q]], py = pyjets_.p[1][ih[q]], pz = pyjets_.p[2][ih[q]];
            const double pt = std::sqrt(px * px + py * py);
            if (!(pt > 0)) continue;
            reta[q] = std::asinh(pz / pt); rphi[q] = std::atan2(py, px); rok[q] = true;
        }
        for (size_t j = 0; j < jets.size(); j++)
        {
            const double jeta = jets[j].eta(), jphi = jets[j].phi_std(), jE = jets[j].E();
            if (jeta < 2.5 || jeta > 4.0) continue;
            int b = -1;
            for (int q = 0; q < kNB; q++) if (jE >= kEdge[q] && jE < kEdge[q + 1]) b = q;
            if (b < 0) continue;
            h_njet.Fill(b);
            double cone[kNC] = {0};
            for (size_t i = 0; i < pE.size(); i++)
            {
                if (std::fabs(pe[i] - jeta) < R) h_band[b]->Fill(DPhi(pp[i], jphi), pE[i]);
                for (int k = 0; k < kNC; k++)
                    if (DR(pe[i], pp[i], jeta, jphi + kRot[k] * M_PI / 180.0) < R) cone[k] += pE[i];
            }
            for (int k = 0; k < kNC; k++) p_cone[b]->Fill(k + 0.5, cone[k]);
            // recoil parton: the outgoing parton the jet is not matched to
            double d0 = rok[0] ? DR(reta[0], rphi[0], jeta, jphi) : 99, d1 = rok[1] ? DR(reta[1], rphi[1], jeta, jphi) : 99;
            if (std::min(d0, d1) >= R) continue;
            const int r = d0 < d1 ? 1 : 0;
            if (!rok[r]) continue;
            h_rwhere[b]->Fill(1);
            h_rdeta[b]->Fill(reta[r] - jeta);
            h_rdphi[b]->Fill(DPhi(rphi[r], jphi));
            if (std::fabs(reta[r] - jeta) < R) h_rwhere[b]->Fill(2);
            if (DR(reta[r], rphi[r], jeta, jphi + M_PI) < R) h_rwhere[b]->Fill(3);
            if (DR(reta[r], rphi[r], jeta, jphi + M_PI / 2) < R || DR(reta[r], rphi[r], jeta, jphi - M_PI / 2) < R) h_rwhere[b]->Fill(4);
        }
    }
    double sigma = pypars_.pari[0];
    double nevd = nev;
    TTree info("info", "run info");
    info.Branch("sigma", &sigma, "sigma/D"); info.Branch("events", &nevd, "events/D");
    info.Fill();
    out.Write();
    out.Close();
    std::printf("recoil_truth: %ld events, sigma = %.4g mb, jets per bin %.0f %.0f %.0f %.0f\n", nev, sigma,
                h_njet.GetBinContent(1), h_njet.GetBinContent(2), h_njet.GetBinContent(3), h_njet.GetBinContent(4));
    return 0;
}
