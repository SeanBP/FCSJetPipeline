#ifndef UNFOLDING_INPUT_H
#define UNFOLDING_INPUT_H

// Matched reco/truth jet pairs for the unfolding playground, read from
// sim_to_jet JetTrees that carry the underlying-event branches of
// underlying_event/UnderlyingEvent.cpp (reco_ue_*, truth_ue_*), with each
// pair's cross-section weight (jet_scale/JetXSecWeights.h).
//
// Two ways to pair jets:
//   kMatchPipeline  the pipeline's own one-to-one matching (match_reco_idx /
//                   match_truth_idx, closest pairs first, dR < R/2)
//   kMatchFms       the STAR FMS EM-jet analysis's
//                   (starSim/src/FmsSimMakeResponseMatrix.cxx in BrightSTAR):
//                   for every reco jet, the nearest truth jet in eta-phi,
//                   accepted if dR < 0.2; not one-to-one (a truth jet can be
//                   used twice), and dphi is wrapped here (FMS did not wrap it)

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "JetParameters.h"
#include "JetXSecWeights.h"

enum MatchMode { kMatchPipeline, kMatchFms };

struct JetPair
{
    // reco: raw energy and pT, the UE estimates, position, photon multiplicity
    float rE, rpt, rueE, ruept, rx, ry, reta, rphi;
    int rnph;
    // truth
    float tE, tpt, tueE, tuept, teta, tphi;
    int tnph;          // truth photons, no energy cut (as FMS; reco needs 1 GeV)
    int em_level;      // highest EM trigger level that fired, -1 none (ProxyEmLevel convention)
    float w;           // cross-section weight
    int file;          // index in the file list (for splitting the sample)
    bool has_ue;       // both UE estimates present
    float rEsub() const { return rE - rueE; }
    float rptsub() const { return rpt - ruept; }
    float tEsub() const { return tE - tueE; }
    float tptsub() const { return tpt - tuept; }
};

inline std::vector<std::string> ReadFileList(const std::string &list)
{
    std::vector<std::string> files;
    std::ifstream in(list.c_str());
    std::string s;
    while (in >> s) files.push_back(s);
    return files;
}

// Reads every matched pair whose reco jet is inside the reco fiducial
// region. Returns false if the weights cannot be computed.
inline bool LoadJetPairs(const std::vector<std::string> &files, MatchMode mode, std::vector<JetPair> &pairs)
{
    std::vector<float> file_weight;
    if (!ComputeXSecWeights(files, "jetTree", file_weight)) return false;

    long n_no_ue_branch = 0;
    for (size_t fidx = 0; fidx < files.size(); fidx++)
    {
        if (file_weight[fidx] < 0) continue;
        TFile *f = TFile::Open(files[fidx].c_str());
        TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("jetTree") : 0;
        if (!t) { if (f) { f->Close(); delete f; } continue; }
        const bool ue = t->GetBranch("reco_ue_E") && t->GetBranch("truth_ue_E");
        if (!ue) n_no_ue_branch++;

        std::vector<float> *rE = 0, *rpt = 0, *rx = 0, *ry = 0, *reta = 0, *rphi = 0, *rueE = 0, *ruept = 0;
        std::vector<float> *tE = 0, *tpt = 0, *teta = 0, *tphi = 0, *tueE = 0, *tuept = 0;
        std::vector<int> *rnph = 0, *tnph = 0, *mr = 0, *mt = 0;
        Int_t em0 = 0, em1 = 0, em2 = 0;
        t->SetBranchStatus("*", 0);
        const char *used[] = {"reco_E", "reco_pt", "reco_x", "reco_y", "reco_eta", "reco_phi", "reco_nphotons", "truth_E", "truth_pt",
                              "truth_eta", "truth_phi", "truth_nphotons", "match_reco_idx", "match_truth_idx", "fcsEM0", "fcsEM1", "fcsEM2"};
        for (int k = 0; k < 17; k++) t->SetBranchStatus(used[k], 1);
        t->SetBranchAddress("reco_E", &rE);
        t->SetBranchAddress("reco_pt", &rpt);
        t->SetBranchAddress("reco_x", &rx);
        t->SetBranchAddress("reco_y", &ry);
        t->SetBranchAddress("reco_eta", &reta);
        t->SetBranchAddress("reco_phi", &rphi);
        t->SetBranchAddress("reco_nphotons", &rnph);
        t->SetBranchAddress("truth_E", &tE);
        t->SetBranchAddress("truth_pt", &tpt);
        t->SetBranchAddress("truth_eta", &teta);
        t->SetBranchAddress("truth_phi", &tphi);
        t->SetBranchAddress("truth_nphotons", &tnph);
        t->SetBranchAddress("match_reco_idx", &mr);
        t->SetBranchAddress("match_truth_idx", &mt);
        t->SetBranchAddress("fcsEM0", &em0);
        t->SetBranchAddress("fcsEM1", &em1);
        t->SetBranchAddress("fcsEM2", &em2);
        if (ue)
        {
            const char *ub[] = {"reco_ue_E", "reco_ue_pt", "truth_ue_E", "truth_ue_pt"};
            for (int k = 0; k < 4; k++) t->SetBranchStatus(ub[k], 1);
            t->SetBranchAddress("reco_ue_E", &rueE);
            t->SetBranchAddress("reco_ue_pt", &ruept);
            t->SetBranchAddress("truth_ue_E", &tueE);
            t->SetBranchAddress("truth_ue_pt", &tuept);
        }

        for (Long64_t ev = 0; ev < t->GetEntries(); ev++)
        {
            t->GetEntry(ev);
            std::vector<std::pair<int, int> > m;
            if (mode == kMatchPipeline)
            {
                for (size_t k = 0; k < mr->size(); k++) m.push_back(std::make_pair(mr->at(k), mt->at(k)));
            }
            else
            {
                for (size_t ri = 0; ri < rE->size(); ri++)
                {
                    int best = -1;
                    double best_dr = 1e9;
                    for (size_t ti = 0; ti < tE->size(); ti++)
                    {
                        double dphi = std::fabs(rphi->at(ri) - tphi->at(ti));
                        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
                        const double deta = reta->at(ri) - teta->at(ti);
                        const double dr = std::sqrt(deta * deta + dphi * dphi);
                        if (dr < best_dr) { best_dr = dr; best = (int)ti; }
                    }
                    if (best >= 0 && best_dr < 0.2) m.push_back(std::make_pair((int)ri, best));
                }
            }
            for (size_t k = 0; k < m.size(); k++)
            {
                const int ri = m[k].first, ti = m[k].second;
                if (!pass_fiducial_cut(rx->at(ri), ry->at(ri), reco_fiducial_buffer_ecal, kFiducialRectEcal)) continue;
                JetPair p;
                p.rE = rE->at(ri); p.rpt = rpt->at(ri); p.rx = rx->at(ri); p.ry = ry->at(ri);
                p.reta = reta->at(ri); p.rphi = rphi->at(ri); p.rnph = rnph->at(ri);
                p.tE = tE->at(ti); p.tpt = tpt->at(ti); p.teta = teta->at(ti); p.tphi = tphi->at(ti); p.tnph = tnph->at(ti);
                p.has_ue = ue && rueE->at(ri) >= 0 && tueE->at(ti) >= 0;
                p.rueE = p.has_ue ? rueE->at(ri) : 0.f; p.ruept = p.has_ue ? ruept->at(ri) : 0.f;
                p.tueE = p.has_ue ? tueE->at(ti) : 0.f; p.tuept = p.has_ue ? tuept->at(ti) : 0.f;
                p.em_level = em2 ? 2 : (em1 ? 1 : (em0 ? 0 : -1));
                p.w = file_weight[fidx];
                p.file = (int)fidx;
                pairs.push_back(p);
            }
        }
        f->Close();
        delete f;
    }
    if (n_no_ue_branch > 0)
        std::cout << "WARNING: " << n_no_ue_branch << " files have no reco/truth UE branches; their pairs have has_ue = false" << std::endl;
    return true;
}

#endif
