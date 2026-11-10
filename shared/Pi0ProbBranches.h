#ifndef PI0_PROB_BRANCHES_H
#define PI0_PROB_BRANCHES_H

// Per-reco-jet pi0 probability branches (shared/Pi0Prob.h), shared by
// sim_to_jet/JetMatcher.cpp and data_to_jet/RecoJets.cpp:
//   reco_pi0_prob     P(>= 1 resolved pi0 | all point-pair masses); 0 for
//                     n < 2, -1 when not computed (status 2 or 3)
//   reco_pi0_logodds  log of the posterior odds (-99: odds 0, 99: prior pi_0 = 0)
//   reco_pi0_npairs   point pairs, n(n-1)/2
//   reco_pi0_best_r   largest pair likelihood ratio S/B (diagnostic) and its
//   reco_pi0_best_i/j two points, as indices among the jet's ECal point
//                     constituents in reco_const_* order
//   reco_pi0_lead_pair  P(the jet's leading point is in a resolved pi0 pair | masses)
//   reco_pi0_lead_prob  P(the leading point is from a pi0): lead_pair plus the
//                     simulated share of merged / partner-lost pi0 (q(n, E));
//                     defined for one-point jets too; -1 if not computed
//   reco_pi0_status   0 OK, 1 n < 2, 2 more than max_points points (not
//                     approximated), 3 not computed (no calibration, or a
//                     hit jet: the probability needs point jets)
// The calibration JSON: $PI0PROB_JSON if set, else ./pi0prob.json; without
// either every jet gets status 3. Data use the S peak refitted on data
// (Pi0ProbCalib::Load(..., data = true)).

#include "Pi0Prob.h"

#include <TTree.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

struct Pi0ProbBranches
{
    std::vector<float> prob, logodds, best_r, lead_pair, lead;
    std::vector<int> npairs, best_i, best_j, status;
    Pi0ProbCalib calib;

    void Book(TTree *t)
    {
        t->Branch("reco_pi0_prob", &prob);
        t->Branch("reco_pi0_logodds", &logodds);
        t->Branch("reco_pi0_npairs", &npairs);
        t->Branch("reco_pi0_best_r", &best_r);
        t->Branch("reco_pi0_best_i", &best_i);
        t->Branch("reco_pi0_best_j", &best_j);
        t->Branch("reco_pi0_status", &status);
        t->Branch("reco_pi0_lead_pair", &lead_pair);
        t->Branch("reco_pi0_lead_prob", &lead);
    }
    // loads the calibration; false (with a message) if there is none
    bool Load(bool data)
    {
        const char *env = std::getenv("PI0PROB_JSON");
        const std::string path = env && *env ? env : "pi0prob.json";
        std::string err;
        if (!calib.Load(path, data, &err))
        {
            std::cerr << "Warning: no pi0 probability calibration (" << err << "); reco_pi0_status = 3 for every jet" << std::endl;
            return false;
        }
        std::cout << "pi0 probability calibration " << path << " (" << calib.version << (data ? ", data S" : ", simulation S") << ")" << std::endl;
        return true;
    }
    void Clear()
    {
        prob.clear(); logodds.clear(); best_r.clear(); lead_pair.clear(); lead.clear();
        npairs.clear(); best_i.clear(); best_j.clear(); status.clear();
    }
    // one jet: its ECal point constituents' 4-vectors and raw energy; compute
    // = false for hit jets (status 3)
    void Push(const std::vector<double> &E, const std::vector<double> &px, const std::vector<double> &py,
              const std::vector<double> &pz, double ejet, bool compute)
    {
        Pi0ProbResult r;
        if (compute) r = Pi0Prob(calib, E, px, py, pz, ejet);
        else { r.status = kPi0ProbNoCalib; r.prob = r.logodds = r.lead_pair = r.lead = -1; r.npairs = (int)E.size() * ((int)E.size() - 1) / 2; }
        prob.push_back(r.prob); logodds.push_back(r.logodds); best_r.push_back(r.best_r);
        npairs.push_back(r.npairs); best_i.push_back(r.best_i); best_j.push_back(r.best_j); status.push_back(r.status);
        lead_pair.push_back(r.lead_pair); lead.push_back(r.lead);
    }
};

#endif
