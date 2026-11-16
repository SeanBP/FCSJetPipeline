// Per-run helpers shared by AnAsymmetry, AnPlots and AnMerge (2026-11-16, the
// merged report flow), so the run selection and the per-run rows cannot drift
// apart between them:
//   RunOf          the run of a yield file (its "runs" entry)
//   ReadPolTable   make_run_polarization.py's table: run -> (P_blue, P_yellow)
//   ReadBadRuns    a "Runnumber | Flag | Notes" bad-run list (D. Kapukchyan's,
//                  or pi0_peak_runs.C's <out>_bad_runs.txt)
//   RunRow, AddRunRow, Write/ReadRunRows  AnPlots' per-run rows (runs.png,
//                  runs.txt), additive over the files of a run
#ifndef AN_RUN_TOOLS_H
#define AN_RUN_TOOLS_H

#include <TDirectory.h>
#include <TH1.h>
#include <TNamed.h>
#include <TString.h>
#include <TTree.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>

#include "AnConfig.h"

// D. Kapukchyan's pi0 QA bad-run list (AnAsymmetry's default)
static const char *const kBadRunListDefault = "/star/u/dkap7827/Fcs2019/FcsAsymRun22/BadRuns.txt";

// the run of an AnHistograms file (its "runs" entry), 0 if it holds none or several
inline int RunOf(TDirectory &f)
{
    TNamed *r = (TNamed *)f.Get("runs");
    if (!r) return 0;
    std::istringstream is(r->GetTitle());
    int run = 0, extra = 0;
    if (!(is >> run) || (is >> extra)) return 0;
    return run;
}

// run -> (P_blue, P_yellow) from "run fill t_mid P_b dP_b P_y dP_y" lines; false if unreadable
inline bool ReadPolTable(const std::string &path, std::map<int, std::pair<double, double> > &table)
{
    std::ifstream t(path.c_str());
    if (!t) return false;
    std::string line;
    while (std::getline(t, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        int run, fill; double tmid, pb, dpb, py, dpy;
        if (is >> run >> fill >> tmid >> pb >> dpb >> py >> dpy) table[run] = std::make_pair(pb, py);
    }
    return true;
}

// Bad-run list: run -> flag ("" when the flag column is empty); false if unreadable
inline bool ReadBadRuns(const std::string &path, std::map<int, std::string> &bad)
{
    std::ifstream in(path.c_str());
    if (!in) return false;
    std::string line;
    while (std::getline(in, line))
    {
        std::istringstream is(line);
        int run = 0;
        if (!(is >> run) || run <= 0) continue;  // header, blank, comment
        std::string rest, flag;
        std::getline(is, rest);
        const size_t a = rest.find('|');
        if (a != std::string::npos)
        {
            const size_t b = rest.find('|', a + 1);
            flag = rest.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1);
            const size_t f0 = flag.find_first_not_of(" \t"), f1 = flag.find_last_not_of(" \t");
            flag = f0 == std::string::npos ? "" : flag.substr(f0, f1 - f0 + 1);
        }
        bad[run] = flag;
    }
    return true;
}

// AnPlots' per-run row; every field is a sum over the run's files
struct RunRow
{
    double ev_read, ev_pass, jets, xpos, xneg, spin[kNSpinStates], lvl[kNTrgLevels], esum, en;
    RunRow() : ev_read(0), ev_pass(0), jets(0), xpos(0), xneg(0), esum(0), en(0)
    {
        for (int p = 0; p < kNSpinStates; p++) spin[p] = 0;
        for (int l = 0; l < kNTrgLevels; l++) lvl[l] = 0;
    }
    void Add(const RunRow &o)
    {
        ev_read += o.ev_read; ev_pass += o.ev_pass; jets += o.jets; xpos += o.xpos; xneg += o.xneg; esum += o.esum; en += o.en;
        for (int p = 0; p < kNSpinStates; p++) spin[p] += o.spin[p];
        for (int l = 0; l < kNTrgLevels; l++) lvl[l] += o.lvl[l];
    }
};

// one yield file's contribution to its run's row: events read and passing,
// default-factor jets per side, spin state and level, corrected-E sum.
// get(name) returns the file's histogram of that name or 0; old_labels: the
// file has no spin_convention (old mapping: state p is physical 3 - p)
template <class Getter>
inline void AddRunRowT(Getter get, bool old_labels, RunRow &r)
{
    TH1 *cut = get("h_cutflow");
    if (cut) r.ev_read += cut->GetBinContent(cut->GetNbinsX());
    for (int l = 0; l < kNTrgLevels; l++)
    {
        TH1 *he = get(TString::Format("h_events_%s", kTrgLevelName[l]).Data());
        if (!he) continue;
        for (int p = 0; p < kNSpinStates; p++)
        {
            r.spin[old_labels ? kNSpinStates - 1 - p : p] += he->GetBinContent(p + 1);
            r.ev_pass += he->GetBinContent(p + 1);
        }
        for (int k = 0; k < kNPhotonClasses; k++)
            for (int s = 0; s < kNSides; s++)
                for (int p = 0; p < kNSpinStates; p++)
                {
                    double n = 0;
                    for (int g = 0; g < kNTrgSlices; g++)
                    {
                        if (kTrgSliceEdge[g] < CutFactor(l) - 1e-9) continue;
                        // JP levels: only jets not firing EM at the default factor (each jet once)
                        for (int tg = 0; tg < (IsJpLevel(l) ? kNEmTags : 1); tg++)
                        {
                            if (IsJpLevel(l) && tg > 0 && kTrgSliceEdge[tg - 1] >= kEmCutFactor - 1e-9) continue;
                            TH1 *hy = get(TString::Format("h_xF_%s_%s_%s_%s%s%s", kTrgLevelName[l], PhotonClassName(k), kSideName[s], kSpinStateName[p],
                                                          kTrgSliceSuffix[g], IsJpLevel(l) ? kEmTagSuffix[tg] : "").Data());
                            if (hy) n += hy->Integral();
                        }
                    }
                    r.jets += n; r.lvl[l] += n;
                    (s ? r.xpos : r.xneg) += n;
                }
    }
    TH1 *ec = get("qa_Ecorr_sel");
    if (ec) { r.esum += ec->GetMean() * ec->GetEntries(); r.en += ec->GetEntries(); }
}

// the same, reading the histograms from a file or directory
inline void AddRunRow(TDirectory &f, RunRow &r)
{
    AddRunRowT([&f](const std::string &n) { return (TH1 *)f.Get(n.c_str()); }, !f.Get("spin_convention"), r);
}

// the rows as TTree "runs" in the current directory
inline void WriteRunRows(const std::map<int, RunRow> &rows)
{
    TTree t("runs", "per-run rows of AnPlots");
    int run; RunRow r;
    t.Branch("run", &run, "run/I");
    t.Branch("ev_read", &r.ev_read, "ev_read/D"); t.Branch("ev_pass", &r.ev_pass, "ev_pass/D");
    t.Branch("jets", &r.jets, "jets/D"); t.Branch("xpos", &r.xpos, "xpos/D"); t.Branch("xneg", &r.xneg, "xneg/D");
    t.Branch("spin", r.spin, TString::Format("spin[%d]/D", kNSpinStates));
    t.Branch("lvl", r.lvl, TString::Format("lvl[%d]/D", kNTrgLevels));
    t.Branch("esum", &r.esum, "esum/D"); t.Branch("en", &r.en, "en/D");
    for (std::map<int, RunRow>::const_iterator it = rows.begin(); it != rows.end(); ++it) { run = it->first; r = it->second; t.Fill(); }
    t.Write();
}

// adds a TTree "runs" (WriteRunRows) to rows; false if there is none
inline bool ReadRunRows(TDirectory &f, std::map<int, RunRow> &rows)
{
    TTree *t = (TTree *)f.Get("runs");
    if (!t || !t->InheritsFrom(TTree::Class())) return false;
    int run; RunRow r;
    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("ev_read", &r.ev_read); t->SetBranchAddress("ev_pass", &r.ev_pass);
    t->SetBranchAddress("jets", &r.jets); t->SetBranchAddress("xpos", &r.xpos); t->SetBranchAddress("xneg", &r.xneg);
    t->SetBranchAddress("spin", r.spin); t->SetBranchAddress("lvl", r.lvl);
    t->SetBranchAddress("esum", &r.esum); t->SetBranchAddress("en", &r.en);
    for (Long64_t i = 0; i < t->GetEntries(); i++) { t->GetEntry(i); rows[run].Add(r); }
    t->ResetBranchAddresses();
    return true;
}

#endif
