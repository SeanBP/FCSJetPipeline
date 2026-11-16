// Merged report flow (Sean 2026-11-16, A_N job efficiency). Each report step
// used to copy and re-read all ~9,000 per-run yield files, at ~1-2.4 s a file
// (~11,000 histograms each). Every step only sums them (A_N with each run's
// polarization as weight), so they are summed once:
//
//   AnMerge -b <yield list> <run polarization table> <out.root>
//       build (one per condor merge job, ~300 files): every histogram summed;
//       the yields (names h_xF*) also weighted by each run's P_blue / P_yellow
//       (directories Pblue / Pyellow: sum over runs of P_run x yield, which
//       is AnAsymmetry's per-bin P weighting); AnPlots' per-run rows (TTree
//       "runs", AnRunTools.h); per-run pi0 mass spectra (directory pi0_runs:
//       m_<run>, title = the pi0 scale factor used, for pi0_peak_runs.C); the
//       files used (TNamed "files": name:run:P_blue:P_yellow). Files whose run
//       is not one run, not in the polarization table, or has no P for a beam
//       are left out, as AnAsymmetry -r does (TNamed "skipped").
//   AnMerge -a <partial list> <out.root>
//       adds -b outputs (or -a outputs) together.
//   AnMerge -x <merged.root> <hists dir> <run polarization table> <out.root>
//          [-B <bad-run list> | -B none] [-F flagged | all] [-Q <pi0-peak bad-run list>]
//       masks: subtracts the yield files (read from <hists dir>) of every run
//       that AnAsymmetry -r would leave out: AnConfig.h kSuppBadRuns, D.
//       Kapukchyan's list (-B, default AnRunTools.h kBadRunListDefault; -F
//       flagged = runs with a flag, all = every listed run) and the pi0-peak
//       list (-Q). All are linear sums of raw counts, so the result equals a
//       sum over the kept runs. Writes TNamed "masked" and "bad_runs" (the
//       note AnAsymmetry prints). Per-run rows and pi0 spectra keep every run.
//
// Then: AnAsymmetry [options] -M <masked.root> <out.root> for every A_N mode
// (seconds each), AnPlots -m <masked.root> <dir>, pi0_peak_runs.C on the
// merged file. The yields h_xF* of the top directory carry no Sumw2 (raw
// counts; AnAsymmetry checks that).
//
// Build: g++ -m64 -O2 AnMerge.cpp -o AnMerge `root-config --cflags --libs` -std=c++11

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH3.h>
#include <TKey.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TDatime.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "AnConfig.h"
#include "AnRunTools.h"

static bool StartsWith(const std::string &a, const char *b) { return a.compare(0, std::strlen(b), b) == 0; }
static bool EndsWith(const std::string &a, const char *b)
{
    const size_t n = std::strlen(b);
    return a.size() >= n && a.compare(a.size() - n, n, b) == 0;
}
static bool IsYield(const std::string &n) { return StartsWith(n, "h_xF"); }
// the pi0 candidates summed for the run-by-run peak (as pi0_peak_runs.C):
// EM levels (no suffix) and JP levels not firing EM (_jponly)
static bool IsPi0PeakInput(const std::string &n)
{
    return (StartsWith(n, "h_xFpi_EM") && !EndsWith(n, "_jponly")) || (StartsWith(n, "h_xFpi_JP") && EndsWith(n, "_jponly"));
}

struct FileEntry { std::string name; int run; double pb, py; };

struct Merged
{
    std::map<std::string, TH1 *> plain, wgt[2];
    std::map<int, RunRow> rows;
    std::map<int, TH1D *> pi0m;
    std::map<int, double> pi0sc;
    std::vector<FileEntry> files;
    std::vector<std::string> skipped;
    std::set<std::string> flags;       // *_missing notes of any input
    std::string spin_conv, trg_rule, pi0_source;
    int n_files = 0;

    static void AddTo(std::map<std::string, TH1 *> &m, const std::string &name, TH1 *h, double c)
    {
        std::map<std::string, TH1 *>::iterator it = m.find(name);
        if (it == m.end())
        {
            TH1 *x = (TH1 *)h->Clone(name.c_str());
            x->SetDirectory(0);
            if (c != 1) { x->Reset(); x->Add(h, c); }
            m[name] = x;
        }
        else it->second->Add(h, c);
    }
    void AddPi0(int run, TH1D *m, double sc)
    {
        std::map<int, TH1D *>::iterator it = pi0m.find(run);
        if (it == pi0m.end()) { TH1D *x = (TH1D *)m->Clone(TString::Format("m_%d", run)); x->SetDirectory(0); pi0m[run] = x; }
        else it->second->Add(m);
        if (sc > 0) pi0sc[run] = sc;
    }
    // the TNamed bookkeeping shared by every file
    bool Notes(TDirectory &f, const std::string &fname)
    {
        TNamed *sc = (TNamed *)f.Get("spin_convention"), *tr = (TNamed *)f.Get("trigger_rule");
        if (!sc || !tr) { std::cerr << "ERROR: " << fname << " has no spin_convention / trigger_rule (yield files before 2026-11-05)" << std::endl; return false; }
        if (spin_conv.empty()) { spin_conv = sc->GetTitle(); trg_rule = tr->GetTitle(); }
        else if (spin_conv != sc->GetTitle() || trg_rule != tr->GetTitle())
        {
            std::cerr << "ERROR: " << fname << " has a different spin_convention or trigger_rule from the earlier files" << std::endl;
            return false;
        }
        const char *fl[3] = {"plead_missing", "collins_missing", "pz_missing"};
        for (int i = 0; i < 3; i++) if (f.Get(fl[i])) flags.insert(fl[i]);
        return true;
    }
    // -b: one yield file
    bool AddYieldFile(const std::string &path, const std::map<int, std::pair<double, double> > &pol)
    {
        TFile f(path.c_str());
        if (f.IsZombie()) { skipped.push_back(path + " (cannot open)"); return true; }
        const int run = RunOf(f);
        std::map<int, std::pair<double, double> >::const_iterator it = pol.find(run);
        std::string base = path.substr(path.rfind('/') + 1);
        if (run == 0 || it == pol.end() || !(it->second.first > 0 && it->second.second > 0))
        {
            skipped.push_back(base + (run == 0 ? " (not one run)" : it == pol.end() ? " (run not in the polarization table)" : " (no polarization for a beam)"));
            return true;
        }
        if (!Notes(f, path)) return false;
        const double P[2] = {it->second.first, it->second.second};
        // the pi0 scale factor of this run (TNamed pi0_scale: "<source>; <run>=<f> ...")
        double sc = 0;
        if (TNamed *ps = (TNamed *)f.Get("pi0_scale"))
        {
            const std::string t = ps->GetTitle(), key = TString::Format(" %d=", run).Data();
            const size_t p = t.find(key);
            if (p != std::string::npos) sc = std::atof(t.c_str() + p + key.size());
            if (pi0_source.empty()) pi0_source = t.substr(0, t.find(';'));
        }
        // every histogram once: summed, yields weighted, pi0 peak input projected
        std::map<std::string, TH1 *> cur;
        TH1D *m = 0;
        TIter next(f.GetListOfKeys());
        while (TKey *k = (TKey *)next())
        {
            if (!TClass::GetClass(k->GetClassName())->InheritsFrom(TH1::Class())) continue;
            TH1 *h = (TH1 *)k->ReadObj();
            h->SetDirectory(0);
            const std::string name = k->GetName();
            cur[name] = h;
            AddTo(plain, name, h, 1);
            if (IsYield(name)) for (int b = 0; b < 2; b++) AddTo(wgt[b], name, h, P[b]);
            if (IsPi0PeakInput(name) && h->InheritsFrom(TH3::Class()))
            {
                TH1D *py = ((TH3 *)h)->ProjectionY("pi0m_tmp");
                py->SetDirectory(0);
                if (!m) { m = (TH1D *)py->Clone("pi0m_file"); m->SetDirectory(0); }
                else m->Add(py);
                delete py;
            }
        }
        AddRunRowT([&cur](const std::string &n) { std::map<std::string, TH1 *>::iterator i = cur.find(n); return i == cur.end() ? (TH1 *)0 : i->second; }, false, rows[run]);
        if (m) { AddPi0(run, m, sc); delete m; }
        for (std::map<std::string, TH1 *>::iterator c = cur.begin(); c != cur.end(); ++c) delete c->second;
        FileEntry e = {base, run, P[0], P[1]};
        files.push_back(e);
        n_files++;
        return true;
    }
    // -a / -x: one merged file
    bool AddMerged(const std::string &path)
    {
        TFile f(path.c_str());
        if (f.IsZombie() || !f.Get("merge_info")) { std::cerr << "ERROR: " << path << " is not an AnMerge output" << std::endl; return false; }
        if (f.Get("masked")) { std::cerr << "ERROR: " << path << " is already masked (-x); add unmasked ones" << std::endl; return false; }
        if (!Notes(f, path)) return false;
        TIter next(f.GetListOfKeys());
        while (TKey *k = (TKey *)next())
            if (TClass::GetClass(k->GetClassName())->InheritsFrom(TH1::Class()))
            {
                TH1 *h = (TH1 *)k->ReadObj(); h->SetDirectory(0);
                AddTo(plain, k->GetName(), h, 1);
                delete h;
            }
        const char *wn[2] = {"Pblue", "Pyellow"};
        for (int b = 0; b < 2; b++)
        {
            TDirectory *d = (TDirectory *)f.Get(wn[b]);
            if (!d) { std::cerr << "ERROR: no " << wn[b] << " in " << path << std::endl; return false; }
            TIter nd(d->GetListOfKeys());
            while (TKey *k = (TKey *)nd())
            {
                TH1 *h = (TH1 *)k->ReadObj(); h->SetDirectory(0);
                AddTo(wgt[b], k->GetName(), h, 1);
                delete h;
            }
        }
        ReadRunRows(f, rows);
        if (TDirectory *d = (TDirectory *)f.Get("pi0_runs"))
        {
            TIter nd(d->GetListOfKeys());
            while (TKey *k = (TKey *)nd())
            {
                TH1D *h = (TH1D *)k->ReadObj(); h->SetDirectory(0);
                AddPi0(std::atoi(k->GetName() + 2), h, std::atof(h->GetTitle()));
                delete h;
            }
        }
        if (TNamed *fl = (TNamed *)f.Get("files"))
        {
            std::istringstream is(fl->GetTitle());
            std::string tok;
            while (is >> tok)
            {
                FileEntry e;
                char nm[512];
                if (std::sscanf(tok.c_str(), "%511[^:]:%d:%lf:%lf", nm, &e.run, &e.pb, &e.py) == 4) { e.name = nm; files.push_back(e); }
            }
        }
        if (TNamed *sk = (TNamed *)f.Get("skipped"))
        {
            std::istringstream is(sk->GetTitle());
            std::string line;
            while (std::getline(is, line)) if (!line.empty()) skipped.push_back(line);
        }
        if (TNamed *nf = (TNamed *)f.Get("n_files")) n_files += std::atoi(nf->GetTitle());
        if (TNamed *ps = (TNamed *)f.Get("pi0_scale_source")) if (pi0_source.empty()) pi0_source = ps->GetTitle();
        return true;
    }
    void Write(const std::string &path, const std::string &info, const std::string &bad_note, bool masked)
    {
        TFile out(path.c_str(), "RECREATE");
        // raw counts: the yields carry no Sumw2 (AnAsymmetry's check)
        for (std::map<std::string, TH1 *>::iterator it = plain.begin(); it != plain.end(); ++it)
        {
            if (IsYield(it->first)) it->second->GetSumw2()->Set(0);
            it->second->Write(it->first.c_str());
        }
        std::set<int> runs;
        std::ostringstream fs;
        for (size_t i = 0; i < files.size(); i++)
        {
            runs.insert(files[i].run);
            fs << (i ? " " : "") << files[i].name << ":" << files[i].run << ":" << TString::Format("%.6f:%.6f", files[i].pb, files[i].py).Data();
        }
        std::ostringstream rs, sk, sc;
        for (std::set<int>::const_iterator r = runs.begin(); r != runs.end(); ++r) rs << (r == runs.begin() ? "" : " ") << *r;
        for (size_t i = 0; i < skipped.size(); i++) sk << skipped[i] << "\n";
        sc << (pi0_source.empty() ? "unknown" : pi0_source) << ";";
        for (std::map<int, double>::const_iterator s = pi0sc.begin(); s != pi0sc.end(); ++s) sc << " " << s->first << "=" << TString::Format("%.6f", s->second).Data();
        TNamed("merge_info", info.c_str()).Write();
        TNamed("n_files", TString::Format("%d", n_files).Data()).Write();
        TNamed("runs_used", rs.str().c_str()).Write();
        TNamed("files", fs.str().c_str()).Write();
        TNamed("skipped", sk.str().c_str()).Write();
        TNamed("spin_convention", spin_conv.c_str()).Write();
        TNamed("trigger_rule", trg_rule.c_str()).Write();
        TNamed("pi0_scale", sc.str().c_str()).Write();
        TNamed("pi0_scale_source", pi0_source.c_str()).Write();
        for (std::set<std::string>::const_iterator f = flags.begin(); f != flags.end(); ++f) TNamed(f->c_str(), "set by at least one input file").Write();
        if (masked) { TNamed("masked", "yes").Write(); TNamed("bad_runs", bad_note.c_str()).Write(); }
        WriteRunRows(rows);
        const char *wn[2] = {"Pblue", "Pyellow"};
        for (int b = 0; b < 2; b++)
        {
            TDirectory *d = out.mkdir(wn[b]);
            d->cd();
            for (std::map<std::string, TH1 *>::iterator it = wgt[b].begin(); it != wgt[b].end(); ++it) it->second->Write(it->first.c_str());
            out.cd();
        }
        TDirectory *d = out.mkdir("pi0_runs");
        d->cd();
        for (std::map<int, TH1D *>::iterator it = pi0m.begin(); it != pi0m.end(); ++it)
        {
            std::map<int, double>::const_iterator s = pi0sc.find(it->first);
            it->second->SetTitle(TString::Format("%.6f", s == pi0sc.end() ? 0.0 : s->second));
            it->second->Write(TString::Format("m_%d", it->first));
        }
        out.Close();
    }
};

int main(int argc, char **argv)
{
    if (argc < 4) { std::cerr << "Usage: see the header of AnMerge.cpp (-b | -a | -x)" << std::endl; return 1; }
    const std::string mode = argv[1];
    Merged M;
    if (mode == "-b" && argc >= 5)
    {
        std::map<int, std::pair<double, double> > pol;
        if (!ReadPolTable(argv[3], pol)) { std::cerr << "ERROR: cannot read " << argv[3] << std::endl; return 1; }
        std::ifstream list(argv[2]);
        std::string fn;
        while (list >> fn) if (!M.AddYieldFile(fn, pol)) return 1;
        if (M.n_files == 0) { std::cerr << "ERROR: no usable yield files in " << argv[2] << std::endl; return 1; }
        M.Write(argv[4], TString::Format("built from %d yield files of %s, P from %s", M.n_files, argv[2], argv[3]).Data(), "", false);
        printf("AnMerge -b: %d files, %zu runs, %zu skipped -> %s\n", M.n_files, M.rows.size(), M.skipped.size(), argv[4]);
        return 0;
    }
    if (mode == "-a")
    {
        std::ifstream list(argv[2]);
        std::string fn;
        int n = 0;
        while (list >> fn) { if (!M.AddMerged(fn)) return 1; n++; }
        if (n == 0) { std::cerr << "ERROR: no inputs in " << argv[2] << std::endl; return 1; }
        M.Write(argv[3], TString::Format("%d partial merges, %d yield files", n, M.n_files).Data(), "", false);
        printf("AnMerge -a: %d partial merges, %d files, %zu runs -> %s\n", n, M.n_files, M.rows.size(), argv[3]);
        return 0;
    }
    if (mode == "-x" && argc >= 6)
    {
        const std::string in = argv[2], hdir = argv[3], poltab = argv[4], out = argv[5];
        std::string bad_path = kBadRunListDefault, bad_mode = "flagged", pi0_path;
        for (int i = 6; i + 1 < argc; i += 2)
        {
            const std::string a = argv[i];
            if (a == "-B") bad_path = argv[i + 1];
            else if (a == "-F") bad_mode = argv[i + 1];
            else if (a == "-Q") pi0_path = argv[i + 1];
            else { std::cerr << "ERROR: unknown option " << a << std::endl; return 1; }
        }
        if (bad_mode != "flagged" && bad_mode != "all") { std::cerr << "ERROR: -F flagged | all" << std::endl; return 1; }
        if (!M.AddMerged(in)) return 1;
        std::map<int, std::pair<double, double> > pol;
        if (!ReadPolTable(poltab, pol)) { std::cerr << "ERROR: cannot read " << poltab << std::endl; return 1; }
        // the bad-run lists, as AnAsymmetry -r reads them
        std::string bad_note;
        std::map<int, std::string> bad, bad_pi0;
        if (bad_path != "none")
        {
            if (!ReadBadRuns(bad_path, bad)) { std::cerr << "ERROR: cannot read bad-run list " << bad_path << " (give -B none to run without one)" << std::endl; return 1; }
            FileStat_t st;
            const bool have_st = gSystem->GetPathInfo(bad_path.c_str(), st) == 0;
            int n_flagged = 0;
            for (std::map<int, std::string>::const_iterator b = bad.begin(); b != bad.end(); ++b) n_flagged += !b->second.empty();
            bad_note = TString::Format("%s (modified %s), mode %s: %zu runs listed, %d with a flag", bad_path.c_str(),
                                       have_st ? TDatime((UInt_t)st.fMtime).AsSQLString() : "unknown", bad_mode.c_str(), bad.size(), n_flagged).Data();
        }
        else bad_note = "no bad-run list";
        bad_note += TString::Format(" + %d runs of AnConfig.h kSuppBadRuns", kNSuppBadRuns).Data();
        if (!pi0_path.empty())
        {
            if (!ReadBadRuns(pi0_path, bad_pi0)) { std::cerr << "ERROR: cannot read the pi0-peak bad-run list " << pi0_path << std::endl; return 1; }
            int nf = 0;
            for (std::map<int, std::string>::const_iterator b = bad_pi0.begin(); b != bad_pi0.end(); ++b) nf += !b->second.empty();
            bad_note += TString::Format(" + pi0-peak list %s: %d flagged runs", pi0_path.c_str(), nf).Data();
        }
        // the runs to remove and why (same order of reasons as AnAsymmetry -r)
        std::map<int, std::string> why;
        for (size_t i = 0; i < M.files.size(); i++)
        {
            const int run = M.files[i].run;
            if (why.count(run)) continue;
            if (const char *w = SuppBadRunReason(run)) { why[run] = TString::Format("%d %s (AnConfig.h)", run, w).Data(); continue; }
            std::map<int, std::string>::const_iterator bi = bad.find(run);
            if (bi != bad.end() && (bad_mode == "all" || !bi->second.empty())) { why[run] = TString::Format("%d%s%s", run, bi->second.empty() ? "" : " ", bi->second.c_str()).Data(); continue; }
            std::map<int, std::string>::const_iterator qi = bad_pi0.find(run);
            if (qi != bad_pi0.end() && !qi->second.empty()) why[run] = TString::Format("%d %s (pi0 peak)", run, qi->second.c_str()).Data();
        }
        // subtract their files, with the P each was merged with
        std::vector<FileEntry> kept;
        int n_sub = 0;
        for (size_t i = 0; i < M.files.size(); i++)
        {
            const FileEntry &e = M.files[i];
            if (!why.count(e.run)) { kept.push_back(e); continue; }
            const std::string p = hdir + "/" + e.name;
            TFile f(p.c_str());
            if (f.IsZombie()) { std::cerr << "ERROR: cannot open " << p << " to subtract it" << std::endl; return 1; }
            const double P[2] = {e.pb, e.py};
            TIter next(f.GetListOfKeys());
            while (TKey *k = (TKey *)next())
            {
                if (!TClass::GetClass(k->GetClassName())->InheritsFrom(TH1::Class())) continue;
                TH1 *h = (TH1 *)k->ReadObj(); h->SetDirectory(0);
                const std::string name = k->GetName();
                std::map<std::string, TH1 *>::iterator it = M.plain.find(name);
                if (it == M.plain.end()) { std::cerr << "ERROR: " << name << " of " << p << " is not in the merged file" << std::endl; return 1; }
                it->second->Add(h, -1);
                if (IsYield(name)) for (int b = 0; b < 2; b++) M.wgt[b][name]->Add(h, -P[b]);
                delete h;
            }
            M.n_files--;
            n_sub++;
        }
        M.files = kept;
        // a negative bin means a subtracted file was not in the sum
        int nneg = 0;
        for (std::map<std::string, TH1 *>::iterator it = M.plain.begin(); it != M.plain.end(); ++it)
        {
            if (!IsYield(it->first)) continue;
            TH1 *h = it->second;
            const int ncell = (h->GetNbinsX() + 2) * (h->GetDimension() > 1 ? h->GetNbinsY() + 2 : 1) * (h->GetDimension() > 2 ? h->GetNbinsZ() + 2 : 1);
            for (int c = 0; c < ncell; c++) if (h->GetBinContent(c) < -1e-6) { nneg++; break; }
        }
        if (nneg) { std::cerr << "ERROR: " << nneg << " yield histograms have negative bins after the subtraction" << std::endl; return 1; }
        bad_note += "; removed:";
        for (std::map<int, std::string>::const_iterator w = why.begin(); w != why.end(); ++w) bad_note += " " + w->second;
        M.Write(out, TString::Format("%s minus %d files of %zu masked runs", in.c_str(), n_sub, why.size()).Data(), bad_note, true);
        printf("AnMerge -x: %zu runs masked (%d files subtracted), %d files kept -> %s\n", why.size(), n_sub, M.n_files, out.c_str());
        for (std::map<int, std::string>::const_iterator w = why.begin(); w != why.end(); ++w) printf("  bad run: %s\n", w->second.c_str());
        return 0;
    }
    std::cerr << "Usage: see the header of AnMerge.cpp (-b | -a | -x)" << std::endl;
    return 1;
}
