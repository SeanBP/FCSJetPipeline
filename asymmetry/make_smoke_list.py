#!/usr/bin/env python3
"""Good-run list and an evenly spread test sample of Run 22 for the A_N chain.

Good runs: the catalog's analysis runs (23005044-23108024 minus
MakeRunList.pl's excluded ranges) that pass the shift-crew flag (RunLog
runStatus, MakeRunList.pl's rule: rtsStatus = 0 and shiftLeaderStatus <= 1),
are not flagged in the pi0 bad-run list (AnAsymmetry's default mode; the
list has a line for every reviewed run, an empty flag meaning it passed)
and have both beams' polarization in run22_polarization.txt. This is the
sample AnAsymmetry -r keeps from a full pass, so its catalog events are the
"full data" an error projection scales to (project_errors.py).

Test sample: --runs runs evenly spaced in the good-run list (sorted by run
number, so across the whole of Run 22), and in each up to --files files
evenly spaced through the run (by file sequence number), so each job's files
span its run.

Usage:
  ./make_smoke_list.py <out prefix> [--runs 60] [--files 40]
      [--catalog catalog_runs_local.txt] [--pol run22_polarization.txt]
      [--bad /star/u/dkap7827/Fcs2019/FcsAsymRun22/BadRuns.txt]
      [--min-files 20]
Writes <prefix>_good_runs.txt (run events files bytes, every good run),
<prefix>_runs.txt (the chosen runs) and <prefix>.list (MakeRunList.pl format,
for submit_an_job.sh -f). --min-files leaves out runs with fewer files than
this from the choice (short runs), not from the good-run list.
"""

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXCLUDED = [(23003063, 23004033), (23007067, 23007076), (23027046, 23027048)]
RUN_MIN, RUN_MAX = 23005044, 23108024


def runlog_bad(first, last):
    q = ("SELECT runNumber, rtsStatus, shiftLeaderStatus FROM runStatus WHERE runNumber BETWEEN %d AND %d "
         "AND deactive=0 AND flavor='ofl'" % (first, last))
    out = subprocess.run(["mysql", "-h", "db04.star.bnl.gov", "--port=3421", "--connect-timeout=10", "-N", "-s",
                          "-e", q, "RunLog"], capture_output=True, text=True, check=True).stdout
    bad = set()
    for line in out.splitlines():
        t = line.split()
        if len(t) == 3 and not (int(t[1]) == 0 and int(t[2]) <= 1):
            bad.add(int(t[0]))
    if not out.strip():
        sys.exit("RunLog.runStatus returned nothing")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--runs", type=int, default=60)
    ap.add_argument("--files", type=int, default=40)
    ap.add_argument("--catalog", default=os.path.join(HERE, "catalog_runs_local.txt"))
    ap.add_argument("--pol", default=os.path.join(HERE, "run22_polarization.txt"))
    ap.add_argument("--bad", default="/star/u/dkap7827/Fcs2019/FcsAsymRun22/BadRuns.txt")
    ap.add_argument("--min-files", type=int, default=20)
    a = ap.parse_args()

    cat = {}
    for line in open(a.catalog):
        f = line.strip().split("::")
        if len(f) == 4:
            r = int(f[0])
            if RUN_MIN <= r <= RUN_MAX and not any(lo <= r <= hi for lo, hi in EXCLUDED):
                cat[r] = (int(float(f[1])), int(f[2]), int(float(f[3])))
    pol = set()
    for line in open(a.pol):
        t = line.split()
        if t and not line.startswith("#") and float(t[3]) > 0 and float(t[5]) > 0:
            pol.add(int(t[0]))
    bad = set()
    for line in open(a.bad):
        t = line.split("|")
        if len(t) > 1 and t[0].strip().isdigit() and t[1].strip():
            bad.add(int(t[0]))
    rl = runlog_bad(RUN_MIN, RUN_MAX)

    good, n_rl, n_bad, n_pol = [], 0, 0, 0
    for r in sorted(cat):
        if r in rl:
            n_rl += 1
        elif r in bad:
            n_bad += 1
        elif r not in pol:
            n_pol += 1
        else:
            good.append(r)
    with open(a.prefix + "_good_runs.txt", "w") as f:
        f.write("# good runs: catalog %s, RunLog shift-crew flag, not in %s, polarization in %s\n" % (a.catalog, a.bad, a.pol))
        f.write("# run events files bytes\n")
        for r in good:
            f.write("%d %d %d %d\n" % ((r,) + cat[r]))
    ev = sum(cat[r][0] for r in good)
    print("catalog analysis runs %d: %d RunLog-flagged, %d flagged in the pi0 bad-run list, %d no polarization -> %d good runs, "
          "%d files, %.3g events" % (len(cat), n_rl, n_bad, n_pol, len(good), sum(cat[r][1] for r in good), ev))

    cand = [r for r in good if cat[r][1] >= a.min_files]
    n = min(a.runs, len(cand))
    chosen = sorted(set(cand[int((i + 0.5) * len(cand) / n)] for i in range(n)))
    with open(a.prefix + "_runs.txt", "w") as f:
        for r in chosen:
            f.write("%d %d %d %d\n" % ((r,) + cat[r]))

    tmp = a.prefix + ".all_files.tmp"
    subprocess.run([os.path.join(HERE, "catalog_query.sh"), "list", tmp] + [str(r) for r in chosen], check=True)
    files = {}
    for line in open(tmp):
        line = line.strip()
        base = line.rsplit("/", 2)[-2]  # <name>.MuDst.root
        run = int(base.split("_")[2])
        seq = int(base.split("_")[-1].split(".")[0])
        files.setdefault(run, []).append((seq, line))
    os.remove(tmp)
    nfiles, nev = 0, 0
    with open(a.prefix + ".list", "w") as f:
        for r in chosen:
            fl = sorted(files.get(r, []))
            m = min(a.files, len(fl))
            for i in range(m):
                line = fl[int((i + 0.5) * len(fl) / m)][1]
                f.write(line + "\n")
                nfiles += 1
                nev += int(line.rsplit("/", 1)[1])
    print("test sample: %d runs (%d-%d), %d files, %.3g events = %.4f of the good-run events -> %s.list"
          % (len(chosen), chosen[0], chosen[-1], nfiles, nev, nev / float(ev), a.prefix))


if __name__ == "__main__":
    main()
