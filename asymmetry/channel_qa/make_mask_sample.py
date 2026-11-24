#!/usr/bin/env python3
"""Mask periods and the event sample for the hot / dead channel masks (Sean
2026-11-24: divide the run into periods, take a sample of events in each,
determine a mask per period, then run the full pass with those masks).

Periods are chunks of consecutive good runs with about equal events, also
split at the ECal gain changes (Run 22 FCS trigger docs: 23027047, 23066056),
so no period mixes gain sets. In each period, --runs runs spread evenly over
the period and, in each, --files MuDsts spread evenly through the run.

Writes, with <prefix> the output prefix:
  <prefix>_periods.txt  period first_run last_run n_runs events
  <prefix>_runs.txt     the sample's runs: run period events files
  <prefix>.list         MakeRunList.pl-format list for submit_an_job.sh -Q 1

  make_mask_sample.py <good_runs.txt> <prefix> [--periods 25] [--runs 6] [--files 4]
"""
import argparse
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
GAIN_CHANGES = (23027047, 23066056)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("good_runs", help="make_smoke_list.py <prefix>_good_runs.txt (run events files bytes)")
    ap.add_argument("prefix")
    ap.add_argument("--periods", type=int, default=25, help="target number of periods (about equal events)")
    ap.add_argument("--runs", type=int, default=6, help="sample runs per period")
    ap.add_argument("--files", type=int, default=4, help="sample MuDsts per sample run")
    ap.add_argument("--min-files", type=int, default=10, help="sample runs need at least this many files")
    a = ap.parse_args()

    runs = []
    for line in open(a.good_runs):
        t = line.split()
        if t and not line.startswith("#"):
            runs.append((int(t[0]), int(t[1]), int(t[2])))
    runs.sort()
    total = float(sum(r[1] for r in runs))
    target = total / a.periods

    # chunks of about equal events, cut at the gain changes
    periods, cur, ev = [], [], 0
    for r in runs:
        crosses = cur and any(cur[-1][0] < g <= r[0] for g in GAIN_CHANGES)
        if cur and (crosses or ev + r[1] / 2.0 > target):
            periods.append(cur)
            cur, ev = [], 0
        cur.append(r)
        ev += r[1]
    if cur:
        periods.append(cur)

    with open(a.prefix + "_periods.txt", "w") as f:
        f.write("# hot / dead channel mask periods (make_mask_sample.py, target %d, gain changes %s)\n" % (a.periods, GAIN_CHANGES))
        f.write("# period first_run last_run n_runs events\n")
        for i, p in enumerate(periods):
            f.write("%d %d %d %d %d\n" % (i, p[0][0], p[-1][0], len(p), sum(r[1] for r in p)))

    chosen = []
    for i, p in enumerate(periods):
        cand = [r for r in p if r[2] >= a.min_files] or p
        n = min(a.runs, len(cand))
        for k in sorted(set(int((j + 0.5) * len(cand) / n) for j in range(n))):
            chosen.append((cand[k][0], i, cand[k][1], cand[k][2]))
    with open(a.prefix + "_runs.txt", "w") as f:
        f.write("# run period events files\n")
        for c in chosen:
            f.write("%d %d %d %d\n" % c)

    tmp = a.prefix + ".all_files.tmp"
    subprocess.run([os.path.join(HERE, "..", "catalog_query.sh"), "list", tmp] + [str(c[0]) for c in chosen], check=True)
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
        for c in chosen:
            fl = sorted(files.get(c[0], []))
            m = min(a.files, len(fl))
            for i in range(m):
                line = fl[int((i + 0.5) * len(fl) / m)][1]
                f.write(line + "\n")
                nfiles += 1
                nev += int(line.rsplit("/", 1)[1])
    print("%d good runs, %.3g events -> %d periods (%d-%d runs each); sample: %d runs, %d files, %.3g events (%.2f%% of the good-run events) -> %s.list"
          % (len(runs), total, len(periods), min(len(p) for p in periods), max(len(p) for p in periods),
             len(chosen), nfiles, nev, 100.0 * nev / total, a.prefix))


if __name__ == "__main__":
    main()
