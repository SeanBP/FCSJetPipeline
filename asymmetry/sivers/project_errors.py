#!/usr/bin/env python3
"""Projects a test pass's statistical A_N errors to the full good-run sample.

Statistical errors scale as 1/sqrt(jets), and jets scale with the events
read, so
    stat_full = stat_test * sqrt(events_test / events_full)
with events_test the catalog events of the MuDsts the test jobs processed
(timing tables: files that made a SimpleTree, in runs with a yield file) and
events_full the catalog events of every good run (make_smoke_list.py's
<prefix>_good_runs.txt: RunLog flag, flagged pi0 bad runs and runs without
polarization removed). This assumes the test runs are representative (spread
evenly over Run 22) in jets per event and polarization. Only test runs in
the good-run list count, the runs AnAsymmetry -r keeps (same cuts); runs it
removed (listed in an.txt's "# bad runs" line) count in neither sample.

Usage:
  ./project_errors.py <test outdir> <an.txt> <good_runs.txt>
Writes <an>_projected.txt (AnAsymmetry's table plus stat_full) and prints the
scale factor and the final points (level all, class all).
"""

import glob
import math
import os
import re
import sys


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    outdir, an_txt, good = sys.argv[1:4]
    good_ev = {}
    for line in open(good):
        t = line.split()
        if t and not line.startswith("#"):
            good_ev[int(t[0])] = int(t[1])
    # runs AnAsymmetry left out (bad-run lists, AnConfig.h kSuppBadRuns): the
    # "# bad runs: ... removed: <run> <reason> ..." header line of an.txt
    removed = set()
    for line in open(an_txt):
        if line.startswith("# bad runs:") and "removed:" in line:
            removed = set(int(r) for r in re.findall(r"\b(2[0-9]{7})\b", line.split("removed:", 1)[1]))
    for r in removed:
        good_ev.pop(r, None)
    hist_runs = set()
    for f in glob.glob(os.path.join(outdir, "hists", "an_hists_*.root")):
        hist_runs.add(int(os.path.basename(f).split("_")[2]))

    ev_test, files_test, runs_test = 0, 0, set()
    for tab in glob.glob(os.path.join(outdir, "timing", "timing_*.txt")):
        mudst, made = {}, set()
        for line in open(tab):
            t = line.split()
            if len(t) < 3 or t[0].startswith("#"):
                continue
            if t[0] == "mudst":
                m = re.search(r"_(\d{8})_(?:raw_)?(\d{7})\.MuDst", t[1])
                mudst[(int(m.group(1)), int(m.group(2)))] = int(t[2])
            elif t[0] == "jets":
                m = re.search(r"_(\d{8})_(\d{7})\.root", t[1])
                made.add((int(m.group(1)), int(m.group(2))))
        for key in made:
            if key[0] in hist_runs and key[0] in good_ev and key in mudst:
                ev_test += mudst[key]
                files_test += 1
                runs_test.add(key[0])
    ev_full = sum(good_ev.values())
    if ev_test == 0:
        sys.exit("no processed files found in %s" % outdir)
    scale = math.sqrt(ev_test / float(ev_full))
    note = ("test: %d runs, %d files, %.4g events; full: %d good runs, %.4g events; fraction %.5f; "
            "stat_full = stat x %.4f; %d runs left out by AnAsymmetry not counted" % (len(runs_test), files_test, ev_test, len(good_ev), ev_full,
                                         ev_test / float(ev_full), scale, len(removed)))
    print(note)

    out = an_txt[:-4] + "_projected.txt" if an_txt.endswith(".txt") else an_txt + "_projected"
    with open(out, "w") as f:
        f.write("# %s\n" % note)
        for line in open(an_txt):
            if line.startswith("# level"):
                f.write(line.rstrip("\n") + " stat_full\n")
            elif line.startswith("#"):
                f.write(line)
            else:
                t = line.split()
                f.write(line.rstrip("\n") + " %.6f\n" % (float(t[7]) * scale))
    print("\nlevel all, class all:  beam  x_F bin  <x_F>  A_N  stat(test)  stat(full, projected)  jets(test)")
    for line in open(out):
        t = line.split()
        if len(t) > 11 and t[0] == "all" and t[1] == "all":
            print("  %-6s %6.3f-%6.3f %7.4f %+8.4f %8.4f %9.5f %9.0f" % (t[2], float(t[3]), float(t[4]), float(t[5]),
                                                                       float(t[6]), float(t[7]), float(t[11]), float(t[8])))
    print("wrote %s" % out)
    print("scale=%.6f" % scale)


if __name__ == "__main__":
    main()
