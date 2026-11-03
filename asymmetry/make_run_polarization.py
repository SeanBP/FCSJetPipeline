#!/usr/bin/env python3
"""Per-run beam polarizations for Run 22, for AnAsymmetry's -r mode.

Each run's polarization is its fill's at the middle of the run, as the FMS
EM-jet analysis (Kabir, Liang analysis note, sec. 5.3):
    P_run = P0 + dP/dt * (t_mid - t_fill_start)
from the RHIC polarimetry fill table (P0 in %, dP/dt in %/h, fill start in
unix time) and the run's start and end from STAR's RunLog (runDescriptor;
fill number from beamInfo). Error: dP0 and d(dP/dt) * (t_mid - t_fill_start)
in quadrature (the fill table's errors; the polarization scale uncertainty
is not included).

Output, one line per run with a fill in the table:
    run fill t_mid P_blue dP_blue P_yellow dP_yellow
as fractions (not %); a beam the table gives no value for (P0 = 0) is -1.
Runs whose fill is not in the table, or whose RunLog times put the run's
middle more than 1 h before its fill's start or more than 24 h after it (a
few runs have start times of 0), are left out and counted.

Usage:
    ./make_run_polarization.py <fill table> <out.txt> [first run = 23005044] [last run = 23108024] [--runlog <cached RunLog dump>]
The fill table is the polarimetry page's text (16 columns: fill energy start
stop P0 +- dP0 dP/dt +- d(dP/dt) for blue, then yellow), e.g.
/star/u/dkap7827/Fcs2019/FcsAsymRun22/Run22PolFromWeb_20250205.txt (copied
2025-02-05; https://wiki.bnl.gov/rhicspin/Run_22_polarization, oldid 6020).
RunLog is read with the mysql client (db04.star.bnl.gov, port 3421), unless
--runlog gives a dump of "run blueFill yellowFill start end" lines.
"""

import math
import subprocess
import sys


def read_fills(path):
    fills = {}
    for line in open(path):
        t = line.split()
        if len(t) != 16 or not t[0].isdigit():
            continue
        fill, start = int(t[0]), int(t[2])
        # tokens: fill E start stop P0 +- dP0 dPdt +- ddPdt  (blue)  P0 +- dP0 dPdt +- ddPdt  (yellow)
        beams = []
        for o in (4, 10):
            p0, dp0, slope, dslope = float(t[o]), float(t[o + 2]), float(t[o + 3]), float(t[o + 5])
            beams.append(None if p0 == 0 else (p0, dp0, slope, dslope))
        fills[fill] = (start, beams)
    return fills


def read_runlog(first, last, dump):
    if dump:
        lines = open(dump).read().splitlines()
    else:
        query = ("SELECT r.runNumber, MIN(b.blueFillNumber), MIN(b.yellowFillNumber), MIN(r.startRunTime), "
                 "MAX(r.endRunTime) FROM RunLog.runDescriptor r JOIN RunLog.beamInfo b ON r.runNumber=b.runNumber "
                 "WHERE r.runNumber BETWEEN %d AND %d GROUP BY r.runNumber" % (first, last))
        out = subprocess.run(["mysql", "-h", "db04.star.bnl.gov", "--port=3421", "--connect-timeout=10", "-N", "-s",
                              "-e", query], capture_output=True, text=True, check=True)
        lines = out.stdout.splitlines()
    runs = []
    for line in lines:
        t = line.split()
        if len(t) == 5:
            runs.append((int(t[0]), int(float(t[1])), int(float(t[2])), int(t[3]), int(t[4])))
    return runs


def main():
    args = [a for a in sys.argv[1:]]
    dump = None
    if "--runlog" in args:
        i = args.index("--runlog")
        dump = args[i + 1]
        del args[i:i + 2]
    if len(args) < 2:
        sys.exit(__doc__)
    fills = read_fills(args[0])
    first = int(args[2]) if len(args) > 2 else 23005044
    last = int(args[3]) if len(args) > 3 else 23108024
    runs = read_runlog(first, last, dump)

    n_out, n_nofill, n_mismatch, n_outside, n_badtime = 0, 0, 0, 0, 0
    with open(args[1], "w") as f:
        f.write("# run fill t_mid P_blue dP_blue P_yellow dP_yellow   (fractions; -1 = no value)\n")
        f.write("# P = P0 + dP/dt (t_mid - fill start), fill table %s\n" % args[0])
        for run, bfill, yfill, start, end in runs:
            if bfill != yfill:
                n_mismatch += 1
            if bfill not in fills:
                n_nofill += 1
                continue
            fstart, beams = fills[bfill]
            t_mid = 0.5 * (start + end)
            dt_h = (t_mid - fstart) / 3600.0
            if dt_h < -1 or dt_h > 24 or end < start:
                n_badtime += 1
                continue
            if dt_h < 0:
                n_outside += 1
            vals = []
            for b in beams:
                if b is None:
                    vals += [-1.0, -1.0]
                else:
                    p0, dp0, slope, dslope = b
                    vals += [(p0 + slope * dt_h) / 100.0, math.hypot(dp0, dslope * dt_h) / 100.0]
            f.write("%d %d %.0f %.4f %.4f %.4f %.4f\n" % (run, bfill, t_mid, *vals))
            n_out += 1
    print("make_run_polarization: %d runs written, %d with no fill in the table, %d with bad times, %d with different "
          "blue/yellow fill numbers, %d with their middle up to 1 h before the fill's start time"
          % (n_out, n_nofill, n_badtime, n_mismatch, n_outside))


if __name__ == "__main__":
    main()
