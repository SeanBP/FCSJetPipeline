#!/usr/bin/env python3
"""Splits a pass's per-run yield files by time in fill (2026-11-17).

Time in fill of a run = its middle (t_mid, asymmetry/run22_polarization.txt,
from the RunLog) minus its fill's start (the polarimetry fill table that
make_run_polarization.py used, column 3). Runs are split into bins of time in
fill, each bin's yield files written as a list for AnAsymmetry -r; the A_N of
the bins are then compared with sys_slices.py --disjoint (a consistency check:
the FMS EM-jet analysis compared runs < 2 h and > 2 h into the fill,
BrightSTAR scripts/ChiSquareBkg.C).

Usage:
  time_in_fill.py <yields.list> <out prefix> [--edges 2] [--quantiles N]
                  [--runs ../run22_polarization.txt] [--fills <fill table>]
  --edges h1,h2,...: bin edges in hours (default 2: < 2 h and >= 2 h)
  --quantiles N: N bins with equal numbers of yield files instead
Writes <out prefix>_bin<i>.list and <out prefix>_runs.txt (run, fill, hours).
Yield files are an_hists_<run>_<tag>.root; a run with no time is left out.
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FILL_TABLE = '/star/u/dkap7827/Fcs2019/FcsAsymRun22/Run22PolFromWeb_20250205.txt'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('yields')
    ap.add_argument('out')
    ap.add_argument('--edges', default='2')
    ap.add_argument('--quantiles', type=int, default=0)
    ap.add_argument('--runs', default=os.path.join(HERE, '..', 'run22_polarization.txt'))
    ap.add_argument('--fills', default=FILL_TABLE)
    a = ap.parse_args()

    fstart = {}
    for line in open(a.fills):
        t = line.split()
        if len(t) > 3 and t[0].isdigit():
            fstart[int(t[0])] = int(t[2])
    hours = {}
    fill_of = {}
    for line in open(a.runs):
        if line.startswith('#'):
            continue
        t = line.split()
        run, fill, tmid = int(t[0]), int(t[1]), float(t[2])
        if fill in fstart:
            hours[run] = (tmid - fstart[fill]) / 3600.0
            fill_of[run] = fill

    files = [l.strip() for l in open(a.yields) if l.strip()]
    rx = re.compile(r'an_hists_(\d+)_')
    entries, missing = [], 0
    for f in files:
        m = rx.search(os.path.basename(f))
        if not m or int(m.group(1)) not in hours:
            missing += 1
            continue
        entries.append((hours[int(m.group(1))], int(m.group(1)), f))
    entries.sort()
    if not entries:
        print('time_in_fill: no yield file with a known run time')
        return 1

    if a.quantiles > 1:
        n = len(entries)
        cuts = [entries[min(n - 1, (i * n) // a.quantiles)][0] for i in range(1, a.quantiles)]
    else:
        cuts = [float(x) for x in a.edges.split(',') if x]
    edges = [-1e9] + cuts + [1e9]
    bins = [[] for _ in range(len(edges) - 1)]
    for h, run, f in entries:
        for i in range(len(edges) - 1):
            if edges[i] <= h < edges[i + 1]:
                bins[i].append((h, run, f))
                break

    with open(a.out + '_runs.txt', 'w') as o:
        o.write('# run fill hours_in_fill (t_mid - fill start); bin edges (h): %s\n' % ' '.join('%.3f' % c for c in cuts))
        for run in sorted(set(r for _, r, _ in entries)):
            o.write('%d %d %.3f\n' % (run, fill_of[run], hours[run]))
    for i, b in enumerate(bins):
        with open('%s_bin%d.list' % (a.out, i), 'w') as o:
            for _, _, f in b:
                o.write(f + '\n')
        lo = '%.2f' % edges[i] if i else '-inf'
        hi = '%.2f' % edges[i + 1] if i < len(bins) - 1 else 'inf'
        hs = [h for h, _, _ in b]
        print('time_in_fill: bin %d [%s, %s) h: %d yield files, %d runs, median %.2f h' % (
            i, lo, hi, len(b), len(set(r for _, r, _ in b)), sorted(hs)[len(hs) // 2] if hs else float('nan')))
    if missing:
        print('time_in_fill: %d yield files left out (run not in the time table)' % missing)
    return 0


if __name__ == '__main__':
    sys.exit(main())
