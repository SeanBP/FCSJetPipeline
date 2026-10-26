#!/usr/bin/env python3
"""Extrapolates benchmark an_job timing tables to a full Run 22 A_N pass.

Usage:
  ./estimate_resources.py <timing dir> [--catalog catalog_runs_local.txt]
      [--files-per-job 40] [--slots 1000,3000,5000] [--run-min 23005044]
      [--run-max 23108024]

Inputs: the timing_<run>_<job>.txt tables of benchmark jobs (an_job.xml) and
the per-run catalog counts (catalog_query.sh counts). Runs are those in
[--run-min, --run-max] minus MakeRunList.pl's excluded ranges; its RunLog
bad-run filter is not applied, so the totals are an upper bound.

Model, per job of n files of one run:
  per file   fetch (wall, per GB) + MuDst -> SimpleTree (wall = a + b * events,
             fitted) + jet finding (mean per file)
  per job    compiles + calibration + histograms (mean)
             + UE: measured per-job time scaled as (n / n_measured)^2, since
             every jet searches the job's whole proxy pool
  memory     largest peak RSS of any stage; scratch: largest scratch use,
             scaled by n / n_measured for the SimpleTree part
"""

import argparse
import glob
import math
import os
import statistics
import sys

EXCLUDED = [(23003063, 23004033), (23007067, 23007076), (23027046, 23027048)]


def read_timing(path):
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        if len(f) < 9:
            continue
        rows.append(dict(stage=f[0], file=f[1], events=float(f[2]), in_kb=float(f[3]), wall=float(f[4]),
                         user=float(f[5]), sys=float(f[6]), rss_kb=float(f[7]), scratch_kb=float(f[8])))
    return rows


def fit_line(x, y):
    n = len(x)
    if n < 2:
        return 0.0, (y[0] if y else 0.0)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    if sxx == 0:
        return my, 0.0
    b = sum((a - mx) * (c - my) for a, c in zip(x, y)) / sxx
    return my - b * mx, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('timing_dir')
    ap.add_argument('--catalog', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'catalog_runs_local.txt'))
    ap.add_argument('--files-per-job', type=int, default=40)
    ap.add_argument('--slots', default='1000,3000,5000')
    ap.add_argument('--run-min', type=int, default=23005044)
    ap.add_argument('--run-max', type=int, default=23108024)
    a = ap.parse_args()

    tables = sorted(glob.glob(os.path.join(a.timing_dir, 'timing_*.txt')))
    if not tables:
        sys.exit('no timing_*.txt in ' + a.timing_dir)
    jobs = [read_timing(t) for t in tables]
    jobs = [j for j in jobs if any(r['stage'] == 'hists' for r in j)]
    print('benchmark: %d complete jobs (%d tables)' % (len(jobs), len(tables)))
    if not jobs:
        sys.exit('no complete job')

    rows = [r for j in jobs for r in j]
    fetch = [r for r in rows if r['stage'] == 'fetch']
    mudst = [r for r in rows if r['stage'] == 'mudst']
    jets = [r for r in rows if r['stage'] == 'jets']
    per_job = {s: [sum(r['wall'] for r in j if r['stage'] == s) for j in jobs] for s in ('ue', 'calib', 'hists')}
    per_job['compile'] = [sum(r['wall'] for r in j if r['stage'] == 'compile') for j in jobs]
    n_meas = statistics.mean(sum(1 for r in j if r['stage'] == 'jets') for j in jobs)

    a0, b0 = fit_line([r['events'] for r in mudst], [r['wall'] for r in mudst])
    fetch_per_gb = sum(r['wall'] for r in fetch) / max(1e-9, sum(r['in_kb'] for r in fetch) / 1e6)
    jets_per_file = statistics.mean(r['wall'] for r in jets)
    cpu_frac = sum(r['user'] + r['sys'] for r in mudst) / max(1e-9, sum(r['wall'] for r in mudst))
    ev_per_file = statistics.mean(r['events'] for r in mudst)
    print('\nper file (%d files, mean %.0f events, %.2f GB):' % (len(mudst), ev_per_file, statistics.mean(r['in_kb'] for r in fetch) / 1e6))
    print('  fetch           %6.1f s  (%.1f s/GB)' % (statistics.mean(r['wall'] for r in fetch), fetch_per_gb))
    print('  MuDst->Simple   %6.1f s  (fit %.1f s + %.2f ms/event; CPU/wall %.2f)' % (
        statistics.mean(r['wall'] for r in mudst), a0, 1e3 * b0, cpu_frac))
    print('  jet finding     %6.1f s' % jets_per_file)
    print('per job (measured n = %.1f files):' % n_meas)
    for s in ('compile', 'ue', 'calib', 'hists'):
        print('  %-15s %6.1f s  (max %.1f)' % (s, statistics.mean(per_job[s]), max(per_job[s])))
    rss = max(r['rss_kb'] for r in rows) / 1024.
    rss_stage = max(rows, key=lambda r: r['rss_kb'])['stage']
    scratch = max(r['scratch_kb'] for r in rows) / 1e6
    st_mb = statistics.mean(r['in_kb'] for r in jets) / 1024.
    out_kb = statistics.mean(r['in_kb'] for j in jobs for r in j if r['stage'] == 'output')
    walls = [sum(r['wall'] for r in j) for j in jobs]
    print('  job wall        %6.1f min (min %.1f, max %.1f)' % (statistics.mean(walls) / 60, min(walls) / 60, max(walls) / 60))
    print('  peak RSS        %6.0f MB (%s)' % (rss, rss_stage))
    print('  peak scratch    %6.2f GB; SimpleTree %.0f MB/file' % (scratch, st_mb))
    print('  output          %6.0f kB histograms per job' % out_kb)

    runs = []
    for line in open(a.catalog):
        f = line.strip().split('::')
        if len(f) < 4:
            continue
        r = int(f[0])
        if r < a.run_min or r > a.run_max or any(lo <= r <= hi for lo, hi in EXCLUDED):
            continue
        runs.append((r, float(f[1]), int(f[2]), float(f[3])))
    n = a.files_per_job
    tot_files = sum(x[2] for x in runs)
    tot_ev = sum(x[1] for x in runs)
    tot_gb = sum(x[3] for x in runs) / 1e9
    njobs = 0
    total_s = 0.0
    job_walls = []
    for r, ev, nf, size in runs:
        c = int(math.ceil(nf / float(n)))
        njobs += c
        evf, gbf = ev / nf, size / nf / 1e9
        for k in range(c):
            m = nf // c + (1 if k < nf % c else 0)
            w = m * (fetch_per_gb * gbf + a0 + b0 * evf + jets_per_file) \
                + statistics.mean(per_job['compile']) + statistics.mean(per_job['calib']) * m / n_meas \
                + statistics.mean(per_job['hists']) * m / n_meas + statistics.mean(per_job['ue']) * (m / n_meas) ** 2
            total_s += w
            job_walls.append(w)
    job_walls.sort()
    print('\nfull pass: runs %d-%d, %d runs, %d files, %.2e events, %.0f TB read' % (
        a.run_min, a.run_max, len(runs), tot_files, tot_ev, tot_gb / 1e3))
    print('  jobs            %d of <= %d files' % (njobs, n))
    print('  job wall        median %.0f min, 99%% %.0f min, max %.0f min' % (
        job_walls[len(job_walls) // 2] / 60, job_walls[int(0.99 * (len(job_walls) - 1))] / 60, job_walls[-1] / 60))
    print('  total           %.0f core-hours (%.0f core-days)' % (total_s / 3600, total_s / 86400))
    for s in [int(x) for x in a.slots.split(',')]:
        print('  on %5d slots   %.1f days' % (s, total_s / 3600 / s / 24))
    print('  memory request  %.0f MB (1.5 x peak RSS, rounded up to 500 MB)' % (500 * math.ceil(1.5 * rss / 500)))
    print('  scratch         %.1f GB per job (one MuDst + %d SimpleTrees)' % (
        statistics.mean(r['in_kb'] for r in fetch) / 1e6 + n * st_mb / 1024., n))
    print('  output          %.1f GB of histograms (%d files; %.0f MB after hadd per run)' % (
        njobs * out_kb / 1e6, njobs, len(runs) * out_kb / 1e3))


if __name__ == '__main__':
    main()
