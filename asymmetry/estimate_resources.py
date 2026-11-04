#!/usr/bin/env python3
"""Extrapolates benchmark an_job timing tables to a full Run 22 A_N pass.

Usage:
  ./estimate_resources.py <timing dir> [--catalog catalog_runs_local.txt]
      [--files-per-job 40] [--slots 1000,3000,5000] [--run-min 23005044]
      [--run-max 23108024] [--runs <good_runs.txt>]

Inputs: the timing_<run>_<job>.txt tables of benchmark jobs (an_job.xml) and
the per-run catalog counts (catalog_query.sh counts). Runs are those in
[--run-min, --run-max] minus MakeRunList.pl's excluded ranges. The RunLog
shift-crew flags and the pi0 bad-run list (applied later, by MakeRunList.pl
and AnAsymmetry -r) are not taken off here, so the totals are an upper bound,
unless --runs gives a run list (first column; make_smoke_list.py's
<prefix>_good_runs.txt) to restrict the pass to.

Model, per job of n files of one run:
  per file   fetch (wall, per GB) + MuDst -> SimpleTree (wall = a + b * events,
             fitted) + jet finding (mean per file)
  per job    compiles (mean, fixed) + UE + calibration + histograms (mean,
             each scaled as n / n_measured: UE is in-event (weighted wedges), so
             each file is independent)
  failures   failed fetches (retried and skipped) cost their measured wall
             time per processed file, reported separately and with the total
  memory     largest peak RSS of any stage
  scratch    an_job.xml deletes each MuDst once its SimpleTree is made and
             each SimpleTree once its jets and UE are done, so the peak is
             one MuDst + its SimpleTree (paired by run and file sequence) +
             the job's jet files (~n x the mean jet file) + the sandbox.
             Reported as percentiles over the benchmark's files and as a
             disk request (1.25 x the largest, rounded up to 1 GB). Benchmarks
             from the older layout (all SimpleTrees kept to the end) give the
             same per-file sizes; only their measured peak differs
"""

import argparse
import glob
import math
import os
import re
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
        return (y[0] if y else 0.0), 0.0
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
    ap.add_argument('--runs', default=None)
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
    failed = [r for r in rows if r['stage'] in ('fetchfail', 'fetchretry')]
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
    fail_s_per_file = sum(r['wall'] for r in failed) / max(1, len(fetch))
    # per-file scratch peak: the MuDst + its SimpleTree, paired by (run, file sequence)
    def key(name):
        m = re.search(r'_(\d{8})_(?:raw_)?(\d{7})', name)
        return (m.group(1), m.group(2)) if m else None
    st_by_key = {}
    for r in jets:
        if key(r['file']):
            st_by_key[key(r['file'])] = r['in_kb']
    pair_gb = sorted((r['in_kb'] + st_by_key[key(r['file'])]) / 1e6 for r in fetch if key(r['file']) in st_by_key)
    jf_kb = [r['in_kb'] / max(1, sum(1 for x in j if x['stage'] == 'jets')) for j in jobs for r in j if r['stage'] == 'jetfiles']
    jf_mb = (statistics.mean(jf_kb) / 1e3) if jf_kb else 0.3
    jf_max_mb = (max(jf_kb) / 1e3) if jf_kb else 3.0
    sandbox_gb = 0.05
    if failed:
        print('  failed fetches  %d retried, %d skipped, of %d files (%.0f s lost in all)' % (
            sum(1 for r in failed if r['stage'] == 'fetchretry'), sum(1 for r in failed if r['stage'] == 'fetchfail'),
            len(fetch) + sum(1 for r in failed if r['stage'] == 'fetchfail'), sum(r['wall'] for r in failed)))
    slow = sorted(r['wall'] for r in fetch)
    print('  fetch wall      median %.0f s, 95%% %.0f s, max %.0f s' % (slow[len(slow) // 2], slow[int(0.95 * (len(slow) - 1))], slow[-1]))
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
    print('  peak scratch    %6.2f GB as measured (this benchmark\'s job layout); SimpleTree %.0f MB/file (max %.0f)' % (
        scratch, st_mb, max(r['in_kb'] for r in jets) / 1024.))
    if pair_gb:
        q = lambda f: pair_gb[min(len(pair_gb) - 1, int(f * (len(pair_gb) - 1)))]
        print('  MuDst + its SimpleTree per file: median %.2f GB, 99%% %.2f GB, max %.2f GB (%d files)' % (q(0.5), q(0.99), pair_gb[-1], len(pair_gb)))
    print('  jet files       %.2f MB per MuDst (max %.1f)' % (jf_mb, jf_max_mb))
    print('  output          %6.0f kB histograms per job' % out_kb)

    keep = None
    if a.runs:
        keep = set(int(l.split()[0]) for l in open(a.runs) if l.strip() and not l.startswith('#'))
    runs = []
    for line in open(a.catalog):
        f = line.strip().split('::')
        if len(f) < 4:
            continue
        r = int(f[0])
        if r < a.run_min or r > a.run_max or any(lo <= r <= hi for lo, hi in EXCLUDED):
            continue
        if keep is not None and r not in keep:
            continue
        runs.append((r, float(f[1]), int(f[2]), float(f[3])))
    n = a.files_per_job
    tot_files = sum(x[2] for x in runs)
    tot_ev = sum(x[1] for x in runs)
    tot_gb = sum(x[3] for x in runs) / 1e9
    njobs = 0
    total_s = 0.0
    job_walls = []
    job_walls_f = []
    fail_total_s = 0.0
    for r, ev, nf, size in runs:
        c = int(math.ceil(nf / float(n)))
        njobs += c
        evf, gbf = ev / nf, size / nf / 1e9
        for k in range(c):
            m = nf // c + (1 if k < nf % c else 0)
            w = m * (fetch_per_gb * gbf + a0 + b0 * evf + jets_per_file) + statistics.mean(per_job['compile']) \
                + (statistics.mean(per_job['ue']) + statistics.mean(per_job['calib']) + statistics.mean(per_job['hists'])) * m / n_meas
            total_s += w
            job_walls.append(w)
            fail_total_s += m * fail_s_per_file
            job_walls_f.append(w + m * fail_s_per_file)
    job_walls.sort()
    job_walls_f.sort()
    print('\nfull pass%s: runs %d-%d, %d runs, %d files, %.2e events, %.0f TB read' % (
        ' (runs in ' + a.runs + ')' if a.runs else '', a.run_min, a.run_max, len(runs), tot_files, tot_ev, tot_gb / 1e3))
    print('  jobs            %d of <= %d files' % (njobs, n))
    print('  job wall        median %.0f min, 99%% %.0f min, max %.0f min' % (
        job_walls[len(job_walls) // 2] / 60, job_walls[int(0.99 * (len(job_walls) - 1))] / 60, job_walls[-1] / 60))
    print('  job wall with failed fetches: median %.0f min, 99%% %.0f min, max %.0f min (%.0f s lost per processed file)' % (
        job_walls_f[len(job_walls_f) // 2] / 60, job_walls_f[int(0.99 * (len(job_walls_f) - 1))] / 60, job_walls_f[-1] / 60, fail_s_per_file))
    print('  total           %.0f core-hours (%.0f core-days); with failed fetches %.0f slot-hours' % (
        total_s / 3600, total_s / 86400, (total_s + fail_total_s) / 3600))
    for s in [int(x) for x in a.slots.split(',')]:
        print('  on %5d slots   %.1f days (%.1f with failed fetches)' % (s, total_s / 3600 / s / 24, (total_s + fail_total_s) / 3600 / s / 24))
    print('  memory request  %.0f MB (1.5 x peak RSS, rounded up to 500 MB)' % (500 * math.ceil(1.5 * rss / 500)))
    if pair_gb:
        q = lambda f: pair_gb[min(len(pair_gb) - 1, int(f * (len(pair_gb) - 1)))]
        typ = q(0.5) + n * jf_mb / 1e3 + sandbox_gb
        worst = pair_gb[-1] + n * jf_max_mb / 1e3 + sandbox_gb
        print('  scratch         typical %.1f GB, 99%% %.1f GB, worst %.1f GB per job (one MuDst + its SimpleTree + %d jet files + sandbox)' % (
            typ, q(0.99) + n * jf_mb / 1e3 + sandbox_gb, worst, n))
        print('  disk request    %d GB (1.25 x worst, rounded up)' % int(math.ceil(1.25 * worst)))
    print('  output          %.1f GB of histograms (%d files; %.0f MB after hadd per run)' % (
        njobs * out_kb / 1e6, njobs, len(runs) * out_kb / 1e3))


if __name__ == '__main__':
    main()
