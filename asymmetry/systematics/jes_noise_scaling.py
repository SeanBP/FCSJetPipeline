#!/usr/bin/env python3
"""Statistical noise in the JES-shifted A_N from data subsets (Sean 2026-11-17).

groups <groups.txt>   (jes_data_subsets.sh: 'k nominal_table jes_table' per
    data group, k = 1 the full sample)
    Delta_g = A_jes - A_nom per bin of each group (1/k of the data). Fit per
    bin, over all groups,
        Delta_g^2 = delta^2 + k s^2
    (the per-k means, weighted by n_k / <Delta^2>_k^2, their Gaussian variance): delta is the JES effect, s the noise of the
    full-sample shift. delta^2 can come out negative (all noise); then
    delta = 0. Also pooled over the bins of a row, in units of each bin's
    full-sample statistical error (rows fitted separately).

Bins are matched on (level, class, beam, x_F low, high); --rows keeps
level:class rows (default all:all, the union). --min-jets drops bins with
fewer jets in a group's table.

Usage:
  jes_noise_scaling.py groups <groups.txt> --out <prefix> [--rows all:all] [--min-jets 20]
"""
import argparse
import math
import sys


def read_table(path):
    rows = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        t = line.split()
        if len(t) < 9:
            continue
        rows[(t[0], t[1], t[2], round(float(t[3]), 4), round(float(t[4]), 4))] = (float(t[6]), float(t[7]), float(t[8]))
    return rows


def lsq(xs_list, ys, ws=None):
    """Weighted least squares y = sum_i c_i x_i (normal equations, up to 2 terms)."""
    ws = ws or [1.0] * len(ys)
    n = len(xs_list[0])
    if n == 1:
        sxx = sum(w * x[0] * x[0] for x, w in zip(xs_list, ws))
        return [sum(w * x[0] * y for x, y, w in zip(xs_list, ys, ws)) / sxx] if sxx > 0 else [float('nan')]
    a = [[sum(w * x[i] * x[j] for x, w in zip(xs_list, ws)) for j in range(2)] for i in range(2)]
    b = [sum(w * x[i] * y for x, y, w in zip(xs_list, ys, ws)) for i in range(2)]
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    if abs(det) < 1e-300:
        return [float('nan'), float('nan')]
    return [(b[0] * a[1][1] - b[1] * a[0][1]) / det, (a[0][0] * b[1] - a[1][0] * b[0]) / det]


def fit_means(pts, terms):
    """Fit the per-k means of y (Delta^2 or d^2) with weights n_k / mean_k^2: a
    squared Gaussian variable has variance 2 <y>^2, so the noisy large-k points
    must not pull the intercept. terms(k) -> regressors. Returns the coefficients."""
    by = {}
    for k, y in pts:
        by.setdefault(k, []).append(y)
    ks = sorted(by)
    means = [sum(by[k]) / len(by[k]) for k in ks]
    # weights from a first unweighted fit (the model, not the noisy means)
    c0 = lsq([terms(k) for k in ks], means)
    model = [sum(ci * ti for ci, ti in zip(c0, terms(k))) for k in ks]
    ws = [len(by[k]) / max(m, 1e-3 * max(means) if max(means) > 0 else 1e-30) ** 2 for k, m in zip(ks, model)]
    return lsq([terms(k) for k in ks], means, ws)


def keep_row(key, rows):
    return (key[0], key[1]) in rows


def cmd_groups(a, rows):
    data = {}  # bin -> [(k, Delta^2)]
    full_stat = {}
    for line in open(a.list):
        if line.startswith('#') or not line.strip():
            continue
        k, pn, pj = line.split()[:3]
        k = int(k)
        tn, tj = read_table(pn), read_table(pj)
        for key in tn:
            if not keep_row(key, rows) or key not in tj:
                continue
            if tn[key][2] < a.min_jets or tj[key][2] < a.min_jets:
                continue
            d = tj[key][0] - tn[key][0]
            data.setdefault(key, []).append((k, d * d))
            if k == 1:
                full_stat[key] = tn[key][1]
    out = open(a.out + '.txt', 'w')
    out.write('# JES shift noise from data subsets (jes_noise_scaling.py groups %s): Delta^2 = delta^2 + k s^2 per bin\n' % a.list)
    out.write('# stat: full-sample statistical error of A_N (nominal); delta = sqrt(max(delta^2, 0)); Delta_full = the full-sample shift\n')
    out.write('# level class beam xF_lo xF_hi groups k_values Delta_full delta2 delta s stat delta/stat s/stat\n')
    pooled = {}
    for key in sorted(data):
        v = data[key]
        ks = sorted(set(k for k, _ in v))
        if len(ks) < 2:
            continue
        d2, s2 = fit_means(v, lambda k: (1.0, float(k)))
        dfull = [math.sqrt(y) for k, y in v if k == 1]
        st = full_stat.get(key, float('nan'))
        delta, s = math.sqrt(max(d2, 0.0)), math.sqrt(max(s2, 0.0))
        out.write('%s %s %s %.3f %.3f %d %s %.6f %+.3e %.6f %.6f %.6f %.3f %.3f\n' % (
            key[0], key[1], key[2], key[3], key[4], len(v), ','.join(str(k) for k in ks), dfull[0] if dfull else float('nan'),
            d2, delta, s, st, delta / st if st > 0 else float('nan'), s / st if st > 0 else float('nan')))
        if st > 0:
            r = pooled.setdefault((key[0], key[1], key[2]), [])
            r += [(k, y / (st * st)) for k, y in v]
    out.write('# pooled per row, in units of the full-sample stat error squared: Delta^2/stat^2 = D + k S\n')
    out.write('# level class beam points D S sqrt(D) sqrt(S)\n')
    for r in sorted(pooled):
        v = pooled[r]
        D, S = fit_means(v, lambda k: (1.0, float(k)))
        out.write('# %s %s %s %d %+.4f %+.4f %.3f %.3f\n' % (r + (len(v), D, S, math.sqrt(max(D, 0)), math.sqrt(max(S, 0)))))
    out.close()
    print(open(a.out + '.txt').read())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd')
    p = sub.add_parser('groups')
    p.add_argument('list')
    p.add_argument('--out', required=True)
    p.add_argument('--rows', default='all:all')
    p.add_argument('--min-jets', type=int, default=20)
    a = ap.parse_args()
    if not a.cmd:
        ap.print_help()
        return 1
    rows = set(tuple(r.split(':')) for r in a.rows.split(','))
    return cmd_groups(a, rows)


if __name__ == '__main__':
    sys.exit(main())
