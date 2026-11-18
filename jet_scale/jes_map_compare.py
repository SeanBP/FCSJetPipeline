#!/usr/bin/env python3
"""JES map statistical noise and map-to-map comparisons (Sean 2026-11-17).

Maps are the log-cubic JSONs of CalibrationMap_em_only: per grid node (i, j)
scale F(E) = A + B l + C l^2 + D l^3, l = ln E (C, D face-wide or per node).
Everything here compares maps at their grid nodes (no interpolation), at a
few fixed jet energies, as the relative difference d = (F_a - F_b) / F_mean.

noise <manifest.txt>   (jes_subset_maps.sh output)
    For each division k, every pair of its disjoint subset maps (1/k of the
    sample each) is independent, so per node and energy
        <d^2>_pairs = 2 k sigma_1^2
    with sigma_1 the relative noise of a map on the whole sample. The
    estimate <d^2> / (2k) must not depend on k if the noise scales as
    1/sqrt(n) and the method does not change with statistics; the share of
    Gaussian-core cells per map (from its log) is printed as the second check.
    Per node the estimate is the variance of the division's k maps about their
    mean (k - 1 degrees of freedom; the same as <d^2>_pairs / 2 without the
    pairs' correlations), / k; pooled over the divisions with weights k - 1.
    Writes <out>.txt (per k and energy: RMS over nodes of sigma_1, its ratio to
    the largest-k value, core share) and <out>_sigma.csv (per node and energy:
    sigma_1 pooled, and sigma1_smooth, the mean of sigma_1^2 over the node and
    its grid neighbours, about 9x the degrees of freedom: the input of "test").
    sigma_1_sm / sigma1_sm (Sean 2026-11-18): the same after smoothing each
    subset map's relative deviation (Gaussian, --smooth grid bins, NaN-aware):
    the noise that survives the smoothing, measured, so the correlations
    between neighbouring nodes (the map's own smoothing of A, B) are included.
    "error" subtracts it from a difference smoothed with the same kernel.

test <map_a.json> <map_b.json> --sigma <out>_sigma.csv [--na N --nb N --nref N]
    Two maps from independent samples (e.g. jesmap_pixring_v1 and the new
    full map at the same settings). Predicted noise per node and energy:
        var_pred = sigma_1^2 (nref / na + nref / nb)
    (nref: events of the sample sigma_1 refers to; na, nb: events of the two
    maps' samples; default all equal). Per energy: RMS of d, RMS predicted and
    their ratio (1 if the difference is all noise: the main number), RMS pull
    d / sqrt(var_pred) with the smoothed sigma_1 (a few % above 1 from its
    remaining estimation noise), share of
    |pull| > 3, mean d (a coherent scale shift; nodes are smoothed together,
    so its error is larger than sqrt(var_pred / n_nodes)), and the excess
    sqrt(<d^2> - <var_pred>) beyond the noise. Writes <out>.txt and
    <out>_nodes.csv.

error --base <map.json> --sigma <noise>_sigma.csv --nbase N [--nref N]
      --alt NAME=<map.json>:N [--alt ...] [--combine rms|max] [--write-map <out.json>]
    The JES map's own error vs energy (Sean 2026-11-18): one standard map, and
    its uncertainty in place of a second map.
      stat(E): the standard map's statistical noise, sigma_1 x sqrt(nref /
               nbase) per node (sigma_1 from "noise" on subsets of the
               standard map's sample, nref its events; default nref = nbase),
               RMS over the nodes.
      sys(E):  per alternative-generator / tune map (a smaller sample is
               enough): its difference to the standard map, d = (F_alt -
               F_base) / F_mean, smoothed over the grid (--smooth bins, as in
               "noise"), with the noise of both maps after the same smoothing
               (sigma1_sm, measured from the smoothed subset maps) removed:
                   sys_alt^2 = <d_sm^2> - <sigma1_sm^2 (nref/nbase + nref/nalt)>
               over the nodes (0 if negative). sys_node is the same without
               smoothing (per-node d and noise), a cross-check: both remove the
               noise at the level it enters. The alternatives combine as the
               RMS (default) or the maximum.
      total(E) = stat (+) sys.
    Diagnostics per alternative: the mean signed difference (a coherent scale
    offset) and, per node, sqrt(max(0, d_sm^2 - noise_sm^2)) as a spatial mean
    and spread (the band of plot_jes_error.C). Writes <out>.txt (per energy),
    <out>_nodes.csv (per node: the heat maps) and, with --write-map, a copy of
    the standard map with the arrays jes_error_E / _stat / _sys / _total
    appended (fractions; AnHistograms --jes-error reads _E and _total for the
    _vjesup / _vjesdn copies). The map's text is kept byte for byte; the arrays
    go before its final brace.
    Note: the alternative maps should use the standard map's settings
    (CalibrationMap_em_only -F, grid); a small sample moves more cells from the
    Gaussian core to the kernel peak, which is a method change, not noise:
    compare the core shares (--alt-dir NAME=<map dir> reads its log).

Usage:
  jes_map_compare.py noise <manifest.txt> --out <prefix> [--json NAME] [--energies ...]
  jes_map_compare.py test <a.json> <b.json> --sigma <csv> --out <prefix> [--na N --nb N --nref N] [--energies ...]
  jes_map_compare.py error --base <map.json> --sigma <csv> --nbase N --alt NAME=<map.json>:N ... --out <prefix>
      [--nref N] [--combine rms|max] [--smooth 5] [--write-map <out.json>] [--energies ...]
Default energies (all subcommands; noise and error must use the same):
  3,4,5,6,8,10,12,15,20,25,30,35,40,50,60,70,80,90,100,120,140,160,180,200,225,250 GeV
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import sys

E_EVAL_MIN = 2.0  # JetParameters.h kJesEvalMinGeV: ApplyCorrections holds E at >= 2 GeV
DEFAULT_ENERGIES = '3,4,5,6,8,10,12,15,20,25,30,35,40,50,60,70,80,90,100,120,140,160,180,200,225,250'


def load_map(path):
    """{(x, y): (A, B, C, D)} at the grid nodes, plus the JSON header."""
    d = json.load(open(path))
    if d.get('model') != 'logcubic':
        sys.exit('jes_map_compare: %s is not a log-cubic map' % path)
    xg, yg = d['x_grid'], d['y_grid']
    nodes = {}
    for b in d['bins']:
        A, B = b.get('A'), b.get('B')
        if A is None or B is None or not (A == A and B == B):
            continue
        key = (round(xg[b['i']], 3), round(yg[b['j']], 3))
        nodes[key] = (A, B, b.get('C', d['C']), b.get('D', d['D']))
    return nodes, d


def scale(p, E):
    l = math.log(max(E, E_EVAL_MIN))
    return p[0] + l * (p[1] + l * (p[2] + l * p[3]))


def core_share(map_dir):
    """Gaussian-core share of the fine-grid points, from the map job's log (-k 1)."""
    for f in sorted(glob.glob(os.path.join(map_dir, 'log', '*.out'))):
        for line in open(f, errors='replace'):
            m = re.search(r'Gaussian core used for (\d+) points .*kernel peak for (\d+)', line)
            if m:
                c, p = int(m.group(1)), int(m.group(2))
                return c / float(c + p) if c + p else float('nan')
    return float('nan')


def rel_diff(pa, pb, E):
    fa, fb = scale(pa, E), scale(pb, E)
    m = 0.5 * (fa + fb)
    return (fa - fb) / m if m > 0 else float('nan')


def cmd_noise(a, energies):
    rows = []
    for line in open(a.manifest):
        if line.startswith('#') or not line.strip():
            continue
        k, g, n, sdir, mdir = line.split()[:5]
        rows.append((int(k), int(g), int(n), mdir))
    maps = {}
    for k, g, n, mdir in rows:
        p = os.path.join(mdir, a.json)
        if not os.path.exists(p):
            print('jes_map_compare: missing %s (map job not finished?), left out' % p)
            continue
        maps.setdefault(k, []).append((g, load_map(p)[0], core_share(mdir), n))
    if not maps:
        sys.exit('jes_map_compare: no maps found')
    # per k and energy: per node <d^2> over the pairs, / (2k)
    est = {}  # (k, E) -> {node: sigma_1^2 estimate}
    est_s = {}  # the same after smoothing each subset map (Sean 2026-11-18: the noise that survives the smoothing)
    for k in sorted(maps):
        ml = maps[k]
        if len(ml) < 2:
            print('jes_map_compare: k=%d has %d map(s), no pair' % (k, len(ml)))
            continue
        common = set.intersection(*[set(m[1]) for m in ml])
        xs = sorted(set(nd[0] for nd in common))
        ys = sorted(set(nd[1] for nd in common))
        for E in energies:
            e, rel = {}, [dict() for _ in ml]
            for node in common:
                f = [scale(m[1][node], E) for m in ml]
                mu = sum(f) / len(f)
                if not mu > 0:
                    continue
                var = sum((x - mu) ** 2 for x in f) / (len(f) - 1) / (mu * mu)  # one subset map's relative variance
                e[node] = var / k
                for g in range(len(ml)):
                    rel[g][node] = ((f[g] - mu) / mu, 0.0)
            est[(k, E)] = e
            if a.smooth > 0:
                sm = [smooth_grid(r, None, xs, ys, a.smooth) for r in rel]
                es = {}
                for node in e:
                    v = [smg[node][0] for smg in sm if node in smg]
                    if len(v) == len(ml):
                        m_ = sum(v) / len(v)
                        es[node] = sum((x - m_) ** 2 for x in v) / (len(v) - 1) / k
                est_s[(k, E)] = es
    ks = sorted(set(k for k, _ in est))
    if not ks:
        sys.exit('jes_map_compare: no division with two maps')
    kref = ks[-1]
    with open(a.out + '.txt', 'w') as f:
        f.write('# JES map noise from disjoint subset maps (jes_map_compare.py noise %s, map %s)\n' % (a.manifest, a.json))
        f.write('# sigma_1 = RMS over nodes of sqrt(var_k / k), var_k the variance of the division\'s k maps: the relative noise of a full-sample map,\n')
        f.write('# estimated from the 1/k subsets; flat in k if the noise scales as 1/sqrt(n). ratio = sigma_1(k) / sigma_1(k=%d)\n' % kref)
        f.write('# core share per map: Gaussian-core points / all fine-grid points (log of each map job)\n')
        for k in sorted(maps):
            f.write('# k=%d maps %s core shares %s files %s\n' % (k, [m[0] for m in maps[k]],
                    ' '.join('%.3f' % m[2] for m in maps[k]), [m[3] for m in maps[k]]))
        f.write('# sigma_1_sm: the same after smoothing each subset map (Gaussian, %.1f grid bins, NaN-aware): the noise that survives the smoothing\n' % a.smooth)
        f.write('# sub_mean / sub_sd: one subset map\'s own relative error (1/k of the sample), sqrt(var_k) per node, spatial mean and spread (plot_jes_noise.C);\n')
        f.write('#   corrected for the small number of maps per division: mean / c4(k), spread less the estimate\'s own scatter\n')
        f.write('# k E_GeV nodes dof sigma_1 ratio sigma_1_sm ratio_sm sub_mean sub_sd\n')
        rms = lambda e: math.sqrt(sum(e.values()) / len(e)) if e else float('nan')
        for E in energies:
            sref = rms(est.get((kref, E), {}))
            sref_s = rms(est_s.get((kref, E), {}))
            for k in ks:
                s1, s1s = rms(est[(k, E)]), rms(est_s.get((k, E), {}))
                # one subset map's error per node, sqrt of a k-map sample variance:
                # biased low by c4(k) (0.80 for 2 maps, 0.92 for 4, 0.97 for 8) and
                # scattered by its own estimation noise (relative variance 1 - c4^2);
                # both are removed from the spatial mean and spread
                c4 = math.sqrt(2.0 / (k - 1)) * math.exp(math.lgamma(k / 2.0) - math.lgamma((k - 1) / 2.0))
                sub = [math.sqrt(v * k) for v in est[(k, E)].values()]
                sm_ = sum(sub) / len(sub) / c4 if sub else float('nan')
                ssd2 = sum((x / c4 - sm_) ** 2 for x in sub) / len(sub) if sub else float('nan')
                ssd = math.sqrt(max(0.0, ssd2 - (1 - c4 * c4) / (c4 * c4) * sm_ * sm_)) if sub else float('nan')
                f.write('%d %g %d %d %.5f %.3f %.5f %.3f %.5f %.5f\n' % (k, E, len(est[(k, E)]), len(maps[k]) - 1, s1, s1 / sref if sref > 0 else float('nan'),
                                                                         s1s, s1s / sref_s if sref_s > 0 else float('nan'), sm_, ssd))
    # pooled per node: the divisions weighted by their degrees of freedom (number of maps - 1)
    dof = {k: len(maps[k]) - 1 for k in ks}
    with open(a.out + '_sigma.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(['x', 'y', 'E', 'sigma1', 'sigma1_smooth', 'dof', 'sigma1_sm', 'smooth_width'])
        for E in energies:
            pooled = {}
            for nd in set().union(*[set(est[(k, E)]) for k in ks]):
                kk = [k for k in ks if nd in est[(k, E)]]
                pooled[nd] = (sum(dof[k] * est[(k, E)][nd] for k in kk) / sum(dof[k] for k in kk), sum(dof[k] for k in kk))
            xs = sorted(set(nd[0] for nd in pooled))
            ys = sorted(set(nd[1] for nd in pooled))
            ix = {x: i for i, x in enumerate(xs)}
            iy = {y: i for i, y in enumerate(ys)}
            at = {(ix[nd[0]], iy[nd[1]]): v[0] for nd, v in pooled.items()}
            for nd in sorted(pooled):
                i, j = ix[nd[0]], iy[nd[1]]
                nb = [at[(i + di, j + dj)] for di in (-1, 0, 1) for dj in (-1, 0, 1) if (i + di, j + dj) in at]
                kk = [k for k in ks if nd in est_s.get((k, E), {})]
                ssm = math.sqrt(sum(dof[k] * est_s[(k, E)][nd] for k in kk) / sum(dof[k] for k in kk)) if kk else float('nan')
                w.writerow([nd[0], nd[1], E, '%.6g' % math.sqrt(pooled[nd][0]), '%.6g' % math.sqrt(sum(nb) / len(nb)), pooled[nd][1],
                            '%.6g' % ssm, a.smooth])
    print(open(a.out + '.txt').read())
    print('jes_map_compare: wrote %s.txt, %s_sigma.csv' % (a.out, a.out))


def cmd_test(a, energies):
    ma, ha = load_map(a.map_a)
    mb, hb = load_map(a.map_b)
    sig = read_sigma(a.sigma, a.raw_sigma)
    fac = a.nref / a.na + a.nref / a.nb
    common = sorted(set(ma) & set(mb))
    out = open(a.out + '.txt', 'w')
    out.write('# map a %s (sample %s, tune %s)\n# map b %s (sample %s, tune %s)\n' % (
        a.map_a, ha.get('sample', '?'), ha.get('tune', '?'), a.map_b, hb.get('sample', '?'), hb.get('tune', '?')))
    out.write('# noise %s; var_pred = sigma_1^2 x (nref/na + nref/nb) = sigma_1^2 x %.3f; %d common nodes (a %d, b %d)\n' % (
        a.sigma, fac, len(common), len(ma), len(mb)))
    out.write('# sigma_1: %s; ratio = rms_d / rms_pred (1: all noise)\n' % ('per node' if a.raw_sigma else 'smoothed over grid neighbours'))
    out.write('# E_GeV nodes rms_d rms_pred ratio rms_pull frac_|pull|>3 mean_d excess\n')
    nodes_out = open(a.out + '_nodes.csv', 'w')
    w = csv.writer(nodes_out)
    w.writerow(['x', 'y', 'E', 'Fa', 'Fb', 'd', 'sigma_pred', 'pull'])
    for E in energies:
        ds, ps, pulls = [], [], []
        for nd in common:
            s1 = sig.get((nd[0], nd[1], float(E)))
            d = rel_diff(ma[nd], mb[nd], E)
            if s1 is None or not (d == d) or not s1 > 0:
                continue
            vp = s1 * s1 * fac
            ds.append(d)
            ps.append(vp)
            pulls.append(d / math.sqrt(vp))
            w.writerow([nd[0], nd[1], E, '%.6f' % scale(ma[nd], E), '%.6f' % scale(mb[nd], E), '%.6f' % d, '%.6f' % math.sqrt(vp), '%.3f' % pulls[-1]])
        n = len(ds)
        if not n:
            out.write('%g 0\n' % E)
            continue
        md2, mp = sum(x * x for x in ds) / n, sum(ps) / n
        ex = math.sqrt(md2 - mp) if md2 > mp else 0.0
        out.write('%g %d %.5f %.5f %.3f %.3f %.4f %+.5f %.5f\n' % (
            E, n, math.sqrt(md2), math.sqrt(mp), math.sqrt(md2 / mp), math.sqrt(sum(p * p for p in pulls) / n),
            sum(1 for p in pulls if abs(p) > 3) / float(n), sum(ds) / n, ex))
    out.close()
    nodes_out.close()
    print(open(a.out + '.txt').read())
    print('jes_map_compare: wrote %s.txt, %s_nodes.csv' % (a.out, a.out))


def read_sigma(path, raw=False, column=None):
    """sigma_1 per (x, y, E): column sigma1 (raw), sigma1_smooth (neighbour-averaged
    per-node noise, default) or sigma1_sm (the noise of the smoothed maps)."""
    sig, width = {}, None
    for r in csv.DictReader(open(path)):
        col = column or ('sigma1' if raw else 'sigma1_smooth')
        if col not in r:
            sys.exit('jes_map_compare: %s has no column %s (rerun noise)' % (path, col))
        sig[(round(float(r['x']), 3), round(float(r['y']), 3), float(r['E']))] = float(r[col])
        width = float(r.get('smooth_width') or 0)
    read_sigma.width = width
    return sig


def smooth_grid(vals, wts, xs, ys, width):
    """NaN-aware Gaussian smoothing on the node grid (separable): returns
    {node: (sum w v / sum w, sum w^2 var / (sum w)^2)} from vals {node: (v, var)}.
    width in grid bins; nodes missing from vals carry no weight."""
    ix = {x: i for i, x in enumerate(xs)}
    iy = {y: i for i, y in enumerate(ys)}
    nx, ny = len(xs), len(ys)
    r = int(math.ceil(3 * width))
    ker = [math.exp(-0.5 * (k / width) ** 2) for k in range(-r, r + 1)]
    # grids: w, w*v, w^2*var with w = 1 at a node with a value
    g = [[[0.0, 0.0, 0.0] for _ in range(ny)] for _ in range(nx)]
    for nd, (v, var) in vals.items():
        g[ix[nd[0]]][iy[nd[1]]] = [1.0, v, var]

    def pass1(src, along_x, power):
        out = [[[0.0, 0.0, 0.0] for _ in range(ny)] for _ in range(nx)]
        for i in range(nx):
            for j in range(ny):
                acc = [0.0, 0.0, 0.0]
                for k in range(-r, r + 1):
                    ii, jj = (i + k, j) if along_x else (i, j + k)
                    if 0 <= ii < nx and 0 <= jj < ny:
                        c = src[ii][jj]
                        kw = ker[k + r]
                        acc[0] += kw * c[0]
                        acc[1] += kw * c[1]
                        acc[2] += kw ** power * c[2]
                out[i][j] = acc
        return out
    # weights / weighted values with the kernel, variances with the kernel squared
    sm = pass1(pass1(g, True, 2), False, 2)
    res = {}
    for nd in vals:
        a = sm[ix[nd[0]]][iy[nd[1]]]
        if a[0] > 0:
            res[nd] = (a[1] / a[0], a[2] / a[0] ** 2)
    return res


def cmd_error(a, energies):
    base, hb = load_map(a.base)
    sig = read_sigma(a.sigma, a.raw_sigma)
    sig_s = read_sigma(a.sigma, column='sigma1_sm') if a.smooth > 0 else {}
    if a.smooth > 0 and abs(read_sigma.width - a.smooth) > 1e-6:
        sys.exit('jes_map_compare error: --smooth %g but %s was made with %g (the noise must be smoothed alike)' % (a.smooth, a.sigma, read_sigma.width))
    nref = a.nref if a.nref > 0 else a.nbase
    alts = []
    for spec in a.alt:
        name, rest = spec.split('=', 1)
        path, n = rest.rsplit(':', 1)
        m, h = load_map(path)
        alts.append((name, m, h, float(n), path))
    if not alts:
        sys.exit('jes_map_compare error: give at least one --alt NAME=<map.json>:N')
    adirs = dict(x.split('=', 1) for x in a.alt_dir)
    xs = sorted(set(nd[0] for nd in base))
    ys = sorted(set(nd[1] for nd in base))
    missing = [E for E in energies if not any(k[2] == float(E) for k in sig)]
    if missing:
        sys.exit('jes_map_compare error: %s has no sigma_1 at E = %s (run noise with the same --energies)' % (a.sigma, missing))
    rows, node_rows = [], []
    for E in energies:
        sb = []  # the standard map's own noise per node
        for nd in base:
            s1 = sig.get((nd[0], nd[1], float(E)))
            if s1 is not None and s1 > 0:
                sb.append(s1 * s1 * nref / a.nbase)
        stat = math.sqrt(sum(sb) / len(sb)) if sb else float('nan')
        stat_sp = math.sqrt(sum((math.sqrt(v) - sum(math.sqrt(u) for u in sb) / len(sb)) ** 2 for v in sb) / len(sb)) if sb else float('nan')
        per = []
        for name, m, h, n, path in alts:
            vals, d2, vp, ds = {}, [], [], []
            for nd in set(base) & set(m):
                s1 = sig.get((nd[0], nd[1], float(E)))
                d = rel_diff(m[nd], base[nd], E)
                if s1 is None or not s1 > 0 or not (d == d):
                    continue
                var = s1 * s1 * (nref / a.nbase + nref / n)
                vals[nd] = (d, var)
                d2.append(d * d)
                vp.append(var)
                ds.append(d)
            if not ds:
                per.append((name,) + (float('nan'),) * 8)
                continue
            md2, mvp = sum(d2) / len(d2), sum(vp) / len(vp)
            sys_node = math.sqrt(md2 - mvp) if md2 > mvp else 0.0
            # smoothed difference (the same kernel as the subset maps in "noise") less
            # the measured noise of smoothed maps, both maps' samples
            sysv, lm, lsd, sm, vs = sys_node, float('nan'), float('nan'), {}, {}
            if a.smooth > 0:
                sm = smooth_grid(dict((nd, (v, 0.0)) for nd, (v, var) in vals.items()), None, xs, ys, a.smooth)
                ds2, vps, loc = [], [], []
                for nd, (v, _) in sm.items():
                    s1 = sig_s.get((nd[0], nd[1], float(E)))
                    if s1 is None or not s1 == s1:
                        continue
                    vs[nd] = s1 * s1 * (nref / a.nbase + nref / n)
                    ds2.append(v * v)
                    vps.append(vs[nd])
                    loc.append(math.sqrt(max(0.0, v * v - vs[nd])))
                if ds2:
                    m2, mp2 = sum(ds2) / len(ds2), sum(vps) / len(vps)
                    sysv = math.sqrt(m2 - mp2) if m2 > mp2 else 0.0
                    lm = sum(loc) / len(loc)
                    lsd = math.sqrt(sum((x - lm) ** 2 for x in loc) / len(loc))
            per.append((name, sum(ds) / len(ds), math.sqrt(md2), math.sqrt(mvp), sysv, lm, lsd, sys_node, len(ds)))
            for nd, (v, var) in vals.items():
                svv = sm.get(nd, (float('nan'), 0.0))[0]
                node_rows.append([name, nd[0], nd[1], E, '%.6f' % v, '%.6f' % math.sqrt(var), '%.6f' % svv,
                                  '%.6f' % math.sqrt(max(0.0, svv ** 2 - vs[nd])) if nd in vs else 'nan'])
        good = [p[4] for p in per if p[4] == p[4]]
        comb = (math.sqrt(sum(x * x for x in good) / len(good)) if a.combine == 'rms' else max(good)) if good else float('nan')
        tot = math.hypot(stat, comb) if comb == comb else stat
        rows.append((E, stat, stat_sp, comb, tot, per))
    names = [x[0] for x in alts]
    with open(a.out + '.txt', 'w') as f:
        f.write('# JES map error vs energy (jes_map_compare.py error): standard map %s (sample %s, tune %s, %g events)\n' % (
            a.base, hb.get('sample', '?'), hb.get('tune', '?'), a.nbase))
        f.write('# noise %s (sigma_1 for %g events); alternatives (combined as %s): %s\n' % (
            a.sigma, nref, a.combine, '; '.join('%s = %s (tune %s, %g events%s)' % (
                x[0], x[4], x[2].get('tune', '?'), x[3],
                ', core share %.3f' % core_share(adirs[x[0]]) if x[0] in adirs else '') for x in alts)))
        f.write('# all relative (fractions). stat: RMS over nodes of the standard map\'s noise; stat_sp: its spatial spread\n')
        f.write('# per alternative: mean_d (signed, alt - base), rms_d and rms_noise per node (both maps); sys: %s; '
                'loc_mean / loc_sd: per node sqrt(max(0, d_sm^2 - noise_sm^2)), spatial mean and spread; sys_node: sqrt(<d^2> - <noise^2>) '
                'without smoothing (cross-check)\n' % (
                    'the smoothed difference (%.1f grid bins) less the measured noise of maps smoothed alike (noise sigma1_sm), RMS over nodes' % a.smooth
                    if a.smooth > 0 else 'sys_node (no smoothing)'))
        f.write('# E_GeV stat stat_sp sys total ' + ' '.join('%s_mean_d %s_rms_d %s_rms_noise %s_sys %s_loc_mean %s_loc_sd %s_sys_node' % ((n,) * 7) for n in names) + '\n')
        for E, stat, stat_sp, comb, tot, per in rows:
            f.write('%g %.5f %.5f %.5f %.5f %s\n' % (E, stat, stat_sp, comb, tot, ' '.join(
                '%+.5f %.5f %.5f %.5f %.5f %.5f %.5f' % tuple(p[1:8]) for p in per)))
    with open(a.out + '_nodes.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(['alt', 'x', 'y', 'E', 'd', 'noise', 'd_smooth', 'sys_smooth'])  # sys_smooth: sqrt(max(0, d_smooth^2 - noise_sm^2))
        w.writerows(node_rows)
    print(open(a.out + '.txt').read())
    if a.write_map:
        txt = open(a.base).read().rstrip()
        if not txt.endswith('}') or '"jes_error_' in txt:
            sys.exit('jes_map_compare error: %s does not end with } or already has jes_error arrays' % a.base)
        fmt = lambda v: ', '.join('%.6g' % x for x in v)
        es = [r[0] for r in rows]
        block = (',\n  "jes_error_note": "%s",\n  "jes_error_E": [%s],\n  "jes_error_stat": [%s],\n  "jes_error_sys": [%s],\n  "jes_error_total": [%s]\n}\n' % (
            'relative error of F(E) vs jet energy (GeV): stat = this map\'s noise (subset maps, %s), sys = %s over %s of the noise-subtracted RMS difference to %s; total = stat (+) sys; jes_map_compare.py error' % (
                os.path.basename(a.sigma), a.combine, len(alts), ', '.join(names)),
            fmt(es), fmt([r[1] for r in rows]), fmt([r[3] for r in rows]), fmt([r[4] for r in rows])))
        with open(a.write_map, 'w') as f:
            f.write(txt[:-1].rstrip() + block)
        json.load(open(a.write_map))  # still valid JSON
        print('jes_map_compare: wrote %s (the standard map with its error arrays)' % a.write_map)
    print('jes_map_compare: wrote %s.txt, %s_nodes.csv' % (a.out, a.out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd')
    pn = sub.add_parser('noise')
    pn.add_argument('manifest')
    pn.add_argument('--json', default='JetEnergyScale_lookup_em_only.json', help='map file name in each map folder')
    pn.add_argument('--smooth', type=float, default=5.0, help='Gaussian width (grid bins) for sigma1_sm; error must use the same')
    pt = sub.add_parser('test')
    pt.add_argument('map_a')
    pt.add_argument('map_b')
    pt.add_argument('--sigma', required=True)
    pt.add_argument('--na', type=float, default=1.0)
    pt.add_argument('--nb', type=float, default=1.0)
    pt.add_argument('--nref', type=float, default=1.0)
    pt.add_argument('--raw-sigma', action='store_true', help='per-node sigma_1, not the neighbour-smoothed one')
    pe = sub.add_parser('error')
    pe.add_argument('--base', required=True)
    pe.add_argument('--sigma', required=True)
    pe.add_argument('--nbase', type=float, required=True)
    pe.add_argument('--nref', type=float, default=0, help='events sigma_1 refers to (default --nbase)')
    pe.add_argument('--alt', action='append', default=[], help='NAME=<map.json>:<events>')
    pe.add_argument('--alt-dir', action='append', default=[], help='NAME=<map dir> (its log: core share)')
    pe.add_argument('--combine', choices=('rms', 'max'), default='rms')
    pe.add_argument('--smooth', type=float, default=5.0, help='Gaussian width in grid bins (as in noise); 0 = per node, no smoothing')
    pe.add_argument('--raw-sigma', action='store_true')
    pe.add_argument('--write-map', default='')
    for p in (pn, pt, pe):
        p.add_argument('--out', required=True)
        p.add_argument('--energies', default=DEFAULT_ENERGIES)
    a = ap.parse_args()
    if not a.cmd:
        ap.print_help()
        return 1
    energies = [float(e) for e in a.energies.split(',')]
    return {'noise': cmd_noise, 'test': cmd_test, 'error': cmd_error}[a.cmd](a, energies)


if __name__ == '__main__':
    sys.exit(main())
