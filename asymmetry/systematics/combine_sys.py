#!/usr/bin/env python3
"""Combine step: one per-bin table with the statistical error, every systematic
component and their quadrature total (Sean 2026-11-17).

Inputs
  --nominal TABLE   the nominal A_N: an AnAsymmetry .txt
                    ('level class beam xF_lo xF_hi <xF> A_N stat jets ...') or a
                    pi0_bgsub.C table ('beam xF_lo xF_hi <xF> A_pi0 stat f ...
                    df sys_f', read as level "all", class "pi0")
  --comp NAME=FILE  a systematic component from a sys_slices.py output (the
                    per-bin Barlow "sys" column; bins it lacks, e.g. below its
                    --min-jets, count 0 and are listed in the header). Repeat
                    per component, e.g. trigger=sys_trig.txt fiducial=sys_fid25_jets.txt
                    vertex=sys_vz60_jets.txt novertex=sys_vtxmeas_jets.txt
                    isolation=pi0_iso_sys/sys_iso_tight.txt
  --comp-max NAME=FILE1,FILE2,...
                    one component as the per-bin maximum of several variations
                    (e.g. fiducial=sys_fid25_jets.txt,sys_fid30_jets.txt)
  --sysf NAME       (pi0 tables) take the nominal table's own sys_f column (the
                    signal-fraction fit, pi0_bgsub.C) as component NAME
  --rel NAME=X      a flat relative component, X x |A_N| (e.g. a scale)
  --xf-rel X        horizontal (x_F) uncertainty, X x |<x_F>| (FMS: flat 3%,
                    default 0.03); written as its own column, not in the total
The polarization scale is not combined: quoted separately (as FMS).

Output <prefix>.txt, per bin:
  level class beam xF_lo xF_hi <xF> A_N stat sys_total stat+sys xF_err sys_<comp> ...
where sys_total = quadrature sum of the components, stat+sys = quadrature sum
of stat and sys_total. plot_an_final.C draws it (stat only / sys only, total
and each component / combined).

Usage:
  combine_sys.py --nominal an_nominal.txt --comp trigger=... [--comp ...]
                 [--comp-max NAME=F1,F2] [--sysf NAME] [--rel NAME=X]
                 [--xf-rel 0.03] [--rows level:class[,...]] --out final_jets
"""
import argparse
import math
import sys


def key_of(t, pi0):
    if pi0:
        return ('all', 'pi0', t[0], round(float(t[1]), 4), round(float(t[2]), 4))
    return (t[0], t[1], t[2], round(float(t[3]), 4), round(float(t[4]), 4))


def read_nominal(path):
    rows = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        t = line.split()
        if len(t) >= 15 and t[0] in ('blue', 'yellow'):  # pi0_bgsub.C
            rows[key_of(t, True)] = dict(xf=float(t[3]), an=float(t[4]), stat=float(t[5]), sysf=float(t[14]))
        elif len(t) >= 9 and t[2] in ('blue', 'yellow'):
            rows[key_of(t, False)] = dict(xf=float(t[5]), an=float(t[6]), stat=float(t[7]), sysf=float('nan'))
    return rows


def read_sys(path):
    """sys_slices.py per-bin output: the Barlow sys column (index 18)."""
    d = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        t = line.split()
        if len(t) < 20:
            continue
        d[key_of(t, False)] = float(t[18])
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--nominal', required=True)
    ap.add_argument('--comp', action='append', default=[])
    ap.add_argument('--comp-max', action='append', default=[])
    ap.add_argument('--sysf', default='')
    ap.add_argument('--rel', action='append', default=[])
    ap.add_argument('--xf-rel', type=float, default=0.03)
    ap.add_argument('--rows', default='')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    nom = read_nominal(a.nominal)
    if not nom:
        sys.exit('combine_sys: no rows in %s' % a.nominal)
    keep = set(tuple(r.split(':')) for r in a.rows.split(',')) if a.rows else None
    comps = []  # (name, dict key -> sys, source)
    for c in a.comp:
        name, path = c.split('=', 1)
        comps.append((name, read_sys(path), path))
    for c in a.comp_max:
        name, paths = c.split('=', 1)
        tabs = [read_sys(p) for p in paths.split(',')]
        merged = {}
        for tb in tabs:
            for k, v in tb.items():
                merged[k] = max(merged.get(k, 0.0), v)
        comps.append((name, merged, 'max of ' + paths))
    rels = [(r.split('=')[0], float(r.split('=')[1])) for r in a.rel]
    names = [c[0] for c in comps] + ([a.sysf] if a.sysf else []) + [r[0] for r in rels]
    if len(set(names)) != len(names):
        sys.exit('combine_sys: component names must differ')

    missing = {c[0]: 0 for c in comps}
    out = open(a.out + '.txt', 'w')
    out.write('# combined A_N table (combine_sys.py): nominal %s\n' % a.nominal)
    for c in comps:
        out.write('# component %s: %s\n' % (c[0], c[2]))
    if a.sysf:
        out.write('# component %s: the nominal pi0 table\'s sys_f (signal-fraction fit)\n' % a.sysf)
    for r in rels:
        out.write('# component %s: flat %.4f x |A_N|\n' % r)
    out.write('# xF_err: %.3f x <x_F> (horizontal; not in the totals). Polarization scale not included (quoted separately).\n' % a.xf_rel)
    out.write('# level class beam xF_lo xF_hi <xF> A_N stat sys_total stat+sys xF_err ' + ' '.join('sys_' + n for n in names) + '\n')
    nrow = 0
    for key in sorted(nom):
        if keep is not None and (key[0], key[1]) not in keep:
            continue
        r = nom[key]
        vals = []
        for name, tab, _ in comps:
            if key not in tab:
                missing[name] += 1
            vals.append(tab.get(key, 0.0))
        if a.sysf:
            vals.append(r['sysf'] if r['sysf'] == r['sysf'] else 0.0)
        for _, x in rels:
            vals.append(x * abs(r['an']))
        tot = math.sqrt(sum(v * v for v in vals))
        out.write('%s %s %s %.3f %.3f %.4f %+.6f %.6f %.6f %.6f %.5f %s\n' % (
            key[0], key[1], key[2], key[3], key[4], r['xf'], r['an'], r['stat'], tot, math.hypot(r['stat'], tot),
            a.xf_rel * abs(r['xf']), ' '.join('%.6f' % v for v in vals)))
        nrow += 1
    out.close()
    print('combine_sys: %d bins, components %s -> %s.txt' % (nrow, ', '.join(names) if names else '(none)', a.out))
    for name, n in missing.items():
        if n:
            print('combine_sys: component %s has no value in %d bins (counted 0; below its --min-jets or absent)' % (name, n))
    return 0


if __name__ == '__main__':
    sys.exit(main())
