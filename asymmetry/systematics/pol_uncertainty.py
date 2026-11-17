#!/usr/bin/env python3
"""Polarization uncertainty of the A_N data set, per beam (2026-11-17).

Follows the RHIC polarimetry group's recipe: Z. Chang, "Example calculation
of fill-to-fill polarization uncertainties" (2018, wiki.bnl.gov/rhicspin,
ExampleFillToFill.pdf) and W. B. Schmidke, BNL-209057-2018-TECH Sec. 4.3:

  P_run  = P0 + dP/dt (t_mid - t0)                          (run at its middle)
  P_fill = sum_run w P_run / sum_run w,   t_fill = sum w t_mid / sum w
  s^2(P_fill) = s^2(P0) + s^2(dP/dt) (t_fill - t0)^2 + (f2f x P_fill)^2
  P_set  = sum_fill W P_fill / sum W,     W = sum of the fill's run weights
  s(f2f) = sqrt(1 - M/N) x sqrt(sum W^2 s^2(P_fill)) / sum W   (0 if M >= N)
  s(profile) = P_set x profile / sqrt(M)
  s(P_set) = P_set x scale (+) s(f2f) (+) s(profile)

w = the run's A_N jet count (our A_N weights P by jets, so this replaces the
luminosity of the recipe). M = fills in the data set, N = fills the group
used for the analyzing power. scale, f2f (relative fill-to-fill systematic),
N and profile are per beam and per running period; Run 22 values are not on
the polarimetry pages (2026-11-17), so the defaults are the Run 17 (255 GeV)
ones of BNL-209057-2018-TECH Tables 6-8 and Sec. 4.3.2 and are PLACEHOLDERS:
scale 1.5 / 1.5 % (conservative placeholder, Sean 2026-11-17; Run 17 prelim 1.1 %), f2f 0 / 2.6 %, N 190 / 191, profile 2.2 % per fill.

Usage:
  pol_uncertainty.py <run_jets.txt> [--runs ../run22_polarization.txt]
      [--fills <fill table>] [--exclude <AnMerge -x log or run list>]
      [--scale-blue 1.5] [--scale-yellow 1.5] [--f2f-blue 0] [--f2f-yellow 2.6]
      [--n-blue 190] [--n-yellow 191] [--profile 2.2] [--out pol_unc.txt]
run_jets.txt: "run jets" lines (dump_run_jets.C on an AnMerge output).
"""
import argparse
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FILL_TABLE = '/star/u/dkap7827/Fcs2019/FcsAsymRun22/Run22PolFromWeb_20250205.txt'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('weights')
    ap.add_argument('--runs', default=os.path.join(HERE, '..', 'run22_polarization.txt'))
    ap.add_argument('--fills', default=FILL_TABLE)
    ap.add_argument('--exclude', default='', help='AnMerge -x log ("bad run: <run>") or a list of runs')
    ap.add_argument('--scale-blue', type=float, default=1.5)
    ap.add_argument('--scale-yellow', type=float, default=1.5)
    ap.add_argument('--f2f-blue', type=float, default=0.0)
    ap.add_argument('--f2f-yellow', type=float, default=2.6)
    ap.add_argument('--n-blue', type=int, default=190)
    ap.add_argument('--n-yellow', type=int, default=191)
    ap.add_argument('--profile', type=float, default=2.2, help='relative profile uncertainty per fill, %%')
    ap.add_argument('--out', default='pol_unc.txt')
    a = ap.parse_args()

    # fill table: fill E start stop  P0 +- dP0 dPdt +- ddPdt (blue)  P0 +- dP0 dPdt +- ddPdt (yellow); %, %/h
    fills = {}
    for line in open(a.fills):
        t = line.split()
        if len(t) < 16 or not t[0].isdigit():
            continue
        try:
            beams = []
            for o in (4, 10):
                beams.append((float(t[o]), float(t[o + 2]), float(t[o + 3]), float(t[o + 5])))
            fills[int(t[0])] = (int(t[2]), beams)
        except ValueError:
            continue
    runinfo = {}
    for line in open(a.runs):
        if line.startswith('#'):
            continue
        t = line.split()
        runinfo[int(t[0])] = (int(t[1]), float(t[2]))
    excl = set()
    if a.exclude:
        for line in open(a.exclude):
            m = re.search(r'bad run:\s*(\d+)', line) or re.match(r'\s*(\d{8})\b', line)
            if m:
                excl.add(int(m.group(1)))
    w = {}
    for line in open(a.weights):
        if line.startswith('#'):
            continue
        t = line.split()
        if len(t) >= 2 and float(t[1]) > 0:
            w[int(t[0])] = float(t[1])

    # group the weighted runs by fill
    byfill, nskip, wskip = {}, 0, 0.0
    for run, wt in w.items():
        if run in excl:
            continue
        if run not in runinfo or runinfo[run][0] not in fills:
            nskip += 1
            wskip += wt
            continue
        fill, tmid = runinfo[run]
        byfill.setdefault(fill, []).append((tmid, wt))

    out = open(a.out, 'w')
    out.write('# polarization uncertainty of the data set (pol_uncertainty.py): weights %s, %d runs excluded, %d runs without P (%.2f%% of the weight)\n'
              % (a.weights, len(excl), nskip, 100 * wskip / max(1e-9, sum(w.values()))))
    out.write('# inputs per beam: scale %%, fill-to-fill sys %%, N fills, profile %% per fill (defaults = Run 17 PLACEHOLDERS)\n')
    out.write('# beam M N P_set%% stat_f2f%% (rel) f2f_corr%% (rel) profile%% (rel) scale%% (rel) total%% (rel)\n')
    params = {0: ('blue', a.scale_blue, a.f2f_blue, a.n_blue), 1: ('yellow', a.scale_yellow, a.f2f_yellow, a.n_yellow)}
    for b in (0, 1):
        name, scale, f2f, nfill = params[b]
        sw = swp = sw2s2 = 0.0
        m = 0
        for fill, rows in byfill.items():
            t0, beams = fills[fill]
            p0, dp0, slope, dslope = beams[b]
            if p0 <= 0:
                continue
            W = sum(r[1] for r in rows)
            tbar = sum(r[0] * r[1] for r in rows) / W
            dt_h = (tbar - t0) / 3600.0
            pf = p0 + slope * dt_h  # = the weighted mean of the runs' P (linear in t)
            s2 = dp0 ** 2 + (dslope * dt_h) ** 2 + (f2f / 100.0 * pf) ** 2
            sw += W
            swp += W * pf
            sw2s2 += W * W * s2
            m += 1
        if sw <= 0:
            continue
        pset = swp / sw
        s_stat = math.sqrt(sw2s2) / sw
        corr = math.sqrt(1 - m / nfill) if m < nfill else 0.0
        s_f2f = corr * s_stat
        s_prof = pset * a.profile / 100.0 / math.sqrt(m)
        s_scale = pset * scale / 100.0
        tot = math.sqrt(s_scale ** 2 + s_f2f ** 2 + s_prof ** 2)
        line = '%s %d %d %.3f %.3f %.3f %.3f %.3f %.3f' % (name, m, nfill, pset, 100 * s_stat / pset, 100 * s_f2f / pset,
                                                              100 * s_prof / pset, 100 * s_scale / pset, 100 * tot / pset)
        out.write(line + '\n')
        print('pol_uncertainty: %-6s M = %d fills (N = %d): P_set = %.2f%%; relative: fill-to-fill %.2f%% (uncorrected %.2f%%), '
              'profile %.2f%%, scale %.2f%% -> total %.2f%%' % (name, m, nfill, pset, 100 * s_f2f / pset, 100 * s_stat / pset,
                                                                100 * s_prof / pset, 100 * s_scale / pset, 100 * tot / pset))
    out.close()
    if nskip:
        print('pol_uncertainty: %d weighted runs have no polarization (%.2f%% of the weight); left out' % (nskip, 100 * wskip / max(1e-9, sum(w.values()))))
    return 0


if __name__ == '__main__':
    sys.exit(main())
