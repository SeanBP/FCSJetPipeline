#!/usr/bin/env python3
"""Photon-class misidentification systematic for the merged classes (Sean
2026-11-17: reported for np12 = 1-2 photons and np3plus = >= 3).

The simulated response (misid_response.cpp, <prefix>_merged.txt) gives per
x_F bin the fraction of each reco merged class whose true class is the other
one. The measured A_N of a class is diluted by the other class's A_N:
  dA_np12    = f(true np3plus | reco np12) x (A_np3plus - A_np12)
  dA_np3plus = f(true np12 | reco np3plus) x (A_np12 - A_np3plus)
with the measured A_N of both (same level, beam and x_F bin of the nominal
AnAsymmetry table). The measured difference carries statistical noise
(stat_diff = the two errors in quadrature), so, as the Barlow test of the cut
variations, only the part beyond it counts:
  sys = f_other x sqrt(max(0, (A_other - A_N)^2 - stat_diff^2))
(dA, the plain estimate, is written too). The response bin is the one holding
the point's |<x_F>| (the response is binned in truth x_F, both beams alike).

Output <prefix>.txt: 'level class beam xF_lo xF_hi <xF> f_other A_N A_other
stat_diff dA sys' per bin; combine_sys.py reads it with --comp-simple.

Usage: misid_sys.py <AnAsymmetry nominal .txt> <response _merged.txt> <out prefix> [--levels all,allJP]
"""
import argparse
import math
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('nominal')
    ap.add_argument('response')
    ap.add_argument('out')
    ap.add_argument('--levels', default='all')
    a = ap.parse_args()
    levels = set(a.levels.split(','))

    resp = []  # (lo, hi, reco class, f_true_np12, f_true_np3plus)
    for line in open(a.response):
        if line.startswith('#'):
            continue
        t = line.split()
        if len(t) >= 8:
            resp.append((float(t[0]), float(t[1]), t[2], float(t[5]), float(t[6])))
    if not resp:
        sys.exit('misid_sys: no rows in %s' % a.response)

    def f_other(cls, axf):
        for lo, hi, rc, f12, f3 in resp:
            if rc == cls and lo <= axf < hi:
                return f3 if cls == 'np12' else f12
        return None

    an = {}
    for line in open(a.nominal):
        if line.startswith('#'):
            continue
        t = line.split()
        if len(t) < 9 or t[1] not in ('np12', 'np3plus') or t[0] not in levels:
            continue
        an[(t[0], t[1], t[2], round(float(t[3]), 4), round(float(t[4]), 4))] = (float(t[5]), float(t[6]), float(t[7]))

    out = open(a.out + '.txt', 'w')
    out.write('# photon-class misidentification (misid_sys.py): nominal %s, response %s\n' % (a.nominal, a.response))
    out.write('# dA = f_other x (A_other - A_N); sys = f_other x sqrt(max(0, (A_other - A_N)^2 - stat_diff^2)) (Barlow)\n')
    out.write('# level class beam xF_lo xF_hi <xF> f_other A_N A_other stat_diff dA sys\n')
    n = 0
    for k in sorted(an):
        lvl, cls, beam, lo, hi = k
        other = 'np3plus' if cls == 'np12' else 'np12'
        ko = (lvl, other, beam, lo, hi)
        if ko not in an:
            continue
        xf, A, s = an[k]
        _, Ao, so = an[ko]
        f = f_other(cls, abs(xf))
        if f is None:
            continue
        d = f * (Ao - A)
        sd = math.hypot(s, so)
        sysv = f * math.sqrt(max(0.0, (Ao - A) ** 2 - sd ** 2))
        out.write('%s %s %s %.3f %.3f %.4f %.4f %+.6f %+.6f %.6f %+.6f %.6f\n' % (lvl, cls, beam, lo, hi, xf, f, A, Ao, sd, d, sysv))
        n += 1
    out.close()
    print('misid_sys: %d bins -> %s.txt' % (n, a.out))
    return 0


if __name__ == '__main__':
    sys.exit(main())
