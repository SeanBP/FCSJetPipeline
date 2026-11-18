#!/usr/bin/env python3
"""Nested-cut systematic from two AnAsymmetry tables (2026-11-17).

The loose sample contains the tight one (for example AnAsymmetry -f 1.15 and
-f 1.3: trigger value >= 1.15 vs >= 1.3 x threshold), every other cut being
the same. Per x_F bin of each (level, class, beam) row:

  delta    = A_loose - A_tight          (the shift of the central value)
  sig_diff = sqrt(|sig_tight^2 - sig_loose^2|)

For a loose sample = slice + tight with inverse-variance weights, sig_diff
equals w sqrt(sig_slice^2 + sig_tight^2) and delta = w (A_slice - A_tight),
w = 1 - sig_loose^2 / sig_tight^2 (the slice's weight in the loose sample).
The slice's own A_N and error are also written:
  A_slice = (A_loose - (1 - w) A_tight) / w,  sig_slice = sig_loose / sqrt(w)

Barlow test: sys = sqrt(delta^2 - sig_diff^2) where |delta| > sig_diff, else 0
(the shift is consistent with a statistical fluctuation). Per row also the
chi^2 = sum (delta / sig_diff)^2 over the bins, its ndf and p-value.

--disjoint: the two samples share no events (e.g. runs early / late in the
fill, time_in_fill.py): delta = A_1 - A_2, sig_diff = sqrt(sig_1^2 + sig_2^2),
no slice; a consistency check (chi^2), sys still printed for reference.

--shifted: the same events with different jet kinematics (the JES systematic,
AnAsymmetry -v jesup / jesdn: jets move between x_F bins and across the cuts), neither
nested nor independent: delta = A_first - A_second and sys = |delta| (the full
shift, conservative: its statistical part from the migrating jets is not
subtracted, the two samples being too correlated for either formula above);
sig_diff, w, pull and chi^2 are not defined (written as -1 / nan).

Bins are matched on (level, class, beam, x_F low edge, high edge), so both
tables must use the same x_F binning (AnAsymmetry rebin / -x).

Usage:
  sys_slices.py <loose.txt> <tight.txt> <out prefix> [--label NAME]
                [--min-jets N] [--rows level:class[,level:class...]] [--disjoint | --shifted]
Writes <out prefix>.txt (per bin) and <out prefix>_summary.txt (per row);
plot_sys_slices.C draws them.
"""
import argparse
import math
import sys


def read_table(path):
    """AnAsymmetry .txt: '# level class beam xF_lo xF_hi <xF> A_N stat jets <|cos|> <P>',
    or a pi0_bgsub.C table: '# beam xF_lo xF_hi <xF> A_pi0 stat f A_sig stat_sig A_sb stat_sb n_sig n_sb ...'
    (read as level "all", class "pi0"; jets = n_sig, the signal-window candidates)."""
    rows, header = {}, []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
                continue
            t = line.split()
            if len(t) >= 13 and t[0] in ('blue', 'yellow'):
                key = ('all', 'pi0', t[0], round(float(t[1]), 4), round(float(t[2]), 4))
                rows[key] = dict(xf=float(t[3]), an=float(t[4]), err=float(t[5]), jets=int(float(t[11])))
                continue
            if len(t) < 9:
                continue
            key = (t[0], t[1], t[2], round(float(t[3]), 4), round(float(t[4]), 4))
            rows[key] = dict(xf=float(t[5]), an=float(t[6]), err=float(t[7]), jets=int(float(t[8])))
    return rows, header


def chi2_pvalue(chi2, ndf):
    """Upper-tail p-value of chi^2 with ndf degrees of freedom (regularized gamma)."""
    if ndf <= 0:
        return float('nan')
    a, x = ndf / 2.0, chi2 / 2.0
    if x <= 0:
        return 1.0
    if x < a + 1:  # series for the lower incomplete gamma
        s, term, n = 1.0 / a, 1.0 / a, a
        while abs(term) > 1e-14 * abs(s):
            n += 1
            term *= x / n
            s += term
        return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - math.lgamma(a)))
    # continued fraction for the upper one
    b, c, d = x + 1 - a, 1e300, 1.0 / (x + 1 - a)
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-14:
            break
    return min(1.0, math.exp(-x + a * math.log(x) - math.lgamma(a)) * h)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('loose')
    ap.add_argument('tight')
    ap.add_argument('out')
    ap.add_argument('--label', default='', help='name of the variation, written in the headers')
    ap.add_argument('--min-jets', type=int, default=100, help='skip bins with fewer tight-sample jets (default 100)')
    ap.add_argument('--rows', default='', help='level:class list to keep (default: all rows)')
    ap.add_argument('--disjoint', action='store_true', help='independent samples: sig_diff = sqrt(sig_1^2 + sig_2^2)')
    ap.add_argument('--shifted', action='store_true', help='same events, shifted kinematics (JES): sys = |delta|')
    a = ap.parse_args()
    if a.disjoint and a.shifted:
        sys.exit('sys_slices: --disjoint and --shifted exclude each other')
    kind_name = 'disjoint-sample check' if a.disjoint else 'shifted-kinematics systematic' if a.shifted else 'nested-cut systematic'

    loose, hl = read_table(a.loose)
    tight, ht = read_table(a.tight)
    keep = set(tuple(r.split(':')) for r in a.rows.split(',') if r)

    out = open(a.out + '.txt', 'w')
    out.write('# %s%s: %s %s, %s %s\n' % (kind_name, ' (' + a.label + ')' if a.label else '',
              'first' if a.disjoint or a.shifted else 'loose', a.loose, 'second' if a.disjoint or a.shifted else 'tight', a.tight))
    for h in hl[:3]:
        out.write('# loose ' + h.lstrip('# ') + '\n')
    for h in ht[:3]:
        out.write('# tight ' + h.lstrip('# ') + '\n')
    out.write('# delta = A_loose - A_tight; sig_diff = sqrt(|sig_t^2 - sig_l^2|); w = slice weight; '
              'sys = sqrt(delta^2 - sig_diff^2) if |delta| > sig_diff else 0\n')
    out.write('# level class beam xF_lo xF_hi <xF>_tight A_loose stat_loose A_tight stat_tight jets_loose jets_tight '
              'w A_slice stat_slice delta sig_diff pull sys flag\n')

    summ = {}
    nflag = nbins = 0
    for key in sorted(tight):
        lvl, cls, beam, lo, hi = key
        if keep and (lvl, cls) not in keep:
            continue
        if key not in loose:
            continue
        L, T = loose[key], tight[key]
        if T['jets'] < a.min_jets or not (T['err'] > 0 and L['err'] > 0):
            continue
        delta = L['an'] - T['an']
        if a.shifted:
            d2, flag = float('nan'), 'shifted'
            w = float('nan')
        elif a.disjoint:
            d2, flag = T['err'] ** 2 + L['err'] ** 2, 'disjoint'
            w = float('nan')
        else:
            d2 = T['err'] ** 2 - L['err'] ** 2
            # loose must be the larger sample: its error is the smaller one
            if L['jets'] == T['jets'] and abs(d2) < 1e-12:
                flag = 'identical'  # no slice jets in this bin: delta = sig_diff = 0
            else:
                flag = 'ok' if d2 > 0 and L['jets'] >= T['jets'] else 'NOT_NESTED'
            if flag == 'NOT_NESTED':
                nflag += 1
            w = 1 - L['err'] ** 2 / T['err'] ** 2 if d2 > 0 else float('nan')
        sd = math.sqrt(abs(d2)) if d2 == d2 else float('nan')
        if w == w and w > 1e-6:
            a_s = (L['an'] - (1 - w) * T['an']) / w
            e_s = L['err'] / math.sqrt(w)
        else:
            a_s = e_s = float('nan')
        pull = delta / sd if sd > 0 else float('nan')
        if a.shifted:
            sysv = abs(delta)
        else:
            sysv = math.sqrt(delta ** 2 - sd ** 2) if sd > 0 and abs(delta) > sd else 0.0
        out.write('%s %s %s %.3f %.3f %.4f %+.6f %.6f %+.6f %.6f %d %d %.4f %+.6f %.6f %+.6f %.6f %+.3f %.6f %s\n' % (
            lvl, cls, beam, lo, hi, T['xf'], L['an'], L['err'], T['an'], T['err'], L['jets'], T['jets'],
            w if w == w else -1, a_s if a_s == a_s else -1, e_s if e_s == e_s else -1, delta, sd if sd == sd else -1, pull, sysv, flag))
        nbins += 1
        s = summ.setdefault((lvl, cls, beam), dict(chi2=0.0, ndf=0, nsig=0, maxsys=0.0, maxrel=0.0))
        if sd > 0:
            s['chi2'] += pull ** 2
            s['ndf'] += 1
        if sysv > 0:
            s['nsig'] += 1
            s['maxsys'] = max(s['maxsys'], sysv)
            s['maxrel'] = max(s['maxrel'], sysv / T['err'])
    out.close()

    with open(a.out + '_summary.txt', 'w') as f:
        kind = ('disjoint-sample consistency check (sys columns for reference only)' if a.disjoint else
                'shifted-kinematics systematic (sys = |delta|; no chi^2)' if a.shifted else 'nested-cut systematic')
        f.write('# %s%s, per row: chi2 = sum (delta/sig_diff)^2 over x_F bins, p-value; '
                'bins with |delta| > sig_diff (expected ~32%% of ndf by chance); largest sys, and largest sys / stat of the second table\n'
                % (kind, ' (' + a.label + ')' if a.label else ''))
        f.write('# level class beam ndf chi2 chi2/ndf p n_sig max_sys max_sys/stat\n')
        for k in sorted(summ):
            s = summ[k]
            p = chi2_pvalue(s['chi2'], s['ndf'])
            f.write('%s %s %s %d %.2f %.2f %.3f %d %.6f %.2f\n' % (k + (s['ndf'], s['chi2'], s['chi2'] / s['ndf'] if s['ndf'] else float('nan'),
                                                                  p, s['nsig'], s['maxsys'], s['maxrel'])))
    print('sys_slices: %d rows, %d bins written to %s.txt; summary %s_summary.txt%s' % (
        len(summ), nbins, a.out, a.out,
        '; WARNING: %d bins not nested (loose error >= tight error)' % nflag if nflag else ''))


if __name__ == '__main__':
    sys.exit(main())
