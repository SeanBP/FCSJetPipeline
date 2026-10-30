#!/usr/bin/env python3
"""Compare ECal gain corrections tower by tower: an anchor file from X. Liang's
pi0 calibration against the conditions DB (dump_ecal_gaincorr.C output).

Either anchor format is accepted and recognised by its column count:
  fcsgaincorr_<day>_<N>.txt         "ehp ns dep ch gaincorr" (mapped to tower ids
                                    with liang_depch_id.txt, "ehp ns dep ch id")
  fcsgaincorrOffline_<day>_<N>.txt  "ns id gaincorr"
Prints the median anchor / DB ratio, its 16-84% spread, the share of towers
within 1% and 5%, and the largest deviations.

Usage: ./compare_gaincorr.py <anchor file> [<DB offline dump> [<depch_id.txt>]]
"""
import statistics
import sys


def read_offline(path):
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) >= 3:
            out[(int(f[0]), int(f[1]))] = float(f[2])
    return out


def read_anchor(path, depch):
    rows = [l.split() for l in open(path) if l.strip()]
    if len(rows[0]) == 5:  # ehp ns dep ch value
        out = {}
        for f in rows:
            ns, dep, ch = int(f[1]), int(f[2]), int(f[3])
            tid = depch.get((ns, dep, ch))
            if tid is not None:
                out[(ns, tid)] = float(f[4])
        return out, "online (ehp ns dep ch)"
    return read_offline(path), "offline (ns id)"


def main():
    anchor = sys.argv[1]
    dbfile = sys.argv[2] if len(sys.argv) > 2 else "db_fcsgaincorrOffline_23007048.txt"
    depfile = sys.argv[3] if len(sys.argv) > 3 else "liang_depch_id.txt"
    depch = {}
    for line in open(depfile):
        e, ns, dep, ch, tid = map(int, line.split())
        depch[(ns, dep, ch)] = tid
    a, fmt = read_anchor(anchor, depch)
    db = read_offline(dbfile)
    keys = sorted(set(a) & set(db))
    r = sorted(a[k] / db[k] for k in keys if db[k] > 0)
    n = len(r)
    print(f"{anchor}: {fmt} format, {len(a)} towers; DB {dbfile}: {len(db)}; compared {n}")
    print(f"  anchor / DB: median {statistics.median(r):.4f}, 16-84% {r[n // 6]:.4f}-{r[5 * n // 6]:.4f}, "
          f"mean {statistics.mean(r):.4f}, rms {statistics.pstdev(r):.4f}")
    print(f"  within 1%: {sum(abs(x - 1) < 0.01 for x in r) / n:.1%}, within 5%: {sum(abs(x - 1) < 0.05 for x in r) / n:.1%}")
    worst = sorted(keys, key=lambda k: -abs(a[k] / db[k] - 1))[:10]
    print("  largest deviations (ns, id: anchor, DB):", ", ".join(f"({k[0]},{k[1]}: {a[k]:.3f}, {db[k]:.3f})" for k in worst))


if __name__ == "__main__":
    main()
