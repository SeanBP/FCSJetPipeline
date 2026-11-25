#!/usr/bin/env python3
"""Point-threshold JES maps by the same-events difference (Sean 2026-11-24).

Only ~1M events (the jes_full_v1 jobs that kept their SimpleTree) can be
rebuilt at another point threshold, and a map from 1M events differs from the
10M map by its noise and by the estimator change at low statistics (more
kernel-peak points). Both largely cancel between two maps of the SAME events,
so the threshold map for the analysis is

    F_thr(E) = F_base(E) + [F_sub,thr(E) - F_sub,ref(E)]

with F_base the standard (10M) map and F_sub,thr / F_sub,ref the subset maps at
the new and the nominal threshold. All three are log-cubic, so the sum is
exact coefficient by coefficient: A, B per grid node, and C, D (face-wide, or
per node where a map has them). Nodes missing from a subset map keep the base
map's values (difference 0). Emin / Emax are the base map's: above the subset
maps' Emax the difference follows their cubic in ln E.

--smooth s: Gaussian smoothing of the per-node A, B differences over the
grid (s grid bins, NaN-aware), as the threshold dependence is expected to be
smooth across the face; 0 (default) = none.

  jes_map_combine.py <base.json> <sub_thr.json> <sub_ref.json> <out.json> [--smooth s] [--label text]
"""
import argparse
import json
import math


def node_map(m):
    return {(b["i"], b["j"]): b for b in m["bins"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("sub_thr")
    ap.add_argument("sub_ref")
    ap.add_argument("out")
    ap.add_argument("--smooth", type=float, default=0.0)
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    base, thr, ref = (json.load(open(f)) for f in (a.base, a.sub_thr, a.sub_ref))
    for m, f in ((thr, a.sub_thr), (ref, a.sub_ref)):
        if m.get("model") != "logcubic" or m["x_grid"] != base["x_grid"] or m["y_grid"] != base["y_grid"]:
            raise SystemExit("jes_map_combine: %s is not a log-cubic map on the base map's grid" % f)
    nb, nt, nr = node_map(base), node_map(thr), node_map(ref)
    # per-node differences of A, B, C, D (C, D per node where a map has them, else face-wide)
    diff = {}
    for key in nb:
        if key in nt and key in nr:
            diff[key] = [nt[key][p] - nr[key][p] for p in ("A", "B")] + \
                        [nt[key].get(p, thr[p]) - nr[key].get(p, ref[p]) for p in ("C", "D")]
    if a.smooth > 0:
        sm = {}
        r = int(math.ceil(3 * a.smooth))
        for (i, j) in nb:
            w = s0 = s1 = 0.0
            for di in range(-r, r + 1):
                for dj in range(-r, r + 1):
                    d = diff.get((i + di, j + dj))
                    if d is None:
                        continue
                    g = math.exp(-0.5 * (di * di + dj * dj) / (a.smooth * a.smooth))
                    w += g; s0 += g * d[0]; s1 += g * d[1]
            if w > 0 and (i, j) in diff:
                sm[(i, j)] = [s0 / w, s1 / w] + diff[(i, j)][2:]
        diff = sm
    per_node_cd = any("C" in b or "D" in b for m in (base, thr, ref) for b in m["bins"])
    out = dict(base)
    out["bins"] = []
    for key, b in nb.items():
        d = diff.get(key, [0.0, 0.0, thr["C"] - ref["C"], thr["D"] - ref["D"]])
        nbn = {"i": b["i"], "j": b["j"], "A": b["A"] + d[0], "B": b["B"] + d[1]}
        if per_node_cd:
            nbn["C"] = b.get("C", base["C"]) + d[2]
            nbn["D"] = b.get("D", base["D"]) + d[3]
        out["bins"].append(nbn)
    if not per_node_cd:
        out["C"] = base["C"] + thr["C"] - ref["C"]
        out["D"] = base["D"] + thr["D"] - ref["D"]
    out["sample"] = "%s + (%s - %s)" % (base.get("sample", "?"), thr.get("sample", "?"), ref.get("sample", "?"))
    out["combined_from"] = {"base": a.base, "sub_thr": a.sub_thr, "sub_ref": a.sub_ref, "smooth": a.smooth, "label": a.label}
    json.dump(out, open(a.out, "w"), indent=1)
    print("jes_map_combine: %d nodes (%d with a subset difference), smooth %g -> %s" % (len(nb), len(diff), a.smooth, a.out))


if __name__ == "__main__":
    main()
