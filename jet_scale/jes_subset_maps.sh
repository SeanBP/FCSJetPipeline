#!/bin/bash
# JES map statistical noise (Sean 2026-11-17): maps built on disjoint subsets
# of one simulation sample, with the same settings as its full map. The pairs
# of independent maps at each subset size give the map's statistical noise and
# its scaling (jes_map_compare.py noise): for maps on 1/k of the sample,
#   Var(map_a - map_b) = 2 k sigma_1^2   (sigma_1: noise of the full-sample map)
# so Var / (2k) must stay constant over k if the noise scales as 1/sqrt(n)
# and the method itself does not change with statistics (the share of
# Gaussian-core cells, from each map's log, must stay put as well).
#
# Each subset is a folder of symlinks to the sample's jet files. Files are
# dealt out per ptHat bin (the "_pt<N>_" in the file name), shuffled with a
# fixed seed, so every subset keeps every ptHat bin (the cross-section
# weights, JetXSecWeights.h, need each bin) with ~1/k of its files. The
# sample's run_params file is copied in for the map provenance. Then one
# submit_jet_scale_em_only.sh job per subset (-k 1: the logs carry the core /
# kernel-peak counts; -s 0).
#
# Usage: jes_subset_maps.sh -i <flat jet dir> -o <outdir> [-k "2 4 8"] [-S <seed>]
#                           [-t em|full] [-F "2"] [-N]
#   -i  folder of the full sample's jet files (the one the full map was built
#       from, truncated jobs already left out)
#   -o  output: <outdir>/subsets/k<k>_g<g>/ (symlinks, files.txt) and
#       <outdir>/maps/k<k>_g<g>/ (each map), manifest.txt (k, group, files,
#       map folder)
#   -k  subset divisions, default "2 4 8" (2 + 4 + 8 = 14 maps)
#   -S  shuffle seed, default 1
#   -t  jet type for submit_jet_scale_em_only.sh, default em
#   -F  free-parameter variants (submit_jet_scale_em_only.sh -F), default "2"
#       (A, B per position, C, D face-wide: the settings of jesmap_pixring_v1)
#   -N  make the subset folders and the manifest only, submit nothing
set -e
IN=""; OUT=""; KLIST="2 4 8"; SEED=1; JTYPE="em"; FREE="2"; DRY=0
while getopts "i:o:k:S:t:F:Nh" opt; do
    case $opt in
        i) IN="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        k) KLIST="$OPTARG" ;;
        S) SEED="$OPTARG" ;;
        t) JTYPE="$OPTARG" ;;
        F) FREE="$OPTARG" ;;
        N) DRY=1 ;;
        h|*) sed -n '2,36p' "$0"; exit 0 ;;
    esac
done
[ -n "$IN" ] && [ -n "$OUT" ] || { sed -n '2,36p' "$0"; exit 1; }
IN="$(cd "$IN" && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="jet_output_"; [ "$JTYPE" = full ] && PREFIX="jet_full_"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
[ ! -e "${OUT}/manifest.txt" ] || { echo "Error: ${OUT}/manifest.txt exists (subsets already made); use a new -o"; exit 1; }
for k in $KLIST; do [[ "$k" =~ ^[2-9]$|^[1-9][0-9]+$ ]] || { echo "Error: -k entries must be integers >= 2 (got '$k')"; exit 1; }; done

# the sample's files, grouped by ptHat bin
ALL="${OUT}/all_files.txt"
find "$IN" -maxdepth 1 -name "${PREFIX}*.root" | sort > "$ALL"
N=$(wc -l < "$ALL")
[ "$N" -gt 0 ] || { echo "Error: no ${PREFIX}*.root in $IN"; exit 1; }
RUNPARAMS=$(ls "${IN}"/run_params_*.txt "${IN}"/../run_params_*.txt 2>/dev/null | head -1 || true)
echo "jes_subset_maps: ${N} files in ${IN}; divisions ${KLIST}; seed ${SEED}; run_params ${RUNPARAMS:-none}"

# assignment: per ptHat bin, files shuffled (seeded), file i of the bin to group i % k
python3 -I - "$ALL" "$SEED" "$OUT" $KLIST <<'EOF'
import os, random, re, sys
path, seed, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
ks = [int(k) for k in sys.argv[4:]]
files = [l.strip() for l in open(path) if l.strip()]
bins = {}
for f in files:
    m = re.search(r'_pt(\d+)_', os.path.basename(f))
    if not m:
        sys.exit('jes_subset_maps: no _pt<N>_ in %s (ptHat bin unknown)' % f)
    bins.setdefault(int(m.group(1)), []).append(f)
man = open(os.path.join(out, 'manifest.txt'), 'w')
man.write('# k group files subset_dir map_dir   (jes_subset_maps.sh; per ptHat bin, file i of the shuffled bin -> group i %% k; seed %d)\n' % seed)
for k in ks:
    groups = [[] for _ in range(k)]
    for b in sorted(bins):
        fl = sorted(bins[b])
        random.Random('%d-%d-%d' % (seed, k, b)).shuffle(fl)
        for i, f in enumerate(fl):
            groups[i % k].append(f)
    for g in range(k):
        name = 'k%d_g%d' % (k, g + 1)
        d = os.path.join(out, 'subsets', name)
        os.makedirs(d)
        with open(os.path.join(d, 'files.txt'), 'w') as fo:
            for f in sorted(groups[g]):
                os.symlink(f, os.path.join(d, os.path.basename(f)))
                fo.write(f + '\n')
        man.write('%d %d %d %s %s\n' % (k, g + 1, len(groups[g]), d, os.path.join(out, 'maps', name)))
    sizes = [len(x) for x in groups]
    print('k=%d: files per group %s; every ptHat bin in every group: %s' % (k, sizes, all(len(bins[b]) >= k for b in bins)))
man.close()
EOF
# disjointness within each division and complete coverage
for k in $KLIST; do
    tot=$(cat "${OUT}"/subsets/k${k}_g*/files.txt | wc -l)
    uniq=$(cat "${OUT}"/subsets/k${k}_g*/files.txt | sort -u | wc -l)
    [ "$tot" = "$N" ] && [ "$uniq" = "$N" ] || { echo "Error: k=${k} subsets are not a partition ($tot listed, $uniq distinct, $N files)"; exit 1; }
done
# provenance: the sample's run_params in each subset folder (the label gets the subset name)
if [ -n "$RUNPARAMS" ]; then
    while read -r k g n sdir mdir; do
        [[ "$k" == \#* ]] && continue
        sed "s/^label=\(.*\)$/label=\1_$(basename "$sdir")/" "$RUNPARAMS" > "${sdir}/$(basename "$RUNPARAMS")"
    done < "${OUT}/manifest.txt"
fi
echo "jes_subset_maps: subsets in ${OUT}/subsets, manifest ${OUT}/manifest.txt"
[ "$DRY" = 1 ] && { echo "jes_subset_maps: -N, nothing submitted"; exit 0; }

while read -r k g n sdir mdir; do
    [[ "$k" == \#* ]] && continue
    name=$(basename "$sdir")
    "${HERE}/submit_jet_scale_em_only.sh" -i "$sdir" -o "$mdir" -l "sub_${name}" -s 0 -k 1 -t "$JTYPE" -F "$FREE"
done < "${OUT}/manifest.txt"
echo "jes_subset_maps: submitted $(grep -vc '^#' "${OUT}/manifest.txt") map jobs; then jes_map_compare.py noise ${OUT}/manifest.txt"
