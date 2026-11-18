#!/bin/bash
# Data-side statistical noise in the JES shift (Sean 2026-11-17). The _vjesup /
# _vjesdn copies use the same events as the nominal yields, with each jet's
# x_F scaled by 1 +- delta(E), so A_N(jes) - A_N(nominal) mixes the true JES effect delta
# with noise from the jets that change x_F bins. The yield files of a pass
# are independent slices of the data: split into k random groups, each group
# gives Delta_g = A_jes - A_nom, and
#     <Delta^2>(k) = delta^2 + k s^2
# (s: the noise of the full sample) separates the two
# (jes_noise_scaling.py groups). k = 1 is the full sample.
# Per group: AnAsymmetry -r (the run selection of smoke_report.sh: per-run
# polarization, the default bad-run list; extra options with -A, e.g. the
# pi0-peak list -Q and variable bins -x), -J, nominal, -v jesup and -v jesdn.
#
# Usage: jes_data_subsets.sh -y <yield file list> -o <outdir> [-k "2 4 8"] [-S <seed>] [-A "<AnAsymmetry options>"]
#   -y  yield files (e.g. a report's yields.list)
#   -o  output: groups/k<k>_g<g>.list, tables/an_{nom,jesup,jesdn}_k<k>_g<g>.txt,
#       groups_jesup.txt and groups_jesdn.txt ("k nominal_table jes_table" per
#       group, k = 1 the full sample)
set -e
YL=""; OUT=""; KLIST="2 4 8"; SEED=1; AOPT=""
while getopts "y:o:k:S:A:h" opt; do
    case $opt in
        y) YL="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        k) KLIST="$OPTARG" ;;
        S) SEED="$OPTARG" ;;
        A) AOPT="$OPTARG" ;;
        h|*) sed -n '2,20p' "$0"; exit 0 ;;
    esac
done
[ -f "$YL" ] && [ -n "$OUT" ] || { sed -n '2,20p' "$0"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASYM="$(dirname "$HERE")"
[ ! -e "$OUT" ] || { echo "Error: $OUT exists; use a new folder"; exit 1; }
mkdir -p "$OUT"/{build,groups,tables}
OUT="$(cd "$OUT" && pwd)"
cp "$ASYM/sivers/AnAsymmetry.cpp" "$ASYM/AnConfig.h" "$ASYM/AnRunTools.h" "$OUT/build/"
container() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 -B /tmp \
        /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd ${OUT}; $1"
}
container "cd build; g++ -m64 -O2 AnAsymmetry.cpp -o AnAsymmetry \`root-config --cflags --libs\` -std=c++11" > "$OUT/build.log" 2>&1 || { tail "$OUT/build.log"; exit 1; }

# the groups: shuffled files (seeded), file i -> group i % k
python3 -I - "$YL" "$SEED" "$OUT/groups" $KLIST <<'EOF'
import os, random, sys
files = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
seed, out = int(sys.argv[2]), sys.argv[3]
with open(os.path.join(out, 'k1_g1.list'), 'w') as f:
    f.write('\n'.join(files) + '\n')
for k in [int(x) for x in sys.argv[4:]]:
    fl = list(files)
    random.Random('%d-%d' % (seed, k)).shuffle(fl)
    for g in range(k):
        with open(os.path.join(out, 'k%d_g%d.list' % (k, g + 1)), 'w') as f:
            f.write('\n'.join(sorted(fl[g::k])) + '\n')
print('jes_data_subsets: %d yield files; divisions 1 %s' % (len(files), ' '.join(sys.argv[4:])))
EOF
: > "$OUT/groups_jesup.txt"; : > "$OUT/groups_jesdn.txt"
for gl in $(ls "$OUT"/groups/k*_g*.list | sort -V); do
    tag=$(basename "$gl" .list); k=${tag#k}; k=${k%%_*}
    for v in nom jesup jesdn; do
        vopt=""; [ "$v" != nom ] && vopt="-v $v"
        container "./build/AnAsymmetry -J ${vopt} ${AOPT} -r groups/${tag}.list ${ASYM}/run22_polarization.txt tables/an_${v}_${tag}.root" > "$OUT/tables/an_${v}_${tag}.log" 2>&1 \
            || { tail -3 "$OUT/tables/an_${v}_${tag}.log"; exit 1; }
    done
    for v in jesup jesdn; do echo "$k $OUT/tables/an_nom_${tag}.txt $OUT/tables/an_${v}_${tag}.txt" >> "$OUT/groups_${v}.txt"; done
    echo "jes_data_subsets: ${tag} ($(wc -l < "$gl") files) done"
done
echo "jes_data_subsets: $OUT/groups_jesup.txt, groups_jesdn.txt; then jes_noise_scaling.py groups <file> --out <prefix>"
