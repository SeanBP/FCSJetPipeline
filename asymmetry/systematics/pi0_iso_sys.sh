#!/bin/bash
# Isolation-threshold systematic of the FMS-style pi0 A_N (Sean 2026-11-17).
# On a finished report folder (smoke_report.sh output: build/AnAsymmetry,
# build/pi0_bgsub.C, and masked.root or yields.list), for each z_em threshold
#   iso:    0.95, 0.98 (nominal), 0.99
#   noniso: 0.85, 0.90 (nominal), 0.95
# it runs AnAsymmetry -J -I <cls>=<z>:sig and :sb and pi0_bgsub.C (the
# background-subtracted A_N, with its own signal fraction), then compares the
# nested pairs with sys_slices.py (loose contains tight, other cuts nominal):
#   iso loose:    z > 0.95 vs z > 0.98      iso tight:    z > 0.98 vs z > 0.99
#   noniso loose: z < 0.95 vs z < 0.90      noniso tight: z < 0.90 vs z < 0.85
# and draws each with plot_sys_slices.C (blue beam). Outputs go to
# <report>/pi0_iso_sys/. Needs yield files made from 2026-11-17 on (z_em bin
# edges 0.85 ... 0.99; AnAsymmetry stops otherwise).
# Usage: pi0_iso_sys.sh <report dir> [AnAsymmetry input, default "-M masked.root"]
set -e
REP="$(cd "$1" && pwd)"
ANIN="${2:--M masked.root}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${REP}/pi0_iso_sys"
[ -x "${REP}/build/AnAsymmetry" ] || { sed -n '2,17p' "$0"; echo "no ${REP}/build/AnAsymmetry"; exit 1; }
mkdir -p "$OUT"
cp "${HERE}/plot_sys_slices.C" "$OUT/"
container() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
        /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd ${REP}; $1"
}
for v in iso=0.95 iso=0.98 iso=0.99 noniso=0.85 noniso=0.9 noniso=0.95; do
    tag="${v/=/_}"
    for w in sig sb; do
        container "./build/AnAsymmetry -J -I ${v}:${w} ${ANIN} pi0_iso_sys/an_${tag}_${w}.root" > "${OUT}/an_${tag}_${w}.log" 2>&1 \
            || { tail -3 "${OUT}/an_${tag}_${w}.log"; exit 1; }
    done
    container "root -b -q -l 'build/pi0_bgsub.C+(\"pi0_iso_sys/an_${tag}_sig.root\",\"pi0_iso_sys/an_${tag}_sb.root\",\"pi0_iso_sys/pi0_${tag}\",\"all\",\"${v}\")'" \
        > "${OUT}/pi0_bgsub_${tag}.log" 2>&1
    echo "pi0_iso_sys: ${v} done"
done
cmp() {  # name loose tight label
    python3 -I "${HERE}/sys_slices.py" "${OUT}/pi0_$2.txt" "${OUT}/pi0_$3.txt" "${OUT}/sys_$1" --label "$4" --min-jets 50
    container "cd pi0_iso_sys; root -b -q -l 'plot_sys_slices.C+(\"sys_$1.txt\",\"sys_$1.png\",\"all\",\"pi0\",\"blue\",\"$4\")'" > "${OUT}/plot_$1.log" 2>&1
}
cmp iso_loose iso_0.95 iso_0.98 "isolated: z_{em} > 0.95 vs > 0.98"
cmp iso_tight iso_0.98 iso_0.99 "isolated: z_{em} > 0.98 vs > 0.99"
cmp noniso_loose noniso_0.95 noniso_0.9 "non-isolated: z_{em} < 0.95 vs < 0.9"
cmp noniso_tight noniso_0.9 noniso_0.85 "non-isolated: z_{em} < 0.9 vs < 0.85"
grep -h "^all" "${OUT}"/sys_*_summary.txt | sed 's/^/  /'
echo "pi0_iso_sys: outputs in ${OUT}"
