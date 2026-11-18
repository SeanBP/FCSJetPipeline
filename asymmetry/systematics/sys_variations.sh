#!/bin/bash
# Cut-variation systematics from AnHistograms' variation copies (Sean
# 2026-11-17; AnConfig.h kVarSuffix, AnAsymmetry -v). On a finished report
# folder (smoke_report.sh output: build/AnAsymmetry, build/pi0_bgsub.C, and
# masked.root or yields.list) it runs, for the nominal selection and each
# variation,
#   photon classes: AnAsymmetry -J [-v <var>]           -> an_<var>.txt
#   pi0 (FMS iso / noniso, nominal z_em thresholds): AnAsymmetry -J -I
#   <cls>:sig / :sb [-v <var>] and pi0_bgsub.C           -> pi0_<cls>_<var>.txt
# and compares them with sys_slices.py:
#   nested (loose = nominal, tight = variation; Barlow systematic):
#     fid25, fid30    jet axis >= 25 / 30 cm inside the ECal edge (nominal 20)
#     vz60            |vz| < 60 cm (nominal 100)
#     vtxmeas         measured vertex only (the no-vertex systematic)
#   shifted (same events, other jet kinematics; sys = |A_N shift|):
#     jesup, jesdn    jet E, pT, x_F x (1 +- delta(E)), the JES map's error
#                     (= the calibration rerun with the map shifted; photon
#                     classes only: the pi0 candidates do not use the JES);
#                     symmetric, the larger side: combine_sys.py --comp-max jes=<up>,<dn>
#   disjoint (a check, not a box):
#     bbcveto         the events the collision cut removes vs the nominal
#                     sample; also their share of the jets per row
#                     (bbcveto_share.txt)
# Each comparison is drawn with plot_sys_slices.C (blue beam; union "all" row
# and pi0). Outputs: <report>/sys_variations/sys_<var>_<sample>.txt,
# _summary.txt, .png; combine_sys.py collects them. Needs yield files made
# from 2026-11-17 on (AnAsymmetry stops otherwise).
# Usage: sys_variations.sh <report dir> [AnAsymmetry input, default "-M masked.root"] [variations, default "fid25 fid30 vz60 vtxmeas bbcveto jesup jesdn"]
set -e
REP="$(cd "$1" && pwd)"
ANIN="${2:--M masked.root}"
VARS="${3:-fid25 fid30 vz60 vtxmeas bbcveto jesup jesdn}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${REP}/sys_variations"
[ -x "${REP}/build/AnAsymmetry" ] || { sed -n '2,26p' "$0"; echo "no ${REP}/build/AnAsymmetry"; exit 1; }
mkdir -p "$OUT"
cp "${HERE}/plot_sys_slices.C" "$OUT/"
container() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 -B /tmp \
        /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd ${REP}; $1"
}
# one selection: photon classes and the two pi0 classes
run_sel() {  # tag, AnAsymmetry variation option ("" = nominal)
    local tag="$1" vopt="$2"
    container "./build/AnAsymmetry -J ${vopt} ${ANIN} sys_variations/an_${tag}.root" > "${OUT}/an_${tag}.log" 2>&1 \
        || { tail -3 "${OUT}/an_${tag}.log"; exit 1; }
    [[ "$tag" == jes* ]] && { echo "sys_variations: ${tag} done (photon classes only)"; return; }
    for cls in iso noniso; do
        for w in sig sb; do
            container "./build/AnAsymmetry -J -I ${cls}:${w} ${vopt} ${ANIN} sys_variations/an_${cls}_${tag}_${w}.root" > "${OUT}/an_${cls}_${tag}_${w}.log" 2>&1 \
                || { tail -3 "${OUT}/an_${cls}_${tag}_${w}.log"; exit 1; }
        done
        container "root -b -q -l 'build/pi0_bgsub.C+(\"sys_variations/an_${cls}_${tag}_sig.root\",\"sys_variations/an_${cls}_${tag}_sb.root\",\"sys_variations/pi0_${cls}_${tag}\",\"all\",\"${cls} ${tag}\")'" \
            > "${OUT}/pi0_bgsub_${cls}_${tag}.log" 2>&1
    done
    echo "sys_variations: ${tag} done"
}
run_sel nominal ""
for v in $VARS; do run_sel "$v" "-v $v"; done

cmp() {  # name, loose table, tight table, label, extra sys_slices options, plot row (level class)
    python3 -I "${HERE}/sys_slices.py" "$2" "$3" "${OUT}/sys_$1" --label "$4" --min-jets 50 $5
    container "cd sys_variations; root -b -q -l 'plot_sys_slices.C+(\"sys_$1.txt\",\"sys_$1.png\",\"$6\",\"$7\",\"blue\",\"$4\")'" > "${OUT}/plot_$1.log" 2>&1 || true
}
declare -A LABEL=([fid25]="fiducial 20 vs 25 cm" [fid30]="fiducial 20 vs 30 cm" [vz60]="|v_{z}| < 100 vs < 60 cm"
                  [vtxmeas]="all vertices vs measured only" [bbcveto]="collision cut: kept vs removed events"
                  [jesup]="JES scale +delta(E)" [jesdn]="JES scale -delta(E)")
for v in $VARS; do
    opt=""; [ "$v" = bbcveto ] && opt="--disjoint"; [[ "$v" == jes* ]] && opt="--shifted"
    cmp "${v}_jets" "${OUT}/an_nominal.txt" "${OUT}/an_${v}.txt" "${LABEL[$v]}" "$opt" all all
    [[ "$v" == jes* ]] && continue
    for cls in iso noniso; do
        cmp "${v}_pi0${cls}" "${OUT}/pi0_${cls}_nominal.txt" "${OUT}/pi0_${cls}_${v}.txt" "${LABEL[$v]} (${cls} #pi^{0})" "$opt" all pi0
    done
done
# share of the jets the collision cut removes, per row (sum over x_F bins and both beams)
if [[ " $VARS " == *" bbcveto "* ]]; then
    python3 -I - "${OUT}/an_nominal.txt" "${OUT}/an_bbcveto.txt" > "${OUT}/bbcveto_share.txt" <<'EOF'
import sys
def rows(p):
    d = {}
    for line in open(p):
        if line.startswith('#'): continue
        t = line.split()
        if len(t) < 9: continue
        d[(t[0], t[1])] = d.get((t[0], t[1]), 0) + float(t[8])
    return d
a, b = rows(sys.argv[1]), rows(sys.argv[2])
print('# level class jets_kept jets_removed removed_share (collision cut east BBC >= 1, TOF > 2; both beams summed)')
for k in sorted(a):
    n, m = a[k], b.get(k, 0)
    print('%s %s %.0f %.0f %.4f' % (k[0], k[1], n, m, m / (n + m) if n + m > 0 else 0))
EOF
fi
grep -h "^all" "${OUT}"/sys_*_summary.txt | sed 's/^/  /'
echo "sys_variations: outputs in ${OUT}"
