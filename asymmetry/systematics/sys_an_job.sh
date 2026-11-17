#!/bin/bash
# One AnAsymmetry -r pass over a list of per-run yield files, as a condor job
# (2026-11-17; like sivers/an_pass_job.sh, which takes a whole hists dir): the
# listed yields are copied to the job's scratch disk, AnAsymmetry is compiled
# from <code dir> (AnAsymmetry.cpp, AnConfig.h, AnRunTools.h,
# run22_polarization.txt side by side) and <out>.root, <out>.txt and <out>.log
# go to <out dir>. Usage (condor arguments):
#   sys_an_job.sh <yield list> <code dir> <out dir> <out> [AnAsymmetry options...]
# The options go before -r. submit_sys_an.sh writes the condor file.
set -e
LIST="$1"; CODE="$2"; OUTDIR="$3"; OUT="$4"; shift 4
SCR="${_CONDOR_SCRATCH_DIR:-$(pwd)}/sysan"
mkdir -p "$SCR/hists"; t0=$(date +%s)
echo "start $(date) on $(hostname): AnAsymmetry $* over $(wc -l < "$LIST") files -> ${OUT}"
xargs -a "$LIST" cp -pL -t "$SCR/hists/"
ls "$SCR"/hists/*.root > "$SCR/yields.list"
echo "copied $(wc -l < "$SCR/yields.list") files in $(( $(date +%s) - t0 )) s"
cp "$CODE"/AnAsymmetry.cpp "$CODE"/AnConfig.h "$CODE"/AnRunTools.h "$CODE"/run22_polarization.txt "$SCR/"
ct() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B "$SCR" /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd $SCR; $1"
}
ct "g++ -m64 -O2 AnAsymmetry.cpp -o AnAsymmetry \`root-config --cflags --libs\` -std=c++11"
ct "./AnAsymmetry $* -r yields.list run22_polarization.txt ${OUT}.root" > "$SCR/${OUT}.log"
cp "$SCR/${OUT}.root" "$SCR/${OUT}.txt" "$SCR/${OUT}.log" "$OUTDIR/"
rm -rf "${SCR:?}/hists"
echo "done $(date), $(( $(date +%s) - t0 )) s"
