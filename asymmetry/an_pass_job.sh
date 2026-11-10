#!/bin/bash
# One AnAsymmetry pass over a pass's yield files as a condor job: copies the
# yields to the job's scratch disk, compiles AnAsymmetry from <code dir>, runs
# it and puts <out>.root, <out>.txt and <log name> into <report dir>. Several
# of these run side by side to redo a report's A_N steps (each pass is ~40
# min over ~2900 per-run files); smoke_report.sh with SMOKE_PLOTS_ONLY=1 then
# makes the plots. Usage (condor arguments):
#   an_pass_job.sh <hists dir> <code dir> <report dir> <out> <log name> [AnAsymmetry options...]
# The options go before -r; the yield list, the polarization table
# (<code dir>/run22_polarization.txt) and <out>.root are added.
set -e
HISTS="$1"; CODE="$2"; REP="$3"; OUT="$4"; LOG="$5"; shift 5
SCR="${_CONDOR_SCRATCH_DIR:-$(pwd)}/anp"
mkdir -p "$SCR/hists"; t0=$(date +%s)
echo "start $(date) on $(hostname): AnAsymmetry $* -> ${OUT}"
find -L "$HISTS" -maxdepth 1 -name 'an_hists_*.root' -mmin +10 > "$SCR/src.list"
xargs -a "$SCR/src.list" cp -p -t "$SCR/hists/"
ls "$SCR"/hists/*.root > "$SCR/yields.list"
echo "copied $(wc -l < "$SCR/yields.list") files in $(( $(date +%s) - t0 )) s"
cp "$CODE"/AnAsymmetry.cpp "$CODE"/AnConfig.h "$CODE"/run22_polarization.txt "$SCR/"
ct() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B "$SCR" /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd $SCR; $1"
}
ct "g++ -m64 -O2 AnAsymmetry.cpp -o AnAsymmetry \`root-config --cflags --libs\` -std=c++11"
ct "./AnAsymmetry $* -r yields.list run22_polarization.txt ${OUT}.root" > "$SCR/$LOG"
cp "$SCR/${OUT}.root" "$SCR/${OUT}.txt" "$SCR/$LOG" "$REP/"
echo "done $(date), $(( $(date +%s) - t0 )) s"
