#!/bin/bash
# AnPlots over a pass's yield files as a condor job: copies the yields to the
# job's scratch disk, compiles AnPlots from <code dir>, writes the pages and
# plots/qa_sums.root (for AnPlots -s / replot_report.sh) into <report dir>/plots.
# Usage (condor arguments): anplots_job.sh <hists dir> <code dir> <report dir>
set -e
HISTS="$1"; CODE="$2"; REP="$3"
SCR="${_CONDOR_SCRATCH_DIR:-$(pwd)}/anpl"
mkdir -p "$SCR/hists"; t0=$(date +%s)
echo "start $(date) on $(hostname)"
find -L "$HISTS" -maxdepth 1 -name 'an_hists_*.root' -mmin +10 > "$SCR/src.list"
xargs -a "$SCR/src.list" cp -p -t "$SCR/hists/"
ls "$SCR"/hists/*.root > "$SCR/yields.list"
echo "copied $(wc -l < "$SCR/yields.list") files in $(( $(date +%s) - t0 )) s"
cp "$CODE"/AnPlots.cpp "$CODE"/AnConfig.h "$SCR/"
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B "$SCR" /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
    "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd $SCR; g++ -m64 -O2 AnPlots.cpp -o AnPlots \`root-config --cflags --libs\` -std=c++11 && ./AnPlots yields.list plots" > "$SCR/anplots.log" 2>&1
mkdir -p "$REP/plots"
cp -r "$SCR/plots/." "$REP/plots/"
cp "$SCR/anplots.log" "$REP/anplots.log"
grep '^AnPlots:' "$SCR/anplots.log"
echo "done $(date), $(( $(date +%s) - t0 )) s"
