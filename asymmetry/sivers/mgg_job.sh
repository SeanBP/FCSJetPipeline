#!/bin/bash
# The two-photon pair-mass fit of a report (mgg_peak.C: the pi0 peak used by
# AnAsymmetry -m) as a condor job: copies the yields to the job's scratch
# disk and puts mgg_peak.png/.pdf and mgg_peak.log into <report dir>.
# Usage (condor arguments): mgg_job.sh <hists dir> <code dir> <report dir>
set -e
HISTS="$1"; CODE="$2"; REP="$3"
SCR="${_CONDOR_SCRATCH_DIR:-$(pwd)}/mgg"
mkdir -p "$SCR/hists"; t0=$(date +%s)
echo "start $(date) on $(hostname)"
find -L "$HISTS" -maxdepth 1 -name 'an_hists_*.root' -mmin +10 > "$SCR/src.list"
xargs -a "$SCR/src.list" cp -p -t "$SCR/hists/"
ls "$SCR"/hists/*.root > "$SCR/yields.list"
echo "copied $(wc -l < "$SCR/yields.list") files in $(( $(date +%s) - t0 )) s"
cp "$CODE"/mgg_peak.C "$SCR/"
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B "$SCR" /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
    "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd $SCR; root -b -q -l 'mgg_peak.C(\"yields.list\",\"mgg_peak.png\")'" > "$SCR/mgg_peak.log" 2>&1
cp "$SCR"/mgg_peak.png "$SCR"/mgg_peak.pdf "$SCR"/mgg_peak.log "$REP/" 2>/dev/null || cp "$SCR"/mgg_peak.log "$REP/"
grep -E '^(peak|sigma)=' "$SCR/mgg_peak.log" || true
echo "done $(date), $(( $(date +%s) - t0 )) s"
