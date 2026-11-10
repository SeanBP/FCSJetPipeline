#!/bin/bash
# Runs smoke_report.sh on a condor worker: copies the pass's yield files to the
# job's scratch disk (reading ~10k histograms per file straight from
# /star/data01 takes ~3 s a file, copying them ~0.15 s), runs every report
# step there, and copies the report back to <outdir>/<report name>. The
# condor resource table (condor_jobs.txt) is filled in from the submit node
# afterwards, if wanted: condor_history does not run on the workers.
# Submitted by submit_report_job.sh; arguments:
#   report_job.sh <outdir> <good_runs.txt> <report name> <code dir>
set -e
OUTDIR="$1"; GOOD="$2"; NAME="$3"; CODE="$4"
SCR="${_CONDOR_SCRATCH_DIR:-${SCRATCH:-$(pwd)}}/rep"
mkdir -p "$SCR/pass/hists"
echo "start $(date) on $(hostname), scratch $SCR"
t0=$(date +%s)
# yield files: copy, keeping times (smoke_report skips files younger than 10 min)
find "$OUTDIR/hists" -maxdepth 1 -name 'an_hists_*.root' -mmin +10 > "$SCR/src.list"
xargs -a "$SCR/src.list" cp -p -t "$SCR/pass/hists/"
echo "copied $(ls "$SCR/pass/hists" | wc -l) of $(wc -l < "$SCR/src.list") yield files in $(( $(date +%s) - t0 )) s ($(du -sh "$SCR/pass/hists" | cut -f1))"
ln -s "$OUTDIR/timing" "$SCR/pass/timing"
for f in "$OUTDIR"/run_params_*.txt; do ln -s "$f" "$SCR/pass/"; done
export SMOKE_BIND="$SCR"
# keep what the report made even if a late step fails
set +e
"$CODE/smoke_report.sh" "$SCR/pass" "$GOOD" "$NAME"
RC=$?
set -e
[ "$RC" = "0" ] || echo "WARNING: smoke_report.sh exited ${RC}; copying back what it made"
# the report back next to the pass (yields.list keeps the scratch paths)
rm -rf "$SCR/pass/$NAME/build"
mkdir -p "$OUTDIR/$NAME"
cp -r "$SCR/pass/$NAME/." "$OUTDIR/$NAME/"
echo "done $(date), $(( $(date +%s) - t0 )) s, smoke_report exit ${RC}"
