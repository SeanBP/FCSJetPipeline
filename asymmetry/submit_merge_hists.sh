#!/bin/bash
# Merges the per-job yield files of an A_N pass into one file per run
# (merge_hists_job.xml), to save disk. AnAsymmetry and smoke_report.sh read
# the merged files like the per-job ones (the run comes from the file's
# "runs" entry, not its name).
# Usage:
#   ./submit_merge_hists.sh -d <pass outdir> [-r <run list>] [-n <runs per job>] [-l <label>] [-N]
#
#   -d  the pass's output directory (its hists/ is merged in place)
#   -r  runs to merge, one per line (first column); default: every run with
#       more than one file in hists/
#   -n  runs per job, default 50
#   -l  label for the generated XML and logs, default "merge"
#   -N  write the XML only
# Only finished jobs' files should be merged: run it once the pass (or the
# batch holding these runs) has no job left in condor. A part arriving later
# is merged by running it again.

set -e
DIR=""; RUNS=""; NPER=50; LABEL="merge"; DRY=0
while getopts "d:r:n:l:Nh" opt; do
  case $opt in
    d) DIR="$OPTARG" ;;
    r) RUNS="$OPTARG" ;;
    n) NPER="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    N) DRY=1 ;;
    h|*) sed -n '2,18p' "$0"; exit 0 ;;
  esac
done
[ -d "$DIR/hists" ] || { sed -n '2,18p' "$0"; exit 1; }
[[ "$NPER" =~ ^[1-9][0-9]*$ ]] || { echo "Error: -n must be a positive integer"; exit 1; }
DIR="$(cd "$DIR" && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"

mkdir -p "$DIR/merge_logs"
LIST="$DIR/merge_logs/runs_${LABEL}.list"
if [ -n "$RUNS" ]; then
    awk '{print $1}' "$RUNS" | sort -u > "$LIST.all"
else
    find "$DIR/hists" -maxdepth 1 -name 'an_hists_*.root' | sed 's#.*/an_hists_\([0-9]*\)_.*#\1#' | sort -u > "$LIST.all"
fi
# only runs that have more than one file now
find "$DIR/hists" -maxdepth 1 -name 'an_hists_*.root' | sed 's#.*/an_hists_\([0-9]*\)_.*#\1#' | sort | uniq -c | awk '$1 > 1 {print $2}' > "$LIST.multi"
sort "$LIST.all" | join - "$LIST.multi" > "$LIST"
rm -f "$LIST.all" "$LIST.multi"
NR=$(wc -l < "$LIST")
if [ "$NR" -eq 0 ]; then echo "Nothing to merge"; echo "Submitted request none"; exit 0; fi
NPROC=$(( (NR + NPER - 1) / NPER ))
echo "Merge: ${NR} runs in ${NPROC} jobs of up to ${NPER}"

XML="${SCRIPT_DIR}/merge_hists_job_${LABEL}.xml"
sed -e "s#{{RUNLIST}}#${LIST}#g" -e "s/{{NPERJOB}}/${NPER}/g" -e "s/{{NPROCESSES}}/${NPROC}/g" \
    -e "s#{{HISTDIR}}#${DIR}/hists#g" -e "s#{{LOGDIR}}#${DIR}/merge_logs#g" \
    "${SCRIPT_DIR}/merge_hists_job.xml" > "$XML"
echo "Generated ${XML}"
[ "$DRY" = "1" ] && { echo "Dry run (-N): not submitted"; exit 0; }

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "asymmetry/merge_hists_job_${LABEL}.xml"
REPORT=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*.report' -newer "${MARKER}" | head -1)
REQID=$(basename "${REPORT}" .report | sed 's/^sched//')
rm -f "${MARKER}"
echo "${REQID}" > "$DIR/merge_logs/request_${LABEL}.txt"
# ahead of the pass's own idle jobs (same user: condor otherwise runs them in submission order; SUMS jobs here come with JobPrio 10)
condor_qedit -constraint "regexp(\"${REQID}\", Arguments)" JobPrio 20 > /dev/null 2>&1 || true
echo "Submitted request ${REQID}"
