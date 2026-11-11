#!/bin/bash
# Runs smoke_report.sh as one condor job (report_job.sh): the yield files are
# copied to the worker's scratch disk first, which makes the report's six
# passes over them minutes instead of hours each on a full pass.
# Usage: ./submit_report_job.sh <outdir> <good_runs.txt> [report name = report] [memory MB = 1000]
# The report code (smoke_report.sh and what it compiles) is frozen into
# <outdir>/<report name>_code at submission. Log: <outdir>/<report name>_job.log
set -e
OUTDIR="$(cd "$1" && pwd)"
GOOD="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
NAME="${3:-report}"; MEM="${4:-1000}"
[ -d "$OUTDIR/hists" ] && [ -f "$GOOD" ] || { sed -n '2,8p' "$0"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# a file of this script's directory, else of the one above (repo layout:
# asymmetry/sivers/ with the pass-level files in asymmetry/; a frozen report
# code directory has them all side by side)
src() { if [ -e "${HERE}/$1" ]; then echo "${HERE}/$1"; else echo "${HERE}/../$1"; fi; }
CODE="$OUTDIR/${NAME}_code"
rm -rf "$CODE"; mkdir -p "$CODE"
for f in smoke_report.sh report_job.sh AnAsymmetry.cpp AnPlots.cpp AnConfig.h plot_an.C plot_an_means.C mgg_peak.C plot_an_mass.C replot_report.sh project_errors.py estimate_resources.py run22_polarization.txt catalog_runs_local.txt run22_run_events_all.txt; do
    cp "$(src $f)" "$CODE/"
done
# the repository's other files the report code includes or reads
for f in $(grep -h -o '#include "[^"]*"' "$(src AnAsymmetry.cpp)" "$(src AnPlots.cpp)" "$(src AnConfig.h)" | sed 's/#include "\(.*\)"/\1/' | sort -u); do
    [ -f "$(src $f)" ] && cp "$(src $f)" "$CODE/"
done
NFILES=$(find "$OUTDIR/hists" -maxdepth 1 -name 'an_hists_*.root' | wc -l)
# requests that fit a freed 1 GB farm slot: a job asking for more waits behind
# the pass's own 1 GB jobs, whose claims condor hands on to each other. The
# worker's scratch disk is shared and large; the copy needs ~1.5 x hists/.
DISK_KB=1000
cat > "$CODE/report.sub" <<EOF
universe       = vanilla
executable     = $CODE/report_job.sh
arguments      = $OUTDIR $GOOD $NAME $CODE
output         = $OUTDIR/${NAME}_job.log
error          = $OUTDIR/${NAME}_job.log
log            = $CODE/report.condor.log
request_memory = $MEM
request_disk   = $DISK_KB
getenv         = true
priority       = 20
queue
EOF
echo "report job: ${NFILES} yield files, disk ${DISK_KB} kB, code frozen in ${CODE}"
condor_submit "$CODE/report.sub"
