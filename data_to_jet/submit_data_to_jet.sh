#!/bin/bash
# Generates and submits a concrete data_to_jet.xml job.
# Usage:
#   ./submit_data_to_jet.sh -f <listfile> -s <0|1> -o <outdir> \
#       [-j <nprocesses>] [-l <label>] [-k <0|1>] [-t <triggers>] [-E <0|1>]
#
#   -f  absolute path to a file catalog list (MakeRunList.pl's
#       "path,filename,events" format: xrootd URL or local path, event
#       count appended as a trailing number)
#   -s  1 = keep intermediate SimpleTree_mudst files, 0 = discard
#   -o  output directory (created if missing)
#   -j  number of parallel jobs; listfile lines split evenly across them, default 1
#   -l  short label for the generated XML's filename, default "run"
#   -k  1 = save stdout/stderr under <outdir>/log/, 0 = discard, default 0
#   -t  comma-separated FCS trigger flag names (see FcsTriggerDefs.h); an
#       event only gets a jetTree entry if one fired. Empty = keep every event.
#   -E  1 = build reco jets from ECAL hits only (no HCAL), using the
#       ECAL-only fiducial boundary (../shared/JetParameters.h's
#       kFiducialRectEcal/reco_fiducial_buffer_ecal). 0/omit = previous
#       ECAL+HCAL behavior.
#
# Example:
#   ./submit_data_to_jet.sh -f /star/u/seanp/FCSJetPipeline/catalog/run22.list \
#       -s 0 -j 50 -o /star/data01/pwg/seanp/pipeline_output/run22_test \
#       -t "fcsJP2,fcsJPA0,fcsJPA1,fcsJPBC0,fcsJPBC1,fcsJPDE0,fcsJPDE1"

set -e

LISTFILE=""
SAVE_SIMPLETREE=""
OUTDIR=""
NPROC=1
LABEL="run"
KEEPLOGS=0
TRIGGER_FILTER=""
EM_ONLY=0

while getopts "f:s:o:j:l:k:t:E:h" opt; do
  case $opt in
    f) LISTFILE="$OPTARG" ;;
    s) SAVE_SIMPLETREE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    j) NPROC="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    t) TRIGGER_FILTER="$OPTARG" ;;
    E) EM_ONLY="$OPTARG" ;;
    h|*)
      sed -n '2,25p' "$0"
      exit 0
      ;;
  esac
done

if [ -z "$LISTFILE" ] || [ -z "$SAVE_SIMPLETREE" ] || [ -z "$OUTDIR" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ ! -f "$LISTFILE" ]; then
    echo "Error: listfile '$LISTFILE' does not exist"
    exit 1
fi

if [ "$EM_ONLY" != "0" ] && [ "$EM_ONLY" != "1" ]; then
    echo "Error: -E must be 0 or 1 (got '$EM_ONLY')"
    exit 1
fi

# SandBox <File> paths in data_to_jet.xml are pipeline-root-relative, so
# star-submit must run from there, not from here.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/data_to_jet_${LABEL}.xml"

# Fail fast on a typo'd trigger name, rather than inside a job with discarded logs.
if [ -n "$TRIGGER_FILTER" ]; then
    VALID_NAMES=$(awk '/kTrigFlagName\[kNTrigFlags\]/{flag=1; next} flag && /};/{flag=0} flag' "${SCRIPT_DIR}/FcsTriggerDefs.h" | sed -e 's/^[[:space:]]*"//' -e 's/",$//')
    IFS=',' read -ra REQUESTED <<< "$TRIGGER_FILTER"
    for name in "${REQUESTED[@]}"; do
        if ! grep -qxF "$name" <<< "$VALID_NAMES"; then
            echo "Error: unknown trigger filter name '$name' (not in FcsTriggerDefs.h)."
            echo "Valid names:"
            echo "$VALID_NAMES" | tr '\n' ' '
            echo
            exit 1
        fi
    done
fi

if [ "$KEEPLOGS" = "1" ]; then
    mkdir -p "${OUTDIR}/log"
    STDOUT_URL="file:${OUTDIR}/log/datatojet_\$JOBID.out"
    STDERR_URL="file:${OUTDIR}/log/datatojet_\$JOBID.err"
else
    STDOUT_URL="file:/dev/null"
    STDERR_URL="file:/dev/null"
fi

sed -e "s#{{LISTFILE}}#${LISTFILE}#g" \
    -e "s/{{SAVE_SIMPLETREE}}/${SAVE_SIMPLETREE}/g" \
    -e "s#{{OUTDIR}}#${OUTDIR}#g" \
    -e "s/{{NPROCESSES}}/${NPROC}/g" \
    -e "s#{{STDOUT_URL}}#${STDOUT_URL}#g" \
    -e "s#{{STDERR_URL}}#${STDERR_URL}#g" \
    -e "s/{{TRIGGER_FILTER}}/${TRIGGER_FILTER}/g" \
    -e "s/{{EM_ONLY}}/${EM_ONLY}/g" \
    "${SCRIPT_DIR}/data_to_jet.xml" > "${OUTXML}"

echo "Generated ${OUTXML}"
echo "listfile=${LISTFILE} save_simpletree=${SAVE_SIMPLETREE} outdir=${OUTDIR} nProcesses=${NPROC} keep_logs=${KEEPLOGS} trigger_filter=${TRIGGER_FILTER:-<none>} em_only=${EM_ONLY}"

mkdir -p "${OUTDIR}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "data_to_jet/data_to_jet_${LABEL}.xml"
SUBMIT_DIR="${PIPELINE_ROOT}/submit"
mkdir -p "${SUBMIT_DIR}"

# sched*/.csh/.package files are read from PIPELINE_ROOT by their original
# path at job start, not transferred by condor -- sweeping immediately races
# job startup and kills jobs. Defer the sweep until condor_wait confirms
# this submission's own jobs are done, in the background.
REPORT=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*.report' -newer "${MARKER}")
REQID=$(basename "${REPORT}" .report | sed 's/^sched//')
# Read the real Log path from the .condor file rather than assuming
# /tmp/$USER/... -- a wrong guess makes condor_wait fail immediately, so
# also gate the sweep on its exit status.
CONDORFILE=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name "sched${REQID}_*.condor" -newer "${MARKER}" | head -1)
CONDORLOG=$(grep -m1 '^Log' "${CONDORFILE}" | sed -e 's/^Log[[:space:]]*=[[:space:]]*//')
nohup bash -c "if condor_wait '${CONDORLOG}' >/dev/null 2>&1; then find '${PIPELINE_ROOT}' -maxdepth 1 -name 'sched${REQID}*' -exec mv -t '${SUBMIT_DIR}' {} +; else echo \"condor_wait failed for ${REQID} (log='${CONDORLOG}') -- NOT sweeping, leaving sched files in place\" >> '${PIPELINE_ROOT}/submit/.sweep_failures.log'; fi" >/dev/null 2>&1 &
disown
rm -f "${MARKER}"
