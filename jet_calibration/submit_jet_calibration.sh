#!/bin/bash
# Generates and submits a concrete jet_calibration.xml job.
# Usage:
#   ./submit_jet_calibration.sh -i <jettrees_dir> -c <calib_json> -j <nprocesses> [-l <label>] [-k <0|1>] [-E <0|1>]
#
#   -i  absolute path to a folder of jet_output_*.root (JetTrees) files;
#       every file in it is modified in place (no undo)
#   -c  absolute path to the calibration JSON lookup table
#   -j  number of parallel jobs to evenly split the folder's files across
#   -l  short label for the generated XML's filename, default "run"
#   -k  1 = save stdout/stderr under <jettrees_dir>/log/, 0 = discard, default 0
#   -E  1 = the JetTrees are EM-only and -c is an EM-only JSON (from
#       submit_jet_scale_em_only.sh / submit_jet_scale.sh -E 1): hand off to
#       submit_jet_calibration_em_only.sh (same folder) with the same
#       arguments. 0/omit = the ECAL+HCAL correction below, unchanged.
#
# Example:
#   ./submit_jet_calibration.sh -i /star/data01/pwg/seanp/tunes/pythia8_mstw2008lo/JetTrees \
#       -c /star/data01/pwg/seanp/tunes/pythia8_mstw2008lo/JetEnergyScale_lookup_pythia8_mstw2008lo.json \
#       -j 10

set -e

JETTREES_DIR=""
CALIB_JSON=""
NPROCESSES=""
LABEL="run"
KEEPLOGS=0
EM_ONLY=0
# every option except -E, forwarded verbatim if -E 1 hands off to the
# EM-only script
FORWARD=()

while getopts "i:c:j:l:k:E:h" opt; do
  case $opt in
    i) JETTREES_DIR="$OPTARG"; FORWARD+=(-i "$OPTARG") ;;
    c) CALIB_JSON="$OPTARG"; FORWARD+=(-c "$OPTARG") ;;
    j) NPROCESSES="$OPTARG"; FORWARD+=(-j "$OPTARG") ;;
    l) LABEL="$OPTARG"; FORWARD+=(-l "$OPTARG") ;;
    k) KEEPLOGS="$OPTARG"; FORWARD+=(-k "$OPTARG") ;;
    E) EM_ONLY="$OPTARG" ;;
    h|*)
      sed -n '2,21p' "$0"
      exit 0
      ;;
  esac
done

if [ "$EM_ONLY" != "0" ] && [ "$EM_ONLY" != "1" ]; then
    echo "Error: -E must be 0 or 1 (got '$EM_ONLY')"
    exit 1
fi

if [ "$EM_ONLY" = "1" ]; then
    exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_jet_calibration_em_only.sh" "${FORWARD[@]}"
fi

if [ -z "$JETTREES_DIR" ] || [ -z "$CALIB_JSON" ] || [ -z "$NPROCESSES" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ ! -d "$JETTREES_DIR" ]; then
    echo "Error: JetTrees directory '$JETTREES_DIR' does not exist"
    exit 1
fi

if [ ! -f "$CALIB_JSON" ]; then
    echo "Error: calibration JSON '$CALIB_JSON' does not exist"
    exit 1
fi

NROOTFILES=$(ls "$JETTREES_DIR"/*.root 2>/dev/null | wc -l)
if [ "$NROOTFILES" -eq 0 ]; then
    echo "Error: no .root files found in '$JETTREES_DIR'"
    exit 1
fi
echo "Found ${NROOTFILES} .root files in ${JETTREES_DIR}"

if [ "$NPROCESSES" -gt "$NROOTFILES" ]; then
    echo "Warning: -j ${NPROCESSES} exceeds the file count (${NROOTFILES}); some jobs will have nothing to do"
fi

echo "This will modify ${NROOTFILES} files in place under ${JETTREES_DIR}. There is no undo."
read -p "Continue? (y/N): " CONFIRM
if [ "$CONFIRM" != "y" ] && [ "$CONFIRM" != "Y" ]; then
    echo "Aborted."
    exit 1
fi

# SandBox <File> paths in jet_calibration.xml are pipeline-root-relative,
# so star-submit must run from there, not from here.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/jet_calibration_${LABEL}.xml"

if [ "$KEEPLOGS" = "1" ]; then
    mkdir -p "${JETTREES_DIR}/log"
    STDOUT_URL="file:${JETTREES_DIR}/log/jetcal_\$JOBID.out"
    STDERR_URL="file:${JETTREES_DIR}/log/jetcal_\$JOBID.err"
else
    STDOUT_URL="file:/dev/null"
    STDERR_URL="file:/dev/null"
fi

sed -e "s#{{JETTREES_DIR}}#${JETTREES_DIR}#g" \
    -e "s#{{CALIB_JSON}}#${CALIB_JSON}#g" \
    -e "s/{{NPROCESSES}}/${NPROCESSES}/g" \
    -e "s#{{STDOUT_URL}}#${STDOUT_URL}#g" \
    -e "s#{{STDERR_URL}}#${STDERR_URL}#g" \
    "${SCRIPT_DIR}/jet_calibration.xml" > "${OUTXML}"

echo "Generated ${OUTXML}"
echo "jettrees_dir=${JETTREES_DIR} calib_json=${CALIB_JSON} nProcesses=${NPROCESSES} keep_logs=${KEEPLOGS}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "jet_calibration/jet_calibration_${LABEL}.xml"
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
