#!/bin/bash
# Generates and submits a concrete jet_scale_em_only.xml job (the EM-only
# counterpart of submit_jet_scale.sh (same folder); there is no -E flag --
# EM-only is hardwired).
# Usage:
#   ./submit_jet_scale_em_only.sh -i <jettrees_dir> -o <outdir> [-n <json_name>] [-s <0|1>] [-l <label>] [-k <0|1>]
#
#   -i  absolute path to a folder of EM-only jet_output_*.root (JetTrees)
#       files, i.e. built with submit_sim_to_jet.sh -E 1
#   -o  output directory (created if missing)
#   -n  output JSON lookup table filename, default JetEnergyScale_lookup_em_only.json
#   -s  1 = keep intermediate jet_calibration.root, 0 = discard, default 0
#   -l  short label for the generated XML's filename, default "run"
#   -k  1 = save stdout/stderr under <outdir>/log/, 0 = discard, default 0
#
# Example:
#   ./submit_jet_scale_em_only.sh -i /star/data01/pwg/seanp/pipeline_output/em_only_production \
#       -o /star/data01/pwg/seanp/em_only_cal \
#       -n JetEnergyScale_lookup_em_only.json

set -e

JETTREES_DIR=""
OUTDIR=""
JSON_NAME="JetEnergyScale_lookup_em_only.json"
SAVE_CALIBRATION_ROOT="0"
LABEL="run"
KEEPLOGS=0

while getopts "i:o:n:s:l:k:h" opt; do
  case $opt in
    i) JETTREES_DIR="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    n) JSON_NAME="$OPTARG" ;;
    s) SAVE_CALIBRATION_ROOT="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    h|*)
      sed -n '2,19p' "$0"
      exit 0
      ;;
  esac
done

if [ -z "$JETTREES_DIR" ] || [ -z "$OUTDIR" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ ! -d "$JETTREES_DIR" ]; then
    echo "Error: JetTrees directory '$JETTREES_DIR' does not exist"
    exit 1
fi

NROOTFILES=$(ls "$JETTREES_DIR"/*.root 2>/dev/null | wc -l)
if [ "$NROOTFILES" -eq 0 ]; then
    echo "Error: no .root files found in '$JETTREES_DIR'"
    exit 1
fi
echo "Found ${NROOTFILES} .root files in ${JETTREES_DIR}"

# SandBox <File> paths in jet_scale_em_only.xml are pipeline-root-relative, so
# star-submit must run from there, not from here.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/jet_scale_em_only_${LABEL}.xml"

if [ "$KEEPLOGS" = "1" ]; then
    mkdir -p "${OUTDIR}/log"
    STDOUT_URL="file:${OUTDIR}/log/jetscale_em_only_\$JOBID.out"
    STDERR_URL="file:${OUTDIR}/log/jetscale_em_only_\$JOBID.err"
else
    STDOUT_URL="file:/dev/null"
    STDERR_URL="file:/dev/null"
fi

sed -e "s#{{JETTREES_DIR}}#${JETTREES_DIR}#g" \
    -e "s#{{OUTDIR}}#${OUTDIR}#g" \
    -e "s/{{JSON_NAME}}/${JSON_NAME}/g" \
    -e "s/{{SAVE_CALIBRATION_ROOT}}/${SAVE_CALIBRATION_ROOT}/g" \
    -e "s#{{STDOUT_URL}}#${STDOUT_URL}#g" \
    -e "s#{{STDERR_URL}}#${STDERR_URL}#g" \
    "${SCRIPT_DIR}/jet_scale_em_only.xml" > "${OUTXML}"

echo "Generated ${OUTXML}"
echo "jettrees_dir=${JETTREES_DIR} outdir=${OUTDIR} json_name=${JSON_NAME} save_calibration_root=${SAVE_CALIBRATION_ROOT} keep_logs=${KEEPLOGS}"

mkdir -p "${OUTDIR}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "jet_scale/jet_scale_em_only_${LABEL}.xml"
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
