#!/bin/bash
# Generates and submits a concrete jet_calibration.xml job.
# Usage:
#   ./submit_jet_calibration.sh -i <jettrees_dir> -c <calib_json> -j <nprocesses> [-l <label>] [-k <0|1>]
#
#   -i  absolute path to a folder of jet_output_*.root (JetTrees) files;
#       every file in it is modified in place (no undo)
#   -c  absolute path to the calibration JSON lookup table
#   -j  number of parallel jobs to evenly split the folder's files across
#   -l  short label used only in the generated XML's filename, default "run"
#   -k  1 to save each job's stdout/stderr under <jettrees_dir>/log/
#       (auto-created; there's no separate output dir for this in-place
#       pipeline), 0 to discard them, default 0
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

while getopts "i:c:j:l:k:h" opt; do
  case $opt in
    i) JETTREES_DIR="$OPTARG" ;;
    c) CALIB_JSON="$OPTARG" ;;
    j) NPROCESSES="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    h|*)
      sed -n '2,16p' "$0"
      exit 0
      ;;
  esac
done

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

# This script lives in jet_calibration/. jet_calibration.xml's SandBox
# <File> paths are relative to the pipeline root (SUMS packages each
# <File> preserving that relative path -- confirmed directly from a real
# job's .package.zip), so star-submit MUST run from there, not from a
# dedicated submit folder. Instead, the sched*/*.package/*.csh/*.list/
# *.condor files star-submit drops in the CWD it runs from are swept into
# submit/ (a dedicated folder) right after submission -- see the
# marker-file logic below.
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
find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*' -newer "${MARKER}" -exec mv -t "${SUBMIT_DIR}" {} +
rm -f "${MARKER}"
