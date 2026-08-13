#!/bin/bash
# Generates and submits a concrete jet_scale.xml job.
# Usage:
#   ./submit_jet_scale.sh -i <jettrees_dir> -o <outdir> [-n <json_name>] [-s <0|1>] [-l <label>] [-k <0|1>]
#
#   -i  absolute path to a folder of jet_output_*.root (JetTrees) files
#   -o  output directory (created if missing)
#   -n  filename for the output JSON lookup table, default
#       JetEnergyScale_lookup.json
#   -s  1 to keep the intermediate jet_calibration.root, 0 to discard it
#       (default 0)
#   -l  short label used only in the generated XML's filename, default "run"
#   -k  1 to save the job's stdout/stderr under <outdir>/log/ (auto-
#       created), 0 to discard them, default 0
#
# Example:
#   ./submit_jet_scale.sh -i /star/data01/pwg/seanp/tunes/pythia8_mstw2008lo/JetTrees \
#       -o /star/data01/pwg/seanp/tunes/pythia8_mstw2008lo \
#       -n JetEnergyScale_lookup_pythia8_mstw2008lo.json

set -e

JETTREES_DIR=""
OUTDIR=""
JSON_NAME="JetEnergyScale_lookup.json"
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
      sed -n '2,16p' "$0"
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

# This script lives in jet_scale/. jet_scale.xml's SandBox <File> paths
# are relative to the pipeline root (SUMS packages each <File> preserving
# that relative path -- confirmed directly from a real job's .package.zip),
# so star-submit MUST run from there, not from a dedicated submit folder.
# Instead, the sched*/*.package/*.csh/*.list/*.condor files star-submit
# drops in the CWD it runs from are swept into submit/ (a dedicated
# folder) right after submission -- see the marker-file logic below.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/jet_scale_${LABEL}.xml"

if [ "$KEEPLOGS" = "1" ]; then
    mkdir -p "${OUTDIR}/log"
    STDOUT_URL="file:${OUTDIR}/log/jetscale_\$JOBID.out"
    STDERR_URL="file:${OUTDIR}/log/jetscale_\$JOBID.err"
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
    "${SCRIPT_DIR}/jet_scale.xml" > "${OUTXML}"

echo "Generated ${OUTXML}"
echo "jettrees_dir=${JETTREES_DIR} outdir=${OUTDIR} json_name=${JSON_NAME} save_calibration_root=${SAVE_CALIBRATION_ROOT} keep_logs=${KEEPLOGS}"

mkdir -p "${OUTDIR}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "jet_scale/jet_scale_${LABEL}.xml"
SUBMIT_DIR="${PIPELINE_ROOT}/submit"
mkdir -p "${SUBMIT_DIR}"
find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*' -newer "${MARKER}" -exec mv -t "${SUBMIT_DIR}" {} +
rm -f "${MARKER}"
