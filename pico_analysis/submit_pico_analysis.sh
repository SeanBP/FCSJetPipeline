#!/bin/bash
# Generates and submits a concrete pico_analysis.xml job.
# Usage:
#   ./submit_pico_analysis.sh -d <mip|hits> -f <listfile> -o <outdir> \
#       [-l <label>] [-k <0|1>]
#
#   -d  analysis mode:
#         mip   MIP peak analysis (mip_ana.C) -- projects tracks onto the
#               FCS ECal/HCal, writes diagnostic plots to mip_ana.pdf
#         hits  flat hit-level tree (pico_to_root.C) -- FCS hit energy/
#               position/detector ID + vertex position per event, no
#               clustering, no jets
#   -f  absolute path to a text file listing PicoDst files, one per line
#   -o  output directory (created if missing)
#   -l  short label used only in the generated XML's filename, default "run"
#   -k  1 to save stdout/stderr under <outdir>/log/, 0 to discard them,
#       default 0
#
# Example:
#   ./submit_pico_analysis.sh -d mip -f infiles_27082025.list \
#       -o /star/data01/pwg/seanp/pipeline_output/mip_test

set -e

MODE=""
LISTFILE=""
OUTDIR=""
LABEL="run"
KEEPLOGS=0

while getopts "d:f:o:l:k:h" opt; do
  case $opt in
    d) MODE="$OPTARG" ;;
    f) LISTFILE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    h|*)
      sed -n '2,21p' "$0"
      exit 0
      ;;
  esac
done

if [ -z "$MODE" ] || [ -z "$LISTFILE" ] || [ -z "$OUTDIR" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ "$MODE" != "mip" ] && [ "$MODE" != "hits" ]; then
    echo "Error: -d must be 'mip' or 'hits' (got '$MODE')"
    exit 1
fi

if [ ! -f "$LISTFILE" ]; then
    echo "Error: listfile '$LISTFILE' does not exist"
    exit 1
fi

# This script lives in pico_analysis/. pico_analysis.xml's SandBox <File>
# paths are relative to the pipeline root (SUMS packages each <File>
# preserving that relative path -- confirmed directly from a real job's
# .package.zip), so star-submit MUST run from there, not from a dedicated
# submit folder. Instead, the sched*/*.package/*.csh/*.list/*.condor
# files star-submit drops in the CWD it runs from are swept into submit/
# (a dedicated folder) right after submission -- see the marker-file
# logic below.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/pico_analysis_${LABEL}.xml"

if [ "$KEEPLOGS" = "1" ]; then
    mkdir -p "${OUTDIR}/log"
    STDOUT_URL="file:${OUTDIR}/log/pico_\$JOBID.out"
    STDERR_URL="file:${OUTDIR}/log/pico_\$JOBID.err"
else
    STDOUT_URL="file:/dev/null"
    STDERR_URL="file:/dev/null"
fi

sed -e "s/{{MODE}}/${MODE}/g" \
    -e "s#{{LISTFILE}}#${LISTFILE}#g" \
    -e "s#{{OUTDIR}}#${OUTDIR}#g" \
    -e "s#{{STDOUT_URL}}#${STDOUT_URL}#g" \
    -e "s#{{STDERR_URL}}#${STDERR_URL}#g" \
    "${SCRIPT_DIR}/pico_analysis.xml" > "${OUTXML}"

echo "Generated ${OUTXML}"
echo "mode=${MODE} listfile=${LISTFILE} outdir=${OUTDIR} keep_logs=${KEEPLOGS}"

mkdir -p "${OUTDIR}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "pico_analysis/pico_analysis_${LABEL}.xml"
SUBMIT_DIR="${PIPELINE_ROOT}/submit"
mkdir -p "${SUBMIT_DIR}"
find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*' -newer "${MARKER}" -exec mv -t "${SUBMIT_DIR}" {} +
rm -f "${MARKER}"
