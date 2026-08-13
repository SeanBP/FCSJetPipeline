#!/bin/bash
# Generates and submits a concrete data_to_jet.xml job.
# Usage:
#   ./submit_data_to_jet.sh -f <listfile> -s <0|1> -o <outdir> \
#       [-j <nprocesses>] [-l <label>] [-k <0|1>] [-t <triggers>]
#
#   -f  absolute path to a file catalog list (one file per line, in the
#       "path,filename,events" format MakeRunList.pl writes: an xrootd URL
#       or local path, with the event count appended as a trailing
#       event-count number)
#   -s  1 to keep the intermediate SimpleTree_mudst files, 0 to discard them
#   -o  output directory (created if missing)
#   -j  number of parallel jobs (nProcesses); the listfile's lines are
#       split evenly across them, default 1
#   -l  short label used only in the generated XML's filename, default "run"
#   -k  1 to save each job's stdout/stderr under <outdir>/log/ (auto-
#       created), 0 to discard them, default 0
#   -t  comma-separated FCS trigger flag names (see FcsTriggerDefs.h /
#       TriggerIDs.txt, e.g. "fcsJP2,fcsJPA0,fcsJPA1,fcsJPBC0,fcsJPBC1,
#       fcsJPDE0,fcsJPDE1" for jet-patch triggers) -- an event only gets a
#       jetTree entry if at least one of them fired. Omit/empty (default)
#       keeps every event; the SimpleTree itself is never filtered.
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

while getopts "f:s:o:j:l:k:t:h" opt; do
  case $opt in
    f) LISTFILE="$OPTARG" ;;
    s) SAVE_SIMPLETREE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    j) NPROC="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    t) TRIGGER_FILTER="$OPTARG" ;;
    h|*)
      sed -n '2,27p' "$0"
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

# This script lives in data_to_jet/. data_to_jet.xml's SandBox <File>
# paths are relative to the pipeline root (SUMS packages each <File>
# preserving that relative path, e.g. data_to_jet/runMudst.C -- confirmed
# directly from a real job's .package.zip), so star-submit MUST run from
# there, not from a dedicated submit folder (an absolute SandBox path
# would repackage everything under the wrong subpath and break the
# `cp data_to_jet/* .` step in data_to_jet.xml). Instead, the
# sched*/*.package/*.csh/*.list/*.condor files star-submit drops in the
# CWD it runs from are swept into submit/ (a dedicated folder) right
# after submission -- see the marker-file logic below (this used to dump
# ~9000 files/~100MB directly into FCSJetPipeline/ over a handful of
# submissions before that cleanup step existed).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/data_to_jet_${LABEL}.xml"

# Fail fast on a typo'd trigger name here, rather than have it silently
# error out inside a job whose logs are discarded by default (-k 0).
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
    "${SCRIPT_DIR}/data_to_jet.xml" > "${OUTXML}"

echo "Generated ${OUTXML}"
echo "listfile=${LISTFILE} save_simpletree=${SAVE_SIMPLETREE} outdir=${OUTDIR} nProcesses=${NPROC} keep_logs=${KEEPLOGS} trigger_filter=${TRIGGER_FILTER:-<none>}"

mkdir -p "${OUTDIR}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "data_to_jet/data_to_jet_${LABEL}.xml"
SUBMIT_DIR="${PIPELINE_ROOT}/submit"
mkdir -p "${SUBMIT_DIR}"
find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*' -newer "${MARKER}" -exec mv -t "${SUBMIT_DIR}" {} +
rm -f "${MARKER}"
