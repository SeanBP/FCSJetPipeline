#!/bin/bash
# Generates and submits a concrete data_to_jet.xml job.
# Usage:
#   ./submit_data_to_jet.sh -f <listfile> -s <0|1> -o <outdir> \
#       [-n <files_per_job>] [-j <nprocesses>] [-l <label>] [-k <0|1>] \
#       [-t <triggers>] [-E <0|1>] [-F <0|1>]
#
#   -f  absolute path to a file catalog list (MakeRunList.pl's
#       "path,filename,events" format: xrootd URL or local path, event
#       count appended as a trailing number)
#   -s  1 = keep intermediate SimpleTree_mudst files, 0 = discard
#   -o  output directory (created if missing)
#   -n  files per job, default 40. The list is reordered by run and by
#       time within a run, and each run is cut into jobs of at most this
#       many consecutive files (evenly sized within the run). The
#       underlying-event step draws a jet's proxy events only from its own
#       job's files and only from its own run, so a job should hold as long
#       a stretch of one run as its running time allows.
#   -j  number of parallel jobs, with the listfile lines split evenly across
#       them in the order given. Overrides -n; jobs then mix runs and times,
#       which the underlying-event step tolerates but does worse with.
#   -l  short label for the generated XML's filename, default "run"
#   -k  1 = save stdout/stderr under <outdir>/log/, 0 = discard, default 0
#   -t  comma-separated FCS trigger flag names (see FcsTriggerDefs.h); an
#       event is only written (to the SimpleTree and so to the jetTree) if
#       one fired. Empty = keep every event.
#   -E  1 = reco jets from ECAL only, with the ECAL-only fiducial boundary
#       (kFiducialRectEcal). 0/omit = ECAL+HCAL.
#   -F  1 = reco jets from FCS points (StMuFcsPoint) instead of tower hits.
#       HCAL has no points, so this is ECAL-only input. 0/omit = hits.
#
# Example:
#   ./submit_data_to_jet.sh -f /star/u/seanp/FCSJetPipeline/catalog/run22.list \
#       -s 0 -j 50 -o /star/data01/pwg/seanp/pipeline_output/run22_test \
#       -t "fcsJP2,fcsJPA0,fcsJPA1,fcsJPBC0,fcsJPBC1,fcsJPDE0,fcsJPDE1"

set -e

LISTFILE=""
SAVE_SIMPLETREE=""
OUTDIR=""
NPROC=""
FILES_PER_JOB=40
LABEL="run"
KEEPLOGS=0
TRIGGER_FILTER=""
EM_ONLY=0
USE_POINTS=0

while getopts "f:s:o:j:n:l:k:t:E:F:h" opt; do
  case $opt in
    f) LISTFILE="$OPTARG" ;;
    s) SAVE_SIMPLETREE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    j) NPROC="$OPTARG" ;;
    n) FILES_PER_JOB="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    t) TRIGGER_FILTER="$OPTARG" ;;
    E) EM_ONLY="$OPTARG" ;;
    F) USE_POINTS="$OPTARG" ;;
    h|*)
      sed -n '2,35p' "$0"
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

if [ "$USE_POINTS" != "0" ] && [ "$USE_POINTS" != "1" ]; then
    echo "Error: -F must be 0 or 1 (got '$USE_POINTS')"
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

# Job layout. Without -j: reorder the list by run, then by the file's
# sequence number (the last five digits of the 7-digit file number, which
# follows time), and give each job consecutive files of one run.
RANGEFILE=""
if [ -z "$NPROC" ]; then
    if ! [[ "$FILES_PER_JOB" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: -n must be a positive integer (got '$FILES_PER_JOB')"
        exit 1
    fi
    mkdir -p "${OUTDIR}"
    SORTEDLIST="${OUTDIR}/filelist_${LABEL}.sorted.list"
    RANGEFILE="${OUTDIR}/filelist_${LABEL}.ranges"
    NAME_RE='_[0-9]{8}_(raw_)?[0-9]{7}\.MuDst\.root'
    NBADNAME=$(grep -vcE "${NAME_RE}" "${LISTFILE}" || true)
    if [ "$NBADNAME" != "0" ]; then
        echo "Error: ${NBADNAME} line(s) of ${LISTFILE} have no <run>_[raw_]<file>.MuDst.root name; cannot order them by run and time. Use -j for the plain even split."
        exit 1
    fi
    sed -E 's#^(.*_([0-9]{8})_(raw_)?([0-9]{2})([0-9]{5})\.MuDst\.root.*)$#\2 \5 \4 \1#' "${LISTFILE}" \
        | sort -k1,1n -k2,2n -k3,3n | cut -d' ' -f4- > "${SORTEDLIST}"
    awk -v N="${FILES_PER_JOB}" '
        { match($0, /_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_(raw_)?[0-9][0-9][0-9][0-9][0-9][0-9][0-9]\.MuDst\.root/)
          run = substr($0, RSTART + 1, 8)
          if (!(run in first)) { first[run] = NR; order[++k] = run }
          cnt[run]++ }
        END { for (i = 1; i <= k; i++) {
                  r = order[i]; n = cnt[r]; c = int((n + N - 1) / N); base = int(n / c); rem = n % c; s = first[r]
                  for (j = 0; j < c; j++) { m = base + (j < rem ? 1 : 0); print s, m; s += m } } }
    ' "${SORTEDLIST}" > "${RANGEFILE}"
    LISTFILE="${SORTEDLIST}"
    NPROC=$(wc -l < "${RANGEFILE}")
    echo "Job layout: $(wc -l < "${SORTEDLIST}") files in ${NPROC} jobs of at most ${FILES_PER_JOB} consecutive files of one run (${RANGEFILE})"
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
    -e "s#{{RANGEFILE}}#${RANGEFILE}#g" \
    -e "s#{{STDOUT_URL}}#${STDOUT_URL}#g" \
    -e "s#{{STDERR_URL}}#${STDERR_URL}#g" \
    -e "s/{{TRIGGER_FILTER}}/${TRIGGER_FILTER}/g" \
    -e "s/{{EM_ONLY}}/${EM_ONLY}/g" \
    -e "s/{{USE_POINTS}}/${USE_POINTS}/g" \
    "${SCRIPT_DIR}/data_to_jet.xml" > "${OUTXML}"

echo "Generated ${OUTXML}"
echo "listfile=${LISTFILE} save_simpletree=${SAVE_SIMPLETREE} outdir=${OUTDIR} nProcesses=${NPROC} keep_logs=${KEEPLOGS} trigger_filter=${TRIGGER_FILTER:-<none>} em_only=${EM_ONLY} use_points=${USE_POINTS}"

mkdir -p "${OUTDIR}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "data_to_jet/data_to_jet_${LABEL}.xml"
SUBMIT_DIR="${PIPELINE_ROOT}/submit"
mkdir -p "${SUBMIT_DIR}"

# Jobs read the sched* files from PIPELINE_ROOT at start-up, so they are moved
# to submit/ only after condor_wait reports the whole submission finished.
REPORT=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*.report' -newer "${MARKER}")
REQID=$(basename "${REPORT}" .report | sed 's/^sched//')

# Record every resolved setting next to the output
GITHASH=$(git -C "${PIPELINE_ROOT}" rev-parse HEAD 2>/dev/null || echo "unknown")
cat > "${OUTDIR}/run_params_${LABEL}_${REQID}.txt" <<EOF
# data_to_jet run parameters -- see FCSJetPipeline/README.md
submitted_utc=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
label=${LABEL}
request_id=${REQID}
pipeline_commit=${GITHASH}
generated_xml=${OUTXML}
listfile=${LISTFILE}
listfile_lines=$(wc -l < "${LISTFILE}")
save_simpletree=${SAVE_SIMPLETREE}
outdir=${OUTDIR}
nProcesses=${NPROC}
files_per_job=$( [ -n "${RANGEFILE}" ] && echo "${FILES_PER_JOB}" || echo "<even split, -j>" )
keep_logs=${KEEPLOGS}
trigger_filter=${TRIGGER_FILTER:-<none>}
em_only=${EM_ONLY}
use_points=${USE_POINTS}
EOF
echo "Wrote run parameters to ${OUTDIR}/run_params_${LABEL}_${REQID}.txt"

CONDORFILE=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name "sched${REQID}_*.condor" -newer "${MARKER}" | head -1)
CONDORLOG=$(grep -m1 '^Log' "${CONDORFILE}" | sed -e 's/^Log[[:space:]]*=[[:space:]]*//')
nohup bash -c "if condor_wait '${CONDORLOG}' >/dev/null 2>&1; then find '${PIPELINE_ROOT}' -maxdepth 1 -name 'sched${REQID}*' -exec mv -t '${SUBMIT_DIR}' {} +; else echo \"condor_wait failed for ${REQID} (log='${CONDORLOG}') -- NOT sweeping, leaving sched files in place\" >> '${PIPELINE_ROOT}/submit/.sweep_failures.log'; fi" >/dev/null 2>&1 &
disown
rm -f "${MARKER}"
