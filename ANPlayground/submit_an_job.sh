#!/bin/bash
# Generates and submits an_job.xml: MuDst -> A_N yield histograms, with a
# per-stage timing table per job.
# Usage:
#   ./submit_an_job.sh -f <listfile> -o <outdir> [-n <files_per_job>] \
#       [-c <calib json>] [-t <triggers>] [-l <label>] [-k <0|1>] [-N]
#
#   -f  MakeRunList.pl-format list ("<url or path>/<events>" per line)
#   -o  output directory; histograms go to <outdir>/hists/, timing tables
#       to <outdir>/timing/
#   -n  files per job, default 40: each run is cut into jobs of at most this
#       many consecutive files (the underlying-event step draws proxies from
#       its own job's files of one run). Also sets the job's scratch use,
#       ~0.1 GB per file of SimpleTree plus one MuDst (~1.2 GB).
#   -c  EM-only point-jet JES lookup JSON, default
#       /star/data01/pwg/seanp/jes_points/jesmap_v3/JetEnergyScale_lookup_points_ue_v3_lognormal.json
#       (copied into ANPlayground/jes/ for the SandBox). It must be a
#       "logcubic" map; the default is of the earlier log-normal form and is
#       refused until a UE-subtracted log-cubic map replaces it
#   -t  FCS trigger names an event must have fired, default
#       fcsEM0,fcsEM1,fcsEM2 (the EM levels the UE step and the histograms use)
#   -l  label for the generated XML, default "an"
#   -k  1 = save stdout/stderr under <outdir>/log/, default 0
#   -N  write the XML and job layout, do not submit

set -e

LISTFILE=""
OUTDIR=""
FILES_PER_JOB=40
CALIB=/star/data01/pwg/seanp/jes_points/jesmap_v3/JetEnergyScale_lookup_points_ue_v3_lognormal.json
TRIGGER_FILTER="fcsEM0,fcsEM1,fcsEM2"
LABEL="an"
KEEPLOGS=0
DRYRUN=0

while getopts "f:o:n:c:t:l:k:Nh" opt; do
  case $opt in
    f) LISTFILE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    n) FILES_PER_JOB="$OPTARG" ;;
    c) CALIB="$OPTARG" ;;
    t) TRIGGER_FILTER="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    N) DRYRUN=1 ;;
    h|*) sed -n '2,27p' "$0"; exit 0 ;;
  esac
done

if [ -z "$LISTFILE" ] || [ -z "$OUTDIR" ]; then
    echo "Missing -f or -o. Run with -h for usage."; exit 1
fi
[ -f "$LISTFILE" ] || { echo "Error: listfile '$LISTFILE' does not exist"; exit 1; }
[ -f "$CALIB" ] || { echo "Error: calibration JSON '$CALIB' does not exist"; exit 1; }
grep -q '"model": *"logcubic"' "$CALIB" || { echo "Error: $CALIB is not a log-cubic EM-only map (CalibrationMap_em_only); pass one with -c"; exit 1; }
[[ "$FILES_PER_JOB" =~ ^[1-9][0-9]*$ ]] || { echo "Error: -n must be a positive integer"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/an_job_${LABEL}.xml"

# Trigger names must exist in FcsTriggerDefs.h
VALID_NAMES=$(awk '/kTrigFlagName\[kNTrigFlags\]/{flag=1; next} flag && /};/{flag=0} flag' "${PIPELINE_ROOT}/data_to_jet/FcsTriggerDefs.h" | sed -e 's/^[[:space:]]*"//' -e 's/",$//')
IFS=',' read -ra REQUESTED <<< "$TRIGGER_FILTER"
for name in "${REQUESTED[@]}"; do
    grep -qxF "$name" <<< "$VALID_NAMES" || { echo "Error: unknown trigger name '$name'"; exit 1; }
done

mkdir -p "${OUTDIR}/hists" "${OUTDIR}/timing"
CALIB_NAME=$(basename "$CALIB")
cp -p "$CALIB" "${SCRIPT_DIR}/jes/${CALIB_NAME}"

# Job layout as submit_data_to_jet.sh: order by run, then by the file's
# sequence number; cut each run into jobs of at most -n consecutive files
SORTEDLIST="${OUTDIR}/filelist_${LABEL}.sorted.list"
RANGEFILE="${OUTDIR}/filelist_${LABEL}.ranges"
NAME_RE='_[0-9]{8}_(raw_)?[0-9]{7}\.MuDst\.root'
NBADNAME=$(grep -vcE "${NAME_RE}" "${LISTFILE}" || true)
[ "$NBADNAME" = "0" ] || { echo "Error: ${NBADNAME} line(s) have no <run>_[raw_]<file>.MuDst.root name"; exit 1; }
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
NPROC=$(wc -l < "${RANGEFILE}")
echo "Job layout: $(wc -l < "${SORTEDLIST}") files in ${NPROC} jobs of at most ${FILES_PER_JOB} consecutive files of one run"

if [ "$KEEPLOGS" = "1" ]; then
    mkdir -p "${OUTDIR}/log"
    STDOUT_URL="file:${OUTDIR}/log/an_\$JOBID.out"
    STDERR_URL="file:${OUTDIR}/log/an_\$JOBID.err"
else
    STDOUT_URL="file:/dev/null"
    STDERR_URL="file:/dev/null"
fi

sed -e "s#{{LISTFILE}}#${SORTEDLIST}#g" \
    -e "s#{{RANGEFILE}}#${RANGEFILE}#g" \
    -e "s/{{NPROCESSES}}/${NPROC}/g" \
    -e "s/{{TRIGGER_FILTER}}/${TRIGGER_FILTER}/g" \
    -e "s/{{CALIB_NAME}}/${CALIB_NAME}/g" \
    -e "s#{{OUTDIR}}#${OUTDIR}#g" \
    -e "s#{{STDOUT_URL}}#${STDOUT_URL}#g" \
    -e "s#{{STDERR_URL}}#${STDERR_URL}#g" \
    "${SCRIPT_DIR}/an_job.xml" > "${OUTXML}"
echo "Generated ${OUTXML}"

cat > "${OUTDIR}/run_params_${LABEL}.txt" <<EOF
# ANPlayground an_job run parameters
submitted_utc=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
label=${LABEL}
pipeline_commit=$(git -C "${PIPELINE_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)
listfile=${LISTFILE}
files=$(wc -l < "${SORTEDLIST}")
nProcesses=${NPROC}
files_per_job=${FILES_PER_JOB}
calib_json=${CALIB}
trigger_filter=${TRIGGER_FILTER}
keep_logs=${KEEPLOGS}
EOF

[ "$DRYRUN" = "1" ] && { echo "Dry run (-N): not submitted"; exit 0; }

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "ANPlayground/an_job_${LABEL}.xml"
SUBMIT_DIR="${PIPELINE_ROOT}/submit"
mkdir -p "${SUBMIT_DIR}"
REPORT=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*.report' -newer "${MARKER}")
REQID=$(basename "${REPORT}" .report | sed 's/^sched//')
echo "request_id=${REQID}" >> "${OUTDIR}/run_params_${LABEL}.txt"

# Move the sched* files to submit/ once condor_wait reports the request done
CONDORFILE=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name "sched${REQID}_*.condor" -newer "${MARKER}" | head -1)
CONDORLOG=$(grep -m1 '^Log' "${CONDORFILE}" | sed -e 's/^Log[[:space:]]*=[[:space:]]*//')
nohup bash -c "if condor_wait '${CONDORLOG}' >/dev/null 2>&1; then find '${PIPELINE_ROOT}' -maxdepth 1 -name 'sched${REQID}*' -exec mv -t '${SUBMIT_DIR}' {} +; fi" >/dev/null 2>&1 &
disown
rm -f "${MARKER}"
echo "Submitted request ${REQID}; condor log ${CONDORLOG}"
