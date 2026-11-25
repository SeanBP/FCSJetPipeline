#!/bin/bash
# Generates and submits an_job.xml: MuDst -> A_N yield histograms, with a
# per-stage timing table per job.
# Usage:
#   ./submit_an_job.sh -f <listfile> -o <outdir> [-n <files_per_job>] \
#       [-c <calib json>] [-t <triggers>] [-l <label>] [-k <0|1>] [-J <N>]
#       [-S <seed>] [-s <pi0 scale table>] [-L <node fail limit>] [-R <triggers>] [-K] [-N]
#
#   -f  MakeRunList.pl-format list ("<url or path>/<events>" per line)
#   -o  output directory; histograms go to <outdir>/hists/, timing tables
#       to <outdir>/timing/
#   -n  files per job, default 40: each run is cut into jobs of at most this
#       many consecutive files (a job's yields are one run's, for the per-run
#       polarization). Scratch does not grow with it: each MuDst and SimpleTree
#       is deleted once used (one MuDst + its SimpleTree at a time, up to
#       ~4.5 GB for the largest MuDsts), plus ~0.3 MB of jet file per MuDst.
#   -c  EM-only point-jet JES lookup JSON, default
#       jet_calibration/JES_Lookups/JetEnergyScale_lookup_em_only.json
#       (copied into asymmetry/jes/ for the SandBox). It must be a
#       "logcubic" map, with its error vs energy (jes_error_* arrays,
#       jet_scale/jes_map_compare.py error) for the JES uncertainty copies
#   -t  FCS trigger names an event must have fired, default the triggers
#       AnHistograms counts: fcsEM0,fcsEM1,fcsEM2 and their _tpc copies (an
#       event listing only a _tpc copy is a valid EM event) and the jet
#       patches fcsJP2, fcsJPA1, fcsJPBC1, fcsJPDE1, fcsJPA0, fcsJPBC0, fcsJPDE0
#       (EM only before 2026-11-05)
#   -R  trigger reference triggers (trg_nominal), added to the filter for the
#       first processed MuDst of each job only; AnHistograms fills the trigger
#       tallies (h_trgref_*) from that file. Default: EM refs fcsDiJP, fcsDiJPAsy,
#       fcsHad0-2, fcsHHT_N/S, fcsHTOT_N/S, fcs_HcalTot; JP refs fcsDY,
#       fcsDYAsy, fcsJPsi, fcsELE2 (2026-11-17: the EM/JP filter keeps 2.6% / 1.5%
#       of them, biased). "" = none
#   -l  label for the generated XML, default "an"
#   -k  1 = save stdout/stderr under <outdir>/log/, default 0
#   -J  keep the calibrated jet files (UE-subtracted, JES applied) of every
#       N-th job in <outdir>/jets/, merged into one jets_<run>_<job>.root per
#       job (hadd in the job's scratch; the tree has event_run, event_id per
#       event): 1 = all, 0 = none (default). About 0.75 MB per MuDst with the
#       JP events and trigger values
#   -Q  1 = channel QA: also write <outdir>/occ/occ_<run>_<job>.root, the
#       per-channel hit counts and energies of the job's MuDsts
#       (channel_qa/channel_occupancy.cpp), for the hot / dead channel masks
#       and per-run QA; 0 = off (default; use 1 for full passes, Sean 2026-11-24)
#   -T  point-threshold copies (_vthr150 / _vthr200, Sean 2026-11-24): "on"
#       (default) = jets, UE and calibration again at FCS_POINT_EMIN = 1.5 and
#       2.0 GeV per SimpleTree, with the maps
#       jet_calibration/JES_Lookups/JetEnergyScale_lookup_em_only_pointthr150/200.json,
#       filled by AnHistograms --thr-copy into the same yield file; "none" = off
#   -M  hot / dead channel mask (channel_qa/channel_mask.cpp's channel_mask.txt,
#       or a filtered copy): masked hits are dropped before reclustering, from
#       the hit list and from the trigger sums, per run period. Default
#       channel_qa/channel_mask_v1.txt (Sean 2026-11-25: applied to all data
#       automatically); "none" = no mask. Runs of the list outside every mask
#       period are reported (they would run unmasked)
#   -S  shuffle the job order with this seed (condor starts jobs roughly in
#       order, so the jobs finished at any time are then a random spread over
#       the runs, not the earliest runs; 0 = run order, default)
#   -s  pi0 energy scale table (asymmetry/pi0_scale.txt format), default
#       asymmetry/pi0_scale.txt (empty: AnConfig.h kPi0EnergyScale for all runs)
#   -L  failed fetches after which a disk node is flagged and its remaining
#       files are skipped (<outdir>/nodes/flagged/), default 3; they are left for
#       a recovery pass (recovery.sh)
#   -K  no catalog node map (every file's node "unknown": no node flagging)
#   -N  write the XML and job layout, do not submit
#
#   The output directory also gets code_<label>.patch (git diff HEAD of the
#   pipeline at submission, the job's SandBox being the working tree) and
#   code_<label>_status.txt (git status, for untracked files).

set -e

LISTFILE=""
OUTDIR=""
FILES_PER_JOB=40
CALIB="$(cd "$(dirname "${BASH_SOURCE[0]}")/../jet_calibration/JES_Lookups" && pwd)/JetEnergyScale_lookup_em_only.json"
TRIGGER_FILTER="fcsEM0,fcsEM1,fcsEM2,fcsEM0_tpc,fcsEM1_tpc,fcsEM2_tpc,fcsJP2,fcsJPA1,fcsJPBC1,fcsJPDE1,fcsJPA0,fcsJPBC0,fcsJPDE0"
REF_TRIGGERS="fcsDiJP,fcsDiJPAsy,fcsHad0,fcsHad1,fcsHad2,fcsHHT_N,fcsHHT_S,fcsHTOT_N,fcsHTOT_S,fcs_HcalTot,fcsDY,fcsDYAsy,fcsJPsi,fcsELE2"
LABEL="an"
KEEPLOGS=0
KEEP_JETS=0
CHANNEL_QA=0
THR_OPT="on"
CHANNEL_MASK="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/channel_qa/channel_mask_v1.txt"
DRYRUN=0
SHUFFLE_SEED=0
PI0_TABLE=""
NODE_FAIL_LIMIT=3
NO_NODEMAP=0

while getopts "f:o:n:c:t:l:k:J:Q:M:T:S:s:L:R:KNh" opt; do
  case $opt in
    f) LISTFILE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    n) FILES_PER_JOB="$OPTARG" ;;
    c) CALIB="$OPTARG" ;;
    t) TRIGGER_FILTER="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    J) KEEP_JETS="$OPTARG" ;;
    Q) CHANNEL_QA="$OPTARG" ;;
    M) CHANNEL_MASK="$OPTARG" ;;
    T) THR_OPT="$OPTARG" ;;
    S) SHUFFLE_SEED="$OPTARG" ;;
    s) PI0_TABLE="$OPTARG" ;;
    L) NODE_FAIL_LIMIT="$OPTARG" ;;
    R) REF_TRIGGERS="$OPTARG" ;;
    K) NO_NODEMAP=1 ;;
    N) DRYRUN=1 ;;
    h|*) sed -n '2,54p' "$0"; exit 0 ;;
  esac
done

if [ -z "$LISTFILE" ] || [ -z "$OUTDIR" ]; then
    echo "Missing -f or -o. Run with -h for usage."; exit 1
fi
[ -f "$LISTFILE" ] || { echo "Error: listfile '$LISTFILE' does not exist"; exit 1; }
[ -f "$CALIB" ] || { echo "Error: calibration JSON '$CALIB' does not exist"; exit 1; }
grep -q '"model": *"logcubic"' "$CALIB" || { echo "Error: $CALIB is not a log-cubic EM-only map (CalibrationMap_em_only); pass one with -c"; exit 1; }
# the map's error vs energy (jet_scale/jes_map_compare.py error): AnHistograms'
# _vjesup / _vjesdn copies (the JES uncertainty); without it they stay empty
if ! grep -q '"jes_error_total"' "$CALIB"; then
    echo "WARNING: $CALIB has no jes_error_E / jes_error_total arrays (jes_map_compare.py error): the JES uncertainty copies (_vjesup / _vjesdn) will be empty"
fi
[[ "$FILES_PER_JOB" =~ ^[1-9][0-9]*$ ]] || { echo "Error: -n must be a positive integer"; exit 1; }
[[ "$KEEP_JETS" =~ ^[0-9]+$ ]] || { echo "Error: -J must be a non-negative integer"; exit 1; }
[[ "$SHUFFLE_SEED" =~ ^[0-9]+$ ]] || { echo "Error: -S must be a non-negative integer"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/an_job_${LABEL}.xml"

# Trigger names must exist in FcsTriggerDefs.h
VALID_NAMES=$(awk '/kTrigFlagName\[kNTrigFlags\]/{flag=1; next} flag && /};/{flag=0} flag' "${PIPELINE_ROOT}/data_to_jet/FcsTriggerDefs.h" | sed -e 's/^[[:space:]]*"//' -e 's/",$//')
IFS=',' read -ra REQUESTED <<< "${TRIGGER_FILTER}${REF_TRIGGERS:+,${REF_TRIGGERS}}"
for name in "${REQUESTED[@]}"; do
    grep -qxF "$name" <<< "$VALID_NAMES" || { echo "Error: unknown trigger name '$name'"; exit 1; }
done

mkdir -p "${OUTDIR}/hists" "${OUTDIR}/timing"
[ "$KEEP_JETS" != "0" ] && mkdir -p "${OUTDIR}/jets"
[ "$CHANNEL_QA" = "0" ] || [ "$CHANNEL_QA" = "1" ] || { echo "Error: -Q must be 0 or 1"; exit 1; }
[ "$CHANNEL_QA" = "1" ] && mkdir -p "${OUTDIR}/occ"
MASK_NAME=""
[ "$CHANNEL_MASK" = "none" ] && CHANNEL_MASK=""
if [ -n "$CHANNEL_MASK" ]; then
    [ -f "$CHANNEL_MASK" ] || { echo "Error: channel mask $CHANNEL_MASK not found"; exit 1; }
    MASK_NAME="${LABEL}_$(basename "$CHANNEL_MASK")"
    cp -p "$CHANNEL_MASK" "${SCRIPT_DIR}/channel_qa/masks/${MASK_NAME}"
    # runs of the list outside every period of the mask run unmasked: report them
    # (periods from the mask's "#P period first last" lines, else its rows)
    MASK_UNCOVERED=$(sed -nE 's#^.*_([0-9]{8})_(raw_)?[0-9]{7}\.MuDst\.root.*$#\1#p' "${LISTFILE}" | sort -u \
        | awk 'NR == FNR { if ($1 == "#P") { a[++n] = $3; b[n] = $4 } else if ($0 !~ /^#/ && NF >= 3) { a[++n] = $2; b[n] = $3 }; next }
               { for (i = 1; i <= n; i++) if (a[i] <= $1 && $1 <= b[i]) next; print }' "$CHANNEL_MASK" -)
    if [ -n "$MASK_UNCOVERED" ]; then
        echo "WARNING: $(echo "$MASK_UNCOVERED" | wc -l) run(s) of the list are outside every period of $(basename "$CHANNEL_MASK") and run unmasked:"
        echo "$MASK_UNCOVERED" | head -20 | tr '\n' ' '; echo
    fi
    echo "Channel mask: $(basename "$CHANNEL_MASK") ($(grep -vc '^#' "$CHANNEL_MASK") channel-period rows)"
else
    echo "Channel mask: none (-M none)"
fi
CALIB_NAME=$(basename "$CALIB")
cp -p "$CALIB" "${SCRIPT_DIR}/jes/${CALIB_NAME}"
# point-threshold copies and their maps (same-events difference on the standard map)
case "$THR_OPT" in
    on)
        THR_COPIES="150 200"
        for t in 150 200; do
            m="${PIPELINE_ROOT}/jet_calibration/JES_Lookups/JetEnergyScale_lookup_em_only_pointthr${t}.json"
            [ -f "$m" ] || { echo "Error: threshold map $m not found"; exit 1; }
            cp -p "$m" "${SCRIPT_DIR}/jes/"
        done
        THR150_NAME=JetEnergyScale_lookup_em_only_pointthr150.json
        THR200_NAME=JetEnergyScale_lookup_em_only_pointthr200.json
        # the threshold maps must be built on the nominal map's sample
        # (jes_map_combine.py: "sample": "<nominal> + (<thr> - <ref>)")
        NOM_SAMPLE=$(grep -o '"sample": *"[^"]*"' "$CALIB" | head -1 | sed 's/.*: *"//; s/"$//')
        for t in 150 200; do
            TS=$(grep -o '"sample": *"[^"]*"' "${SCRIPT_DIR}/jes/JetEnergyScale_lookup_em_only_pointthr${t}.json" | head -1 | sed 's/.*: *"//; s/"$//')
            case "$TS" in
                "${NOM_SAMPLE} + ("*) ;;
                *) echo "WARNING: threshold map ${t} is built on \"${TS}\", not on the nominal map's sample \"${NOM_SAMPLE}\"" ;;
            esac
        done
        echo "Point-threshold copies: _vthr150 / _vthr200 (${THR150_NAME}, ${THR200_NAME})"
        ;;
    none) THR_COPIES=""; THR150_NAME="${CALIB_NAME}"; THR200_NAME="${CALIB_NAME}"; echo "Point-threshold copies: none (-T none)" ;;
    *) echo "Error: -T must be on or none"; exit 1 ;;
esac
# pi0 energy scale table (-s): copied into asymmetry/pi0_scale/ for the SandBox
PI0SCALE="pi0_scale.txt"
if [ -n "$PI0_TABLE" ]; then
    [ -f "$PI0_TABLE" ] || { echo "Error: -s ${PI0_TABLE} not found"; exit 1; }
    mkdir -p "${SCRIPT_DIR}/pi0_scale"
    cp -p "$PI0_TABLE" "${SCRIPT_DIR}/pi0_scale/${LABEL}.txt"
    PI0SCALE="pi0_scale/${LABEL}.txt"
fi

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
if [ "$SHUFFLE_SEED" != "0" ]; then
    # random job order, reproducible from the seed
    awk -v seed="${SHUFFLE_SEED}" 'BEGIN { srand(seed) } { printf "%.12f %s\n", rand(), $0 }' "${RANGEFILE}" \
        | sort -k1,1 | cut -d' ' -f2- > "${RANGEFILE}.tmp" && mv "${RANGEFILE}.tmp" "${RANGEFILE}"
    echo "Job order shuffled (seed ${SHUFFLE_SEED})"
fi
# each file's disk node, aligned with the sorted list (catalog_query.sh nodes),
# for the in-job bad-node check; flags and per-node failures go to <outdir>/nodes
NODEFILE="${SORTEDLIST%.list}.nodes"
NODEDIR="${OUTDIR}/nodes"
mkdir -p "${NODEDIR}/fail" "${NODEDIR}/flagged"
# nodes already known to be down (bad_nodes.txt) start flagged, so their files
# are recorded as skipped (and recoverable) rather than silently missing
for nd in $(grep -v '^#' "${SCRIPT_DIR}/bad_nodes.txt" 2>/dev/null | awk 'NF { print $1 }'); do
    [ -e "${NODEDIR}/flagged/${nd}" ] || echo "$(date '+%F %T') ${nd}: listed in bad_nodes.txt at submission" > "${NODEDIR}/flagged/${nd}"
done
[[ "$NODE_FAIL_LIMIT" =~ ^[1-9][0-9]*$ ]] || { echo "Error: -L must be a positive integer"; exit 1; }
if [ "$NO_NODEMAP" = "1" ] || ! "${SCRIPT_DIR}/catalog_query.sh" nodes "${SORTEDLIST}" "${NODEFILE}"; then
    echo "No node map (-K or the catalog query failed): node flagging off"
    awk '{ print "unknown" }' "${SORTEDLIST}" > "${NODEFILE}"
fi
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
    -e "s/{{REF_TRIGGERS}}/${REF_TRIGGERS}/g" \
    -e "s/{{CALIB_NAME}}/${CALIB_NAME}/g" \
    -e "s/{{THR_COPIES}}/${THR_COPIES}/g" \
    -e "s/{{THR150_NAME}}/${THR150_NAME}/g" \
    -e "s/{{THR200_NAME}}/${THR200_NAME}/g" \
    -e "s/{{KEEP_JETS}}/${KEEP_JETS}/g" \
    -e "s/{{CHANNEL_QA}}/${CHANNEL_QA}/g" \
    -e "s/{{CHANNEL_MASK}}/${MASK_NAME}/g" \
    -e "s#{{OUTDIR}}#${OUTDIR}#g" \
    -e "s#{{STDOUT_URL}}#${STDOUT_URL}#g" \
    -e "s#{{STDERR_URL}}#${STDERR_URL}#g" \
    -e "s#{{PI0SCALE}}#${PI0SCALE}#g" \
    -e "s#{{NODEFILE}}#${NODEFILE}#g" \
    -e "s#{{NODEDIR}}#${NODEDIR}#g" \
    -e "s/{{NODE_FAIL_LIMIT}}/${NODE_FAIL_LIMIT}/g" \
    "${SCRIPT_DIR}/an_job.xml" > "${OUTXML}"
echo "Generated ${OUTXML}"

cat > "${OUTDIR}/run_params_${LABEL}.txt" <<EOF
# asymmetry an_job run parameters
submitted_utc=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
label=${LABEL}
pipeline_commit=$(git -C "${PIPELINE_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)
listfile=${LISTFILE}
files=$(wc -l < "${SORTEDLIST}")
nProcesses=${NPROC}
files_per_job=${FILES_PER_JOB}
calib_json=${CALIB}
calib_jes_error=$(grep -q '"jes_error_total"' "$CALIB" && echo yes || echo no)
trigger_filter=${TRIGGER_FILTER}
ref_triggers=${REF_TRIGGERS}
keep_logs=${KEEPLOGS}
keep_jets_every=${KEEP_JETS}
channel_qa=${CHANNEL_QA}
channel_mask=${CHANNEL_MASK:-none}
point_threshold_copies=${THR_COPIES:-none}${THR_COPIES:+ (maps ${THR150_NAME} ${THR200_NAME})}
shuffle_seed=${SHUFFLE_SEED}
pi0_scale_table=asymmetry/${PI0SCALE}${PI0_TABLE:+ (from ${PI0_TABLE})}
node_fail_limit=${NODE_FAIL_LIMIT}
node_file=${NODEFILE}
code_uncommitted_files=$(git -C "${PIPELINE_ROOT}" status --porcelain 2>/dev/null | grep -c . || true)
EOF
cp -p "${SCRIPT_DIR}/${PI0SCALE}" "${OUTDIR}/pi0_scale_${LABEL}.txt"
git -C "${PIPELINE_ROOT}" diff HEAD > "${OUTDIR}/code_${LABEL}.patch" 2>/dev/null || true
git -C "${PIPELINE_ROOT}" status --short > "${OUTDIR}/code_${LABEL}_status.txt" 2>/dev/null || true

[ "$DRYRUN" = "1" ] && { echo "Dry run (-N): not submitted"; exit 0; }

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "asymmetry/an_job_${LABEL}.xml"
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
