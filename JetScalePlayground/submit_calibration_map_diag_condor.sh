#!/bin/bash
# Runs the diagnostics fork of CalibrationMap's position-bin fit in
# parallel via TRUE separate condor jobs -- N single-core jobs, each fitting
# a disjoint shard of the ~1800 independent position bins (see
# CalibrationMap.cpp's compute_shard), merged locally afterward. This is
# NOT the old fork()-based in-process parallelism (CalibrationMap.cpp's
# compute_all_regions, n_procs_hint>1) -- that path is hard-disabled and
# stays that way; it silently corrupted results because every worker was a
# fork() of one already-ROOT/Minuit2-initialized process (see its header
# comment). A condor job is a genuinely separate process from a clean start,
# so there's no shared post-init state to corrupt. Verified byte-identical
# JSON output against the single-process path on a synthetic dataset before
# this was written (see session notes) -- the shard/merge split changes
# nothing about the fit itself, each position bin's fit is completely
# independent of every other.
#
# Also deliberately NOT one job requesting many cores (this script's
# previous approach, before the condor-job-per-shard rewrite): a
# request_cpus=8 job sat idle for over an hour on this pool once, since it
# needs one machine with 8 simultaneously-free cores on a pool dominated by
# single-core jobs. N separate single-core jobs (like sim_to_jet's own
# nProcesses=500 array) match what the pool actually handles well.
#
# Blocks until every shard job finishes (condor_wait), verifies each one
# actually produced its shard CSV (no silent failures), then runs the merge
# step locally (fast -- no fitting, just I/O + the existing
# smoothing/decay/clamp/fiducial-mask/save pipeline).
#
# Usage:
#   ./submit_calibration_map_diag_condor.sh -i <jet_calibration.root> -o <output.json> [-E <0|1>] [-n <n_shards>]
#
#   -i  absolute path to a jet_calibration.root (from
#       run_jet_energy_scale_diag.sh or production JetEnergyScaleFineGrid)
#   -o  output JSON path (created if its directory doesn't exist)
#   -E  1 = the input came from an EM-only JetEnergyScaleFineGrid run
#       (narrows the skew-Gaussian position/width domain to match its
#       Emax=43), 0 = default ECAL+HCAL domain (Emax=200). Default 0.
#   -n  number of shards / condor jobs. Default 25 (~70 position bins per
#       shard on the ~1800-bin em_only production grid). More shards means
#       faster wall time up to the point condor queue/startup overhead per
#       job dominates -- 25 is a reasonable starting point, not tuned.
#
# Example:
#   ./submit_calibration_map_diag_condor.sh \
#       -i em_only_production_jet_calibration_r10_emax100.root \
#       -o em_only_production_lookup_r10_emax100.json -E 1 -n 25

set -e

INPUT=""
OUTPUT=""
EM_ONLY=0
NSHARDS=25

while getopts "i:o:E:n:h" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    E) EM_ONLY="$OPTARG" ;;
    n) NSHARDS="$OPTARG" ;;
    h|*)
      sed -n '2,35p' "$0"
      exit 0
      ;;
  esac
done

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: input file '$INPUT' does not exist"
    exit 1
fi

if [ "$EM_ONLY" != "0" ] && [ "$EM_ONLY" != "1" ]; then
    echo "Error: -E must be 0 or 1 (got '$EM_ONLY')"
    exit 1
fi

if ! [[ "$NSHARDS" =~ ^[0-9]+$ ]] || [ "$NSHARDS" -lt 1 ]; then
    echo "Error: -n must be a positive integer (got '$NSHARDS')"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"

mkdir -p "$(dirname "$OUTPUT")" 2>/dev/null || true
OUTPUT="$(cd "$(dirname "$OUTPUT")" && pwd)/$(basename "$OUTPUT")"

STAMP="$(date +%Y%m%d_%H%M%S)_$$"
LOGDIR="${SCRIPT_DIR}/condor_logs/calibmap_${STAMP}"
mkdir -p "$LOGDIR"
BINPATH="${LOGDIR}/.calibration_map_diag_bin"

echo "Building CalibrationMap (once, shared read-only by every shard job)..."
BUILD_CSH="${LOGDIR}/build.csh"
cat > "$BUILD_CSH" << EOF
#!/bin/csh
module purge
module load root-5.34.38
cd ${SCRIPT_DIR}
g++ -m64 -pthread CalibrationMap.cpp -o ${BINPATH} \`root-config --cflags --libs\` -lMinuit2 -std=c++11
EOF
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 -B /tmp /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh "$BUILD_CSH"
if [ ! -x "$BINPATH" ]; then
    echo "ERROR: build failed, no binary at ${BINPATH}"
    exit 1
fi
rm -f "$BUILD_CSH"

WRAPPER="${LOGDIR}/shard.csh"
cat > "$WRAPPER" << EOF
#!/bin/csh
set SHARD = \$1
${BINPATH} "${INPUT}" "" ${EM_ONLY} 0 --shard \${SHARD} ${NSHARDS} ${LOGDIR}/shard_\${SHARD}.csv
EOF
chmod +x "$WRAPPER"

SUBMIT="${LOGDIR}/shard.condor"
cat > "$SUBMIT" << EOF
universe        = vanilla
notification    = never
executable      = /bin/env
arguments       = "singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 -B /tmp /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh ${WRAPPER} \$(Process)"
initialdir      = ${SCRIPT_DIR}
request_cpus    = 1
request_memory  = 2048
kill_sig        = SIGINT
output          = ${LOGDIR}/shard_\$(Process).out
error           = ${LOGDIR}/shard_\$(Process).err
log             = ${LOGDIR}/shard.log
queue ${NSHARDS}
EOF

echo "Submitting ${NSHARDS} single-core shard jobs..."
condor_submit "$SUBMIT"

echo "Waiting for all ${NSHARDS} shards to finish (log: ${LOGDIR}/shard.log)..."
if ! condor_wait "${LOGDIR}/shard.log"; then
    echo "ERROR: condor_wait failed -- check ${LOGDIR}/shard.log"
    exit 1
fi

# condor_wait only confirms every job left the queue, not that each
# succeeded -- check every shard's own exit code and output file, same
# never-fail-silently standard as the fork-based path's shard I/O checks.
FAILED=0
for ((i=0; i<NSHARDS; i++)); do
    EXITCODE=$(grep -a "return value" "${LOGDIR}/shard.log" | sed -n "$((i+1))p" | sed -n 's/.*return value \([0-9]*\).*/\1/p') || true
    if [ "$EXITCODE" != "0" ]; then
        echo "ERROR: shard $i exited with code '${EXITCODE:-unknown}' -- see ${LOGDIR}/shard_${i}.out / .err"
        FAILED=1
    fi
    if [ ! -f "${LOGDIR}/shard_${i}.csv" ]; then
        echo "ERROR: shard $i produced no output file (${LOGDIR}/shard_${i}.csv missing)"
        FAILED=1
    fi
done

if [ "$FAILED" = "1" ]; then
    echo "ERROR: one or more shards failed -- not merging. Logs kept at ${LOGDIR}"
    exit 1
fi

echo "All ${NSHARDS} shards succeeded. Merging..."
SHARD_LIST="${LOGDIR}/shard_list.txt"
> "$SHARD_LIST"
for ((i=0; i<NSHARDS; i++)); do
    echo "${LOGDIR}/shard_${i}.csv" >> "$SHARD_LIST"
done

MERGE_CSH="${LOGDIR}/merge.csh"
cat > "$MERGE_CSH" << EOF
#!/bin/csh
module purge
module load root-5.34.38
${BINPATH} "${INPUT}" "" ${EM_ONLY} 0 --merge ${SHARD_LIST} "${OUTPUT}"
EOF
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 -B /tmp /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh "$MERGE_CSH"
MERGESTATUS=$?
rm -f "$MERGE_CSH"

if [ "$MERGESTATUS" != "0" ]; then
    echo "ERROR: merge step failed (exit $MERGESTATUS) -- shard CSVs kept at ${LOGDIR}"
    exit 1
fi

echo "Done. Lookup table: ${OUTPUT}"
rm -f "$BINPATH" "$WRAPPER" "$SUBMIT"
echo "Shard logs kept at ${LOGDIR} for reference (not auto-deleted)."
