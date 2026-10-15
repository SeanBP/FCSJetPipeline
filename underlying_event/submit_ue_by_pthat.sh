#!/bin/bash
# Submits the cross-job underlying-event estimate for one or more simulation
# samples: one condor job per ptHat bin, pooling that bin's events from every
# sample given (proxy events never cross ptHat bins, see UnderlyingEvent.cpp).
# Each job copies its bin's file pairs to scratch, runs UnderlyingEvent there
# and copies the rewritten jet_output files back over the originals
# (run_ue_bin.sh). Use it for samples whose SimpleTrees were kept (-s 1),
# e.g. a jet-energy-scale calibration sample (sim_to_jet -C 1) with "-m minbias".
#
# Usage:
#   ./submit_ue_by_pthat.sh -i <sample_dir>[,<sample_dir>...] -w <workdir> \
#       [-m minbias|default] [-J "<em_only> <use_points> <jet_pt_cuts>"] [-l <label>] [-k <0|1>] [-n]
#
#   -i  comma-separated sim_to_jet output folders, either split by
#       sim_to_jet/organize_sim_output.sh (jet_output/pt<N>/, SimpleTree/pt<N>/)
#       or flat
#   -w  folder for the per-bin lists, the job logs and a summary (created if
#       missing); must be on a shared disk the jobs can read
#   -m  "minbias" (default): every fiducial jet estimated, untriggered events
#       matched with untriggered ones; "default": the analysis behaviour
#       (EM-triggered jets above the reco pT threshold only)
#   -J  rebuild every jet file from its SimpleTree with JetMatcher first, with
#       these options and no trigger filter (e.g. "1 1 0" for an EM-only,
#       point-based calibration sample)
#   -l  label for the generated XML's filename, default "run"
#   -k  1 = keep each job's stdout/stderr under <workdir>/log/ (default 1)
#   -n  dry run: write the lists and the XML, but do not submit
#
# Jet files without a SimpleTree of the same name are listed in
# <workdir>/unpaired.txt and left out.

set -e
INPUTS=""; WORKDIR=""; MODE="minbias"; LABEL="run"; KEEPLOGS=1; DRYRUN=0; JMARGS=""
while getopts "i:w:m:J:l:k:nh" opt; do
  case $opt in
    i) INPUTS="$OPTARG" ;; w) WORKDIR="$OPTARG" ;; m) MODE="$OPTARG" ;;
    J) JMARGS="$OPTARG" ;; l) LABEL="$OPTARG" ;; k) KEEPLOGS="$OPTARG" ;; n) DRYRUN=1 ;;
    h|*) sed -n '2,30p' "$0"; exit 0 ;;
  esac
done
if [ -z "$INPUTS" ] || [ -z "$WORKDIR" ]; then echo "Missing -i or -w. Run with -h for usage."; exit 1; fi
case "$MODE" in
  minbias) UE_ARGS="0 minbias" ;;
  default) UE_ARGS="" ;;
  *) echo "Error: -m must be minbias or default"; exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
mkdir -p "$WORKDIR/log"
WORKDIR="$(cd "$WORKDIR" && pwd)"

# All pairs, then one list per ptHat bin
ALL="$WORKDIR/all_pairs.txt"
: > "$ALL"; : > "$WORKDIR/unpaired.txt"
IFS=',' read -ra DIRS <<< "$INPUTS"
for d in "${DIRS[@]}"; do
  d="$(cd "$d" && pwd)"
  find "$d" -maxdepth 3 -type f -name 'jet_output_*.root' | sort | while read -r jet; do
    base=$(basename "$jet")
    sname=${base/jet_output_/SimpleTree_}
    bin=$(echo "$base" | sed -n -e 's/^.*_\(pt[0-9]*\)_[0-9]*\.root$/\1/p')
    simple=""
    for cand in "$d/SimpleTree/$bin/$sname" "$d/$sname"; do
      if [ -f "$cand" ]; then simple="$cand"; break; fi
    done
    if [ -z "$simple" ] || [ -z "$bin" ]; then echo "$jet" >> "$WORKDIR/unpaired.txt"; continue; fi
    echo "$bin $jet $simple" >> "$ALL"
  done
done
BINS=$(cut -d' ' -f1 "$ALL" | sort -u | sort -t t -k2 -n)
NBINS=$(echo "$BINS" | grep -c . || true)
if [ "$NBINS" -eq 0 ]; then echo "Error: no file pairs found"; exit 1; fi

: > "$WORKDIR/bins.txt"
i=0; max_mb=0
for b in $BINS; do
  grep "^$b " "$ALL" | cut -d' ' -f2- > "$WORKDIR/bin_$i.list"
  mb=$(cut -d' ' -f1,2 --output-delimiter=$'\n' "$WORKDIR/bin_$i.list" | xargs du -cm | tail -1 | cut -f1)
  echo "$i $b $(wc -l < "$WORKDIR/bin_$i.list") pairs ${mb} MB" >> "$WORKDIR/bins.txt"
  [ "$mb" -gt "$max_mb" ] && max_mb=$mb
  i=$((i + 1))
done
MIN_STORAGE=$(( max_mb * 5 / 2 + 1000 ))
echo "$(wc -l < "$ALL") file pairs in $NBINS ptHat bins (largest ${max_mb} MB, requesting ${MIN_STORAGE} MB scratch);"
echo "$(wc -l < "$WORKDIR/unpaired.txt") unpaired jet files; per-bin summary in $WORKDIR/bins.txt"

if [ "$KEEPLOGS" = "1" ]; then
  STDOUT_URL="file:${WORKDIR}/log/ue_\$JOBINDEX.out"; STDERR_URL="file:${WORKDIR}/log/ue_\$JOBINDEX.err"
else
  STDOUT_URL="file:/dev/null"; STDERR_URL="file:/dev/null"
fi
OUTXML="${SCRIPT_DIR}/ue_by_pthat_${LABEL}.xml"
sed -e "s#{{WORKDIR}}#${WORKDIR}#g" -e "s#{{UE_ARGS}}#${UE_ARGS}#g" -e "s#{{MIN_STORAGE}}#${MIN_STORAGE}#g" \
    -e "s#{{STDOUT_URL}}#${STDOUT_URL}#g" -e "s#{{STDERR_URL}}#${STDERR_URL}#g" -e "s#{{JETMATCHER_ARGS}}#${JMARGS}#g" \
    -e "s#nProcesses=\"1\"#nProcesses=\"${NBINS}\"#" "${SCRIPT_DIR}/ue_by_pthat.xml" > "$OUTXML"
echo "Generated $OUTXML"

GITHASH=$(git -C "${PIPELINE_ROOT}" rev-parse HEAD 2>/dev/null || echo "unknown")
cat > "$WORKDIR/run_params_${LABEL}.txt" <<EOF
# submit_ue_by_pthat.sh run parameters
submitted_utc=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
inputs=${INPUTS}
mode=${MODE}
ue_args=${UE_ARGS}
jetmatcher_args=${JMARGS:-<none>}
n_bins=${NBINS}
min_storage_mb=${MIN_STORAGE}
pipeline_commit=${GITHASH}
generated_xml=${OUTXML}
EOF

if [ "$DRYRUN" = "1" ]; then echo "Dry run: not submitted"; exit 0; fi
cd "${PIPELINE_ROOT}"
star-submit "underlying_event/ue_by_pthat_${LABEL}.xml"
