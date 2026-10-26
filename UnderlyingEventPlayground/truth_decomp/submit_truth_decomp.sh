#!/bin/bash
# Submits truth_decomp.xml: one job per ptHat bin of a sim_to_jet sample.
# Usage: ./submit_truth_decomp.sh -i <ue_by_pthat workdir with bin_<n>.list> -w <workdir> [-s <sample filter>] [-m <max jets>] [-n <n_max>] [-o "<options>"] [-k 0|1]
#   -i  folder of per-bin pairs lists ("<jet_output> <SimpleTree>"), e.g. a
#       submit_ue_by_pthat.sh workdir
#   -s  keep only pairs whose path contains this (e.g. fms_tune_calib_v1)
#   -m  jets sampled per bin (default 20000); -n proxy events per jet (default 500)
#   -o  ue_truth_decomposition options: trig (EM trigger levels), mixed (the
#       bins of one list pooled by cross section; give one combined list)
set -e
FILTER=""; MAXJ=20000; NMAX=500; KEEP=1; OPTS=""
while getopts "i:w:s:m:n:o:k:h" opt; do
  case $opt in
    i) IN=$(readlink -f "$OPTARG") ;; w) WORK="$OPTARG" ;; s) FILTER="$OPTARG" ;;
    m) MAXJ="$OPTARG" ;; n) NMAX="$OPTARG" ;; o) OPTS="$OPTARG" ;; k) KEEP="$OPTARG" ;; *) sed -n '2,10p' "$0"; exit 0 ;;
  esac
done
[ -d "$IN" ] && [ -n "$WORK" ] || { sed -n '2,10p' "$0"; exit 1; }
mkdir -p "$WORK/out" "$WORK/log"; WORK=$(readlink -f "$WORK")
n=0
for f in $(ls "$IN"/bin_*.list | sort -V); do
    grep -F -- "$FILTER" "$f" > "$WORK/bin_${n}.list" || true
    [ -s "$WORK/bin_${n}.list" ] && n=$((n + 1))
done
echo "$n bins"
if [ "$KEEP" = "1" ]; then SO="file:$WORK/log/\$JOBID.out"; SE="file:$WORK/log/\$JOBID.err"; else SO="file:/dev/null"; SE="file:/dev/null"; fi
cd "$(dirname "$(readlink -f "$0")")/../.."   # the FCSJetPipeline root
XML=UnderlyingEventPlayground/truth_decomp/truth_decomp_$(basename "$WORK").xml
sed -e "s#{{WORKDIR}}#$WORK#g" -e "s#{{NPROC}}#$n#" -e "s#{{MAXJETS}}#$MAXJ#" -e "s#{{NMAX}}#$NMAX#" -e "s#{{OPTS}}#$OPTS#" \
    -e "s#{{STDOUT_URL}}#$SO#" -e "s#{{STDERR_URL}}#$SE#" UnderlyingEventPlayground/truth_decomp/truth_decomp.xml > "$XML"
star-submit "$XML"
