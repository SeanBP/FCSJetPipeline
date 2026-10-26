#!/bin/bash
# Submits sep_test.xml: one job over a list of SimpleTrees (data_to_jet -s 1
# output of one run).
# Usage: ./submit_sep_test.sh -i <SimpleTree list> -o <outdir> [-l <label>]
set -e
LABEL=septest
while getopts "i:o:l:h" opt; do
  case $opt in
    i) LIST=$(readlink -f "$OPTARG") ;; o) OUTDIR="$OPTARG" ;; l) LABEL="$OPTARG" ;;
    *) sed -n '2,4p' "$0"; exit 0 ;;
  esac
done
[ -f "$LIST" ] && [ -n "$OUTDIR" ] || { sed -n '2,4p' "$0"; exit 1; }
cd "$(dirname "$(readlink -f "$0")")/../.."   # the FCSJetPipeline root
mkdir -p "$OUTDIR"/{sep_veto,sep,nosep,log}; OUTDIR=$(readlink -f "$OUTDIR")
mb=$(cat "$LIST" | xargs du -cm | tail -1 | cut -f1)
XML=UnderlyingEventPlayground/sep_test/sep_test_${LABEL}.xml
sed -e "s#{{LIST}}#$LIST#" -e "s#{{OUTDIR}}#$OUTDIR#g" -e "s#{{MIN_STORAGE}}#$((mb * 2 + 2000))#" \
    -e "s#{{STDOUT_URL}}#file:$OUTDIR/log/\$JOBID.out#" -e "s#{{STDERR_URL}}#file:$OUTDIR/log/\$JOBID.err#" \
    UnderlyingEventPlayground/sep_test/sep_test.xml > "$XML"
star-submit "$XML"
