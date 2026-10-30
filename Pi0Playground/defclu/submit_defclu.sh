#!/bin/bash
# Submits defclu.xml: MuDst -> SimpleTree with default FCS clustering, for every
# file of a MakeRunList.pl-format list.
# Usage: ./submit_defclu.sh -w <workdir> -l <list> [-n <files per job, default 5>] [-t <trigger filter>]
#   -t  comma-separated trigger flags; an event is written only if one fired.
#       Default fcsEM0,fcsEM1,fcsEM2; "none" = every event
set -e
NPER=5; TRIG=fcsEM0,fcsEM1,fcsEM2
while getopts "w:l:n:t:h" opt; do
  case $opt in w) WORK="$OPTARG" ;; l) LIST="$OPTARG" ;; n) NPER="$OPTARG" ;; t) TRIG="$OPTARG" ;; *) sed -n '2,6p' "$0"; exit 0 ;; esac
done
[ -n "$WORK" ] && [ -n "$LIST" ] || { sed -n '2,6p' "$0"; exit 1; }
mkdir -p "$WORK/simpletrees" "$WORK/log"; WORK=$(readlink -f "$WORK"); LIST=$(readlink -f "$LIST")
N=$(grep -c . "$LIST"); NPROC=$(( (N + NPER - 1) / NPER ))
echo "$N files, $NPER per job: $NPROC jobs"
cd "$(dirname "$(readlink -f "$0")")/../.."   # the FCSJetPipeline root
sed -e "s#{{WORKDIR}}#$WORK#g" -e "s#{{LIST}}#$LIST#g" -e "s#{{NPERJOB}}#$NPER#g" -e "s#{{NPROC}}#$NPROC#" \
    -e "s#{{TRIGFILTER}}#$([ "$TRIG" = none ] && echo "" || echo "$TRIG")#g" \
    Pi0Playground/defclu/defclu.xml > "$WORK/defclu_$(basename "$LIST" .list).xml"
star-submit "$WORK/defclu_$(basename "$LIST" .list).xml"
