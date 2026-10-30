#!/bin/bash
# Submits in_event.xml: one in_event_ue job per pairs list.
# Usage: ./submit_in_event.sh -w <workdir> -l "<list> <list> ..." [-m <max jets, 2000000>] [-k 0|1] [-e <emulation dir>]
#   lists of "<jet_output> <SimpleTree>" pairs: one data run or one ptHat bin per list;
#   job i writes <workdir>/out/in_event_<i>.root
set -e
MAXJ=2000000; KEEP=1; LISTS=""; EMU="-"
while getopts "w:l:m:k:e:h" opt; do
  case $opt in
    w) WORK="$OPTARG" ;; l) LISTS="$OPTARG" ;; m) MAXJ="$OPTARG" ;; k) KEEP="$OPTARG" ;; e) EMU="$OPTARG" ;; *) sed -n '2,5p' "$0"; exit 0 ;;
  esac
done
[ -n "$WORK" ] && [ -n "$LISTS" ] || { sed -n '2,5p' "$0"; exit 1; }
mkdir -p "$WORK/out" "$WORK/log"; WORK=$(readlink -f "$WORK")
n=0
for l in $LISTS; do cp "$l" "$WORK/list_${n}.list"; echo "list_${n}: $l"; n=$((n + 1)); done
if [ "$KEEP" = "1" ]; then SO="file:$WORK/log/\$JOBID.out"; SE="file:$WORK/log/\$JOBID.err"; else SO="file:/dev/null"; SE="file:/dev/null"; fi
cd "$(dirname "$(readlink -f "$0")")/../.."   # the FCSJetPipeline root
XML=$WORK/in_event.xml
sed -e "s#{{WORKDIR}}#$WORK#g" -e "s#{{NPROC}}#$n#" -e "s#{{MAXJETS}}#$MAXJ#" -e "s#{{EMUDIR}}#$EMU#" \
    -e "s#{{STDOUT_URL}}#$SO#" -e "s#{{STDERR_URL}}#$SE#" UnderlyingEventPlayground/in_event/in_event.xml > "$XML"
star-submit "$XML"
