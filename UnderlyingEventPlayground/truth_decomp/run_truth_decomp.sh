#!/bin/bash
# One job of submit_truth_decomp.sh: copies one ptHat bin's file pairs to
# scratch and runs both builds of ue_truth_decomposition on them
# (separation: sep_veto and sep; no separation: nosep).
#
# Usage (inside the job, td_sep and td_nosep built here):
#   ./run_truth_decomp.sh <pairs list> <tag> [max_jets] [n_max] [options: trig mixed]
set -e
LIST="$1"; TAG="$2"; MAXJ="${3:-20000}"; NMAX="${4:-500}"; shift 4 || shift $#; OPTS="$*"
mkdir -p in; : > local.list
while read -r jet simple; do
    [ -z "$jet" ] && continue
    cp "$jet" "$simple" in/
    echo "$PWD/in/$(basename "$jet") $PWD/in/$(basename "$simple")" >> local.list
done < "$LIST"
echo "run_truth_decomp.sh: $(wc -l < local.list) pairs ($(du -sh in | cut -f1))"
./td_sep local.list td_sep_${TAG}.root "$MAXJ" "$NMAX" $OPTS | grep -v '^#'
./td_nosep local.list td_nosep_${TAG}.root "$MAXJ" "$NMAX" $OPTS | grep -v '^#'
rm -rf in local.list
