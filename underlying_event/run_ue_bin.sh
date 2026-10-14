#!/bin/bash
# One job of submit_ue_by_pthat.sh: the underlying-event estimate for one
# group of file pairs (one ptHat bin of a simulation sample), run in the job's
# scratch folder. Shared disks are slow to read and many open files there slow
# them for everyone, so the pairs are first copied here; UnderlyingEvent runs
# on the copies; then only the rewritten jet files are copied back, each to
# <file>.ue_new beside the original and renamed over it, so a file on disk is
# never left half written. Nothing is copied back if UnderlyingEvent fails.
#
# Usage (inside the job, UnderlyingEvent already built here):
#   ./run_ue_bin.sh <pairs list> [UnderlyingEvent options, e.g. "0 minbias"]
# The list has one line per pair: "<jet_output file> <SimpleTree file>".

set -e
LIST="$1"; shift
if [ ! -f "$LIST" ]; then echo "run_ue_bin.sh: no pairs list '$LIST'"; exit 1; fi

# Each pair in its own folder: two samples can hold files of the same name
mkdir -p ue_in
: > ue_local.list
n=0
while read -r jet simple; do
    [ -z "$jet" ] && continue
    n=$((n + 1))
    mkdir -p "ue_in/$n"
    cp "$jet" "ue_in/$n/"
    cp "$simple" "ue_in/$n/"
    echo "$PWD/ue_in/$n/$(basename "$jet") $PWD/ue_in/$n/$(basename "$simple")" >> ue_local.list
done < "$LIST"
echo "run_ue_bin.sh: copied $n file pairs to scratch ($(du -sh ue_in | cut -f1))"

./UnderlyingEvent ue_local.list "$@"

n=0
while read -r jet simple; do
    [ -z "$jet" ] && continue
    n=$((n + 1))
    cp "ue_in/$n/$(basename "$jet")" "$jet.ue_new"
    mv "$jet.ue_new" "$jet"
done < "$LIST"
echo "run_ue_bin.sh: copied $n rewritten jet files back"
rm -rf ue_in ue_local.list
