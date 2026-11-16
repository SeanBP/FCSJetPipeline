#!/bin/bash
# Merged report flow, step 1 (Sean 2026-11-16; AnMerge.cpp): sums a pass's
# per-run yield files in parallel condor jobs, so the report reads them once
# instead of once per A_N mode. Each job copies its chunk of yield files to its
# scratch disk, runs AnMerge -b, and writes <outdir>/merge/partial_<i>.root.
# smoke_report.sh / submit_report_job.sh then find the partials and use the
# merged flow (AnMerge -a, pi0_peak_runs, AnMerge -x, AnAsymmetry -M,
# AnPlots -m).
# Usage: ./submit_an_merge.sh <outdir> [files per job = 300] [memory MB = 1000]
#   ~0.4 s per file, so 300 files is ~2 min of merging plus the copy: about
#   30 jobs for a full pass (~9,000 files). Incremental: run it again after
#   new yield files arrive (recovery.py link) and only those are merged, as
#   further partials (chunks numbered on; merge/merged_yields.list holds every
#   file already in a partial). The code of each round is frozen into
#   <outdir>/merge/code_<round>. Logs: <outdir>/merge/logs/merge_<i>.log
set -e
OUTDIR="$(cd "$1" && pwd)"
NPER="${2:-300}"; MEM="${3:-1000}"
[ -d "$OUTDIR/hists" ] || { sed -n '2,15p' "$0"; exit 1; }
[[ "$NPER" =~ ^[1-9][0-9]*$ ]] || { echo "Error: files per job must be a positive integer"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src() { if [ -e "${HERE}/$1" ]; then echo "${HERE}/$1"; else echo "${HERE}/../$1"; fi; }
M="$OUTDIR/merge"
mkdir -p "$M/chunks" "$M/logs"
touch "$M/merged_yields.list"
# yield files older than 10 min (condor may still be copying the newest), not yet in a partial
find -L "$OUTDIR/hists" -maxdepth 1 -name 'an_hists_*.root' -mmin +10 | sort > "$M/all_yields.list"
sort "$M/merged_yields.list" | comm -23 "$M/all_yields.list" - > "$M/new_yields.list"
NF=$(wc -l < "$M/new_yields.list")
[ "$NF" -gt 0 ] || { echo "nothing new to merge in $OUTDIR/hists ($(wc -l < "$M/merged_yields.list") files already in partials)"; exit 0; }
# chunks numbered after the existing ones
OFFSET=$(ls "$M"/chunks/chunk_*.list 2>/dev/null | wc -l)
ROUND=$(ls -d "$M"/code_* 2>/dev/null | wc -l)
CODE="$M/code_${ROUND}"; mkdir -p "$CODE"
for f in AnMerge.cpp AnRunTools.h AnConfig.h run22_polarization.txt an_merge_job.sh; do cp "$(src $f)" "$CODE/"; done
TMPC=$(mktemp -d)
split -l "$NPER" -d -a 4 "$M/new_yields.list" "$TMPC/c_"
NJ=0
for c in $(ls "$TMPC"/c_* | sort); do
    mv "$c" "$M/chunks/chunk_$(printf "%04d" $(( OFFSET + NJ ))).list"; NJ=$(( NJ + 1 ))
done
rmdir "$TMPC"
cat "$M/new_yields.list" >> "$M/merged_yields.list"
cat > "$CODE/merge.sub" <<SUB
universe       = vanilla
executable     = $CODE/an_merge_job.sh
arguments      = $OUTDIR $CODE \$(Process) $OFFSET
output         = $M/logs/merge_r${ROUND}_\$(Process).log
error          = $M/logs/merge_r${ROUND}_\$(Process).log
log            = $CODE/merge.condor.log
request_memory = $MEM
request_disk   = 1000
getenv         = true
queue $NJ
SUB
echo "merge round ${ROUND}: ${NF} new yield files in ${NJ} jobs of up to ${NPER} (chunks from $(printf "%04d" ${OFFSET})); code frozen in ${CODE}"
condor_submit "$CODE/merge.sub"
