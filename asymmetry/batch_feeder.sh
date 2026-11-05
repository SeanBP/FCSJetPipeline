#!/bin/bash
# Submits the next batch of an A_N pass once the earlier batches have drained
# enough, so no more than MAXJOBS of the pass's jobs are in condor at once.
# Runs as a condor scheduler-universe job on the submit node (see
# batch_feeder.sub), so it keeps going after logout; condor_q, star-submit
# and the shared disks all work from there (tested on starsub04 2026-11-05).
#
# Usage: batch_feeder.sh <state dir> <max jobs> [check interval s = 600] [give up after days = 14]
#   <state dir> holds:
#     requests.txt   SUMS request ids already submitted for the pass, one per
#                    line (their idle + running + held jobs are counted)
#     batches.txt    the batches still to submit, in order, one per line:
#                    <njobs> <command>   (run with bash -c, from <state dir>)
#   Each submitted batch is removed from batches.txt and its request id
#   (from the command's "request_id=" or "Submitted request" line) appended
#   to requests.txt; everything is logged to <state dir>/feeder.log.
#   Edit batches.txt to change or stop what is still to come.

STATE="$1"; MAXJOBS="$2"; INTERVAL="${3:-600}"; DAYS="${4:-14}"
[ -d "$STATE" ] && [[ "$MAXJOBS" =~ ^[0-9]+$ ]] || { sed -n '2,20p' "$0"; exit 1; }
LOG="${STATE}/feeder.log"
log() { echo "$(date '+%F %T') $*" >> "$LOG"; }
END=$(( $(date +%s) + DAYS * 86400 ))
log "feeder start on $(hostname), max ${MAXJOBS} jobs, every ${INTERVAL} s"

queued() {
    local n=0 r c
    while read -r r; do
        [ -n "$r" ] || continue
        c=$(condor_q seanp -constraint "regexp(\"${r}\", Args) || regexp(\"${r}\", Arguments) || regexp(\"${r}\", Cmd)" -af ProcId 2>/dev/null | grep -c .)
        n=$(( n + c ))
    done < "${STATE}/requests.txt"
    echo "$n"
}

while [ "$(date +%s)" -lt "$END" ]; do
    NEXT=$(grep -m1 -v '^[[:space:]]*\(#\|$\)' "${STATE}/batches.txt" 2>/dev/null)
    if [ -z "$NEXT" ]; then log "no batches left; feeder done"; exit 0; fi
    NB=$(echo "$NEXT" | awk '{print $1}')
    CMD=$(echo "$NEXT" | cut -d' ' -f2-)
    Q=$(queued)
    if [ $(( Q + NB )) -le "$MAXJOBS" ]; then
        log "queued ${Q} + next ${NB} <= ${MAXJOBS}: submitting: ${CMD}"
        OUT=$(cd "$STATE" && bash -c "$CMD" 2>&1)
        echo "$OUT" >> "$LOG"
        REQ=$(echo "$OUT" | sed -n 's/^Submitted request \([0-9A-F]*\).*/\1/p' | tail -1)
        if [ -z "$REQ" ]; then log "ERROR: no request id from the submission; stopping (batches.txt unchanged)"; exit 1; fi
        echo "$REQ" >> "${STATE}/requests.txt"
        # drop the batch just submitted
        awk -v line="$NEXT" 'done || $0 != line { print; next } { done = 1 }' "${STATE}/batches.txt" > "${STATE}/batches.txt.tmp" \
            && mv "${STATE}/batches.txt.tmp" "${STATE}/batches.txt"
        log "submitted request ${REQ}"
        continue
    fi
    sleep "$INTERVAL"
done
log "gave up after ${DAYS} days with batches left"
exit 1
