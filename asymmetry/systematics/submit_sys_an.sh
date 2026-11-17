#!/bin/bash
# Submits AnAsymmetry passes for systematic comparisons (2026-11-17): one condor
# job per "name|yield list|AnAsymmetry options" line of a job file. The code
# (sivers/AnAsymmetry.cpp, AnConfig.h, AnRunTools.h, run22_polarization.txt)
# is frozen into <out dir>/code first, so later edits do not change the jobs.
# Usage: submit_sys_an.sh <out dir> <job file> [-N (write the condor file only)]
# Then: sys_slices.py <loose>.txt <tight>.txt <prefix> (nested cuts) or
#       sys_slices.py --disjoint <a>.txt <b>.txt <prefix> (independent samples).
set -e
OUTDIR="$(mkdir -p "$1" && cd "$1" && pwd)"; JOBS="$2"; NOSUB="$3"
[ -f "$JOBS" ] || { sed -n '2,9p' "$0"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE="${OUTDIR}/code"
mkdir -p "$CODE" "${OUTDIR}/logs"
cp "${HERE}/../sivers/AnAsymmetry.cpp" "${HERE}/../AnConfig.h" "${HERE}/../AnRunTools.h" "${HERE}/../run22_polarization.txt" \
   "${HERE}/sys_an_job.sh" "$CODE/"
cp "$JOBS" "${OUTDIR}/jobs.txt"
( cd "${HERE}/.." && git rev-parse HEAD 2>/dev/null; git -C "${HERE}/.." diff -- sivers/AnAsymmetry.cpp AnConfig.h ) > "${CODE}/code_info.txt" || true
SUB="${OUTDIR}/sys_an.sub"
cat > "$SUB" <<EOF
universe       = vanilla
executable     = ${CODE}/sys_an_job.sh
log            = ${OUTDIR}/logs/jobs.condor.log
request_memory = 1500
request_disk   = 6000000
getenv         = true
priority       = 20
EOF
while IFS='|' read -r name list opts; do
    [ -z "$name" ] || [ "${name:0:1}" = "#" ] && continue
    cp "$list" "${OUTDIR}/${name}.list"
    cat >> "$SUB" <<EOF

arguments      = ${OUTDIR}/${name}.list ${CODE} ${OUTDIR} ${name} ${opts}
output         = ${OUTDIR}/logs/${name}.out
error          = ${OUTDIR}/logs/${name}.err
queue
EOF
done < "$JOBS"
echo "wrote $SUB ($(grep -c '^queue' "$SUB") jobs)"
[ "$NOSUB" = "-N" ] || condor_submit "$SUB"
