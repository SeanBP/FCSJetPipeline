#!/bin/bash
# Re-runs AnHistograms on the kept jet files of an A_N pass (an_job -J 1),
# for a changed yield selection (AnConfig.h), without new reconstruction.
# Usage:
#   ./submit_rehist.sh -i <pass outdir> -o <new outdir> [-n <files per job>] [-l <label>] [-s <pi0 scale table>] [-N]
#
#   -i  the pass's output directory: its jets/jets_*.root are the input
#   -o  new output directory: yields go to <new outdir>/hists/ under the
#       pass's names (an_hists_<run>_<job>.root). timing/ and run_params_*
#       are linked from the pass, so smoke_report.sh <new outdir> ... works
#       (event counts, resources) as on the pass itself
#   -n  jet files per job, default 100 (each file is ~20-30 MB, a few s of
#       AnHistograms; files are copied to scratch one at a time)
#   -l  label for the generated XML, default "rehist"
#   -s  pi0 energy scale table (asymmetry/pi0_scale.txt format; e.g. a
#       pi0_peak_runs.C <out>_scale_per_run.txt or _per_period.txt). Default:
#       asymmetry/pi0_scale.txt
#   -N  write the XML only
# Writes <new outdir>/rehist_params_<label>.txt (inputs, AnConfig.h trigger
# rule, pi0 scale table; its copy pi0_scale_<label>.txt), code diff and the
# request id.

set -e
IN=""; OUT=""; NPER=100; LABEL="rehist"; DRY=0; PI0_TABLE=""
while getopts "i:o:n:l:s:Nh" opt; do
  case $opt in
    i) IN="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    n) NPER="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    s) PI0_TABLE="$OPTARG" ;;
    N) DRY=1 ;;
    h|*) sed -n '2,22p' "$0"; exit 0 ;;
  esac
done
[ -d "$IN/jets" ] && [ -n "$OUT" ] || { sed -n '2,22p' "$0"; exit 1; }
[[ "$NPER" =~ ^[1-9][0-9]*$ ]] || { echo "Error: -n must be a positive integer"; exit 1; }
IN="$(cd "$IN" && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"

mkdir -p "$OUT/hists"
OUT="$(cd "$OUT" && pwd)"
LIST="$OUT/jets_${LABEL}.list"
find "$IN/jets" -maxdepth 1 -name 'jets_*.root' | sort > "$LIST"
NF=$(wc -l < "$LIST")
[ "$NF" -gt 0 ] || { echo "Error: no jets_*.root in $IN/jets"; exit 1; }
NPROC=$(( (NF + NPER - 1) / NPER ))
echo "Re-histogram: ${NF} jet files in ${NPROC} jobs of up to ${NPER}"

# the pass's timing tables and parameters, for smoke_report.sh
[ -e "$OUT/timing" ] || ln -s "$IN/timing" "$OUT/timing"
for f in "$IN"/run_params_*.txt; do [ -e "$OUT/$(basename "$f")" ] || ln -s "$f" "$OUT/$(basename "$f")"; done

# pi0 energy scale table (-s): copied into asymmetry/pi0_scale/ for the SandBox
PI0SCALE="pi0_scale.txt"
if [ -n "$PI0_TABLE" ]; then
    [ -f "$PI0_TABLE" ] || { echo "Error: -s ${PI0_TABLE} not found"; exit 1; }
    mkdir -p "${SCRIPT_DIR}/pi0_scale"
    cp -p "$PI0_TABLE" "${SCRIPT_DIR}/pi0_scale/${LABEL}.txt"
    PI0SCALE="pi0_scale/${LABEL}.txt"
fi
XML="${SCRIPT_DIR}/rehist_job_${LABEL}.xml"
sed -e "s#{{LISTFILE}}#${LIST}#g" -e "s/{{NPERJOB}}/${NPER}/g" -e "s/{{NPROCESSES}}/${NPROC}/g" -e "s#{{OUTDIR}}#${OUT}#g" \
    -e "s#{{PI0SCALE}}#${PI0SCALE}#g" "${SCRIPT_DIR}/rehist_job.xml" > "$XML"
{
    echo "# rehist parameters"
    echo "submitted_utc=$(date -u +%FT%TZ)"
    echo "input=${IN}"
    echo "files=${NF}"
    echo "nProcesses=${NPROC}"
    echo "files_per_job=${NPER}"
    echo "pipeline_root=${PIPELINE_ROOT}"
    echo "pipeline_commit=$(git -C "${PIPELINE_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)"
    grep -E '^static const double k(EmCutFactor|JpCutFactor|EmVetoMargin|Pi0EnergyScale) ' "${SCRIPT_DIR}/AnConfig.h"
    echo "pi0_scale_table=asymmetry/${PI0SCALE}${PI0_TABLE:+ (from ${PI0_TABLE})}"
} > "$OUT/rehist_params_${LABEL}.txt"
cp -p "${SCRIPT_DIR}/${PI0SCALE}" "$OUT/pi0_scale_${LABEL}.txt"
git -C "${PIPELINE_ROOT}" diff HEAD > "$OUT/code_${LABEL}.patch" 2>/dev/null || true
echo "Generated ${XML}"
[ "$DRY" = "1" ] && { echo "Dry run (-N): not submitted"; exit 0; }

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "asymmetry/rehist_job_${LABEL}.xml"
REPORT=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*.report' -newer "${MARKER}" | head -1)
REQID=$(basename "${REPORT}" .report | sed 's/^sched//')
rm -f "${MARKER}"
echo "request_id=${REQID}" >> "$OUT/rehist_params_${LABEL}.txt"
echo "Submitted request ${REQID}"
