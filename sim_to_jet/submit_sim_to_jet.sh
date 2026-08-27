#!/bin/bash
# Generates and submits a concrete sim_to_jet.xml job.
# Usage:
#   ./submit_sim_to_jet.sh -g <pythia8|pythia6> -t <tune_param> -p <ptcut> \
#       -n <nevents> -s <0|1> -o <outdir> [-j <nprocesses>] [-l <label>] \
#       [-k <0|1>] [-f <triggers>] [-c <run_number>] [-e <filter_ethr>]
#
#   -g  generator: pythia8 or pythia6
#   -t  tune/PDF: pythia8 PDF:pSet (8=CTEQ6L1, 5=MSTW2008LO, 3=MRST LO*,
#       21=NNPDF3.1sx, <=0=default); pythia6 PyTune (325=Perugia STAR, 0=default)
#   -p  ptHatMin cut in GeV, flat for every job (ignored if -w given)
#   -n  events per job
#   -s  1 = keep intermediate SimpleTree, 0 = discard
#   -o  output directory (created if missing)
#   -j  number of parallel jobs (nProcesses), default 1
#   -w  absolute path to a JSON of relative ptHatMin proportions to sample
#       across the -j jobs (see pthat_distribution_optimized.json); overrides -p
#   -l  short label for the generated XML's filename, default "run"
#   -k  1 = save stdout/stderr under <outdir>/log/, 0 = discard, default 0
#   -f  comma-separated FCS trigger flag names (see FcsTriggerDefs.h); an
#       event only gets a jetTree entry if one fired. Empty = keep every event.
#   -c  real Run 22 run number to pin the FCS gain/gainCorrection calibration
#       era to (see RunNumberToDate() in runSimBfc.C). 0/omit = BFC's own
#       nominal date.
#   -e  FcsJetFilter's generator-level forward-flux accept threshold, GeV.
#       An event is kept only if projected flux onto one FCS arm exceeds
#       this (StRoot/StarGenerator/FILT/FcsJetFilter.cxx). Default 50.0,
#       the original hardcoded value before this flag existed. Pass 0 to
#       disable the filter (accept every generated event).
#
# Example:
#   ./submit_sim_to_jet.sh -g pythia8 -t 8 -p 10 -n 500 -s 0 \
#       -o /star/data01/pwg/seanp/pipeline_output/cteq6l1_pt10_test \
#       -f "fcsJP2,fcsJPA0,fcsJPA1,fcsJPBC0,fcsJPBC1,fcsJPDE0,fcsJPDE1" \
#       -c 23101043 -e 50.0

set -e

GENERATOR=""
TUNE_PARAM=""
PTCUT=""
NEVENTS=""
SAVE_SIMPLETREE=""
OUTDIR=""
NPROC=1
PTCUT_JSON=""
LABEL="run"
KEEPLOGS=0
TRIGGER_FILTER=""
CALIB_RUN=0
FILTER_ETHR=50.0

while getopts "g:t:p:n:s:o:j:w:l:k:f:c:e:h" opt; do
  case $opt in
    g) GENERATOR="$OPTARG" ;;
    t) TUNE_PARAM="$OPTARG" ;;
    p) PTCUT="$OPTARG" ;;
    n) NEVENTS="$OPTARG" ;;
    s) SAVE_SIMPLETREE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    j) NPROC="$OPTARG" ;;
    w) PTCUT_JSON="$OPTARG" ;;
    l) LABEL="$OPTARG" ;;
    k) KEEPLOGS="$OPTARG" ;;
    f) TRIGGER_FILTER="$OPTARG" ;;
    c) CALIB_RUN="$OPTARG" ;;
    e) FILTER_ETHR="$OPTARG" ;;
    h|*)
      sed -n '2,34p' "$0"
      exit 0
      ;;
  esac
done

if [ -z "$GENERATOR" ] || [ -z "$TUNE_PARAM" ] || [ -z "$NEVENTS" ] || [ -z "$SAVE_SIMPLETREE" ] || [ -z "$OUTDIR" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ -z "$PTCUT" ] && [ -z "$PTCUT_JSON" ]; then
    echo "Missing required argument: either -p <ptcut> or -w <json> is required. Run with -h for usage."
    exit 1
fi

if [ "$GENERATOR" != "pythia8" ] && [ "$GENERATOR" != "pythia6" ]; then
    echo "Error: -g must be 'pythia8' or 'pythia6' (got '$GENERATOR')"
    exit 1
fi

if [ -n "$PTCUT_JSON" ] && [ ! -f "$PTCUT_JSON" ]; then
    echo "Error: -w json file '$PTCUT_JSON' does not exist"
    exit 1
fi

if ! [[ "$CALIB_RUN" =~ ^[0-9]+$ ]]; then
    echo "Error: -c must be a plain run number (got '$CALIB_RUN')"
    exit 1
fi

# SandBox <File> paths in sim_to_jet.xml are pipeline-root-relative, so
# star-submit must run from there, not from here.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/sim_to_jet_${LABEL}.xml"

# Fail fast on a typo'd trigger name, rather than inside a job with discarded logs.
if [ -n "$TRIGGER_FILTER" ]; then
    VALID_NAMES=$(awk '/kTrigFlagName\[kNTrigFlags\]/{flag=1; next} flag && /};/{flag=0} flag' "${SCRIPT_DIR}/FcsTriggerDefs.h" | sed -e 's/^[[:space:]]*"//' -e 's/",$//')
    IFS=',' read -ra REQUESTED <<< "$TRIGGER_FILTER"
    for name in "${REQUESTED[@]}"; do
        if ! grep -qxF "$name" <<< "$VALID_NAMES"; then
            echo "Error: unknown trigger filter name '$name' (not in FcsTriggerDefs.h)."
            echo "Valid names:"
            echo "$VALID_NAMES" | tr '\n' ' '
            echo
            exit 1
        fi
    done
fi

if [ "$KEEPLOGS" = "1" ]; then
    mkdir -p "${OUTDIR}/log"
    STDOUT_URL="file:${OUTDIR}/log/simtojet_\$JOBID.out"
    STDERR_URL="file:${OUTDIR}/log/simtojet_\$JOBID.err"
else
    STDOUT_URL="file:/dev/null"
    STDERR_URL="file:/dev/null"
fi

export XML_IN="${SCRIPT_DIR}/sim_to_jet.xml"
export XML_OUT="${OUTXML}"
export GENERATOR TUNE_PARAM PTCUT NEVENTS SAVE_SIMPLETREE OUTDIR NPROC PTCUT_JSON STDOUT_URL STDERR_URL TRIGGER_FILTER CALIB_RUN FILTER_ETHR

python3 << 'PYEOF'
import json
import os

with open(os.environ["XML_IN"]) as f:
    xml = f.read()

ptcut_json = os.environ.get("PTCUT_JSON", "")
nproc = os.environ["NPROC"]

if ptcut_json:
    with open(ptcut_json) as f:
        data = json.load(f)
    ref_total = data["reference_total"]
    bins = data["bins"]

    lines = []
    lines.append("        set round_jobs = {}".format(nproc))
    lines.append("        set reference_jobs = {}".format(ref_total))
    lines.append("        @ scaled_index = ${JOBINDEX} * ${reference_jobs} / ${round_jobs}")
    lines.append("")

    cum = 0
    for i, b in enumerate(bins[:-1]):
        cum += b["count"]
        kw = "if" if i == 0 else "else if"
        lines.append("        {} (${{scaled_index}} &lt; {}) then".format(kw, cum))
        lines.append("            set ptcut = {}".format(b["ptcut"]))
    lines.append("        else")
    lines.append("            set ptcut = {}".format(bins[-1]["ptcut"]))
    lines.append("        endif")

    ptcut_block = "\n".join(lines)
else:
    ptcut_block = "        set ptcut = {}".format(os.environ["PTCUT"])

xml = xml.replace("{{GENERATOR}}", os.environ["GENERATOR"])
xml = xml.replace("{{TUNE_PARAM}}", os.environ["TUNE_PARAM"])
xml = xml.replace("{{PTCUT_BLOCK}}", ptcut_block)
xml = xml.replace("{{NEVENTS}}", os.environ["NEVENTS"])
xml = xml.replace("{{SAVE_SIMPLETREE}}", os.environ["SAVE_SIMPLETREE"])
xml = xml.replace("{{OUTDIR}}", os.environ["OUTDIR"])
xml = xml.replace('nProcesses="1"', 'nProcesses="{}"'.format(nproc))
xml = xml.replace("{{STDOUT_URL}}", os.environ["STDOUT_URL"])
xml = xml.replace("{{STDERR_URL}}", os.environ["STDERR_URL"])
xml = xml.replace("{{TRIGGER_FILTER}}", os.environ.get("TRIGGER_FILTER", ""))
xml = xml.replace("{{CALIB_RUN}}", os.environ.get("CALIB_RUN", "0"))
xml = xml.replace("{{FILTER_ETHR}}", os.environ.get("FILTER_ETHR", "50.0"))

with open(os.environ["XML_OUT"], "w") as f:
    f.write(xml)
PYEOF

echo "Generated ${OUTXML}"
if [ -n "$PTCUT_JSON" ]; then
    echo "Generator=${GENERATOR} tune/pdf=${TUNE_PARAM} ptcut_json=${PTCUT_JSON} nevents=${NEVENTS} save_simpletree=${SAVE_SIMPLETREE} outdir=${OUTDIR} nProcesses=${NPROC} keep_logs=${KEEPLOGS} trigger_filter=${TRIGGER_FILTER:-<none>} calib_run=${CALIB_RUN} filter_ethr=${FILTER_ETHR}"
else
    echo "Generator=${GENERATOR} tune/pdf=${TUNE_PARAM} ptcut=${PTCUT} nevents=${NEVENTS} save_simpletree=${SAVE_SIMPLETREE} outdir=${OUTDIR} nProcesses=${NPROC} keep_logs=${KEEPLOGS} trigger_filter=${TRIGGER_FILTER:-<none>} calib_run=${CALIB_RUN} filter_ethr=${FILTER_ETHR}"
fi

mkdir -p "${OUTDIR}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "sim_to_jet/sim_to_jet_${LABEL}.xml"
SUBMIT_DIR="${PIPELINE_ROOT}/submit"
mkdir -p "${SUBMIT_DIR}"

# sched*/.csh/.package files are read from PIPELINE_ROOT by their original
# path at job start, not transferred by condor -- sweeping immediately races
# job startup and kills jobs. Defer the sweep until condor_wait confirms
# this submission's own jobs are done, in the background.
REPORT=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*.report' -newer "${MARKER}")
REQID=$(basename "${REPORT}" .report | sed 's/^sched//')
# Read the real Log path from the .condor file rather than assuming
# /tmp/$USER/... -- a wrong guess makes condor_wait fail immediately, so
# also gate the sweep on its exit status.
CONDORFILE=$(find "${PIPELINE_ROOT}" -maxdepth 1 -name "sched${REQID}_*.condor" -newer "${MARKER}" | head -1)
CONDORLOG=$(grep -m1 '^Log' "${CONDORFILE}" | sed -e 's/^Log[[:space:]]*=[[:space:]]*//')
nohup bash -c "if condor_wait '${CONDORLOG}' >/dev/null 2>&1; then find '${PIPELINE_ROOT}' -maxdepth 1 -name 'sched${REQID}*' -exec mv -t '${SUBMIT_DIR}' {} +; else echo \"condor_wait failed for ${REQID} (log='${CONDORLOG}') -- NOT sweeping, leaving sched files in place\" >> '${PIPELINE_ROOT}/submit/.sweep_failures.log'; fi" >/dev/null 2>&1 &
disown
rm -f "${MARKER}"
