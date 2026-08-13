#!/bin/bash
# Generates and submits a concrete sim_to_jet.xml job.
# Usage:
#   ./submit_sim_to_jet.sh -g <pythia8|pythia6> -t <tune_param> -p <ptcut> \
#       -n <nevents> -s <0|1> -o <outdir> [-j <nprocesses>] [-l <label>] \
#       [-k <0|1>] [-f <triggers>] [-c <run_number>]
#
#   -g  generator: pythia8 or pythia6
#   -t  tune/PDF param: pythia8 PDF:pSet int (8=CTEQ6L1, 5=MSTW2008LO,
#       3=MRST LO*, 21=NNPDF3.1sx, <=0=default) or pythia6 PyTune int
#       (325=Perugia STAR, 0=Pythia6 default)
#   -p  ptHatMin cut in GeV (single flat value for every job; ignored if
#       -w is given)
#   -n  events per job
#   -s  1 to keep the intermediate SimpleTree, 0 to discard it
#   -o  output directory (created if missing)
#   -j  number of parallel jobs (nProcesses), default 1
#   -w  optional: absolute path to a JSON file giving the relative
#       proportions of ptHatMin cuts to sample across the -j jobs (see
#       pthat_distribution_optimized.json for the format/an example);
#       overrides -p
#   -l  short label used only in the generated XML's filename, default "run"
#   -k  1 to save each job's stdout/stderr under <outdir>/log/ (auto-
#       created), 0 to discard them, default 0
#   -f  comma-separated FCS trigger flag names (see FcsTriggerDefs.h /
#       TriggerIDs.txt, e.g. "fcsJP2,fcsJPA0,fcsJPA1,fcsJPBC0,fcsJPBC1,
#       fcsJPDE0,fcsJPDE1" for jet-patch triggers) -- an event only gets a
#       jetTree entry if at least one of them fired. Omit/empty (default)
#       keeps every event; the SimpleTree itself is never filtered.
#       (Called -f here, not -t like data_to_jet's equivalent flag,
#       since -t is already taken by the tune/PDF param above.)
#   -c  a real Run 22 run number (e.g. 23101043) whose FCS gain/
#       gainCorrection calibration era should be used for the simulated
#       detector response, instead of BFC's own nominal "y2023" date.
#       Omit/0 (default) = don't override. See RunNumberToDate() in
#       runSimBfc.C for the run-number -> calendar-date conversion
#       (verified empirically -- the run number's year code is NOT the
#       real calendar year, e.g. 22359013 is really 2021-12-25).
#
# Example:
#   ./submit_sim_to_jet.sh -g pythia8 -t 8 -p 10 -n 500 -s 0 \
#       -o /star/data01/pwg/seanp/pipeline_output/cteq6l1_pt10_test \
#       -f "fcsJP2,fcsJPA0,fcsJPA1,fcsJPBC0,fcsJPBC1,fcsJPDE0,fcsJPDE1" \
#       -c 23101043

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

while getopts "g:t:p:n:s:o:j:w:l:k:f:c:h" opt; do
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
    h|*)
      sed -n '2,44p' "$0"
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

# This script lives in sim_to_jet/. sim_to_jet.xml's SandBox <File> paths
# are relative to the pipeline root (SUMS packages each <File> preserving
# that relative path, e.g. sim_to_jet/runSimBfc.C -- confirmed directly
# from a real job's .package.zip), so star-submit MUST run from there,
# not from a dedicated submit folder. Instead, the sched*/*.package/*.csh/
# *.list/*.condor files star-submit drops in the CWD it runs from are
# swept into submit/ (a dedicated folder) right after submission -- see
# the marker-file logic below.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "${SCRIPT_DIR}")"
OUTXML="${SCRIPT_DIR}/sim_to_jet_${LABEL}.xml"

# Fail fast on a typo'd trigger name here, rather than have it silently
# error out inside a job whose logs are discarded by default (-k 0).
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
export GENERATOR TUNE_PARAM PTCUT NEVENTS SAVE_SIMPLETREE OUTDIR NPROC PTCUT_JSON STDOUT_URL STDERR_URL TRIGGER_FILTER CALIB_RUN

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

with open(os.environ["XML_OUT"], "w") as f:
    f.write(xml)
PYEOF

echo "Generated ${OUTXML}"
if [ -n "$PTCUT_JSON" ]; then
    echo "Generator=${GENERATOR} tune/pdf=${TUNE_PARAM} ptcut_json=${PTCUT_JSON} nevents=${NEVENTS} save_simpletree=${SAVE_SIMPLETREE} outdir=${OUTDIR} nProcesses=${NPROC} keep_logs=${KEEPLOGS} trigger_filter=${TRIGGER_FILTER:-<none>} calib_run=${CALIB_RUN}"
else
    echo "Generator=${GENERATOR} tune/pdf=${TUNE_PARAM} ptcut=${PTCUT} nevents=${NEVENTS} save_simpletree=${SAVE_SIMPLETREE} outdir=${OUTDIR} nProcesses=${NPROC} keep_logs=${KEEPLOGS} trigger_filter=${TRIGGER_FILTER:-<none>} calib_run=${CALIB_RUN}"
fi

mkdir -p "${OUTDIR}"

cd "${PIPELINE_ROOT}"
MARKER=$(mktemp)
star-submit "sim_to_jet/sim_to_jet_${LABEL}.xml"
SUBMIT_DIR="${PIPELINE_ROOT}/submit"
mkdir -p "${SUBMIT_DIR}"
find "${PIPELINE_ROOT}" -maxdepth 1 -name 'sched*' -newer "${MARKER}" -exec mv -t "${SUBMIT_DIR}" {} +
rm -f "${MARKER}"
