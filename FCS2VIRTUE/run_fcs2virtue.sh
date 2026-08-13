#!/bin/bash
# Compiles and runs FCS2VIRTUE, converting a SimpleTree file into a
# VIRTUE event-display JSON. Always rebuilds fresh (fast, single-file
# compile) so it never runs against a stale binary.
# Usage:
#   ./run_fcs2virtue.sh -i <simpletree.root> -n <n_events> -o <output.json> \
#       [-t] [-T] [-r] [-j <jettree.root>]
#
#   -i  absolute path to the input SimpleTree file
#   -n  number of events to convert
#   -o  output JSON path (created if its directory doesn't exist)
#   -t  include mcparticle tracks
#   -T  include truth jets (green)
#   -r  include reco jets (orange)
#   -j  jetTree file to read jets from; required if -T and/or -r is given
#       (truth and reco jets always come from the same file)
#
# Examples:
#   Hits/blocks only:
#     ./run_fcs2virtue.sh -i SimpleTree_pt10_0.root -n 10 -o out.json
#
#   Hits + tracks + both jet types:
#     ./run_fcs2virtue.sh -i SimpleTree_pt10_0.root -n 10 -o out.json \
#         -t -T -r -j jet_output_0.root

set -e

INPUT=""
NEVENTS=""
OUTPUT=""
DO_TRACKS=0
DO_TRUTH=0
DO_RECO=0
JETSFILE=""

while getopts "i:n:o:tTrj:h" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    n) NEVENTS="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    t) DO_TRACKS=1 ;;
    T) DO_TRUTH=1 ;;
    r) DO_RECO=1 ;;
    j) JETSFILE="$OPTARG" ;;
    h|*)
      sed -n '2,24p' "$0"
      exit 0
      ;;
  esac
done

if [ -z "$INPUT" ] || [ -z "$NEVENTS" ] || [ -z "$OUTPUT" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: input SimpleTree file '$INPUT' does not exist"
    exit 1
fi

if { [ "$DO_TRUTH" = "1" ] || [ "$DO_RECO" = "1" ]; } && [ -z "$JETSFILE" ]; then
    echo "Error: -T/-r requires -j <jettree.root>. Run with -h for usage."
    exit 1
fi

if [ -n "$JETSFILE" ] && [ ! -f "$JETSFILE" ]; then
    echo "Error: jets file '$JETSFILE' does not exist"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FLAGS=""
[ "$DO_TRACKS" = "1" ] && FLAGS="$FLAGS --tracks"
[ "$DO_TRUTH" = "1" ] && FLAGS="$FLAGS --truth-jets"
[ "$DO_RECO" = "1" ] && FLAGS="$FLAGS --reco-jets"
[ -n "$JETSFILE" ] && FLAGS="$FLAGS --jets-file ${JETSFILE}"

mkdir -p "$(dirname "$OUTPUT")"

TMPCSH="$(mktemp /tmp/run_fcs2virtue.XXXXXX.csh)"
trap 'rm -f "$TMPCSH"' EXIT

cat > "$TMPCSH" << EOF
#!/bin/csh
stardev
module purge
module load fastjet-3.3.4
module load root-5.34.38
cd ${SCRIPT_DIR}
g++ -m64 FCS2VIRTUE.cpp ../shared/VectorDict.cxx -o .fcs2virtue_bin \`root-config --cflags --libs\` \`fastjet-config --cxxflags --libs\` -std=c++11
if ( \$status != 0 ) then
    echo "ERROR: compilation failed"
    exit 1
endif
./.fcs2virtue_bin "${INPUT}" ${NEVENTS} "${OUTPUT}" ${FLAGS}
set RUNSTATUS = \$status
rm -f .fcs2virtue_bin
exit \$RUNSTATUS
EOF

singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 -B /tmp /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh "$TMPCSH"
