#!/bin/bash
# Compiles and runs this folder's CalibrationMap on a jet_calibration.root
# "bins" tree, single process, producing the JES lookup JSON. For a faster
# run use submit_calibration_map_diag_condor.sh (condor shards).
#
# Usage:
#   ./run_calibration_map_diag.sh -i <jet_calibration.root> -o <output.json> [-E <0|1>]
#
#   -i  absolute path to a jet_calibration.root (from
#       run_jet_energy_scale_diag.sh or production JetEnergyScaleFineGrid)
#   -o  output JSON path (created if its directory doesn't exist)
#   -E  1 = the input came from an EM-only JetEnergyScaleFineGrid run
#       (narrows the skew-Gaussian position/width domain to match its Emax=43), 0 =
#       default ECAL+HCAL domain (Emax=200). Default 0.
#
# Example:
#   ./run_calibration_map_diag.sh \
#       -i em_only_jet_calibration_diag.root -o em_only_lookup_diag.json -E 1

set -e

INPUT=""
OUTPUT=""
EM_ONLY=0

while getopts "i:o:E:h" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    E) EM_ONLY="$OPTARG" ;;
    h|*)
      sed -n '2,18p' "$0"
      exit 0
      ;;
  esac
done

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: input file '$INPUT' does not exist"
    exit 1
fi

if [ "$EM_ONLY" != "0" ] && [ "$EM_ONLY" != "1" ]; then
    echo "Error: -E must be 0 or 1 (got '$EM_ONLY')"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"

mkdir -p "$(dirname "$OUTPUT")" 2>/dev/null || true
OUTPUT="$(cd "$(dirname "$OUTPUT")" && pwd)/$(basename "$OUTPUT")"

TMPCSH="$(mktemp /tmp/run_calibration_map_diag.XXXXXX.csh)"
trap 'rm -f "$TMPCSH"' EXIT

cat > "$TMPCSH" << EOF
#!/bin/csh
module purge
module load root-5.34.38
cd ${SCRIPT_DIR}
g++ -m64 CalibrationMap.cpp -o .calibration_map_diag_bin \`root-config --cflags --libs\` -lMinuit2 -std=c++11
if ( \$status != 0 ) then
    echo "ERROR: compilation failed"
    exit 1
endif
./.calibration_map_diag_bin "${INPUT}" "${OUTPUT}" ${EM_ONLY} 1
set RUNSTATUS = \$status
rm -f .calibration_map_diag_bin
exit \$RUNSTATUS
EOF

singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 -B /tmp /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh "$TMPCSH"
