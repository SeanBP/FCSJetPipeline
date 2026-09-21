#!/bin/bash
# Compiles and runs the diagnostics fork of JetEnergyScaleFineGrid (this
# folder's copy -- same algorithm as production, plus an EM-only-
# conditional Emax, see its own header comment) against a folder of
# sim_to_jet JetTrees, producing a jet_calibration.root "bins" tree.
# Always rebuilds fresh.
#
# Usage:
#   ./run_jet_energy_scale_diag.sh -i <jettrees_dir> -o <output.root> [-E <0|1>] [-r <match_radius>] [-f <E0_floor>] [-x <Emax_ceiling>]
#
#   -i  absolute path to a folder of jet_output_*.root (JetTrees) files
#   -o  output ROOT path (created if its directory doesn't exist)
#   -E  1 = EM-only JetTrees (ECAL-only fiducial rectangle/buffer), 0 = default
#       ECAL+HCAL fiducial rectangle. Default 0.
#   -r  override match_radius (cm) -- the pooling radius around each grid
#       point (see JetEnergyScaleFineGrid.cpp's own declaration). Default
#       5.0 (the file's hardcoded value) if omitted.
#   -f  override E0 (GeV) -- the truth-energy binning floor (see
#       JetEnergyScaleFineGrid.cpp's own declaration). Default 10 (the
#       file's hardcoded value) if omitted.
#   -x  override Emax (GeV) -- the truth-energy binning ceiling (see
#       JetEnergyScaleFineGrid.cpp's own declaration). Default 70 for
#       EM-only / 200 otherwise (the file's hardcoded values) if omitted.
#
# Example:
#   ./run_jet_energy_scale_diag.sh \
#       -i /star/data01/pwg/seanp/pipeline_output/em_only_production \
#       -o em_only_jet_calibration_diag.root -E 1

set -e

INPUT=""
OUTPUT=""
EM_ONLY=0
MATCH_RADIUS=""
E0_FLOOR=""
EMAX_CEIL=""

while getopts "i:o:E:r:f:x:h" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    E) EM_ONLY="$OPTARG" ;;
    r) MATCH_RADIUS="$OPTARG" ;;
    f) E0_FLOOR="$OPTARG" ;;
    x) EMAX_CEIL="$OPTARG" ;;
    h|*)
      sed -n '2,26p' "$0"
      exit 0
      ;;
  esac
done

# argv positions matter (JetEnergyScaleFineGrid.cpp reads them in order) --
# fill in the file's own defaults for any earlier slot left blank so a later
# override doesn't land in the wrong argv position.
if [ -n "$E0_FLOOR" ] || [ -n "$EMAX_CEIL" ]; then
    if [ -z "$MATCH_RADIUS" ]; then
        MATCH_RADIUS="5.0"
    fi
fi
if [ -n "$EMAX_CEIL" ] && [ -z "$E0_FLOOR" ]; then
    E0_FLOOR="10"
fi

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Missing required argument. Run with -h for usage."
    exit 1
fi

if [ ! -d "$INPUT" ]; then
    echo "Error: input directory '$INPUT' does not exist"
    exit 1
fi

if [ "$EM_ONLY" != "0" ] && [ "$EM_ONLY" != "1" ]; then
    echo "Error: -E must be 0 or 1 (got '$EM_ONLY')"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_DIR="$(cd "${SCRIPT_DIR}/../shared" && pwd)"
INPUT="$(cd "$INPUT" && pwd)"

mkdir -p "$(dirname "$OUTPUT")" 2>/dev/null || true
OUTPUT="$(cd "$(dirname "$OUTPUT")" && pwd)/$(basename "$OUTPUT")"

TMPCSH="$(mktemp /tmp/run_jet_energy_scale_diag.XXXXXX.csh)"
trap 'rm -f "$TMPCSH"' EXIT

cat > "$TMPCSH" << EOF
#!/bin/csh
module purge
module load fastjet-3.3.4
module load root-5.34.38
cd ${SCRIPT_DIR}
g++ -m64 JetEnergyScaleFineGrid.cpp ${SHARED_DIR}/VectorDict.cxx -o .jet_energy_scale_diag_bin \\
    -I${SHARED_DIR} \`root-config --cflags --libs\` \`fastjet-config --cxxflags --libs\` -std=c++11
if ( \$status != 0 ) then
    echo "ERROR: compilation failed"
    exit 1
endif
./.jet_energy_scale_diag_bin "${INPUT}/" ${EM_ONLY} ${MATCH_RADIUS} ${E0_FLOOR} ${EMAX_CEIL}
set RUNSTATUS = \$status
mv jet_calibration.root "${OUTPUT}"
rm -f .jet_energy_scale_diag_bin
exit \$RUNSTATUS
EOF

singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 -B /tmp /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh "$TMPCSH"

echo "Output written to ${OUTPUT}"
