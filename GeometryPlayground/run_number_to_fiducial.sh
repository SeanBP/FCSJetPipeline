#!/bin/bash
# Pure run-number -> FCS fiducial-boundary tool.
#
# Usage: ./run_number_to_fiducial.sh <run_number>
#
# Given a real Run 22 run number, prints the north/south ECAL/HCAL
# fiducial rectangles (JetParameters.h-ready cut_x_inner/outer_north/south,
# cut_y_min/max_north/south) for that run's actual calibration era, pulled
# directly from the STAR conditions DB (Geometry/fcs/fcsDetectorPosition) --
# no simulation or input file needed.
#
# Two stages under the hood:
#   1. get_fcs_geometry.C (ROOT/CINT, needs the STAR/singularity
#      environment) -- reads StFcsDb for the given run's date and dumps
#      the 4 corners each of ECAL north/south, HCAL north/south to
#      fcs_corners.txt. Also flags (does not block) runs before the FCS
#      detector-position survey was finalized in the DB (2021-12-21).
#   2. fiducial_from_corners (plain C++, no ROOT/STAR dependency) -- turns
#      those corners into the actual fiducial rectangles, per side.
#
# NOTE (Aug 2026 investigation): the detector position table was found to
# be constant across essentially all of Run 22 (checked from 2021-12-21
# through the y2023 nominal date) -- so in practice this only needs to be
# run once per run *period* that includes a real detector move, not per
# run. Re-run it if you have reason to think the geometry changed (a new
# run period, a documented detector move) and want to confirm before
# trusting cached values.

set -e

if [ -z "$1" ]; then
    echo "Usage: $0 <run_number>"
    exit 1
fi
RUN_NUMBER="$1"

if ! [[ "$RUN_NUMBER" =~ ^[0-9]+$ ]]; then
    echo "Error: run_number must be a plain number (got '$RUN_NUMBER')"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ ! -x ./fiducial_from_corners ]; then
    echo "Building fiducial_from_corners..."
    g++ -std=c++11 -O2 -o fiducial_from_corners fiducial_from_corners.cpp
fi

cat > .run_number_to_fiducial_stage1.csh << EOF
#!/bin/csh
cd ${SCRIPT_DIR}
source /star/nfs4/AFS/star/group/star_cshrc.csh
stardev
root4star -b -q get_fcs_geometry.C\\(${RUN_NUMBER}\\)
EOF

echo "=== Stage 1: querying StFcsDb for run ${RUN_NUMBER} ==="
STAGE1_LOG=".run_number_to_fiducial_stage1.log"
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
    /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif \
    tcsh .run_number_to_fiducial_stage1.csh > "$STAGE1_LOG" 2>&1

if ! grep -qa "Wrote fcs_corners.txt" "$STAGE1_LOG"; then
    echo "Error: stage 1 (get_fcs_geometry.C) did not complete successfully."
    echo "See $STAGE1_LOG"
    exit 1
fi

if grep -qa "WARNING: date=" "$STAGE1_LOG"; then
    echo ""
    grep -a "WARNING: date=" "$STAGE1_LOG"
    echo ""
fi

echo "=== Stage 2: computing fiducial rectangles ==="
./fiducial_from_corners fcs_corners.txt

rm -f .run_number_to_fiducial_stage1.csh
