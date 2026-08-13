# GeometryPlayground

Minimal, standalone, run-number -> FCS detector-geometry tool. Given a real
Run 22 run number, pulls the FCS ECAL/HCAL detector-position geometry for
that run's calibration era directly from the STAR conditions DB
(`Geometry/fcs/fcsDetectorPosition`) and computes the fiducial-region
rectangle used to cut reconstructed jets near the detector edge -- no
simulation, no input event file, nothing but a run number.

This is where the `cut_x_inner`/`cut_x_outer`/`cut_y_min`/`cut_y_max` values
baked into `FCSJetPipeline/shared/JetParameters.h` came from, and the tool
to re-derive/verify them if the detector position ever changes.

## Files

- **`get_fcs_geometry.C`** -- Stage 1. A minimal ROOT/STAR BFC chain (DB +
  `StFcsDbMaker` only) that resolves a run number to its real calendar date
  (via the `RunLog.runDescriptor` DB, same mechanism
  `StFcsDbMaker/macro/fcsTimedepGainCorr_db.C` uses), points `St_db_Maker`
  at that date, and dumps the 4 corner `StarXYZ` positions of each of ECAL
  north/south and HCAL north/south to `fcs_corners.txt`. Flags (does not
  block) dates before 2021-12-21, when the FCS detector-position survey was
  finalized in the DB -- see `kFcsSurveyValidFromDate` in the file.
- **`fiducial_from_corners.cpp`** -- Stage 2. Plain C++ (no ROOT/STAR
  dependency, compiles standalone with `g++`). Reads `fcs_corners.txt`,
  projects the HCAL corners onto the ECAL plane, and reports each side's
  fiducial rectangle (HCAL's own projected footprint -- see the file's
  header comment for why that, rather than the true ECAL/HCAL overlap
  polygon, is the right definition, and for the built-in cross-check that
  would flag it if that assumption ever stopped holding).
- **`run_number_to_fiducial.sh`** -- Runs both stages end to end.
- **`fcs_corners_run22.txt`** -- saved stage-1 output (corner positions) for
  the one geometry that covers essentially all of Run 22, produced from run
  22359013 (an arbitrary representative run inside the stable window -- see
  the file's own header). Feed straight into stage 2:
  `./fiducial_from_corners fcs_corners_run22.txt`.
- **`fcs_detector_offsets_run22.txt`** -- saved `getDetectorOffset()` values
  (full double precision) for the same run/geometry, one line per detector.
  This is the raw data behind the north/south symmetry finding below (x
  mirrored, y identical between north/south, for both ECAL and HCAL).

## Usage

```
./run_number_to_fiducial.sh <run_number>
```

e.g. `./run_number_to_fiducial.sh 22359013`. Requires the STAR/singularity
environment for stage 1 only (`fiducial_from_corners` is compiled on first
run with a plain `g++`, no container needed for stage 2). The script
launches its own `singularity exec ... stardev ...` for stage 1, so it can
be run directly from a bare shell.

## Findings this tool established (Aug 2026)

- The detector-position table is constant across essentially all of Run 22
  (checked from the 2021-12-21 survey date through the nominal `y2023` BFC
  date) -- geometry only needs re-deriving if a documented detector move
  happens in a later run period, not per run.
- North and south were checked for asymmetry (full double-precision DB
  offset comparison, plus extremal hit positions in a real muDst) and found
  exactly symmetric in both the idealized DB geometry and real reconstructed
  hit extents (hit *counts* differ ~6%, but that's occupancy, not geometry).
  `JetParameters.h` therefore uses a single shared rectangle applied to both
  sides via `fabs(x)`, not separate north/south values.

## Requesting a run before 2021-12-21

The DB returns a pre-survey placeholder position for these dates (round
numbers, not the real calibrated position) -- `get_fcs_geometry.C` prints a
`WARNING:` line but does not block. Do not trust the resulting fiducial
numbers for such a run.
