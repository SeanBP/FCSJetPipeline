# UnderlyingEventPlayground

Workshop area for an underlying-event subtraction for FCS jets. Not part of
the pipeline -- nothing here is used by `sim_to_jet`, `data_to_jet`,
`jet_scale` or `jet_calibration`.

- `fms_emjet_reference.txt` -- parameters of the STAR FMS EM-jet analysis
  (BrightSTAR, branch `emjet-gpc`) that this study is meant to complement:
  jet definition, truth-level cuts, the off-axis-cone UE subtraction, analysis
  cuts, energy correction, simulation setup, and where each differs from this
  pipeline.
- `min_opposite_angle.cpp` -- for every point in the EM-only reco fiducial
  region, the smallest azimuthal separation a point on the opposite half can
  have from it (180 deg = back-to-back), with and without requiring the two
  points to be at least 2R apart in eta-phi. Build/run line is at the top of
  the file; writes `min_opposite_angle.pdf/.png/.root`.
- `max_dr_jet_placement.cpp` -- for every point P in the fiducial region, how
  far a jet axis J (either half) can be kept, in eta-phi, from both P and the
  point 180 deg in azimuth from P, while having zero overlap with a jet at P
  (dR(J, P) >= 2R). Draws the map and the worst case, with the three R-radius
  cones projected on the surface and the jet-axis regions allowed by a 0.1,
  0.2 and 0.3 R tolerance. Writes `max_dr_jet_placement.pdf/.png/.root`.
- `UESelection.h` -- the UE event-selection geometry (P', jet distance, best
  achievable distance, the 2R and 0.3 R constants), shared by
  `max_dr_jet_placement.cpp` and `ue_map.cpp` so both use one definition.
- `ue_map.cpp` -- builds the UE energy map on half the FCS ECAL (x > 0, the
  JES fine grid) from `sim_to_jet -E 1 -F 1 -s 1` output: per grid point,
  events with no reco jet within 2R and a tag jet in the 0.3 R region, the
  summed point energy within R (mean, sigma), and the truth-jet
  contamination within R (fraction of events, mean, sigma). Every event is
  used as is and left-right mirrored, and weighted by cross section. Means
  and sigmas come from a Gaussian fit where it can be trusted, otherwise
  numerically; the `*_method` grids flag which. Build/run line is at the top
  of the file.

## Working choices

- The UE correction is an energy subtraction (not pT), from an averaged map
  in position, to be built separately for data, simulated reco and truth.
  No pT binning yet.
- Event selection for the UE cone at a point P: no reco jet within 2R of P,
  and a reco jet inside the region found by `max_dr_jet_placement.cpp` with
  a 0.3 R tolerance (distance to the nearer of P and P' within 0.3 R of the
  best achievable).
- The cone sums the same points the jets are built from (ECAL, 1 GeV to
  sqrt(s)/2), so sub-threshold deposits are not counted.
- Truth-jet contamination of the cone is recorded as a diagnostic.
