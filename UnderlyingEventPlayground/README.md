# UnderlyingEventPlayground

Studies behind the underlying-event (UE) subtraction in `../underlying_event/`
(the in-event pixel ring, `PixelRingUE.h`, `UnderlyingEvent.cpp`); nothing
here is used by the pipeline. Each source has its build and run line at the
top.

The earlier methods (the near off-axis proxy-event search of
`underlying_event/ProxySelection.h` and `ProxySearch.h`, the separation and
veto tests, the truth decomposition of the proxy cone, the position and ptHat
studies) were retired on 2026-10-30. Their sources are in the git history;
sources, plots and binaries are also in
`/star/data01/pwg/seanp/archive/FCSJetPipeline_superseded_2026-10-30.tar.gz`.

## Method

Per jet, from the same event: the jet's eta band |eta - eta_jet| < R, at
least 2R from the jet axis and within 150 deg of it in azimuth, other jets not
excluded, weighted per point so that its area-per-eta matches the jet cone's.
The UE density times the jet's live cone area (reco) or pi R^2 (truth) is
subtracted in energy, per jet, at both levels, before the JES calibration.
Details in `../underlying_event/PixelRingUE.h` and the top-level README.

## Studies

- `fms_emjet_reference.txt`: the STAR FMS EM-jet analysis (BrightSTAR
  `emjet-gpc`) and how this pipeline compares.
- `in_event/`: the in-event estimate, from same-event cones to the pixel
  ring. `in_event_ue.cpp` (condor: `submit_in_event.sh`) writes per-jet
  estimates for every variant; `in_event_geo.h` holds the geometry (FCS ECal
  pixels, the ring variants, the area matching) and `in_event_util.h` small
  helpers. Then:
  - `in_event_compare.cpp`, `leak_ring.cpp`, `trig_compare.cpp`,
    `emu_compare.cpp`: data against simulation, jet leakage, trigger levels;
  - `ann_compare.cpp`, `ann_edge.cpp`, `ann_phi_share.cpp`, `sep_scan.cpp`,
    `annulus_variants.cpp`: the eta-weighted annulus, its edge cases and the
    separation from the jet;
  - `px_check.cpp`, `px_weight_check.cpp`: the pixel ring, per-pixel against
    per-point weights;
  - `ring_cut_check.cpp`, `ring_cap_check.cpp`: the cut around the recoil
    side (ring area kept, eta coverage, density for caps 180/150/120 deg);
  - `truth_scan.cpp`, `truth_scan_plot.cpp`: the azimuthal wedge and cone
    scans at truth level (FCS simulation), as the FMS test
    (`FMSJetPipeline/compare/wedge_density.cpp`);
  - `annulus_diagram.cpp`: diagram of the ring for one jet.
- `recoil_truth/`: generator-level Pythia6 (FMS tune) recoil test: where
  the recoil parton and its energy fall relative to the jet's eta band and
  the off-axis cones, with and without multi-parton interactions. Needs
  `pythia6428_nopdf.f` (the STAR Pythia 6.4.28 source without its dummy
  PDFSET/STRUCTM; not in the repository).
- `energy_flow/cal_chain_compare.cpp`: the ECal reconstruction chain (hits,
  clusters, points), data against simulation.

## Results

The method and its numbers (2026-10): azimuthal scans on FMS Run 15 data
and FCS simulation put the end of the jet's own energy at about 2R and the
start of the recoil rise at about 120-130 deg. On the FMS the ring on FCS-like
geometry is within about 10% of the FMS two-cone (+-90 deg) estimate (full
FMS 4%), taken as the UE systematic; at truth level the ring is within 2% of
the two cones. A 120 deg cap loses eta coverage (32% of FCS data jets below
0.9) without lowering the density, so 150 deg is kept. Per-point and per-pixel
weights agree within 0.5% on Run 22.
