# FCSJetPipeline

SUMS pipelines for STAR FCS jets (Run 22 pp, sqrt(s) = 508 GeV), one folder
per stage. Each folder holds its XML template, a `submit_*.sh` wrapper (run it
with `-h` for its options) and its sources.

| Stage | Input -> output |
|---|---|
| `sim_to_jet/` | PYTHIA + starsim -> BFC/MuDst -> SimpleTree -> matched truth/reco jets (+ in-job underlying event) |
| `data_to_jet/` | official MuDst -> SimpleTree -> reco jets (+ in-job underlying event) |
| `underlying_event/` | underlying-event (UE) estimate per jet (in-event pixel ring); runs inside the two jobs above |
| `jet_scale/` | simulated JetTrees -> jet-energy-scale (JES) lookup JSON and its closure test |
| `jet_calibration/` | applies a JES JSON to JetTrees, in place |
| `pico_analysis/` | PicoDst diagnostics, not part of the jet chain |

Shared resources at the root:

- `shared/JetParameters.h`: jet radius, fiducial rectangles, thresholds, the
  EM-only truth-particle selection, sqrt(s) per run, x_F.
- `shared/VectorDict.*`: the `std::vector` dictionary every program reading
  jetTree vector branches links against.
- `StRoot/`: `StSimpleReaderMaker`, `StSimpleTreeMaker` and a local
  `FcsJetFilter` (settable threshold), and a local `StFcsFastSimulatorMaker`
  (the dev version with an ECal energy scale, `sim_to_jet -K`). `.sl73_x8664_gcc485/` is their
  prebuilt cache; jobs do not rebuild it, so rebuild (`stardev`, `cons`)
  after editing `StRoot/`.

The `*Playground/` folders and `FCS2VIRTUE/` (SimpleTree -> VIRTUE event
display JSON) hold studies and tools that the pipelines do not use. Analysis
plots live outside this repository, in `FCSJetAnalysis/`.

## Conventions

- Run every `submit_*.sh` from anywhere; it calls `star-submit` from the
  pipeline root, since the XML SandBox paths are relative to it. Generated
  XMLs (`<template>_<label>.xml`) and the `sched*` files are job-specific
  and are not tracked. `sched*` files are moved to `submit/` once
  `condor_wait` reports the submission finished.
- Logs are off by default; `-k 1` writes each job's stdout/stderr to
  `<outdir>/log/` while it runs. Keep it off for large batches.
- `sim_to_jet` and `data_to_jet` write `<outdir>/run_params_<label>_<id>.txt`
  with every resolved option, the submission time and the git commit.
- `-1` (or `-9999` for positions) marks a value that does not apply, e.g.
  `sigma_job` in data or `Spin_config` in simulation.
- FCS clusters and points are remade from the MuDst hits in both
  `runMudst.C` (data) and `readMudst.C` (simulation), with the same
  thresholds whatever library made the MuDst: `StFcsMuToEventMaker` puts the
  hits into an StEvent, `StFcsClusterMaker` clusters Ecal towers above
  0.4 GeV (6 sigma above the worst-case Run 22 pedestal noise, with margin)
  with seeds at 1 GeV (P24ia), and `StFcsPointMaker` fits the points.
  `StSimpleReaderMaker::SetFcsFromStEvent()` writes those, not the MuDst's;
  `Cal_pt_parentClusterId` is then the cluster's index within its detector
  (P24ia MuDsts leave it unset). The DEV library the jobs run seeds at
  0.05 GeV and P24ia keeps towers above ~0.045 GeV, so production clusters
  differ between data and simulation
  (`UnderlyingEventPlayground/energy_flow/cal_chain_compare.cpp`).

## sim_to_jet

```
./submit_sim_to_jet.sh -g <pythia8|pythia6> -t <tune> -p <ptHatMin> -n <nevents> -s <N> -o <outdir> \
    [-P <ptHatMax>] [-j <njobs>] [-w <json>] [-l <label>] [-k <0|1>] [-f <triggers>] [-c <run>] \
    [-e <ethr>] [-E <0|1>] [-F <0|1>] [-T <tune_pp>] [-M <0|1>] [-C <0|1>] [-m <0|1>] [-K <k>] [-S <seed>]
```

One job generates `-n` events and writes
`jet_output_<gen>_<tune>_pt<ptHatMin>_<job>.root`, and the SimpleTree for every
`-s`-th job (`-s 1` every job, `-s 0` none; a large sample keeps a fraction for
checks). Main options:

| Flag | Meaning |
|---|---|
| `-w` | JSON of relative ptHat-range shares across the `-j` jobs (`pthat_distribution_optimized.json`: disjoint 5 GeV bins, exact weights); replaces `-p`/`-P` |
| `-f` | keep only events that fired one of these FCS triggers (`FcsTriggerDefs.h` names, emulated with `fcs_trg_base`) |
| `-c` | run number fixing the FCS gain calibration era (default: BFC's nominal date) |
| `-e` | `FcsJetFilter` forward-flux threshold, GeV (default 0 = off: any filter biases low-energy EM-only jets) |
| `-E 1` | EM-only reco jets (ECAL only) with the ECAL fiducial rectangle |
| `-F 1` | reco jets from FCS points (`StMuFcsPoint`) instead of tower hits |
| `-M 1` | the FMS EM-jet generator settings: pythia6, PyTune 370, PARP(90)=0.213, PARP(91)=1, MSEL=1 |
| `-C 1` | JES calibration sample: no jet pT minimums; no `-f` |
| `-K` | ECal photon energy scale k applied to the deposited energy before digitization (local `StFcsFastSimulatorMaker`), default 0.947: the simulation's pi0 peak then sits at the PDG mass under the gain calibration's own selection, as in the calibrated data (`Pi0Playground/pi0_vs_pdg`) |
| `-S` | fixed generator seeds (seed base + job index); 0 = random |

Folder layout: after a submission finishes,
`organize_sim_output.sh <outdir>` moves the files into
`jet_output/pt<N>/` and `SimpleTree/pt<N>/`, so no folder holds thousands of
files.

**Jets.** Anti-kt, R = 0.7, with jet areas (`reco_area`, `truth_area`) as the
FMS analysis: fastjet active area, ghost area 0.04, |ghost rapidity| < 5
(`JetAreaDefinition()` in `JetParameters.h`). Reco jets must be inside the fiducial rectangle
minus a margin (10 cm, or 20 cm for EM-only). Matching is one-to-one within
dR < R/2. ECAL+HCAL truth jets may extend R/2 outside the reco rectangle.
EM-only truth jets are clustered from photons created within 100 um of the
primary vertex (prompt and pi0/eta decay photons, not GEANT secondaries).
Only truth jets that could be matched are kept: axis within R/2 in eta-phi of
the reco fiducial region. Photons are kept in a full-azimuth eta band 1.5R
beyond that region's eta range (2.34-3.57, so 1.29-4.62), which holds every
particle of a stored jet and of its UE cones. EM-only jet pT minimums: reco 2 GeV, truth 1 GeV (as the FMS
analysis), unless `-C 1`. FCS points must have 1 GeV < E <= sqrt(s)/2.

**Cross sections.** Every tree carries `sigma_job` (PYTHIA cross section,
mb), `N_job` (trials) and the job's ptHat range (`mc_pthatmin_gev`,
`mc_pthatmax_gev`). `jet_scale/JetXSecWeights.h` turns them into per-ptHat-bin
event weights. Pythia 6's cross section arrives byte-shifted from the STAR
interface and is decoded in `readMudst.C` (and again on read). Jobs where
PYTHIA 6 crashed carry no valid `sigma_job` and cannot be weighted.

**Triggers in simulation.** `fcs_trg_base` (via `StFcsTriggerSimMaker`,
whose own `Make()` is bypassed) computes 20 of the 64 trigger flags; the
rest stay 0.

## data_to_jet

```
perl MakeRunList.pl -f mylist.list -l <max_files>      # file catalog query (wave 1, good runs)
./submit_data_to_jet.sh -f <listfile> -s <0|1> -o <outdir> [-n <files_per_job>] [-j <njobs>] \
    [-l <label>] [-k <0|1>] [-t <triggers>] [-E <0|1>] [-F <0|1>]
```

The list is cut into jobs of `-n` (default 10) consecutive files; every file
is processed on its own (the UE comes from each jet's own event). `-j`
instead splits the list evenly across that many jobs. Each
SimpleTree gives one `jet_output_<run>_<file>.root` with the same suffix.

Per event the trees carry the fired FCS trigger flags (resolved from
trigger id and run number, `FcsTriggerDefs.h`, generated from
`TriggerIDs.txt` by `generate_trigger_defs.py`; three copies, regenerate
together and rebuild StRoot), `Spin_config` (spin DB, covers 2021-12-14 to
2022-04-18 only), the event vertex (TPC, VPD, then BBC) and sqrt(s).
`MakeRunList.pl` needs the Perl library in
`/star/u/dkap7827/Tools2/Tools/PerlScripts`.

## underlying_event

The UE of a jet comes from its own event, by the pixel ring
(`PixelRingUE.h`): the jet's eta band |eta - eta_jet| < R from 2R off the
jet axis out to 150 deg in azimuth (the 30 deg around 180, the recoil side,
are left out), other jets kept, over the whole ECal (the fiducial region is
only for jets). The ring is weighted so that its area-per-eta matches the jet
cone's: in eta bins of 0.05, a point weighs h_jet(b) / h_ring(b), the ratio
of the jet cone's to the ring's (live) area in its bin; bins the ring covers
under 10% of the jet's area carry no weight. The density is
sum(E r) / sum(r h_ring), so empty parts of the ring count as area, and the
UE is the density times the jet cone's area: its live area on the ECal at
reco level (3 cm pixels, 3 x 3 sub-pixels for the areas), pi R^2 at truth
level (EM-only truth photons, full azimuth, exact areas).

`UnderlyingEvent` runs inside the sim_to_jet and data_to_jet jobs, one file
at a time (no pooling of events, no trigger requirement), and rewrites each
jet file with, per jet: `reco_ue_rho`, `reco_ue_area`, `reco_ue_E` (= rho x
area, the UE to subtract), `reco_ue_pt`, `reco_ue_cov` (share of the jet's
live area in eta bins the ring reaches) and, in simulation, `truth_ue_rho`,
`truth_ue_E`, `truth_ue_pt`. The energies are -1 where the ring reaches none
of the jet's eta bins; such jets have no UE estimate and are dropped by the
JES and calibration steps.

The method and its numbers (2026-10, `UnderlyingEventPlayground/in_event/`):
azimuthal scans on FMS Run 15 data and FCS simulation put the end of the
jet's own energy at about 2R and the start of the recoil rise at about
120-130 deg; on the FMS, the ring on FCS-like geometry is within about 10% of
the FMS two-cone (+-90 deg) estimate (full FMS 4%), which is taken as the UE
systematic. Weighting per point and per pixel agree within 0.5% on Run 22.

## jet_scale

```
./submit_jet_scale_em_only.sh -i <jettrees_dir/> -o <outdir> [-n <json>] [-s <0|1>] [-l <label>] [-k <0|1>]
```

EM-only JES from simulated JetTrees (`sim_to_jet -E 1 -C 1`, with the
in-job UE). Reads a flat folder (for an organized sample, a folder of links).
Everything is done on UE-subtracted jets: the energy bins, the scale and its
fit, and the closure.

1. `JetEnergyScaleFineGrid_em_only`: on a 2.5 cm grid (5 cm match radius,
   10 cm from 37 GeV, where forward jets come from the heaviest-weight low-ptHat
   events and points are short of statistics) and 20 truth-energy bins (3-100 GeV), the cross-section-weighted median of
   scale = (reco_E - reco_ue_E) / (truth_E - truth_ue_E) and of the response.
   Medians close exactly under numerical inversion (arXiv:1609.05195).
   `noue` turns the UE subtraction off (for checks only).
2. `CalibrationMap_em_only`: describes the scale vs response as
   `A + B lnE + C ln^2E + D ln^3E`, with one shared shape (C, D from a single
   least-squares fit over every position, A and B free per position) and A, B
   per position (the level and slope); energy bins with n_eff < 3 are left
   out (a median of one or two effective jets is biased). The fits are linear and exact, so the
   step takes about a second. The A, B maps are smoothed (as the scale at
   3 and 30 GeV) and masked to the fiducial region. The evaluators clamp E to
   the range of all fitted responses [Emin, Emax] (the polynomial has no
   plateau). The JSON says `"model": "logcubic"`, with the global `C`, `D`,
   `Emin`, `Emax` and bins `{i, j, A, B}`; the evaluators
   refuse maps of earlier forms (including those with a `K`).
3. `jes_closure`: the map applied in memory exactly as `jet_calibration`
   does, on the same UE-subtracted weighted jets: weighted medians of
   E_corr / E_truth per truth-energy bin and radial band, and the map against
   the medians it was fitted to (`jes_closure.txt`, `.root`, `.png`).

   Known limitation: below about 9 GeV the response varies within an ECal
   tower (about 8% higher at tower boundaries than at centres, 3-8.6 GeV;
   within 0.5% above 8.6 GeV). The map's 2.5 cm grid with 5 cm pooling only
   calibrates the tower average, so single low-energy jets are off by up to
   +-5% and closure near the acceptance edge is off by about 0.6% there. Not
   corrected: A_N jets are above about 20 GeV (the FMS EM-jet analysis,
   PRD 103, 092009, used x_F 0.14-0.37 at 500 GeV, i.e. about 35-90 GeV).

`submit_jet_scale.sh` (`JetEnergyScaleFineGrid`, `CalibrationMap`) is the
ECAL+HCAL version, with a log-quadratic model and Gaussian-fit modes;
`-E 1` forwards to the EM-only script.

## jet_calibration

```
./submit_jet_calibration_em_only.sh -i <jettrees_dir> -c <json> -j <njobs> [-l <label>] [-k <0|1>]
```

Adds `reco_E_corr` and `reco_x_F_corr` to every file **in place** (asks for
confirmation). The EM-only map is evaluated at the UE-subtracted energy,
`reco_E_corr = (reco_E - reco_ue_E) / scale`. Jets without a UE estimate
(`reco_ue_E < 0`) or outside the map get `reco_E_corr = reco_x_F_corr = -1`
and must be dropped by the analysis. x is folded to
positive for the lookup, y is not. `submit_jet_calibration.sh` is the
ECAL+HCAL version (also writes `reco_E_uniform_corr`); `-E 1` forwards to the
EM-only script. The two JSON schemas are not interchangeable.

## pico_analysis

```
./submit_pico_analysis.sh -d <mip|hits> -f <listfile> -o <outdir> [-l <label>] [-k <0|1>]
```

`mip`: MIP peaks from tracks projected onto the FCS (`mip_ana.C`). `hits`:
flat FCS hit tree (`pico_to_root.C`). Single job, ACLiC-compiled at run time.
