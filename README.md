# FCSJetPipeline

SUMS pipelines for STAR FCS jets (Run 22 pp, sqrt(s) = 508 GeV), one folder
per stage. Each folder holds its XML template, a `submit_*.sh` wrapper (run it
with `-h` for its options) and its sources.

| Stage | Input -> output |
|---|---|
| `sim_to_jet/` | PYTHIA + starsim -> BFC/MuDst -> SimpleTree -> matched truth/reco jets (+ in-job underlying event) |
| `data_to_jet/` | official MuDst -> SimpleTree -> reco jets (+ in-job underlying event) |
| `underlying_event/` | underlying-event (UE) estimate per jet; runs inside the two jobs above, or over a whole simulation sample |
| `jet_scale/` | simulated JetTrees -> jet-energy-scale (JES) lookup JSON |
| `jet_calibration/` | applies a JES JSON to JetTrees, in place |
| `pico_analysis/` | PicoDst diagnostics, not part of the jet chain |

Shared resources at the root:

- `shared/JetParameters.h`: jet radius, fiducial rectangles, thresholds, the
  EM-only truth-particle selection, sqrt(s) per run, x_F.
- `shared/VectorDict.*`: the `std::vector` dictionary every program reading
  jetTree vector branches links against.
- `StRoot/`: `StSimpleReaderMaker`, `StSimpleTreeMaker` and a local
  `FcsJetFilter` (settable threshold). `.sl73_x8664_gcc485/` is their
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

## sim_to_jet

```
./submit_sim_to_jet.sh -g <pythia8|pythia6> -t <tune> -p <ptHatMin> -n <nevents> -s <0|1> -o <outdir> \
    [-P <ptHatMax>] [-j <njobs>] [-w <json>] [-l <label>] [-k <0|1>] [-f <triggers>] [-c <run>] \
    [-e <ethr>] [-E <0|1>] [-F <0|1>] [-T <tune_pp>] [-M <0|1>] [-C <0|1>]
```

One job generates `-n` events and writes
`jet_output_<gen>_<tune>_pt<ptHatMin>_<job>.root` (and the SimpleTree with
`-s 1`). Main options:

| Flag | Meaning |
|---|---|
| `-w` | JSON of relative ptHat-range shares across the `-j` jobs (`pthat_distribution_optimized.json`: disjoint 5 GeV bins, exact weights); replaces `-p`/`-P` |
| `-f` | keep only events that fired one of these FCS triggers (`FcsTriggerDefs.h` names, emulated with `fcs_trg_base`) |
| `-c` | run number fixing the FCS gain calibration era (default: BFC's nominal date) |
| `-e` | `FcsJetFilter` forward-flux threshold, GeV (default 50, 0 = off) |
| `-E 1` | EM-only reco jets (ECAL only) with the ECAL fiducial rectangle |
| `-F 1` | reco jets from FCS points (`StMuFcsPoint`) instead of tower hits |
| `-M 1` | the FMS EM-jet generator settings: pythia6, PyTune 370, PARP(90)=0.213, PARP(91)=1, MSEL=1 |
| `-C 1` | JES calibration sample: no jet pT minimums and a minimum-bias UE step; no `-f` |

Folder layout: after a submission finishes,
`organize_sim_output.sh <outdir>` moves the files into
`jet_output/pt<N>/` and `SimpleTree/pt<N>/`, so no folder holds thousands of
files.

**Jets.** Anti-kt, R = 0.7. Reco jets must be inside the fiducial rectangle
minus a margin (10 cm, or 20 cm for EM-only); truth jets may extend outside
it by R/2 (R for EM-only) in eta-phi. Matching is one-to-one within
dR < R/2. EM-only truth jets are clustered from photons created within
100 um of the primary vertex (prompt and pi0/eta decay photons, not GEANT
secondaries). EM-only jet pT minimums: reco 2 GeV, truth 1 GeV (as the FMS
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

The list is ordered by run and time and cut into jobs of at most `-n`
(default 40) consecutive files of one run, so the in-job UE step finds proxy
events from the same run. `-j` instead splits the list evenly. Each
SimpleTree gives one `jet_output_<run>_<file>.root` with the same suffix.

Per event the trees carry the fired FCS trigger flags (resolved from
trigger id and run number, `FcsTriggerDefs.h`, generated from
`TriggerIDs.txt` by `generate_trigger_defs.py`; three copies, regenerate
together and rebuild StRoot), `Spin_config` (spin DB, covers 2021-12-14 to
2022-04-18 only), the event vertex (TPC, VPD, then BBC) and sqrt(s).
`MakeRunList.pl` needs the Perl library in
`/star/u/dkap7827/Tools2/Tools/PerlScripts`.

## underlying_event

The UE under a jet is the mean ECal point energy in an R cone at the jet's
position, averaged over other events (proxy events) of the same group (run in
data, ptHat bin in simulation) and EM trigger level whose own jet has a
similar pT (within 5%) and lies near 90 deg in azimuth from the cone, at least
2.5R away (`ProxySelection.h`, `ProxySearch.h`). Every jet_output file gets
`reco_ue_E`, `reco_ue_pt`, `reco_ue_nproxy`, `reco_ue_pt_window`, and in
simulation the truth-level equivalents `truth_ue_*` (from truth jets and
EM-only truth photons).

The in-job step only sees one job's events. For a simulation sample, rerun it
over all files of each ptHat bin, one condor job per bin, inputs copied to
scratch:

```
./submit_ue_by_pthat.sh -i <sample>[,<sample>...] -w <workdir> [-m minbias|default] [-n]
```

## jet_scale

```
./submit_jet_scale_em_only.sh -i <jettrees_dir/> -o <outdir> [-n <json>] [-s <0|1>] [-l <label>] [-k <0|1>]
```

EM-only JES from simulated JetTrees (`sim_to_jet -E 1`; for the calibration
use a `-C 1` sample with its UE rerun in minimum-bias mode). Reads a flat
folder (for an organized sample, a folder of links).

1. `JetEnergyScaleFineGrid_em_only`: on a 2.5 cm grid (5 cm match radius) and
   20 truth-energy bins (3-100 GeV), the cross-section-weighted median of
   scale = (reco_E - reco_ue_E) / (truth_E - truth_ue_E) and of the response.
   Medians close exactly under numerical inversion (arXiv:1609.05195).
   `noue` turns the UE subtraction off.
2. `CalibrationMap_em_only`: per position, fits
   `A + B/E + C/E^2 + D exp(-(E-E0)^2 / 2 sigma^2)` (sigmaL below E0, sigmaR
   above), smooths the maps and masks them to the fiducial region. For speed
   it can run as condor shards (`--shard` / `--merge`, see
   `JetScalePlayground/submit_calibration_map_diag_condor.sh`).

`submit_jet_scale.sh` (`JetEnergyScaleFineGrid`, `CalibrationMap`) is the
ECAL+HCAL version, with a log-quadratic model and Gaussian-fit modes;
`-E 1` forwards to the EM-only script.

## jet_calibration

```
./submit_jet_calibration_em_only.sh -i <jettrees_dir> -c <json> -j <njobs> [-l <label>] [-k <0|1>]
```

Adds `reco_E_corr` and `reco_x_F_corr` to every file **in place** (asks for
confirmation). The EM-only map is evaluated at the UE-subtracted energy,
`reco_E_corr = (reco_E - reco_ue_E) / scale`. Jets without a UE estimate keep
`reco_E`; jets outside the map keep the subtracted energy. x is folded to
positive for the lookup, y is not. `submit_jet_calibration.sh` is the
ECAL+HCAL version (also writes `reco_E_uniform_corr`); `-E 1` forwards to the
EM-only script. The two JSON schemas are not interchangeable.

## pico_analysis

```
./submit_pico_analysis.sh -d <mip|hits> -f <listfile> -o <outdir> [-l <label>] [-k <0|1>]
```

`mip`: MIP peaks from tracks projected onto the FCS (`mip_ana.C`). `hits`:
flat FCS hit tree (`pico_to_root.C`). Single job, ACLiC-compiled at run time.
