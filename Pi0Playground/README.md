# Pi0Playground

Two-photon (π⁰, η) invariant mass from FCS ECal points and clusters, data
against simulation, to check the ECal energy scale that the jet energy scale
inherits from the simulation.

## Programs

- **`pi0_pairs.cpp`**: reads SimpleTrees, data or simulation, and fills mass
  histograms. Header comment lists every histogram.
  - Events: no FCS LED, an EM trigger fired (simulation: emulated flags).
  - Points: ECal, E > 1 GeV, inside the ECal rectangle less half a tower.
  - Pairs: same detector half, Z < 0.7 (FMS π⁰ cut), mass from the origin.
  - Binned in pair energy, point distance at the ECal, Z, and cluster class
    (both points from one cluster, or from two).
  - Variants: with the event vertex (data); with the shower tails that the
    0.4 GeV reclustering threshold dropped added back (hits 0.05–0.4 GeV
    within 8 cm); the ECal gain calibration's own selection
    (`StFcsPi0FinderForEcal`: highest-energy cluster pair, clusters > 1 GeV,
    Z < 0.7, ≤ 40 hits above 1 GeV).
  - Simulation: points matched to GEANT photons; pairs from one π⁰ are
    tagged, and their mass is also computed with true energies (position
    effect only) and true positions (energy effect only).
  - Weights σ / events per ptHat bin, absolute; a job must hold every file of
    each ptHat bin it reads.
  ```
  ./pi0_pairs <SimpleTree list> <out.root> <data|sim> [trig|all]
  ```
- **`pi0_compare.cpp`**: fits the peaks (Gaussian + quadratic; the
  calibration selection with A. Ogawa's `massfit.C` model) and prints the
  data/simulation tables; writes `<base>.pdf`.
  ```
  ./pi0_compare <data.root> <sim.root> <base>
  ```
- **`pi0_pairs.xml`, `submit_pi0_pairs.sh`**: condor, one job per list.
  ```
  ./submit_pi0_pairs.sh -w <workdir> "<list>:data" "<list>:sim" ...
  ```
- **`reference/`**: the ECal gain calibration code, copied 2026-10-28 from
  `/star/data01/pwg/liangxl/mywork/fcs_OO_Ecal/fcs_fastoffline_run25/`
  (X. Liang, github.com/xilinliang/FCS-ECal-pi0-reconstruction):
  `StFcsPi0FinderForEcal.{cxx,h}`, `runMudst_liang.C`, the READMEs and the
  per-tower fit macro; and A. Ogawa's `massfit_akio.C` (also in
  `FCSJetAnalysis/reference/massfit.C`). The calibration accepts a tower when
  its peak is within 10% of m_π0.

Build from this folder:
```
g++ -m64 -O2 pi0_pairs.cpp ../shared/VectorDict.cxx -o pi0_pairs -I../shared -I../jet_scale `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
g++ -m64 -O2 pi0_compare.cpp -o pi0_compare `root-config --cflags --libs` -std=c++11
```

## First result (2026-10-28, run 23007048 reclustered vs triggered simulation)

Outputs in `/star/data01/pwg/seanp/pi0_study/` (`run_v3`, `pi0_data_vs_sim_v3.{txt,pdf}`).
Simulation: `jes_nofilt_v1` SimpleTrees (10% of files) and `ue_trig_5to20_v1`.

- Points, Z < 0.7: data peak 113–120 MeV at every pair energy; simulation
  123–136 MeV, rising with energy. Data / simulation − 1 = −8 to −9% at
  10–30 GeV, −15% at 30–45 GeV.
- The same for pairs from two separate clusters (−9 to −11% at 5–15 GeV,
  −6% at 20–30 GeV) and in every Z bin: not a two-photon splitting effect.
- With the dropped shower tails added back: simulation 135–137 MeV (the π⁰
  mass), data 124–127 MeV; −7 to −9%. The tail fraction is the same in
  data and simulation (6.5% vs 5.9% of 5–10 GeV points): not the 0.4 GeV
  threshold.
- The calibration's own selection: data 104–107 MeV (−21 to −23% of m_π0),
  simulation 125–128 MeV; data / simulation − 1 = −16 to −17%.
- η (fits less stable): data 425–440 MeV against 526–535 MeV in simulation
  at 30–70 GeV, −18%.
- The data hits carry the current DB gain correction (tower 50: implied
  1.32 from E / ADC, DB 1.30 for 2022-01-07).
- Simulation truth: tagged π0 peak 10% low at 5–10 GeV, 1% at 30–45 GeV,
  from the point energy response (point E / photon E 0.88–0.96); positions
  put the photons 2–4% too far apart.

So the data ECal energies of this run sit 8–17% below the simulation's for the
same π⁰ selection. One run only; the gain correction varies through Run 22.

## Gain corrections: DB against the calibration's anchor (`gaincorr/`)

`gaincorr/dump_ecal_gaincorr.C` dumps the DB ECal gain corrections for a run
in both of X. Liang's formats; `gaincorr/compare_gaincorr.py` compares them
with an anchor file, mapping (ns, dep, ch) to tower ids with his
`depch_id.txt` (copied as `liang_depch_id.txt`; the DB's own mapping agrees
for all 1496 towers).

Run 23007048 (DB date 2022-01-07 13:35) against `fcsgaincorr_007_final.txt`
(received 2026-10-28): anchor / DB median 0.9983, 16–84% 0.9969–0.9995,
98% of towers within 1%, all within 5%. The DB holds the day-007 π⁰
calibration, with about 0.2% added (the time-dependent correction). 102
towers sit at the 1.2 starting value in the anchor (not refit).

## Calibration selection, data against simulation, triggered (2026-10-28)

`defclu/`: `runMudst_defclu.C` (data_to_jet/runMudst.C with StFcsClusterMaker's
default thresholds, seed 0.05 / towers 0.01 GeV, as the calibration and the
`fms_tune_calib_v2` simulation), `defclu.xml`, `submit_defclu.sh`. Data:
run 23007048 (120 files, 155k EM-triggered events), 23020046, 23036016,
23048051 (40 files each, 52–67k). Simulation: `fms_tune_calib_v2`, 150 files
per ptHat bin 0–30, unweighted (`pi0_pairs ... simuw`). `pi0_calib_compare.cpp`
fits A. Ogawa's model (mean 90–190 MeV, width 6–45 MeV) and a Gaussian +
quadratic cross-check. Output: `/star/data01/pwg/seanp/pi0_study/calib_v3.{txt,pdf}`.

Data / simulation, any EM trigger, all rings (Gaussian + quadratic; Ogawa's
model agrees within 1% where its width is sensible):
4–6 GeV 0.90, 6–8 0.94, 8–10 0.90, 10–14 0.89, 14–20 0.92, 20–30 0.95,
30–40 0.95. Rings: inner 0.88–0.97, middle 0.87–0.94, outer 0.88–0.95.
Runs of days 20, 36, 48 sit 0–5% above day 7. The data peaks are 1.5–3
times wider than the simulation's.

## EM trigger QA (2026-10-28): "empty" EM2 events, hot towers, emulation

Tools in `trigemu/` (`StFcsTrgReEmu` + `reemulate.C`: re-runs the FCS trigger
emulator on SimpleTree hits for a chosen run number, DB date, trigger version
and ADC scale; `emu_agreement.cpp`, `empty_em2_events.cpp`,
`empty_em2_radius.cpp`) and `hot_tower_scan.cpp`. Outputs in
`/star/data01/pwg/seanp/pi0_study/trigger_qa/`.

**"Empty" EM2 events.** Real fcsEM2 events whose highest offline ECal cluster
has E_T < 1.5 GeV, far below the 5 GeV threshold:

| Run (day) | Real EM2 events | Empty share | HCal E, empty / normal (GeV) |
|---|---|---|---|
| 23007048 (7) | 53,558 | 31.0% | 24 / 26 |
| 23020046 (20) | 22,729 | 32.0% | 25 / 27 |
| 23036016 (36) | 15,759 | 53.0% | 53 / 57 |
| 23048051 (48) | 20,393 | 51.7% | 57 / 61 |

They carry ~42 GeV of ECal energy (largest hit ~5 GeV) against ~78 GeV (17
GeV) in normal EM2 events. Their largest hit sits within 50 cm of the beam
(about half) or beyond 110 cm (a third), where a 5 GeV deposit has an offline
E_T of ~1 GeV, and concentrates on a few DEP boards: (ns 0, dep 4) ~30%, then
(ns 0/1, dep 22) and (ns 0/1, dep 5) ~8-12% each, against an even spread over
48 boards for normal events. The online E_T weighting of those edge boards
appears not to follow sin(theta): an online trigger configuration issue,
worth raising with A. Ogawa / D. Kapukchyan. Raw EM2 rates per run track it,
so run QA should use an offline-matched trigger rate instead.

**Hot towers** (largest hit in > 2% of triggered events): only det 0 id 540
(x -87.5, y -47.1 cm) in run 23048051, 14.1% of events, mean 6.7 GeV. Busiest
tower elsewhere: 0.5-0.7%.

**Prescales** (RunLog.l0TriggerSet): fcsEM2 3.32 / 2.82 / 1 / 1.36 for runs
23007048 / 23020046 / 23036016 / 23048051; fcsEM1 17.9 / 14.3 / 2.6 / 5.8;
fcsEM0 221 / 179 / 39 / 72. The recorded bits are post-prescale.

**Emulator against hardware** (run 23036016, fcsEM2 unprescaled, data hits
with their ADC sums, run settings, trgSelect 202209): P(real | emulated) =
0.80-1.0 for fiducial top clusters above 3 GeV, but P(emulated | real) is
only 0.50-0.66 there (0.04 below 3 GeV). Scaling the ADC by 1.22 or 1.5
barely raises it (0.75, 0.78 above 7 GeV) while precision drops: the
mismatch is not one scale factor. The simulation's own emulation (run 1,
DB 2023-04-10) is reproduced exactly by `reemulate.C` (0 of 6059 events
differ).

## Energy scale or overlap? (2026-10-28, `pi0_overlap_study.cpp`, outputs `pi0_study/v7_overlap.txt`)

Calibration selection, default clustering, common 45 MeV tower floor; data all
events (day 7, 209k) and EM-triggered (days 7 + 20) against simulation
(`fms_tune_calib_v2`, unweighted).

- eta (photons ~4x further apart, little overlap): data/sim 0.925 +- 0.025,
  0.916 +- 0.018, 0.920 +- 0.012 for lower photon energy 3-5, 5-8, 8-20 GeV
  (all events); 0.92 / 0.87 with the EM bit. The deficit is there without
  overlap: an energy scale effect.
- pi0 against single-photon energy, all events: 0.91-0.93 (1-5 GeV), 0.94
  (5-6), 0.96-0.97 (6-20); EM bit 0.90-0.93 throughout.
- Cluster width (energy-weighted RMS radius): data 1-2% wider than simulation
  below 4 GeV, 3-5% above. Too small to make a 5-10% mass shift by overlap.
- Separation bins are not diagnostic at fixed pair energy (the mass is then
  nearly fixed by kinematics).
- Inclusive over all pair energies: data 136.7 MeV, simulation 142.4 MeV
  (0.960), matching the 139 MeV of the calibration-chain macro: the inclusive
  ratio mixes pair energies whose peaks shift, and does not measure the scale.
- The isolated-tower (MIP) spectrum shows no clean MIP peak in either sample:
  inconclusive.

Data photons are ~7% low relative to the simulation (0.92-0.93 from the eta
and the pi0 below 5 GeV), rising towards ~0.96 for single photons above 6 GeV
in all-event data.
