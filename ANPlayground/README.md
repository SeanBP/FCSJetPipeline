# ANPlayground

Sizing study for a full Run 22 A_N(x_F) pass over the data. It measures wall time,
memory, scratch and output per job and extrapolates them to the full sample.
Unfolding is not included.

Each job runs the whole chain on up to 40 MuDsts of one run. It keeps only
the A_N yield histograms, so no SimpleTrees or jet trees are stored.

| Step | Program | Note |
|---|---|---|
| 1. MuDst → SimpleTree | `data_to_jet/runMudst.C` | events that fired fcsEM0, fcsEM1 or fcsEM2 only |
| 2. Jet finding | `data_to_jet/RecoJets` (`1 1`) | EM-only jets from FCS points, with areas |
| 3. UE subtraction | `underlying_event/UnderlyingEvent` | proxy events from the job's own files (one run) |
| 4. JES calibration | `jet_calibration/ApplyCorrections_em_only` | calibrates the UE-subtracted energy: `reco_E_corr = (reco_E - reco_ue_E) / scale` |
| 5. A_N histograms | `AnHistograms` | yields vs x_F (`reco_x_F_corr`) |

SimpleTrees are deleted after step 3, and jet files after step 5.

## Files

- **`AnConfig.h`**: selection and binning.
  - EM levels and thresholds. `kEmThresholdGeV` is a PLACEHOLDER: about
    1.5 / 2.6 / 3.7 GeV, from where EM1 and EM2 turn on in run 23007048,
    scaled by the stage-2 thresholds 63 / 111 / 158.
  - The 1.1 × threshold factor.
  - Point-multiplicity classes, spin states and x_F bins.
- **`AnHistograms.cpp`**: fills, per jet, `h_xF_<EM level>_<np class>_<side>_<spin>`.
  - EM level: the highest of fcsEM0/1/2 fired. The `_tpc` copies and
    fcsEM3 are not used.
  - np class: `reco_nphotons`, in classes 0–4 and ≥ 5.
  - Side: the sign of x.
  - Spin: `Spin_config` 5/6/9/10; events with other values are dropped.
  - Cuts: calibrated (`reco_E_corr > 0`), `pass_physical_jet_cut`, and
    corrected pT = `reco_E_corr / cosh(eta)` > 1.1 × that event's threshold.
  - Also written: a cut flow, events per spin state per EM level, corrected
    pT before the threshold cut, and the run numbers.
  - The classes are kept apart so they can be combined when A_N is formed.
    Files of one run or fill add with `hadd`.
- **`an_job.xml`, `submit_an_job.sh`**: the condor job. Each stage runs
  under `/usr/bin/time`. Every job writes
  `timing/timing_<run>_<job>.txt`, one row per stage call: wall, user and
  sys time, peak RSS, and scratch use afterwards. Histograms go to `hists/`.
  ```
  ./submit_an_job.sh -f <list> -o <outdir> [-n 40] [-c <JES json>] [-k 1] [-N]
  ```
  The default JES map is
  `jes_points/jesmap_v3/JetEnergyScale_lookup_points_ue_v3_lognormal.json`,
  of the earlier log-normal form, which `ApplyCorrections_em_only` no longer
  reads: pass a log-cubic map with `-c` until a UE-subtracted one replaces it.
- **`catalog_query.sh`**:
  - `counts <out>` writes per-run events, files and bytes. The P24ia st_fwd
    copy of 2026-10-20 is `catalog_runs_local.txt`.
  - `list <out> <runs...> [-l N]` writes a MuDst list.
- **`estimate_resources.py <outdir>/timing`**: fits the per-file costs from
  the benchmark (MuDst step = a + b·events, fetch per GB, jet finding per
  file) and the per-job costs (UE scaled as files², calibration and
  histograms linear in files). It applies them to every run in the catalog
  and prints:
  - number of jobs and the job wall-time distribution
  - total core-hours and campaign length for a given number of slots
  - memory and scratch to request, and output volume
- **`benchmark_12runs.list`**: the first 40 files of 12 runs spread over
  23005044–23108024. Gives 12 jobs and 480 files.

## Dataset (catalog, 2026-10-20)

P24ia st_fwd of production_pp500_2022, on disk:

| Selection | Runs | Files | Events | Size |
|---|---|---|---|---|
| All | 3368 | 954,613 | 8.81e9 | 1199 TB |
| Analysis: 23005044–23108024, MakeRunList.pl exclusions | 2993 | 881,486 | 7.91e9 | 1079 TB |

- Files per run: median 337, 90% below 470, maximum 902.
- Average file: about 9,000 events, 1.22 GB.
- The RunLog bad-run filter is not applied, so these totals are an upper bound.
- The spin DB covers runs up to 2022-04-18 (run 23108), with between-fill
  gaps of about 2%.

## Earlier timing reference

`data_to_jet` on run 23007048 (request 4107497D…, 2026-10-13): 3 jobs of 40
files, EM filter, point jets, in-job UE.

- Wall time 38–61 min per job; CPU 34–57 min.
- Peak memory 238–253 MB.
- MuDst → SimpleTree is about 85% of the wall time (0.8–1.2 min per
  file), and xrdcp most of the rest. Jet finding and UE together take
  about 1 min per job.
