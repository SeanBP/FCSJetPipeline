# asymmetry

The Run 22 EM-jet A_N(x_F) analysis: the condor job that turns MuDsts into
per-run yield histograms, the A_N step (bad runs, per-run polarization,
statistical errors, final plot), the yield and sanity plots, the projection
of the errors to the full sample, and the sizing study that extrapolates the
job's wall time, memory, scratch and output to the full sample. Unfolding is
not included.

## Smoke test (a sample spread over Run 22)

```
./make_smoke_list.py <outdir>/<prefix> --runs 60 --files 40
./submit_an_job.sh -f <outdir>/<prefix>.list -o <outdir> -l <prefix> -k 1 -J 1
./smoke_report.sh <outdir> <outdir>/<prefix>_good_runs.txt     # once the jobs end
```

- `make_smoke_list.py`: the good runs (catalog analysis runs that pass the
  RunLog shift-crew flag, are not flagged in the pi0 bad-run list and have
  polarization) and, from them, runs evenly spaced over Run 22 with files
  evenly spaced through each run.
- `smoke_report.sh` writes `<outdir>/report/`: A_N (`an.root`, `an.txt`,
  `an_vs_xf*.png`), the errors projected to all good runs
  (`an_projected.txt`, `an_vs_xf_projected*.png`), yield and sanity plots
  and per-run trends (`plots/`), and wall time, memory and scratch per
  stage, per job (also condor's record) and extrapolated to the full pass
  (`resources.txt`, `condor_jobs.txt`).
- `smoke_v1` (2026-11-03): 60 runs 23008012-23107043, 2357 files, 2.08e7
  events, 0.27% of the 7.85e9 good-run events, in
  `/star/data01/pwg/seanp/asymmetry/smoke_v1`.

Each job runs the whole chain on up to 40 MuDsts of one run. It keeps the
A_N yield and sanity histograms; no SimpleTrees are stored, and the
calibrated jet files (about 0.3 MB per MuDst) only for every N-th job with
`-J N` (`-J 1` all, default none), in `<outdir>/jets/`.

| Step | Program | Note |
|---|---|---|
| 1. MuDst → SimpleTree | `data_to_jet/runMudst.C` | events listing fcsEM0, fcsEM1, fcsEM2 or a `_tpc` copy only |
| 2. Jet finding | `data_to_jet/RecoJets` (`1 1`) | EM-only jets from FCS points, with areas |
| 3. UE subtraction | `underlying_event/UnderlyingEvent` | in-event eta-weighted side wedges: each jet's own event, its eta band from 2R out to 150° in azimuth |
| 4. JES calibration | `jet_calibration/ApplyCorrections_em_only` | calibrates the UE-subtracted energy: `reco_E_corr = (reco_E - reco_ue_E) / scale` |
| 5. A_N histograms | `AnHistograms` | yields vs x_F (`reco_x_F_corr`) and sanity distributions |
| 6. A_N (offline) | `AnAsymmetry -r`, `plot_an.C`, `AnPlots`, `project_errors.py` | per-run yield files: flagged bad runs dropped (D. Kapukchyan's list), per-run polarization, A_N with statistical errors; yield and sanity plots; projected errors |

SimpleTrees are deleted after step 3, and jet files after step 5 unless
the job keeps them (`-J`). Step 6 runs
on the job outputs, not in the job.

## Files

- **`AnConfig.h`**: selection and binning.
  - EM levels and thresholds. `kEmThresholdGeV` = 2.0 / 3.5 / 5.0 GeV, the
    nominal online values (64 / 112 / 160 counts at 0.0316 GeV), compared
    with uncorrected jet pT; taken as given until the turn-ons are measured.
  - The 1.15 × threshold factor (as the FMS EM-jet analysis) and the 2 GeV corrected pT minimum.
  - Point-multiplicity classes, spin states and x_F bins.
- **`AnHistograms.cpp`**: fills, per jet, `h_xF_<EM level>_<np class>_<side>_<spin>`.
  - EM level: the lowest listed EM trigger (fcsEM0/1/2, `_tpc` copies
    included) whose threshold the jet passes. fcsEM3 is not used.
  - np class: `reco_nphotons`, in classes 0–4 and ≥ 5.
  - Side: the sign of x.
  - Spin: `Spin_config` 5/6/9/10; events with other values are dropped.
  - Cuts: calibrated (`reco_E_corr > 0`), `pass_physical_jet_cut`, and
    corrected pT = `reco_E_corr / cosh(eta)` > 2 GeV and, as FMS, uncorrected
    pT (`reco_pt`) > 1.15 × the threshold of any EM trigger listed in the event
    (`fcsEMk` or `fcsEMk_tpc`); events with `fcs_led` set are dropped.
  - Listed triggers are those that recorded the event (after prescale); the
    jet's EM level is the lowest listed level whose threshold it passes.
  - `h_abscos_*`, `h_xfsum_*`: x_F weighted by |cos phi| and by x_F, for each
    bin's mean |cos phi| and mean corrected x_F.
  - Also written: a cut flow, events per spin state per EM level, corrected
    pT before the threshold cut, and the run numbers.
  - Sanity distributions `qa_<variable>_evt` (every jet of the events
    passing the event cuts) and `qa_<variable>_sel` (the jets in the
    yields): eta, phi, uncorrected, UE and corrected energy, pT, x_F, points,
    jets per event, and the x-y map (`qa_xy_*`, also per side); per event the
    vertex z and source, the listed EM triggers and `Spin_config`.
  - The classes are kept apart so they can be combined when A_N is formed.
    Files of one run or fill add with `hadd`.
- **`AnAsymmetry.cpp`**: A_N from the yields (not part of the job):
  ```
  ./AnAsymmetry <yields.root> <out.root> <P_blue> <P_yellow> [rebin = 2]
  ./AnAsymmetry -r <yield file list> run22_polarization.txt <out.root> [rebin = 2] [list | none] [flagged | all]
  ```
  - `-r`: one yield file per run, summed; each bin uses the mean polarization
    of its jets (runs weighted by their jets, as FMS weighted by events).
  - `-r` also drops bad runs: by default the runs with a flag in D.
    Kapukchyan's pi0 QA list (`/star/u/dkap7827/Fcs2019/FcsAsymRun22/BadRuns.txt`,
    read at each call). As of 2026-02-07 the list has a line for every
    analysis run up to 23024003 (413) and 59 have a flag (LowStat, EnShift,
    DoublePeak, ...), so an empty flag is read as a run that passed; mode
    `all` would drop the whole reviewed period. Runs after 23024003 are not
    reviewed yet. The output's TNamed `bad_runs` records the list's path,
    date and the runs removed.
  - Points at the bin's mean corrected x_F; rebin 2 = x_F bins of 0.05.
  - Statistical errors: Poisson from the raw counts in each cross ratio,
    combined over the other beam's spin states, divided by <P><|cos phi|>;
    weighted yields are refused. EM levels add (a jet is in one level), so
    level "all" is the final result.
  - Also writes `<out>.txt`: every point with its statistical error, jets,
    <|cos phi|> and <P>.
- **`plot_an.C`**: the final A_N vs x_F (blue at x_F > 0, yellow at x_F < 0,
  statistical errors only), `.png` and `.pdf`; optional error scale and
  points at zero for a projection:
  ```
  root -b -q 'plot_an.C("an.root","an_vs_xf.png","all","all"[,err_scale[,at_zero]])'
  ```
  - Left-right cross ratio per beam, the other beam's spin held fixed (one
    cross ratio per state, averaged); A_N = eps / (P <|cos phi|>).
  - Blue: left = x > 0, at x_F > 0. Yellow: left = x < 0, at x_F < 0 (STAR
    +x is south; FCS north points have x < 0).
- **`make_run_polarization.py`, `run22_polarization.txt`**: polarization per
  run, the fill's P0 + dP/dt at the middle of the run (polarimetry fill table,
  copy in /star/u/dkap7827/Fcs2019/FcsAsymRun22/Run22PolFromWeb_20250205.txt;
  run times and fill from RunLog). 4662 runs; 2983 of the 2996 catalog
  analysis runs (13 left out: fill not in the table, or bad RunLog times).
  ```
  ./make_run_polarization.py <fill table> run22_polarization.txt
  ```
- **`AnPlots.cpp`**: `./AnPlots <yield file list> <out dir>`: cut flow,
  jet distributions (both stages), x-y maps, event quantities, pT by EM
  level with the threshold cuts, yields vs x_F per side and spin state (per
  EM level), per photon class and per level, and per-run trends (events,
  jets per event, x > 0 / x < 0, spin-state fractions, EM level mix, mean
  energy; `runs.txt`). PNGs and one `plots.pdf`.
- **`project_errors.py`**: `./project_errors.py <outdir> <an.txt> <good_runs.txt>`:
  stat_full = stat x sqrt(events processed / good-run events), from the
  catalog events of the files the jobs processed; writes
  `<an>_projected.txt` and prints the factor for `plot_an.C`'s `err_scale`.
- **`an_toy_closure.cpp`**: toy jets with injected blue and yellow A_N,
  defined as left = spin x beam direction, unequal spin-state luminosities
  and a side-dependent efficiency; AnAsymmetry recovers both signs and sizes
  (0.05 bins: chi2 37.4 / 40 on 2e7 jets; two runs with different P through
  `-r`: 33.7 / 40).
- **`an_job.xml`, `submit_an_job.sh`**: the condor job. Each stage runs
  under `/usr/bin/time`. An xrootd fetch that fails or takes over 5 min
  is cut off and tried once more (`fetchretry` row); a second failure is
  recorded as `fetchfail` and the file skipped. Normal fetches take seconds;
  stalls of over an hour were seen on 2026-11-03. Every job writes
  `timing/timing_<run>_<job>.txt`, one row per stage call: wall, user and
  sys time, peak RSS, and scratch use afterwards. Histograms go to `hists/`.
  ```
  ./submit_an_job.sh -f <list> -o <outdir> [-n 40] [-c <JES json>] [-k 1] [-J <N>] [-N]
  ```
  The default JES map is
  `jet_calibration/JES_Lookups/JetEnergyScale_lookup_em_only.json`.
- **`catalog_query.sh`**:
  - `counts <out>` writes per-run events, files and bytes. The P24ia st_fwd
    copy of 2026-10-20 is `catalog_runs_local.txt`.
  - `list <out> <runs...> [-l N]` writes a MuDst list.
- **`estimate_resources.py <outdir>/timing [--runs good_runs.txt]`**: fits the per-file costs from
  the benchmark (MuDst step = a + b·events, fetch per GB, jet finding per
  file) and the per-job costs (UE, calibration and histograms linear in
  files; compiles fixed). It applies them to every run in the catalog
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
- These catalog totals apply neither the RunLog shift-crew flags nor the pi0
  bad-run list (both remove runs later: MakeRunList.pl and `AnAsymmetry -r`),
  so they are an upper bound.
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
