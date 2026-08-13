# FCSJetPipeline

Four distinct SUMS pipelines, each in its own self-contained folder:
simulation-to-jet, real-production-data-to-jet, jet-energy-scale
calibration building, and calibration application. Plus one standalone
utility folder, `pico_analysis/`, for PicoDst-based diagnostics -- not
part of the jet-finding chain (see its own section below). Each folder
holds its own XML template, its own `submit_*.sh` wrapper (with a working
`-h`), and every source file it needs -- nothing pipeline-specific is
shared between folders. Only genuinely cross-pipeline resources live at
the pipeline root, all together under `shared/`: `JetParameters.h`
(physics constants/cuts used by every jet-finding step except
`CalibrationMap`) and the compiled `vector<float>`/`vector<int>`
CollectionProxy dictionary (`VectorDict.cxx/.h`) that any step reading
jetTree's `std::vector` branches must link against -- `JetMatcher`,
`JetEnergyScaleFineGrid`, `RecoJets`, `ApplyCorrections`, and the
standalone `FCS2VIRTUE` tool. Also at the pipeline root: `StRoot/` (the
two locally compiled makers, `StSimpleReaderMaker`/`StSimpleTreeMaker`),
and `.sl73_x8664_gcc485/` (their precompiled binary cache).

Every pipeline shares two more conventions:

- **Logs are opt-in.** By default a job's stdout/stderr are discarded
  (sent to `/dev/null`); pass `-k 1` to any wrapper to save them instead,
  under an auto-created `log/` folder inside that pipeline's output
  directory (or `<jettrees_dir>/log/` for `jet_calibration`, which has no
  separate output directory since it corrects files in place).
- **No stray list files.** Any working list a job builds internally
  (line-range chunks, per-job file lists) is `rm -f`'d at the end of the
  job script; nothing but real output ever leaves scratch.

## Layout

```
FCSJetPipeline/
├── shared/                 pipeline-root resources used by multiple stages: JetParameters.h (physics constants/cuts) and the VectorDict.cxx/.h CollectionProxy dictionary (for jetTree's std::vector branches)
├── StRoot/                 source for the two locally-compiled makers, shared by sim_to_jet and data_to_jet
├── .sl73_x8664_gcc485/      precompiled binaries for that StRoot/ source (rebuild if StRoot/ changes)
├── sim_to_jet/               starsim -> BFC/MuDst -> jet finding (simulation)
│   ├── sim_to_jet.xml, submit_sim_to_jet.sh
│   └── pthat_distribution_optimized.json      relative ptHatMin proportions for -w
├── data_to_jet/              MuDst -> SimpleTree -> jet finding (real production data)
│   └── data_to_jet.xml, submit_data_to_jet.sh, MakeRunList.pl
├── jet_scale/                 builds a JES calibration JSON from a folder of JetTrees
│   └── jet_scale.xml, submit_jet_scale.sh
├── jet_calibration/             applies a calibration JSON to JetTrees, in place
│   └── jet_calibration.xml, submit_jet_calibration.sh
├── pico_analysis/                 PicoDst diagnostics utility, not part of the jet-finding chain
│   └── pico_analysis.xml, submit_pico_analysis.sh, mip_ana.C, pico_to_root.C
└── output/                       default scratch space; jobs are usually pointed at an explicit -o instead
```

## Job-splitting, standardized

- **`sim_to_jet`**: `-n` events generated per job, `-j` total number of
  parallel jobs. Each job is fully independent (its own generated events),
  so there's no fixed input to split -- one job runs the whole chain on
  one file.
- **`data_to_jet`** and **`jet_calibration`**: the input (a file catalog
  list, or a folder of JetTrees) is a fixed set of items. `-j` evenly
  splits that fixed set across jobs -- contiguous, non-overlapping ranges,
  so no two concurrent jobs ever touch the same input file. A job that
  gets several items still processes them one at a time (see "1:1
  SimpleTree/JetTree linkage" below), it never combines them.
- **`jet_scale`**: the one deliberate exception. It reads an entire folder
  as a single whole-population fit (memory/CPU-bound, not parallelizable
  across files the way the others are), so it is always exactly one job;
  there is no `-j`.

## sim_to_jet

starsim -> BFC/MuDst reconstruction -> jet finding, all in one job, on a
single generated SimpleTree file.

```
./submit_sim_to_jet.sh -g <pythia8|pythia6> -t <tune_param> -p <ptcut> \
    -n <nevents> -s <0|1> -o <outdir> [-j <nprocesses>] [-l <label>] \
    [-w <json>] [-k <0|1>] [-f <triggers>]
```

| Flag | Meaning |
|---|---|
| `-g` | generator: `pythia8` or `pythia6` |
| `-t` | pythia8: `PDF:pSet` integer (8=CTEQ6L1, 5=MSTW2008LO, 3=MRST LO*, 21=NNPDF3.1sx, <=0=default). pythia6: `PyTune` integer (325=Perugia STAR, 0=default) |
| `-p` | pTHatMin cut, GeV -- a single flat value used by every job; ignored if `-w` is given |
| `-n` | events to generate per job |
| `-s` | `1` to keep the intermediate `SimpleTree_*.root`, `0` to discard it after jet finding runs |
| `-o` | output directory (created if missing) |
| `-j` | total number of parallel jobs, default 1 |
| `-w` | optional: absolute path to a JSON file of relative ptHatMin proportions to sample across the `-j` jobs, instead of one flat `-p` value -- overrides `-p` |
| `-l` | short label used in the generated XML's filename, default `run` |
| `-k` | `1` to save stdout/stderr under `<outdir>/log/`, `0` to discard them, default `0` |
| `-f` | comma-separated FCS trigger flag names (see FCS trigger flags below); an event only gets a jetTree entry if at least one of them fired. Omit/empty (default) keeps every event. Called `-f` here, not `-t` like `data_to_jet`'s equivalent flag, since `-t` is already taken by the tune/PDF param above |

Example -- 500 pythia8/CTEQ6L1 events at ptHatMin=10 GeV, discard the SimpleTree:

```
./submit_sim_to_jet.sh -g pythia8 -t 8 -p 10 -n 500 -s 0 \
    -o /star/data01/pwg/seanp/pipeline_output/cteq6l1_pt10_test
```

Each run produces `jet_output_<generator>_<tune>_pt<ptcut>_<run_number>.root`
(and `SimpleTree_..._<run_number>.root` if `-s 1`) in the output directory.
`run_number` is `JOBINDEX`, so raising `-j` above 1 runs that many
independent jobs in parallel without filename collisions.

**Weighted pT-cut sampling (`-w`).** Instead of one flat `-p` cut, `-w
<json>` draws each job's ptHatMin from a table of relative proportions,
so low-pT (bulk of the cross section) and high-pT (rare tail) both get
proportionate statistics instead of either flooding or starving the
sample. `pthat_distribution_optimized.json` ships with the table
currently in use (fit to `pythia8_nnpdf23lo`'s reference distribution);
format:

```json
{
  "reference_total": 20000,
  "bins": [
    { "ptcut": 0, "count": 21 },
    { "ptcut": 2, "count": 218 },
    ...
  ]
}
```

`count` is that bin's share out of `reference_total`; `submit_sim_to_jet.sh`
scales the proportion to however many jobs `-j` actually submits (does
not require submitting `reference_total` jobs at once). Since this
substitution is multi-line csh, not a single value, the wrapper generates
it with Python rather than `sed`.

**FCS trigger flags.** Both the SimpleTree (`data`, written by
`StSimpleReaderMaker`) and the JetTree (`jetTree`, written by
`JetMatcher`) carry one `int` branch per known FCS trigger name from
`TriggerIDs.txt` -- 1 if that trigger fired on the event, 0 otherwise,
same convention as `data_to_jet`. Simulated events have no real online
trigger decision to read, so `readMudst.C` constructs an
`StFcsTriggerSimMaker` (the FCS trigger-logic emulator, `fcs_trg_base`)
and passes it to `StSimpleReaderMaker::SetTriggerSim()`, which drives it
directly on each event's FCS hits instead of doing the real-data id
lookup. `StFcsTriggerSimMaker::Make()` itself is never called (it
segfaults reading a muDst with no `StEvent` present) -- it's kept in the
chain only for `Init()`/`InitRun()` (gains, pedestals, `stage_version`)
via `SetActive(kFALSE)`, and `StSimpleReaderMaker` drives its
`fcs_trg_base` engine directly. Only ~20 of the 64 flags are computable
this way (the physics-condition categories `fcs_trg_base` models); the
rest stay 0 for simulated events. `trgSelect=202209` (default, set in
`readMudst.C`) should match whatever trigger-algorithm version was live
for the run period being modeled; see `TriggerPlayground` if this needs
revisiting.

**ECAL/HCAL gain calibration.** `runSimBfc.C` deliberately does *not*
call `StFcsDb::forceUniformGain()`/`forceUniformGainCorrection()` --
`StFcsDb` is left on its default `GAINMODE::DB`, which pulls real
calibrated gain and per-tower gain-correction values from the STAR
conditions DB. ECAL's `gainCorr` is genuinely per-tower (not flat);
HCAL's is currently flat 1.0, correctly reflecting that HCAL has no
per-tower calibration yet, while HCAL's base `gain` still includes the
known 1.3 electronics-gain factor. `data_to_jet`'s `runMudst.C` is
unaffected -- it reads already-officially-produced muDSTs and never
touches `StFcsDb`'s gain settings itself.

**Pinning a specific calibration era (`-c`/`calibRun`).** By default the
BFC chain queries the conditions DB at its own nominal `y2023` chain-flag
date (`20230410`), so every simulated job uses the same single gain/
gain-correction era regardless of which real run period it's meant to
model. `submit_sim_to_jet.sh`'s `-c <run_number>` flag (passed through
to `runSimBfc.C`'s `calibRun` parameter, and from there to
`sim_to_jet.xml`'s `{{CALIB_RUN}}`) overrides this: given a real Run 22
run number (e.g. `23101043`), `runSimBfc.C` converts it to a calendar
date+time via `RunNumberToDateTime()` and calls `St_db_Maker::SetDateTime()`
with it before `chain->Init()`, so `StFcsDb` (still on default `GAINMODE::DB`,
per above) picks up that run's actual calibration era instead. Omit `-c`
(or pass `0`) to keep the previous, unchanged nominal-date behavior.
`RunNumberToDateTime()` queries the real, authoritative `RunLog.runDescriptor`
DB (`db04.star.bnl.gov`) for the run's actual start time, falling back to
`RunNumberToDateArithmetic()` (day-of-year arithmetic, `itim=0`) if that
DB is unreachable -- e.g. a SUMS worker node with more restricted network
access than an interactive session. The arithmetic fallback accounts for
a non-obvious STAR run-number quirk: the run number's 2-digit year code
is *not* the real calendar year, it's one less (e.g. run `22359013`'s
real date is 2021-12-25, not 2022-12-25). Each job's log prints a
`Gain check (calibRun=..., ...)` line after `EventLoop()` (representative
ECAL/HCAL channel gain+gainCorr) so it's obvious at a glance whether the
override took effect as expected.

**Fiducial-region geometry cut (`../shared/JetParameters.h`).** `cut_x_inner`/
`cut_x_outer`/`cut_y_min`/`cut_y_max` (used by both `data_to_jet` and
`sim_to_jet`, and by `jet_scale`/`jet_calibration` downstream) are
reproducible directly from the STAR conditions DB via a standalone tool in
`../../GeometryPlayground/`: `run_number_to_fiducial.sh <run_number>` pulls
`StFcsDb`'s detector-position geometry for that run's era and computes the
ECAL/HCAL fiducial rectangle from it. See the comment above `cut_x_inner`
in `JetParameters.h` for provenance. Detector position is constant across
essentially all of Run 22 (one transition, 2021-12-21), so a single
hardcoded rectangle is valid pipeline-wide -- re-run the tool only if
there's reason to think the geometry changed (a different run period, a
documented detector move).

## data_to_jet

Real STAR MuDst data has already been through official offline production
(BFC, FCS clustering), so this pipeline just reads it and finds jets:
MuDst -> SimpleTree (`runMudst.C`) -> jet finding, for every file in a
catalog list, split evenly across `-j` jobs.

**Step 1 -- get a file list.**

```
perl MakeRunList.pl -f mylist.list -l 50
```

Queries the STAR file catalog (via the `StarFileList` Perl module) for
real `daq_reco_mudst` files and writes one URL/path per line, with the
event count appended as a trailing number (each pipeline job strips it
automatically). `-l` caps how many files it returns; see `perldoc
MakeRunList.pl` for all options (run range, production tag, etc. are
currently fixed for Run 22 pp500 forward production in the script itself).
This depends on another user's shared Perl library path
(`/star/u/dkap7827/Tools2/Tools/PerlScripts`); if that ever moves, this
script needs its `use lib` line updated.

**Step 2 -- run the pipeline.**

```
./submit_data_to_jet.sh -f <listfile> -s <0|1> -o <outdir> \
    [-j <nprocesses>] [-l <label>] [-k <0|1>]
```

| Flag | Meaning |
|---|---|
| `-f` | absolute path to the list file from step 1 |
| `-s` | `1` to keep each intermediate `SimpleTree_mudst_*.root`, `0` to discard it |
| `-o` | output directory (created if missing) |
| `-j` | number of parallel jobs to evenly split the list's lines across, default 1 |
| `-l` | short label used in the generated XML's filename, default `run` |
| `-k` | `1` to save each job's stdout/stderr under `<outdir>/log/`, `0` to discard them, default `0` |

Each job fetches its assigned files (`xrdcp` if xrootd URLs, plain `cp`
otherwise) and runs `runMudst.C` on each to get its SimpleTree. Runs
inside the same singularity container as `sim_to_jet`; `root4star` only
reliably resolves that way under SUMS here.

**1:1 SimpleTree/JetTree linkage.** A job assigned several files does
*not* combine them into one `RecoJets` call -- it loops and calls
`RecoJets` once per SimpleTree file, naming each `jet_output_<suffix>.root`
with the same run/segment suffix as the `SimpleTree_mudst_<suffix>.root`
it came from, so the two can be linked back together later to look up the
same event in both.

**FCS trigger flags.** Both the SimpleTree (`data`, written by
`StSimpleReaderMaker`) and the JetTree (`jetTree`, written by
`RecoJets`) carry one `int` branch per known FCS trigger name from
`TriggerIDs.txt` -- 1 if that trigger fired on the event, 0 otherwise.
Resolved from `StMuEvent::triggerIdCollection().nominal().triggerIds()`
against `FcsTriggerDefs.h`, a generated lookup table (see
`generate_trigger_defs.py`) matching each fired numeric trigger id to a
name using *both* the id and the run number, since the online trigger
system reassigns numeric ids to different named triggers across the run
period. (`sim_to_jet` resolves flags differently -- see its own
section.) `FcsTriggerDefs.h` exists as three copies that must be
regenerated together: one under `StRoot/StSpinPool/StSimpleReaderMaker/`
and one each in `data_to_jet/` and `sim_to_jet/`. If `TriggerIDs.txt`
changes, rerun `python3 generate_trigger_defs.py` from `data_to_jet/`
(its source-of-truth copy) and rebuild StRoot.

**Per-event spin configuration.** Both the SimpleTree and JetTree also
carry a single `Spin_config` `int` branch: `StSpinDbMaker`'s combined
spin4 code (0-10, see `StSpinDbMaker::optimizeTables()` for the encoding
legend), resolved per event from the real bunch-crossing id against the
offline `Calibrations/rhic` conditions DB. **`-1` means "no data"** --
the offline spin DB has real per-fill entries only for **2021-12-14
through 2022-04-18** (see `SpinPlayground/run_number_to_spin_coverage.sh
<run_number>` to check a specific run). Any event outside that window
correctly gets `Spin_config=-1`, not a bug. Resolved in
`StSimpleReaderMaker::Make()` (`StSpinDbMaker` attached via
`SetSpinDb()`, only called by `data_to_jet/runMudst.C`). `sim_to_jet`
carries the same branch for schema consistency but never attaches an
`StSpinDbMaker`, so it's always `-1` there.

## jet_scale

Takes a folder of `jet_output_*.root` (JetTrees) files -- from either
`sim_to_jet` or `data_to_jet` -- and produces the JSON lookup table used
to correct reco jet energies. Two steps, one job (see "Job-splitting"
above for why there's no `-j`): `JetEnergyScaleFineGrid` scans every file
in the folder and fits a scale/response Gaussian on a fine (E, x, y) grid
(references the shared `shared/JetParameters.h` for the fiducial-cut
geometry), then `CalibrationMap` fits a smooth model per position bin and
writes the final JSON (this step does not reference `JetParameters.h`,
only the intermediate fit tree).

```
./submit_jet_scale.sh -i <jettrees_dir> -o <outdir> [-n <json_name>] \
    [-s <0|1>] [-l <label>] [-k <0|1>]
```

| Flag | Meaning |
|---|---|
| `-i` | absolute path to a folder of `jet_output_*.root` files |
| `-o` | output directory (created if missing) |
| `-n` | filename for the output JSON lookup table, default `JetEnergyScale_lookup.json` |
| `-s` | `1` to keep the intermediate `jet_calibration.root`, `0` to discard it (default `0`) |
| `-l` | short label used in the generated XML's filename, default `run` |
| `-k` | `1` to save stdout/stderr under `<outdir>/log/`, `0` to discard them, default `0` |

Example:

```
./submit_jet_scale.sh -i /star/data01/pwg/seanp/tunes/pythia8_mstw2008lo/JetTrees \
    -o /star/data01/pwg/seanp/tunes/pythia8_mstw2008lo \
    -n JetEnergyScale_lookup_pythia8_mstw2008lo.json
```

It can also just be compiled and run directly on the login node instead
of through SUMS, if preferred (see the two `g++`/`./` lines in
`jet_scale.xml`'s `<command>`).

## jet_calibration

Takes a folder of `jet_output_*.root` (JetTrees) files and a calibration
JSON (from `jet_scale`) and adds the corrected branches (`reco_E_corr`,
`reco_x_F_corr`, `reco_E_uniform_corr`) to those files. **This modifies
the files in place** -- `ApplyCorrections` opens each with `TFile
"UPDATE"` and overwrites the tree; there are no separate output files and
no undo. References the shared `shared/JetParameters.h`
(`computeFeynmanX`, for `reco_x_F_corr`).

`ApplyCorrections.cpp` reflects `reco_x` to positive before looking it up
in the JES JSON grid (`fabs(x)`), matching `jet_scale`'s
`JetEnergyScaleFineGrid.cpp`, which folds x the same way when building
that grid (`reflect_to_positive_half`). `reco_y` is looked up unreflected
(`y`, not `fabs(y)`) -- the y-grid was never folded on the build side, so
folding it here would look up the wrong cell for every negative-y jet.

```
./submit_jet_calibration.sh -i <jettrees_dir> -c <calib_json> -j <nprocesses> \
    [-l <label>] [-k <0|1>]
```

| Flag | Meaning |
|---|---|
| `-i` | absolute path to a folder of `jet_output_*.root` files; every file in it is modified in place |
| `-c` | absolute path to the calibration JSON lookup table |
| `-j` | number of parallel jobs to evenly split the folder's files across |
| `-l` | short label used in the generated XML's filename, default `run` |
| `-k` | `1` to save each job's stdout/stderr under `<jettrees_dir>/log/`, `0` to discard them, default `0` |

The wrapper asks for a `y`/`N` confirmation before submitting, since this
is destructive. Example:

```
./submit_jet_calibration.sh -i /star/data01/pwg/seanp/tunes/pythia8_mstw2008lo/JetTrees \
    -c /star/data01/pwg/seanp/tunes/pythia8_mstw2008lo/JetEnergyScale_lookup_pythia8_mstw2008lo.json \
    -j 10
```

## pico_analysis

**Not part of the jet-finding chain.** A standalone diagnostics utility
for Run 22 pp500 forward-tracking production data, which is primarily
distributed as PicoDst rather than MuDst. Neither mode here produces a
JetMatcher/RecoJets-compatible SimpleTree, so its output can't feed
`jet_scale` or `jet_calibration` -- it's kept alongside the jet pipelines
for future PicoDst-based analyses, not folded into them. Always a single
job (both modes process their whole input list as one aggregate; no `-j`).

```
./submit_pico_analysis.sh -d <mip|hits> -f <listfile> -o <outdir> \
    [-l <label>] [-k <0|1>]
```

| Flag | Meaning |
|---|---|
| `-d` | `mip` (MIP peak analysis, `mip_ana.C`: projects tracks onto the FCS ECal/HCal, writes diagnostic plots to `mip_ana.pdf`) or `hits` (flat hit-level tree, `pico_to_root.C`: FCS hit energy/position/detector ID + vertex position per event, no clustering, no jets) |
| `-f` | absolute path to a text file listing PicoDst files, one per line |
| `-o` | output directory (created if missing) |
| `-l` | short label used in the generated XML's filename, default `run` |
| `-k` | `1` to save stdout/stderr under `<outdir>/log/`, `0` to discard them, default `0` |

Example -- build a list of PicoDst files and run the MIP analysis on them:

```
ls -1 /gpfs01/star/pwg_tasks/FwdCalib/PROD/forwardCrossSection_2022/27082025/*.root > mylist.list
./submit_pico_analysis.sh -d mip -f $(pwd)/mylist.list \
    -o /star/data01/pwg/seanp/pipeline_output/mip_test
```

No local compilation needed; the `stardev` environment already includes
the PicoDst/Forward-detector classes these macros use, and `root4star`'s
ACLiC (`+`) compiles each macro inline at runtime.

## Notes

- Every `submit_*.sh` script runs `star-submit` from the pipeline root
  (required: each XML's SandBox `<File>` paths are relative to that
  directory, and SUMS packages them preserving that relative path; an
  absolute path would repackage everything under the wrong subpath and
  break each pipeline's `cp <subdir>/* .` flattening step). Right after
  the `star-submit` call, each script sweeps the newly-created
  `sched*`/`*.package`/`*.csh`/`*.list`/`*.condor` files into `submit/`, a
  dedicated scratch folder, so they don't clutter the pipeline root.
- To hand-edit and submit without a wrapper, see the usage block at the
  top of the relevant `<pipeline>/<pipeline>.xml`.
- If `StRoot/` is ever edited, rebuild `.sl73_x8664_gcc485/` (normal
  `stardev`/cons build flow) before submitting -- jobs use the precompiled
  cache, they do not rebuild from source.
