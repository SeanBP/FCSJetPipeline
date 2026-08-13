# TriggerPlayground

Minimal, standalone FCS trigger-simulation validation tool. Runs STAR's
real `fcs_trg_base` trigger-decision engine (via `StFcsTriggerSimMaker`) on
a real muDst's actual FCS hits, and cross-checks the simulated decision
against that event's actually-recorded fired trigger id -- i.e. answers
"does the trigger simulator agree with what really fired?"

This is the validation that justified using `fcs_trg_base` as
`sim_to_jet`'s per-event trigger source (`StSimpleReaderMaker::SetTriggerSim()`
in the production pipeline, `FCSJetPipeline/StRoot/StSpinPool/StSimpleReaderMaker/`)
for *simulated* events, which have no real fired trigger id of their own.
Real data (`data_to_jet`) never uses the simulator -- it reads the real
fired trigger id straight from the muDst (see `FcsTriggerDefs.h` in the
production pipeline). Nothing here is used at production run time; it's a
one-off correctness check on the simulation engine itself.

## Files

- **`TrgSimTest.C`** -- driver macro. Sets up a minimal chain (DB +
  `StFcsDbMaker` + `StFcsTriggerSimMaker`, no full `StEvent`) and hands the
  whole per-event loop to `TrgSimHelperMaker::Make()`.
- **`StRoot/StSpinPool/TrgSimHelperMaker/`** -- `cons`-built maker that
  drives `fcs_trg_base` directly (bypassing `StFcsTriggerSimMaker::Make()`,
  which segfaults outside a full `StEvent` chain) and accumulates a
  real-vs-simulated confusion matrix per named trigger. See its header for
  why the per-event loop had to live in compiled code, not the interpreted
  macro (a CINT interpreter-corruption issue, not a logic bug).
- **`TriggerIDs.txt`** -- trigger-id-to-name table (name/id/valid run
  range), needed to map the muDst's raw fired trigger id(s) to names for
  the cross-check.
- **`stage_params.txt`** -- `fcs_trg_base` threshold values, since
  `StFcsTriggerSimMaker`'s DB-threshold path (`setThresholdDb()`) is an
  unimplemented stub; this is the real, working mechanism
  (`setThresholdFile()`), copied from `StFcsTriggerSimMaker/files/`.
- **`build_helper.csh`** -- `cons`-builds `StRoot/` above into
  `.sl73_x8664_gcc485/` (not shipped here -- rebuild before first run).
- **`run_trgsim_test.csh`** -- runs `TrgSimTest.C` with its defaults.

## Usage

```
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
    /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh build_helper.csh
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
    /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh run_trgsim_test.csh
```

No local muDst is shipped here -- `TrgSimTest.C`'s default `file` argument
is a `root://...` xrootd URL (a real Run 22 file, run 22359013) that
`StMuDstMaker`/`TFile::Open` stream directly, no download needed. Pass any
other accessible muDst path (local or `root://`) as the first argument to
test a different run.

## Interpreting the output

Real L1 trigger bits are heavily **prescaled** (only a fraction of events
where the physics condition was true get flagged "fired" in the real
data), while the simulator reproduces the raw physics condition with no
prescaling applied. So `sim=1, real=0` is *expected and benign* for
high-rate triggers -- it is not evidence the simulator is wrong.

The real validation question is **P(sim=1 | real=1)**: whenever a trigger
really fired, did the simulator also compute that its physics condition
was met? That should be ~100% if `trgSelect` and the thresholds are
correct; a low value there is a genuine red flag. `TrgSimHelperMaker`
prints this breakdown per named trigger at `Finish()`.
