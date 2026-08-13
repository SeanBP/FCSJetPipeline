# SpinPlayground

Two minimal, standalone tools for the STAR offline spin/polarization
pattern DB (`Calibrations/rhic`'s `spinV124`/`spinStar` tables): one checks
whether a *run* has coverage at all, the other actually retrieves the real
per-*event* spin state.

This is where the `Spin_config` branch (documented in
`FCSJetPipeline/README.md` and `StSimpleReaderMaker.h`) came from, and the
tools to re-check it if STAR's spin group backfills more of the DB later
(the production code re-queries the DB live every run, so it will pick up
new entries automatically -- see the coverage caveats below for what to
re-verify if that happens).

## Files

- **`run_number_to_spin_coverage.sh`** -- run-number -> DB coverage.
  Pure bash + `mysql` client, no ROOT/StRoot/singularity needed: both DB
  servers (`db04.star.bnl.gov` for `RunLog`, `dbx.star.bnl.gov:3316` for
  `Calibrations_rhic`) are reachable directly from this host. Its header
  comment has the full investigation writeup -- read it before trusting a
  coverage verdict blindly, especially the two documented stale-fallback
  traps below.
- **`read_spin_state.C`** -- muDst event -> actual resolved spin state.
  Needs the STAR/ROOT environment (reads real FCS/StEvent data from a
  muDst, not just a DB timestamp check). For each event, prints the run
  number, 7-bit bunch crossing id, `StSpinDbMaker::isValid()`, and the
  resolved `spin4usingBX7()` combined code (0-10) decoded into
  yellow/blue up/down -- exactly the same calls production's
  `StSimpleReaderMaker::Make()` uses to fill the `Spin_config` branch, just
  printed instead of written to a tree. Streams a known-good real Run 22
  muDst via xrootd by default (no local file needed); pass a different
  path/URL as the first argument to check another run.

## Usage

```
./run_number_to_spin_coverage.sh <run_number>
```

e.g. `./run_number_to_spin_coverage.sh 22359013`.

```
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
    /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
    'cd /star/u/seanp/FCSJetPipeline/SpinPlayground; source /star/nfs4/AFS/star/group/star_cshrc.csh; stardev; root4star -b -q read_spin_state.C'
```

or pass a specific file/event count: `read_spin_state.C("<path-or-root://...>", 50)`.

## Why a direct timestamp check, not date-sampling

`spinV124`/`spinStar` entries are narrow, fill-scoped windows (a few hours
each) with real gaps between fills -- unlike e.g. the FCS geometry/gain
tables, there's no "most recent entry persists until superseded" here. A
date-sampled scan can land in a gap and look empty even where the DB is
well populated nearby, so this script resolves the run's *exact* real start
time (via `RunLog.runDescriptor`) and checks it against each entry's
`[beginTime, endTime)` window directly.

## Two verdicts, and why

- **VERDICT 1** mirrors `StSpinDbMaker::isValid()` literally (all of
  `spinV124`/`spinStar`/`spinBXmask` resolve to *some* entry) -- but
  `isValid()` has no concept of a stale/fallback entry, so it can return
  `TRUE` even where there's no real per-fill data.
- **VERDICT 2** is the one that actually matters for per-event spin info:
  `spinV124` + `spinStar` both resolve to a *real* (non-suspect,
  short-duration) entry. `spinBXmask` is excluded -- it has zero genuine
  Run 22 entries and always resolves via a stale 2013 default.

Known stale-fallback traps (both flagged automatically by the script's
`>7 day` duration heuristic): `spinBXmask`'s 2013 open-ended default, and
one `spinStar` entry (`dataID 22639`, inserted 2025, comment
"Run23091017, oleg") backdated to cover 2022-04-01 onward.

## Finding this tool established (Aug 2026)

Real per-fill coverage in `spinV124`/`spinStar` exists only for
**2021-12-14 through 2022-04-18**, with no gap in that window. Zero
coverage either side of it, for the rest of Run 22. Verified: run
22359013 -> covered; run 23101043 -> not covered (falls in the post-04-18
gap).
