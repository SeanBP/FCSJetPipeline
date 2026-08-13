# scan_hcal_gain_db

Automated DB scanner for `data_to_jet`'s retroactive HCAL gain-correction
machinery (`FCSJetPipeline/data_to_jet/hcal_gain_corrections/`). Queries
the real STAR conditions DB across a date range, detects when HCAL's
per-tower gain correction actually changes (a new real calibration
period), and writes the correction-factor file(s) + `manifest.txt` that
`data_to_jet` reads -- automatically, with no manual DB inspection.

## Why this exists

As of this investigation (Aug 2026), HCAL has no real per-tower gain
calibration in the conditions DB yet -- `StFcsDb::getGainCorrection()`
returns a flat `1.0` for every HCAL channel, for all of Run 22 (see
`../README.md`). `data_to_jet/hcal_gain_corrections/` exists so that,
once STAR's FCS calibration group determines real corrections,
already-produced SimpleTree files can be corrected retroactively without
rerunning `runMudst.C` from raw muDst files. This script is the tool that
turns "real corrections now exist in the DB" into "the pipeline's
correction files are updated" -- rerun it, and it does the rest.

## Usage

```
root4star -b -q scan_hcal_gain_db.C(20211201,20220701,7)
```

(those are the defaults -- calling with no arguments does the same
thing). Needs the STAR/singularity environment, same invocation pattern
as `../get_fcs_gain.C` / `GeometryPlayground`'s tools.

- `startDate`, `endDate`: `YYYYMMDD`, the date range to scan (defaults
  generously cover all of Run 22 -- see `../../SpinPlayground/README.md`
  and `../../GeometryPlayground/README.md` for the real, narrower
  data-taking window this comfortably contains).
- `stepDays`: sampling interval. Default 7 (weekly) -- real Run 22 ECAL
  calibration periods (the one detector that *does* currently vary,
  used to validate this script's period-detection logic directly) are
  spaced weeks to ~2 months apart, per `StFcsDbMaker/macro/
  fcsTimedepGainCorr_db.C`'s own period table, so weekly sampling
  comfortably resolves real period boundaries without excessive DB
  round-trips. Shrink it if a faster change is suspected.

**CAUTION**: this overwrites `manifest.txt` and adds/updates correction
files directly in `data_to_jet/hcal_gain_corrections/`. It backs up any
existing `manifest.txt` to `manifest.txt.bak` first (overwriting any
previous backup) -- this repo has no version control, so that single
backup is the only undo available. The correction-factor files themselves
aren't backed up (they're regenerated deterministically by this same
script, so there's nothing to lose by overwriting them).

## How it works

1. Enumerates every real HCAL channel once via `StFcsDb::maxId()` +
   `getDepfromId()` (the official iteration STAR's own DB-fill code
   uses, not hand-guessed loop bounds) -- the channel *map* is stable
   across all of Run 22, so this only needs doing once, not per date.
2. Steps through `[startDate, endDate)` in `stepDays` increments. At each
   sampled date, re-pins `St_db_Maker`'s date (`SetDateTime`) and forces
   `StFcsDb` to refresh (`StFcsDbMaker::InitRun()`), then reads
   `getGainCorrection()` for every enumerated channel -- the same live
   framework call `../get_fcs_gain.C` uses for a single date, just
   repeated across many dates in one process.
3. Whenever *any* channel's value differs from the previous sample, that
   marks a period boundary: the just-finished period (and its full
   per-channel value snapshot) gets written out as a correction file, and
   a new period starts at the changed date.
4. Writes `manifest.txt` with one row per detected period, pointing at
   its correction file.

### Why live API calls, not raw SQL

STAR's conditions DB (`Calibrations_fcs.fcsHcalGainCorr`, same versioned
schema as everywhere else in this project -- `dataID`, `beginTime`,
`endTime`, `flavor`, `deactive`, etc.) can have multiple overlapping
entries for the same date range (an old baseline plus a newer patch, for
instance). `SpinPlayground/README.md` documents a case where this
project's own diagnostic tooling couldn't fully confirm which entry
`StDbLib`'s real conflict-resolution rule would pick in such a case. To
sidestep that question entirely, this script never inspects the raw DB
tables -- it always asks `StFcsDb::getGainCorrection()` directly, the
exact same call the production pipeline itself uses, so the result is
guaranteed to match what `data_to_jet` would see live, by construction.

### Validation (Aug 2026)

- Repeated `SetDateTime()`+`InitRun()` calls within one process correctly
  refresh `StFcsDb`'s cached values -- confirmed directly by scanning two
  dates known to give different ECAL `gainCorrection` (1.21 vs ~1.56) and
  confirming both the change and the reversion back to the first date's
  value.
- Channel enumeration output matches `generate_placeholder.C`'s
  independently-written enumeration exactly (520/520 rows, byte-identical
  `ehp ns dep ch` columns).
- Change-detection was validated against ECAL (temporarily pointed at
  detectors 0/1 instead of 2/3, in a throwaway copy) -- correctly
  discovered **9 distinct real calibration periods** across Run 22,
  matching the period count in STAR's own official reference macro
  (`StFcsDbMaker/macro/fcsTimedepGainCorr_db.C`'s `NPERIOD=9` table).
- Run for real against HCAL (its actual purpose): correctly found zero
  changes across all of Run 22 and produced one period, consistent with
  the flat-`1.0` finding in `../README.md`.
