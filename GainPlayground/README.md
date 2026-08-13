# GainPlayground

How FCS per-channel gain calibration is stored, updated, and actually
consumed by `sim_to_jet` and `data_to_jet` -- with a minimal standalone
example of retrieving it. Unlike `GeometryPlayground`'s finding (detector
position is fixed for essentially all of Run 22), **gain/gainCorrection is
genuinely time-dependent across the run** -- there is no single "Run 22
gain" to record, only a mechanism to look one up for a given run/date.

**Important, action-relevant finding**: as of this investigation (Aug
2026), HCAL's `getGainCorrection()` is a flat `1.0` placeholder for all of
Run 22 -- there is no real per-tower HCAL calibration in the conditions DB
yet (see the two-run comparison below; ECAL, by contrast, already has a
genuine per-tower correction, e.g. gainCorrection=1.21-1.43). **Once real
per-tower HCAL gain corrections are measured and become available, every
`data_to_jet` SimpleTree/jetTree already produced under the flat-1.0
assumption will need to be retroactively corrected** -- the muDst hit
energies they were built from went through official production once,
using whatever HCAL calibration existed at production time, and won't
silently update. The machinery to apply such a correction retroactively
(without needing to rerun `runMudst.C` from raw muDst files) already
exists -- see `data_to_jet/hcal_gain_corrections/README.md` -- currently
loaded with a DB-verified all-`1.0` no-op (confirmed by actually scanning
the conditions DB, see `scan_hcal_gain_db/README.md`, not just assumed).
**When real per-tower HCAL corrections do become available, rerun
`scan_hcal_gain_db/scan_hcal_gain_db.C`** -- it will detect the new real
calibration period(s) and regenerate the pipeline's correction files and
`manifest.txt` automatically, no manual DB inspection needed.

## What's actually in the DB

`StFcsDb` reads two related-but-different quantities, both from
`Calibrations/fcs` (date-indexed, same `St_db_Maker`/`StDbLib` mechanism as
everything else in this project), both via the *same* `StFcsDbMaker` that
also loads `Geometry/fcs` -- one `InitRun()` call pulls both:

- **`getGain(det, id)`** -- base electronics gain per channel (from
  `fcsEcalGain`/`fcsHcalGain`/`fcsPresGain`). Mostly a fixed hardware
  property (e.g. HCAL's known 1.3x electronics-gain factor), though it can
  also shift between run periods (see the two-run comparison below).
- **`getGainCorrection(det, id)`** -- per-tower relative calibration
  correction (from `fcsEcalGainCorr`/`fcsHcalGainCorr`), the thing that
  actually varies channel-to-channel and, importantly, **run-to-run**.
  There's also `getGainOnline()` (threshold-related, `*GainOnline` tables)
  and `getEtGain()`, not covered here -- see `StFcsDb.h` if needed.

**How the DB gets filled**: not by this pipeline. STAR's FCS calibration
group periodically derives per-tower gain corrections from real physics
data (e.g. pi0-mass-fit calibration) at specific "anchor" runs, then
interpolates between anchor points across a run period using BBC/ZDC
luminosity scalers to track gain drift over time -- see
`StFcsDbMaker/macro/fcsTimedepGainCorr_db.C` (an official STAR macro, not
part of this pipeline) for the actual interpolation logic. This project
only ever *reads* the result of that process via `StFcsDb`.

## Files

- **`get_fcs_gain.C`** -- minimal BFC chain (DB + `StFcsDbMaker` only, no
  simulation/input file needed), same pattern as
  `GeometryPlayground/get_fcs_geometry.C` (reuses the identical
  `RunNumberToDateTime()` RunLog-based run→date resolution). Prints and
  saves `getGain()`/`getGainCorrection()` for one representative channel
  each of ECAL (det 0) and HCAL (det 2), for a given `calibRun`.
- **`fcs_gain_example_22359013.txt`** -- saved example output for run
  22359013 (the same representative run used elsewhere in this project).
  A single-run snapshot, not a reusable constant -- see the file's own
  header and the time-dependence demonstration below.
- **`scan_hcal_gain_db/`** -- automated DB scanner that regenerates
  `data_to_jet/hcal_gain_corrections/`'s correction files and
  `manifest.txt` directly, by actually querying the conditions DB across
  Run 22 and detecting real calibration-period boundaries. Rerun this
  whenever HCAL's real per-tower corrections change or first become
  available -- see its own README.md for how it works and how it was
  validated.

## Usage

```
root4star -b -q get_fcs_gain.C(<run_number>)
```

e.g. `root4star -b -q get_fcs_gain.C(22359013)`, via the usual
singularity + `stardev` environment (see `GeometryPlayground/README.md`
for the exact invocation pattern).

## Gain genuinely changes across Run 22 (unlike geometry)

Two runs, same channels, run via `get_fcs_gain.C`:

| Run | Date | ECAL(det0,id50) gain | ECAL gainCorrection | HCAL(det2,id50) gain | HCAL gainCorrection |
|---|---|--:|--:|--:|--:|
| 22359013 | 2021-12-25 | 0.0053 | 1.21 | 0.00689 | 1.0 |
| 23101043 | 2022-04-11 | 0.0053 | 1.4266 | 0.00834 | 1.0 |

ECAL's per-tower correction and HCAL's base gain both shifted noticeably
over ~3.5 months; HCAL's gainCorrection stays flat 1.0 at both dates,
consistent with HCAL not having per-tower calibration (its base gain
already folds in an average ECAL-derived correction instead -- see
`sim_to_jet/runSimBfc.C`'s own comment on this). **Conclusion: always
resolve gain for the specific run/date you care about; never assume a
value found for one run applies to another.**

## How each pipeline actually uses this

### `sim_to_jet` -- gain IS applied, deliberately

`StFcsFastSimulatorMaker` (the maker that turns truth-level GEANT energy
deposits into a realistic simulated detector response) calls
`mFcsDb->getGain(det,id)` and `mFcsDb->getGainCorrection(det,id)` directly
(`StFcsFastSimulatorMaker.cxx:313-314`) while building each simulated hit.
`StFcsDb` is left on its default `GAINMODE::DB` (no
`forceUniformGain()`/`forceUniformGainCorrection()` override) in
`sim_to_jet/runSimBfc.C`, so this pulls the *real* per-channel calibration
from the conditions DB -- and `runSimBfc.C`'s `calibRun` parameter (same
`RunNumberToDateTime()` mechanism as this playground) pins *which* real
run's calibration era gets used, so simulated events reproduce the actual
per-tower nonuniformity of a specific real data-taking period rather than
a flat placeholder. `runSimBfc.C` prints a "Gain check" sanity-check line
after `chain->EventLoop()` confirming which values actually got used --
this playground's `get_fcs_gain.C` is that same check, standalone.

### `data_to_jet` -- live gain is NOT applied here; a retroactive HCAL correction is

No `getGain`/`getGainCorrection` call appears anywhere in
`data_to_jet/runMudst.C` or `StSimpleReaderMaker` -- correct, not an
oversight. `data_to_jet` reads real, already-officially-reconstructed
MuDst files; gain calibration was already applied once, for real, by
STAR's official production chain (`StFcsClusterMaker` etc., run by the
production team at BFC time using the real per-channel calibration for
that real run) before the MuDst ever reached this pipeline.
`data_to_jet` just reads each `StMuFcsHit`'s already-calibrated energy
field as its starting point -- reapplying *live* gain here would
double-correct it.

On top of that starting point, `StSimpleReaderMaker` *does* apply one
thing: an optional retroactive HCAL correction, multiplying each HCAL
hit's energy by a per-channel factor read from a plain text file (not the
DB) -- see the "Important, action-relevant finding" note above and
`data_to_jet/hcal_gain_corrections/README.md` for the full mechanism.
`runMudst.C` calls `SetHcalRetroactiveGainCorr()` unconditionally, but as
long as that folder's correction file stays all-`1.0`, this has
mathematically zero effect (multiplying by 1.0 is exact, not an
approximation) -- it only starts mattering once real per-tower HCAL
corrections replace the placeholder.
