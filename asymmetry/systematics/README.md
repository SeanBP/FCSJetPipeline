# A_N systematics

Tools for the A_N systematic checks (started 2026-11-17). The plan, item by item,
is in FCSJetAnalysis/ANALYSIS_ROADMAP.md ("Systematics plan").

## Method

**Nested cut variations.** A loose cut contains a tight one, for example trigger
value ≥ 1.15 vs ≥ 1.3 × threshold. Vary one cut at a time and keep the others at
nominal. Per x_F bin:

- Δ = A_loose − A_tight, both computed with the same cross-ratio code (AnAsymmetry).
- σ_diff = √|σ_tight² − σ_loose²|. For inverse-variance weights this equals
  w √(σ_slice² + σ_tight²), where the slice is loose minus tight and
  w = 1 − σ_loose²/σ_tight² is the slice's weight.
- Barlow test: sys = √(Δ² − σ_diff²) if |Δ| > σ_diff, else 0.
- Per row: χ² = Σ (Δ/σ_diff)² over bins, as a check for trends.

The significant pieces of the different variations add in quadrature per bin.

**Disjoint samples.** These are consistency checks, such as runs early vs late in
the fill. Here σ_diff = √(σ_1² + σ_2²), and the result is a χ², not a
systematic box.

**Trigger cut.** The nominal cut is 1.15 × the period's turn-on x50 (since
2026-11-17). For EM the threshold and the HCal veto move together. The trigger
systematic compares the nominal sample (loose, 1.15) with 1.3 (tight, a
subset): `sys_slices.py an_f115.txt an_f130.txt ...`.

## Files

| file | what it does |
|---|---|
| `sys_slices.py` | Δ, σ_diff, w, the slice's A_N, Barlow sys and pull per bin from two AnAsymmetry `.txt` tables (`--disjoint` for independent samples). Writes `<prefix>.txt` and `<prefix>_summary.txt` (χ², p-value per row). |
| `plot_sys_slices.C` | One row: A_N of the tight, loose and slice samples on the common A_N axes, and Δ ± σ_diff with the Barlow band. |
| `time_in_fill.py` | Hours in fill per run: t_mid from `../run22_polarization.txt` minus the fill start from the polarimetry fill table. Splits a yield list into bins: `--edges` in hours (default 2, as in FMS) or `--quantiles N`. |
| `sys_an_job.sh`, `submit_sys_an.sh` | AnAsymmetry -r over a list of yield files as condor jobs. This is the per-file flow; the merged flow is faster for whole passes. |

## Merged flow for subsets

The subsets are run lists, for example by time in fill. For each subset:

1. Make a directory whose `hists/` holds symlinks to the subset's yield files.
2. Run `sivers/submit_an_merge.sh <dir>`.
3. Run `AnMerge -a`, then `AnMerge -x`. The masking only subtracts runs that are
   in the merge.
4. Run `AnAsymmetry -M` per sample and option. Adding the partials of all subsets
   gives the full sample.
