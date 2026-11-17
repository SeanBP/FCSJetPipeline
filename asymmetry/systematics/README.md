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
| `pi0_iso_sys.sh` | π⁰ isolation-threshold systematic (iso 0.95 / 0.98 / 0.99, noniso 0.85 / 0.9 / 0.95) on a report folder. |
| `sys_variations.sh` | The cut variations stored as yield copies (AnConfig.h `kVarSuffix`, `AnAsymmetry -v`): fiducial 25 / 30 cm, \|vz\| < 60, measured vertex only (nested, Barlow), and the collision-cut failures (disjoint check, plus their share per row in `bbcveto_share.txt`); photon classes (`-J`) and iso / noniso π⁰. On a report folder; outputs in `<report>/sys_variations/`. |
| `combine_sys.py` | The combine step: nominal table + components (`--comp NAME=sys_slices output`, `--comp-max` for the larger of several variations, `--sysf` for the π⁰ signal-fraction column, `--rel` for flat relative pieces) → per-bin stat, each component, quadrature total, stat ⊕ sys and the x_F uncertainty (`--xf-rel`, FMS 3%). The polarization scale is quoted separately. |
| `plot_an_final.C` | The final plots of one row from a combined table: `_stat`, `_sys` (total boxes), `_sys_components`, `_combined` (stat bars + sys boxes), `_total` (stat ⊕ sys bars), on the common A_N axes. |
| `misid_response.cpp` | Photon-class response from simulation (JES production `jet_output_*`, cross-section weighted): reco class × true class per truth-x_F bin, true photons with E >= 1 GeV (`--truth-emin`), offline trigger-value cut (the production's emulated flags are 0); `<out>.txt` per class and `<out>_merged.txt` for np12 / np3plus. |
| `misid_sys.py` | Misidentification systematic for the merged classes: dA = f(other \| reco) × (A_other − A), Barlow-reduced by the difference's statistical error; `combine_sys.py --comp-simple misid=<out>.txt`. |
| `pol_uncertainty.py`, `dump_run_jets.C` | Fill-to-fill and profile polarization uncertainty (a note; not used: only the overall scale is quoted, as FMS). |

## Merged flow for subsets

The subsets are run lists, for example by time in fill. For each subset:

1. Make a directory whose `hists/` holds symlinks to the subset's yield files.
2. Run `sivers/submit_an_merge.sh <dir>`.
3. Run `AnMerge -a`, then `AnMerge -x`. The masking only subtracts runs that are
   in the merge.
4. Run `AnAsymmetry -M` per sample and option. Adding the partials of all subsets
   gives the full sample.

## Full pass: from yields to the final plots

1. Report (`sivers/submit_report_job.sh`), giving `masked.root` and `build/`.
2. Trigger: `AnAsymmetry -f 1.3` vs the default → `sys_slices.py`.
3. Cuts: `sys_variations.sh <report>`; π⁰ isolation: `pi0_iso_sys.sh <report>`.
4. `combine_sys.py`, per result (e.g. union photon classes; iso / noniso π⁰ with
   `--sysf`), then `plot_an_final.C` per row.

Example (photon classes):

```
combine_sys.py --nominal sys_variations/an_nominal.txt --comp trigger=sys_trig.txt \
    --comp-max fiducial=sys_variations/sys_fid25_jets.txt,sys_variations/sys_fid30_jets.txt \
    --comp vertex=sys_variations/sys_vz60_jets.txt --comp novertex=sys_variations/sys_vtxmeas_jets.txt \
    --out final_jets
root -b -q 'plot_an_final.C+("final_jets.txt","final_all_all","all","all","EM + JP union")'
```

Photon-class misidentification: `misid_response` on the JES simulation, then `misid_sys.py` on the nominal table (merged classes, Sean 2026-11-17).
