# A_N systematics

Tools for the A_N systematic checks (started 2026-11-17). The plan, item by item,
is in FCSJetAnalysis/ANALYSIS_ROADMAP.md ("Systematics plan").

## Method

**Every systematic is symmetric per bin (Sean 2026-11-18).** A component with
two sides (up / down, or several variations) takes the larger |ΔA_N| of its
sides and applies it as ±; components then add in quadrature. Nothing in the
chain carries an asymmetric box (`combine_sys.py`, `plot_an_final.C`).

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
| `sys_slices.py` | Δ, σ_diff, w, the slice's A_N, Barlow sys and pull per bin from two AnAsymmetry `.txt` tables (`--disjoint` for independent samples, `--shifted` for the same events with other jet kinematics, the JES: sys = \|Δ\|). Writes `<prefix>.txt` and `<prefix>_summary.txt` (χ², p-value per row). |
| `plot_sys_slices.C` | One row: A_N of the tight, loose and slice samples on the common A_N axes, and Δ ± σ_diff with the Barlow band. |
| `time_in_fill.py` | Hours in fill per run: t_mid from `../run22_polarization.txt` minus the fill start from the polarimetry fill table. Splits a yield list into bins: `--edges` in hours (default 2, as in FMS) or `--quantiles N`. |
| `sys_an_job.sh`, `submit_sys_an.sh` | AnAsymmetry -r over a list of yield files as condor jobs. This is the per-file flow; the merged flow is faster for whole passes. |
| `pi0_iso_sys.sh` | π⁰ isolation-threshold systematic (iso 0.95 / 0.98 / 0.99, noniso 0.85 / 0.9 / 0.95) on a report folder. |
| `sys_variations.sh` | The cut variations stored as yield copies (AnConfig.h `kVarSuffix`, `AnAsymmetry -v`): fiducial 25 / 30 cm, \|vz\| < 60, measured vertex only (nested, Barlow), the JES ±δ(E) and point-threshold 1.5 / 2.0 GeV copies (shifted: \|A_N shift\|, larger side via `combine_sys.py --comp-max`), and the collision-cut failures (disjoint check, plus their share per row in `bbcveto_share.txt`); photon classes (`-J`) and iso / noniso π⁰. On a report folder; outputs in `<report>/sys_variations/`. |
| `combine_sys.py` | The combine step: nominal table + components (`--comp NAME=sys_slices output`, `--comp-max` for the larger of several variations, `--sysf` for the π⁰ signal-fraction column, `--rel` for flat relative pieces) → per-bin stat, each component, quadrature total, stat ⊕ sys and the x_F uncertainty (`--xf-rel`, FMS 3%). The polarization scale is quoted separately. |
| `plot_an_final.C` | The final plots of one row from a combined table: `_stat`, `_sys` (total boxes), `_sys_components`, `_combined` (stat bars + sys boxes), `_total` (stat ⊕ sys bars), on the common A_N axes. |
| `misid_response.cpp` | Photon-class response from simulation (JES production `jet_output_*`, cross-section weighted): reco class × true class per truth-x_F bin, true photons with E >= 1 GeV (`--truth-emin`), offline trigger-value cut (the production's emulated flags are 0); `<out>.txt` per class and `<out>_merged.txt` for np12 / np3plus. |
| `misid_sys.py` | Misidentification systematic for the merged classes: dA = f(other \| reco) × (A_other − A), Barlow-reduced by the difference's statistical error; `combine_sys.py --comp-simple misid=<out>.txt`. |
| `jes_data_subsets.sh` | The pass's yield files in k random groups (default 2, 4, 8, plus the full sample): per group the union A_N nominal, `-v jesup` and `-v jesdn` (`AnAsymmetry -r`) → `groups_jesup.txt`, `groups_jesdn.txt`. |
| `jes_noise_scaling.py` | `groups`: Δ_g² = δ² + k s² per bin (per-k means, variance-weighted): the JES effect δ and the full-sample noise s of the `_vjesup` / `_vjesdn` shift, also pooled per row in stat units. Optional check only: the JES systematic is the raw \|Δ\| (Sean 2026-11-18: this noise is taken as negligible). |
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
    --comp-max jes=sys_variations/sys_jesup_jets.txt,sys_variations/sys_jesdn_jets.txt \
    --out final_jets
root -b -q 'plot_an_final.C+("final_jets.txt","final_all_all","all","all","EM + JP union")'
```

Photon-class misidentification: `misid_response` on the JES simulation, then `misid_sys.py` on the nominal table (merged classes, Sean 2026-11-17).

## JES uncertainty: the standard map and its error (Sean 2026-11-18)

**The A_N systematic.** The calibration is repeated with the map shifted by its
own error up and down, and the bin-by-bin A_N difference is the systematic,
symmetric: sys = max(|A_up − A|, |A_dn − A|). In the pass this is one job:
`AnHistograms --jes-error=calib.json` fills `_vjesup` / `_vjesdn`, every jet's E,
pT and x_F × (1 ± δ), δ at the map's own lookup energy (UE-subtracted reco E),
same events and raw trigger values. That equals rerunning
`jet_calibration/ApplyCorrections_em_only <list> <map> --jes-shift=+1` (or `-1`,
F → F / (1 ± δ)) and a nominal AnHistograms: checked on one 23087036 file,
540 histograms per side, bin for bin identical (2026-11-18). Then
`sys_variations.sh` (`sys_jesup_jets.txt`, `sys_jesdn_jets.txt`, `--shifted`:
|Δ|) and `combine_sys.py --comp-max jes=<up>,<dn>`. Photon classes only.

One standard map; its error vs energy is part of the map (`jes_error_*` arrays)
and AnHistograms applies it as a coherent scale ±δ(E) (`_vjesup` / `_vjesdn`).

1. **Statistical (the standard map's own noise).** `jet_scale/jes_subset_maps.sh
   -i <flat jet dir of the standard sample> -o <dir>`: maps on disjoint 1/2,
   1/4, 1/8 subsets (stratified per ptHat bin, same settings; `-N` folders only).
   `jet_scale/jes_map_compare.py noise <dir>/manifest.txt`: σ₁ per node and
   energy (flat in k if the noise scales as 1/√n; the core share must stay put).
2. **Systematic (generator / tune).** A smaller sample of each alternative, mapped
   at the standard settings. `jes_map_compare.py error --base <map> --sigma
   <noise>_sigma.csv --nbase N --alt NAME=<map>:N ... --write-map <map with
   error>`: per energy, the difference smoothed over the grid (Gaussian, 5 bins)
   less the noise of both maps after the same smoothing, measured from the
   smoothed subset maps (`noise` writes it; the alternative's scaled by its
   smaller sample): sys² = <d_sm²> − <noise_sm²> over the grid; RMS (or max) over
   the alternatives, total = stat ⊕ sys. The unsmoothed per-node version is
   printed as a cross-check (`_sys_node`). `plot_jes_error.C` draws it. The written map
   replaces the standard copy (`jet_calibration/JES_Lookups/`).
   **Plots (ROOT):** `jet_scale/plot_jes_noise.C("<noise>.txt", <full sample, M
   events>, "<out>", emin, emax)`: the standard map's statistical error vs energy
   (band: spatial std. dev.) and vs sample size (subset maps measured at N/k, the
   1/√N curve, the full map extrapolated). `jet_scale/plot_jes_error.C("<err>.txt",
   "<out>")`: the systematic per alternative vs energy (band: spatial std. dev.)
   and stat / sys / total. The per-node subset error is corrected for the few maps
   per division (mean / c₄(k): 0.80, 0.92, 0.97 for 2, 4, 8 maps) and the band for
   that estimate's own scatter.
3. **Check of the projection (Sean):** `jes_map_compare.py test old.json new.json
   --sigma ...`: two maps from independent samples at the same settings; ratio
   of observed to predicted RMS difference (1 = all noise).
4. **Data-side noise of the shifted A_N** (jets moving between x_F bins): taken
   as negligible (Sean 2026-11-18); the A_N systematic is the raw |Δ|, not
   noise-subtracted. `jes_data_subsets.sh` + `jes_noise_scaling.py groups`
   remain as an optional check.

Tests (2026-11-17/18, synthetic maps and tables with known noise): σ₁ flat in k;
pure noise gives ratio 0.99, RMS pull 1.00; an injected energy-dependent offset
(2.9%, 1.9%, 1.1%, 0.26% at 3, 10, 30, 100 GeV) is recovered with the smoothed
estimate (2.9%, 1.9%, 1.1%, 0.19%; per node 2.8%, 1.8%, 1.0%, 0); zero error reproduces nominal exactly in `_vjesup` / `_vjesdn`; the
groups fit is unbiased over 300 trials. δ is resolved only if it is not much
below the noise.
