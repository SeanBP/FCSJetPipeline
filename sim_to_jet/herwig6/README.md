# herwig6

Herwig6 generator scripts, preserved for future use but **not wired into
the live pipeline** (not referenced by `sim_to_jet.xml` or
`runSimBfc.C`'s `{{GENERATOR}}` switch).

**Known issue**: the last test run (see `BUG_SUMMARY.txt`) found a
reproducible forward-plane azimuthal asymmetry -- events coherently favor
one lab-frame quadrant (ratio ~2.75:1). Root cause unresolved as of
2026-07-31. **Do not trust this generator's output until that's fixed.**

## Files

- `starsim_herwig6_filter.C`, `starsim_herwig6_default.xml`, `runHerwig6.C`
  -- the generator scripts.
- `BUG_SUMMARY.txt` -- writeup of the asymmetry bug: symptom, what was
  ruled out, and where to look next.

`runSimBfc.C`'s `proc.Contains("herwig")` branch (which recognizes a
herwig `.fzd` input file) is unaffected by any of this.
