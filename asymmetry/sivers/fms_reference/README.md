# FMS EM-jet A_N reference points

STAR FMS EM-jet A_N vs x_F (x_F > 0), overlaid on the Run 22 FCS union plots by
`replot_report.sh` (`an_vs_xf_union{,_np12,_np3plus}_fms.png`, via the last
argument of `plot_an.C`). Format: `set <legend label>` starts a set, then
`x_F A_N stat` per line; `#` starts a comment. Statistical errors only; neither
source tabulates systematics here.

| file | overlaid on | sets |
|---|---|---|
| `fms_all.txt` | all jets | Run 11, 500 GeV, published |
| `fms_np12.txt` | 1 or 2 photons | Run 17, 510 GeV, preliminary, N_gamma = 1, 2 |
| `fms_np3plus.txt` | >= 3 photons | Run 17, 510 GeV, preliminary, N_gamma >= 3; Run 11, 500 GeV, published, N_gamma > 2 |
| `fms_run17_multiplicity.txt` | (not overlaid) | Run 17, 510 GeV, preliminary, N_gamma = 1, 2 / 3 / 4, 5 |

Sources:
- Published: STAR, PRD 103, 092009 (2021), arXiv:2012.11428, Fig. 9 (Run 11
  p+p 500 GeV EM-jets; all, and photon multiplicity > 2). Values from HEPData
  record ins1837609, Tables 2 and 4 from Fig 9 (points at the mean x_F). The
  paper has no Run 11 1-2 photon class vs x_F.
- Preliminary: W. Zhang (STAR), PoS(DIS2025)110, Fig. 1 right (N_gamma = 1, 2 /
  3 / 4, 5) and Fig. 3 left (red: Run 17 510 GeV, N_gamma = 1, 2 and >= 3). Not
  tabulated. Read 2026-11-11 from the PDF's vector paths, not from pixels:
  - each axis is calibrated on its tick marks;
  - each marker centre comes from its path, and the error from the error-bar caps;
  - the N_gamma = 1, 2 points of Fig. 1 and Fig. 3 agree to 1e-5;
  - points sit at the bin centres;
  - where ROOT hid the bar inside the marker, the error is below the marker
    radius (0.0006), and that value is stored with a comment.

Differences from our points to keep in mind: the FMS uses 2.8 < eta < 3.8, a
vertex cut, the off-axis-cone UE correction and (Run 17) energy unfolded to
the particle level; the 2017 polarization scale uncertainty is 1.1%.
