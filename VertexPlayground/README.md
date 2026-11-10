# VertexPlayground

Run 22 st_fwd primary-vertex study (2026-11-09). The production SimpleTrees'
`event_vertex_z` came from `StBbcTriggerDetector::zVertex()`, whose hard-coded
pre-2003 offsets make it junk for Run 22 (uncorrelated with BBC TAC or EPD).
StSimpleReaderMaker now also writes uncut vertex branches: `bbc_tac_diff`,
`bbc_vz_tac` ((TAC diff - 4096) * -0.2475 cm), `bbc_vz_det`, `epd_vz` (EPD
average-TAC vertex, after D. Kapukchyan's StFwdAnaEpdQaAndVert), `vpd_vz`.

`vertex_test_v1/` — the condor test that read those branches from real data:

| file | what |
|---|---|
| `job.sh`, `job.sub`, `job2.sh`, `job2.sub`, `inner.csh` | per-run job: xrdcp each MuDst, `runMudst.C` with the production trigger filter, `slim.C` keeps only the vertex branches (`vtx_<run>_<n>.root`) |
| `slim.C` | vertex-branch slimmer |
| `analyze.C`, `extra.C` | first look: vertex availability, z distributions, BBC vs EPD correlation |
| `perrun.C`, `lowruns.C`, `r106.C` | per-run availability checks and low-vertex runs |
| `vertex_summary.C` | 155-run summary plots (per-run fractions, mean and RMS of z) |

Inputs and outputs live in `/star/data01/pwg/seanp/vertex_test_v1` (paths in
the scripts are absolute to it).

Results (12 benchmark runs, then 155 runs / 528k EM+JP events):
- BBC TAC vertex present in 98% of events (79% with |z| < 150 cm), RMS ~112 cm,
  mean ~+12 cm; BBC-EPD core sigma ~81 cm, slope 0.59: without a slewing
  correction too poor for kinematics.
- EPD vertex: 0% for runs up to ~23039034 (no hits pass 1 < nMIP < 15,
  TAC > 50), ~95% from ~23050004 on; RMS ~56 cm, mean ~-12 cm (25 cm offset to
  BBC).
- Low-vertex runs: 23055022, 23082017 (already suppressed), 23106010.
- Effect at the FCS (z = 721.5 cm): vz = +-60 cm moves eta by ~0.08 and pT by
  ~8%, x_F by < 0.3%.

The simulation's vertex width is now configurable (`submit_sim_to_jet.sh -Z`,
default 30 cm); the pi0prob_v2 production used 60 cm, closer to data.
