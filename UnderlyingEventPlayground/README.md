# UnderlyingEventPlayground

Studies of the underlying-event (UE) subtraction implemented in
`../underlying_event/` (`ProxySelection.h`, `ProxySearch.h`,
`UnderlyingEvent.cpp`); nothing here is used by the pipeline. Each source has
its build and run line at the top.

## Method

The near off-axis proxy-event method. The UE under a jet is the mean ECal
point energy (same point selection as the jets: ECAL, 1 GeV to sqrt(s)/2) in
an R cone at the jet's own position, averaged over other events of the same
group (run in data, ptHat bin in simulation) and highest EM trigger level,
whose own reco jet:

- has pT within 5% of the jet's (widened in 5% steps up to 50% only if fewer
  than 4 proxy events are found);
- sits near 90 deg in azimuth from the cone, in a band whose half-width gives
  every cone position the same allowed eta-phi area;
- is at least 2.5R from the cone, and no other jet of that event is closer.

Every event that passes is used. The same is done at truth level in
simulation, with truth jets and EM-only truth photons. The subtraction is in
energy, per jet, at both levels, before the JES calibration.

## Studies

- `fms_emjet_reference.txt`: the STAR FMS EM-jet analysis (BrightSTAR
  `emjet-gpc`) and how this pipeline compares.
- `truth_dphi.cpp`, `recoil_rate.cpp`: how often a second (truth) jet is
  back-to-back with the first, i.e. why the cone should avoid the recoil
  direction.
- `proxy_region_area.cpp`: the band half-width that gives every cone position
  the target allowed area.
- `proxy_pt_window.cpp`: cost and UE bias of the proxy pT window width.
- `proxy_trigger.cpp`: cost of matching the trigger, and whether the UE
  depends on it.
- `proxy_ue.cpp`: the full search jet by jet (cost, proxies found, reach in
  events), in simulation or data.
- `ue_position_map.cpp`: mean UE against jet position, in pT ranges.
- `ue_pthat_dependence.cpp`: effect of drawing simulation proxies from the
  jet's own ptHat bin (below).

## Results

**ptHat grouping in simulation** (`ue_pthat_dependence.cpp`, 2026-10, first
calibration sample `fms_tune_calib_v1`, 966 files). Simulation proxies come
from the jet's own ptHat bin, data proxies from the natural event mix. The
cross-section-weighted RMS shift of the UE estimate caused by per-bin pooling,
as a percentage of the mean jet energy, at fixed jet pT and radial position
(inner / middle / outer thirds of the acceptance):

| Jet pT (GeV) | Reco | Truth |
|---|---|---|
| 0.25-0.5 | 12.7 / 7.7 / 3.4 | 5.0 / 3.7 / 2.2 |
| 0.5-1 | 4.0 / 2.8 / 2.0 | 2.4 / 2.1 / 1.4 |
| 1-2 | 2.3 / 1.9 / 1.4 | 1.1 / 1.2 / 1.2 |
| 2-3 | 0.8 / 0.8 / 0.8 | 0.8 / 0.7 / 0.8 |
| 3-4 | 0.4 / 0.4 / 0.6 | 0.3 / 0.5 / 0.6 |
| 4-6 | 0.3 / 0.4 / 0.4 | 0.1 / 0.2 / 0.3 |
| 6-10 | 0.4 / 0.4 / 0.5 | 0.3 / 0.3 / 0.5 |

Below 1% for analysis jets (pT >= 2 GeV), so the per-bin grouping is kept.
Reco and truth shift the same way (0.15-0.3 GeV above 2 GeV pT), so the effect
on the JES ratio is smaller still. 60-77% of the weight is in the 5-10 GeV
ptHat bin. At fixed jet pT and position the mean UE falls with ptHat (above
ptHat 20 GeV it is about half the 0-10 GeV value), at both levels; not yet
understood.
