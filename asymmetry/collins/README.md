# asymmetry/collins

Collins-type asymmetry of the leading pi0 candidate in FCS EM jets: the
azimuthal modulation sin(phi_S - phi_H) of the pi0 around the jet axis. The
yields are filled in the same pass as the A_N (Sivers-type) yields
(../AnHistograms.cpp calls CollinsHistograms.h), into the same yield files.

## Per-jet inputs (jet files; shared/Pi0Prob.h, Pi0LeadPair)

- leading pi0 candidate pair: the leading ECal point and the partner whose pair
  mass is closest to the reconstructed pi0 peak at that pair energy
  (reco_pi0_lead_i, reco_pi0_lead_partner, reco_pi0_lead_mgg)
- frame: j = jet direction, y' = (z x j)/|z x j| (normal of the plane of the
  blue beam +z and the jet), x' = y' x j
- reco_pi0_lead_phih: phi_H = atan2(p.y', p.x') of the pair momentum p
- reco_pi0_lead_jt: |p - (p.j) j|
- reco_pi0_lead_zpair: pair energy / raw jet energy
- reco_pi0_lead_gup: g = sin(phi_S - phi_H), phi_S the azimuth of spin up (+y)
  in the same frame (+pi/2 for jets at x > 0, -pi/2 at x < 0)
- reco_pi0_lead_pair: P_leadpair, the probability that the leading point is in
  a resolved pi0 pair (the selection; P_lead also counts merged and
  partner-lost pi0, which have no correct pair)

## Yields (CollinsHistograms.h, binning in CollinsConfig.h)

Jets with E <= sqrt(s)/2 at the default trigger factor (no suffix) or JP not
EM at it (_jponly, for the EM + JP union):

    h_xFg_<level>_<side>_<spin state>_<lp0|lp1|lp2>[_jponly]   TH3D x_F x g x z
    h_absg_<level>_<lp>[_jponly], h_xfsumg_<level>_<lp>[_jponly]  TH2D x_F x z, sums of |g| and x_F

x_F is the jet's (0-0.6, 0.05 bins), z = reco_pi0_lead_zpair (0, 0.3, 0.5,
0.7, 1), lp0/1/2 = P_leadpair < 0.5, 0.5-0.8, >= 0.8.

## Extraction (to write)

The FCS covers phi_S only near +-pi/2 (the two sides), but phi_H covers 0-2pi
(the pi0 is near the jet axis, <j_T> ~ 0.2 GeV), so g samples [-1, 1]. As the
A_N left-right method with <|cos phi|>: within each side, the cross ratio of
the g > 0 and g < 0 jets over the spin states (luminosity, acceptance and the
inclusive A_N cancel), divided by P <|g|> (h_absg / yield); sides averaged.
g in 20 bins also allows a binned fit in phi_C as a check.

## Simulation (Pi0Playground/pi0_prob)

pi0_prob_perf page 3: purity / efficiency of "chosen pair = the leading
photon's true pi0 pair" vs a P_leadpair cut (pi0prob_v2: ~0.80 purity at ~25%
efficiency for P_leadpair >= 0.8; two-point jets 0.58). Still to do: the phi_H
resolution (reconstructed pair vs true pi0).
