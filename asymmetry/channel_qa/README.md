# Hot / dead FCS channel masks

Plan (Sean 2026-11-24): divide Run 22 into periods, take a sample of events in
each, find the hot and dead channels of each period, and run the full A_N pass
with those masks. If the masks are large enough to matter for the JES map, build
one more JES map with the masked channels and use it as a systematic (like the
point-threshold one).

## Steps

1. **Periods and sample**: `make_mask_sample.py <good_runs.txt> <prefix>`
   (`--periods 25 --runs 6 --files 4`). Periods are chunks of consecutive good
   runs with about equal events, also split at the ECal gain changes (23027047,
   23066056). Writes `<prefix>_periods.txt`, `<prefix>_runs.txt` and the MuDst
   list `<prefix>.list`.
2. **Channel statistics**: `../submit_an_job.sh -f <prefix>.list -o <dir> -n 4 -Q 1`.
   With `-Q 1` the A_N job runs `channel_occupancy` on every SimpleTree before
   deleting it (`<dir>/occ/occ_<run>_<file>.root`): per channel of ECal N/S and
   HCal N/S, hits with E > 0, hits above the reclustering tower threshold (ECal
   0.4 GeV, HCal 0.1 GeV), their summed E and E², the channel position, and the
   events read. Same EM / JP trigger filter as the A_N sample.
3. **Masks**: `channel_mask <periods.txt> <runs.txt> <dir>/occ <outdir>`. Per
   period, each channel's hit rate and mean hit energy over the median of its
   neighbours (within 1.5 tower pitches; at the corners, its ring). Flags: hot
   (rate > 3× and > 5σ), dead (no hits, or rate < 0.2×), hotE (mean E > 2×),
   coldE (mean E < 0.5×). First thresholds, to be tuned on the maps.
   `channel_mask.txt` (period first_run last_run det id x y flag ...),
   `channel_mask_summary.txt`, `channel_ratio_<det>_p<k>.png`.
4. **Applying a mask**: automatic (Sean 2026-11-25). `../submit_an_job.sh`
   applies `channel_mask_v1.txt` (this folder) unless given `-M <other mask>`,
   or `-M none` for no mask. The mask's `#P period first last` header lines list
   every period, masked channels or not. The submit script warns about runs of
   the list outside every period, since those run unmasked. A new mask version
   replaces the default in `submit_an_job.sh`. `runMudst.C`
   loads `channel_mask.txt` into `StFcsMuToEventMaker` (masked hits get zero
   energy before reclustering, so they enter no cluster or point) and
   `StSimpleReaderMaker` (masked hits are left out of the hit list and the
   offline trigger sums), for the runs of each row's period
   (`StRoot/StSpinPool/StSimpleReaderMaker/FcsChannelMask.h`). Which flags go
   into the applied mask is chosen by filtering `channel_mask.txt`.
