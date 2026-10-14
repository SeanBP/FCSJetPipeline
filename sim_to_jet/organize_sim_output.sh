#!/bin/bash
# Splits a finished sim_to_jet output folder into per-ptHat-bin subfolders
# (file names unchanged); log/ and run_params_*.txt stay in place:
#   <outdir>/jet_output/pt<N>/jet_output_<gen>_<tune>_pt<N>_<job>.root
#   <outdir>/SimpleTree/pt<N>/SimpleTree_<gen>_<tune>_pt<N>_<job>.root
# Never overwrites (mv -n), so it is safe to rerun. Run only after every job
# has finished.
#
# Usage: ./organize_sim_output.sh <outdir>

set -e
OUTDIR="$1"
if [ -z "$OUTDIR" ] || [ ! -d "$OUTDIR" ]; then
    echo "Usage: $0 <sim_to_jet output folder>"
    exit 1
fi

n_moved=0; n_skipped=0
for kind in jet_output SimpleTree; do
    while IFS= read -r f; do
        base=$(basename "$f")
        bin=$(echo "$base" | sed -n -e 's/^.*_\(pt[0-9]*\)_[0-9]*\.root$/\1/p')
        if [ -z "$bin" ]; then
            echo "skipping $base (no _pt<N>_<job>.root in its name)"
            n_skipped=$((n_skipped + 1))
            continue
        fi
        mkdir -p "$OUTDIR/$kind/$bin"
        mv -n "$f" "$OUTDIR/$kind/$bin/"
        n_moved=$((n_moved + 1))
    done < <(find "$OUTDIR" -maxdepth 1 -type f -name "${kind}_*.root")
done
echo "moved $n_moved files, skipped $n_skipped"
for kind in jet_output SimpleTree; do
    [ -d "$OUTDIR/$kind" ] && echo "$kind/: $(find "$OUTDIR/$kind" -mindepth 1 -maxdepth 1 -type d | wc -l) ptHat-bin folders, $(find "$OUTDIR/$kind" -type f | wc -l) files"
done
