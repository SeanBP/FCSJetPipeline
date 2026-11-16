#!/bin/bash
# One merge job of submit_an_merge.sh: copies chunk <process + offset>'s yield files to the
# scratch disk (reading ~11k histograms per file straight from /star/data01 is
# slow), runs AnMerge -b and writes <outdir>/merge/partial_<i>.root.
# Arguments: an_merge_job.sh <outdir> <code dir> <process> [chunk offset = 0]
set -e
OUTDIR="$1"; CODE="$2"; I=$(( $3 + ${4:-0} ))
IDX=$(printf "%04d" "$I")
LIST="$OUTDIR/merge/chunks/chunk_${IDX}.list"
SCR="${_CONDOR_SCRATCH_DIR:-${SCRATCH:-$(pwd)}}/merge_${IDX}"
mkdir -p "$SCR/hists"
t0=$(date +%s)
echo "start $(date) on $(hostname): chunk ${IDX}, $(wc -l < "$LIST") files"
xargs -a "$LIST" cp -p -t "$SCR/hists/"
ls "$SCR"/hists/*.root > "$SCR/yields.list"
echo "copied $(wc -l < "$SCR/yields.list") files in $(( $(date +%s) - t0 )) s"
cp "$CODE"/AnMerge.cpp "$CODE"/AnRunTools.h "$CODE"/AnConfig.h "$CODE"/run22_polarization.txt "$SCR/"
ct() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B "$SCR" /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd $SCR; $1"
}
ct "g++ -m64 -O2 AnMerge.cpp -o AnMerge \`root-config --cflags --libs\` -std=c++11"
ct "./AnMerge -b yields.list run22_polarization.txt partial.root"
# write under a temporary name, then rename, so the report never reads a partial copy
cp "$SCR/partial.root" "$OUTDIR/merge/.partial_${IDX}.root.tmp"
mv "$OUTDIR/merge/.partial_${IDX}.root.tmp" "$OUTDIR/merge/partial_${IDX}.root"
rm -rf "${SCR:?}"
echo "done $(date), $(( $(date +%s) - t0 )) s"
