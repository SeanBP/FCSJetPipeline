#!/bin/bash
# File catalog queries for the A_N resource estimate (P24ia st_fwd MuDsts of
# production_pp500_2022, as data_to_jet/MakeRunList.pl).
#
#   ./catalog_query.sh counts <out.txt> [local|hpss]
#       one line per run: run::events::files::bytes (default storage local)
#   ./catalog_query.sh list <out.list> <run> [<run> ...] [-l <files per run>]
#       MakeRunList.pl-format list ("root://.../<file>.MuDst.root/<events>")
#       of the given runs' files on disk, optionally only the first N of each
#
# Runs in the container (get_file_list.pl needs the STAR environment).

set -e
COND="production=P24ia,trgsetupname=production_pp500_2022,filetype=daq_reco_mudst,sname2=st_fwd"
XRD="root://xrdstar.rcf.bnl.gov:1095/"
star() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
        /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; $1"
}

case "$1" in
    counts)
        OUT="$2"; ST="${3:-local}"
        star "get_file_list.pl -keys 'grp(runnumber),sum(events),count(filename),sum(size)' -cond '${COND},storage=${ST}' -limit 0" > "$OUT"
        echo "$(wc -l < "$OUT") runs -> $OUT"
        ;;
    list)
        OUT="$2"; shift 2
        RUNS=(); NPER=0
        while [ $# -gt 0 ]; do
            if [ "$1" = "-l" ]; then NPER="$2"; shift 2; else RUNS+=("$1"); shift; fi
        done
        : > "$OUT"
        for r in "${RUNS[@]}"; do
            star "get_file_list.pl -keys 'path,filename,events' -cond '${COND},storage=local,runnumber=${r}' -limit 0" \
                | sort -t: -k3 | awk -F'::' -v x="$XRD" -v n="$NPER" 'n == 0 || NR <= n { print x $1 "/" $2 "/" $3 }' >> "$OUT"
        done
        echo "$(wc -l < "$OUT") files -> $OUT"
        ;;
    *)
        sed -n '2,14p' "$0"; exit 1 ;;
esac
