#!/bin/bash
# File catalog queries for the A_N resource estimate (P24ia st_fwd MuDsts of
# production_pp500_2022, as data_to_jet/MakeRunList.pl).
#
#   ./catalog_query.sh counts <out.txt> [local|hpss]
#       one line per run: run::events::files::bytes (default storage local)
#   ./catalog_query.sh list <out.list> <run> [<run> ...] [-l <files per run>]
#       MakeRunList.pl-format list ("root://.../<file>.MuDst.root/<events>")
#       of the given runs' files on disk, optionally only the first N of each,
#       leaving out files on the nodes in bad_nodes.txt
#   ./catalog_query.sh filter <in.list> <out.list>
#       an existing list without the files on the nodes in bad_nodes.txt (one
#       catalog query for the whole production, ~4 min)
#
# bad_nodes.txt (next to this script): distributed-disk nodes whose xrootd
# does not serve data, one per line (# comments). Their files have no other
# disk copy, so they are left out of the job lists until the node is fixed;
# then empty the file. 2026-11-04: qgp017 (86,167 of 954,589 copies, ~9%;
# reads hang with 0 bytes while locate/stat work).
#
# Runs in the container (get_file_list.pl needs the STAR environment).

set -e
COND="production=P24ia,trgsetupname=production_pp500_2022,filetype=daq_reco_mudst,sname2=st_fwd"
XRD="root://xrdstar.rcf.bnl.gov:1095/"
BADNODES="$(dirname "$(readlink -f "$0")")/bad_nodes.txt"
# the bad nodes as an awk regexp ("^$" matches no node when the file is empty)
bad_re() {
    local re
    re=$(grep -v '^#' "$BADNODES" 2>/dev/null | awk 'NF { printf "%s%s", (n++ ? "|" : ""), $1 }')
    echo "${re:-^\$}"
}
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
            star "get_file_list.pl -keys 'path,filename,events,node' -cond '${COND},storage=local,runnumber=${r}' -limit 0" \
                | awk -F'::' -v bad="$(bad_re)" '$4 !~ bad' | sort -t: -k3 \
                | awk -F'::' -v x="$XRD" -v n="$NPER" 'n == 0 || NR <= n { print x $1 "/" $2 "/" $3 }' >> "$OUT"
        done
        echo "$(wc -l < "$OUT") files -> $OUT (bad nodes left out: $(bad_re))"
        ;;
    filter)
        IN="$2"; OUT="$3"
        [ -f "$IN" ] && [ -n "$OUT" ] || { sed -n '2,24p' "$0"; exit 1; }
        MAP=$(mktemp)
        star "get_file_list.pl -keys 'filename,node' -cond '${COND},storage=local' -limit 0" \
            | awk -F'::' -v bad="$(bad_re)" '$2 ~ bad { print $1 }' > "$MAP"
        awk 'NR == FNR { b[$1] = 1; next } { n = split($0, p, "/"); if (!(p[n - 1] in b)) print }' "$MAP" "$IN" > "$OUT"
        echo "$(wc -l < "$IN") -> $(wc -l < "$OUT") files -> $OUT ($(wc -l < "$MAP") catalog files on bad nodes: $(bad_re))"
        rm -f "$MAP"
        ;;
    *)
        sed -n '2,24p' "$0"; exit 1 ;;
esac
