#!/bin/bash
# The separation test on FCS data (one job, in its scratch folder): rebuilds
# the EM-only point jets (with areas) of a list of SimpleTrees, as
# data_to_jet does, then runs the UE step (underlying_event/UnderlyingEvent)
# on copies of them three ways:
#   sep_veto  separation 2.5R and the jet veto (-DPROXY_JET_VETO=1)
#   sep       the pipeline method: separation, no veto, as the FMS test
#   nosep     no separation (-DPROXY_MIN_DR_R=0; the veto, which uses the same
#             distance, does nothing)
# Writes <variant>/jet_output_<suffix>.root.
#
# Usage (inside the job; RecoJets and the three UnderlyingEvent builds here):
#   ./run_sep_test.sh <SimpleTree list>
set -e
LIST="$1"
mkdir -p in jets
: > pairs.list
while read -r st; do
    [ -z "$st" ] && continue
    cp "$st" in/
    b=$(basename "$st"); suf=${b#SimpleTree_mudst_}; suf=${suf%.root}
    echo "in/$b" > one.list
    ./RecoJets one.list "" 1 1 | tail -1
    mv jet_output.root "jets/jet_output_${suf}.root"
    echo "jet_output_${suf}.root $PWD/in/$b" >> pairs.list
done < "$LIST"
echo "run_sep_test.sh: $(wc -l < pairs.list) files, jets made"
for v in sep_veto sep nosep; do
    mkdir -p "$v"; cp jets/*.root "$v/"
    (cd "$v" && sed "s#^#$PWD/#" ../pairs.list > pairs.list && "../UE_$v" pairs.list | tail -4 && rm pairs.list)
done
rm -rf in jets one.list pairs.list
