#!/bin/bash
# Submits pi0_pairs.xml: one pi0_pairs job per "<list>:<data|sim>" pair.
# A simulation list (SimpleTrees) must hold every file of each ptHat bin in it (weights are
# sigma / events per bin within the job). Job i writes <workdir>/out/pi0_pairs_<i>.root.
# Usage: ./submit_pi0_pairs.sh -w <workdir> "<list>:<mode>" ["<list>:<mode>" ...]
#   <mode>: data or sim, optionally "+all" (no EM trigger requirement), e.g. "list:sim+all"
set -e
while getopts "w:h" opt; do
  case $opt in w) WORK="$OPTARG" ;; *) sed -n '2,6p' "$0"; exit 0 ;; esac
done
shift $((OPTIND - 1))
[ -n "$WORK" ] && [ $# -gt 0 ] || { sed -n '2,6p' "$0"; exit 1; }
mkdir -p "$WORK/out" "$WORK/log"; WORK=$(readlink -f "$WORK")
n=0
for a in "$@"; do
  cp "${a%%:*}" "$WORK/list_${n}.list"; echo "${a##*:}" | tr '+' ' ' > "$WORK/mode_${n}.txt"
  echo "job $n: ${a##*:} $(wc -l < "${a%%:*}") files (${a%%:*})"; n=$((n + 1))
done
cd "$(dirname "$(readlink -f "$0")")/.."   # the FCSJetPipeline root
sed -e "s#{{WORKDIR}}#$WORK#g" -e "s#{{NPROC}}#$n#" Pi0Playground/pi0_pairs.xml > "$WORK/pi0_pairs.xml"
star-submit "$WORK/pi0_pairs.xml"
