#!/bin/bash
# Redraws every A_N plot of a report, and the yield and sanity pages in
# plots/ (AnPlots -s, from plots/qa_sums.root), from what the report folder
# already holds, in about a minute, with the current code of this directory:
# for cosmetic changes (axes, text, colours, a new plot of existing graphs)
# without rerunning AnAsymmetry over the yields.
#
# Usage: ./replot_report.sh <report dir>
#   Reads, in <report dir>: an.root, an.txt, an_union.root, an_mass_in.root,
#   an_mass_out.root (if there), projection.log (the error scale and the
#   sample line: made by project_errors.py) and mgg_peak.log (the pi0 peak).
#   Writes the plots (.png and .pdf) and plot_an.log there; the macros used
#   are copied to <report dir>/build first.
#   AN_YLO, AN_YHI, AN_XR (environment): the common axes of every A_N plot,
#   default -0.02, 0.08, 0.6 (Sean 2026-11-06: all A_N plots on one scale).
#   AN_ERRMAX (environment): bins with a larger stat. error are not drawn
#   (default 0.022; 0 = all).
# smoke_report.sh calls this after AnAsymmetry and project_errors.py, so the
# list of plots lives here only.
# Needs the STAR SL7 container (run from the host). SMOKE_BIND=<dir> binds a
# further directory into it (a condor job's scratch disk).

set -e
REP="$(cd "$1" && pwd)"
[ -f "${REP}/an.root" ] && [ -f "${REP}/an.txt" ] && [ -f "${REP}/projection.log" ] || { sed -n '2,19p' "$0"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# a file of this script's directory, else of the one above (repo layout:
# asymmetry/sivers/ with the pass-level files in asymmetry/; a frozen report
# code directory has them all side by side)
src() { if [ -e "${HERE}/$1" ]; then echo "${HERE}/$1"; else echo "${HERE}/../$1"; fi; }
mkdir -p "${REP}/build"
cp "$(src plot_an.C)" "$(src plot_an_means.C)" "$(src plot_an_mass.C)" "${REP}/build/"
# the FMS results to overlay (fms_reference/README.md), if this code has them
FMSREF="$(src fms_reference)"
[ -d "${FMSREF}" ] && { mkdir -p "${REP}/build/fms_reference"; cp "${FMSREF}"/*.txt "${REP}/build/fms_reference/"; }

container() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 ${SMOKE_BIND:+-B ${SMOKE_BIND}} \
        /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd ${REP}; $1"
}

SCALE=$(grep '^scale=' "${REP}/projection.log" | cut -d= -f2)
PEAK=$(sed -n 's/^peak=//p' "${REP}/mgg_peak.log" 2>/dev/null | tail -1)
SIGMA=$(sed -n 's/^sigma=//p' "${REP}/mgg_peak.log" 2>/dev/null | tail -1)
# what the plots show: "<test runs> of <good runs> good runs, <events> events = <percent> of the Run 22 good-run events"
SAMPLE=$(awk '/^test:/ { gsub(/[;,]/, ""); printf "%s of %s good runs, %.2g events = %.1f%% of the Run 22 good-run events", $2, $9, $6, 100 * $15 }' "${REP}/projection.log")
AN_YLO="${AN_YLO:--0.02}"; AN_YHI="${AN_YHI:-0.08}"; AN_XR="${AN_XR:-0.6}"
AX="${AN_YLO},${AN_YHI},${AN_XR}"
# bins with a statistical error above AN_ERRMAX are not drawn: default 0.022 =
# the largest stat. error among the published FMS 500 GeV A_N points (PRD 103,
# 092009; Sean 2026-11-11: drop bins with too few jets); 0 = draw all
AN_ERRMAX="${AN_ERRMAX:-0.022}"

{
    p() { echo "root -b -q -l 'build/plot_an.C(\"$1\",\"$2\",\"$3\",\"$4\",${5:-1},${6:-0},10,\"${SAMPLE}\",${AX},\"${7:-}\",${AN_ERRMAX},${8:-0})' || exit 1"; }
    # EM levels together (all), the union, the JP levels together (allJP)
    p an.root an_vs_xf.png all all
    p an_union.root an_vs_xf_union.png all all
    p an.root an_vs_xf_allJP.png allJP all
    # errors projected to the full good-run sample
    p an.root an_vs_xf_projected.png all all "${SCALE}" 0
    p an.root an_vs_xf_projected_zero.png all all "${SCALE}" 1
    # single trigger levels
    for l in EM0 EM1 EM2 JP0 JP1 JP2; do p an.root "an_vs_xf_${l}.png" "$l" all; done
    # photon classes: EM levels, JP levels, union
    for c in np1 np2 np12 np3plus; do
        p an.root "an_vs_xf_${c}.png" all "$c"
        p an.root "an_vs_xf_allJP_${c}.png" allJP "$c"
        p an_union.root "an_vs_xf_union_${c}.png" all "$c"
    done
    # the union with the FMS results overlaid (Sean 2026-11-11): all jets,
    # 1 or 2 photons, >= 3 photons; *_fms.png beside the plain plots
    if [ -d "${REP}/build/fms_reference" ]; then
        p an_union.root an_vs_xf_union_fms.png all all 1 0 build/fms_reference/fms_all.txt
        p an_union.root an_vs_xf_union_np12_fms.png all np12 1 0 build/fms_reference/fms_np12.txt
        p an_union.root an_vs_xf_union_np3plus_fms.png all np3plus 1 0 build/fms_reference/fms_np3plus.txt
        # two-photon jets in the pi0 mass window vs the FMS isolated pi0
        # (z_em > 0.98; Sean 2026-11-11, until the z_pair binning of the next pass)
        [ -f "${REP}/an_mass_in.root" ] && p an_mass_in.root an_vs_xf_mass_in_fms.png all np2 1 0 build/fms_reference/fms_iso_pi0.txt
    fi
    echo "root -b -q -l 'build/plot_an_means.C(\"an.txt\",\"an_means.png\",\"all\",\"${SAMPLE}\",${AN_XR})' || exit 1"
    if [ -f "${REP}/an_mass_in.root" ] && [ -f "${REP}/an_mass_out.root" ] && [ -n "${PEAK}" ]; then
        echo "root -b -q -l 'build/plot_an_mass.C(\"an_mass_in.root\",\"an_mass_out.root\",\"an_vs_xf_mass.png\",${PEAK},\"${SAMPLE}\",${AX},${AN_ERRMAX},${SIGMA:-0})' || exit 1"
    fi
} > "${REP}/plot_an.csh"
container "tcsh plot_an.csh" > "${REP}/plot_an.log" 2>&1 || { tail "${REP}/plot_an.log"; exit 1; }
# the yield and sanity pages (plots/), from the sums AnPlots saved when it read the yields
if [ -f "${REP}/plots/qa_sums.root" ]; then
    cp "$(src AnPlots.cpp)" "$(src AnConfig.h)" "${REP}/build/"
    container "g++ -m64 -O2 build/AnPlots.cpp -o build/AnPlots \`root-config --cflags --libs\` -std=c++11 && ./build/AnPlots -s plots/qa_sums.root plots" > "${REP}/anplots_replot.log" 2>&1 \
        || { tail "${REP}/anplots_replot.log"; exit 1; }
    echo "plots/ redrawn from plots/qa_sums.root"
else
    echo "no plots/qa_sums.root (report from before 2026-11-06): plots/ left as is"
fi
echo "$(grep -c 'has been created' "${REP}/plot_an.log") files written in ${REP} (axes A_N ${AN_YLO} to ${AN_YHI}, x_F +-${AN_XR})"
