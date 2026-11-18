#!/bin/bash
# Offline A_N steps on a finished an_job pass (step 6 and the checks), into
# <outdir>/report/:
#   an.root, an.txt            AnAsymmetry -r (bad runs, per-run polarization;
#                              "all" = EM levels, "allJP" = JP levels)
#   an_union.root, .txt        the same with -J: "all" = EM + JP union, each
#                              jet once
#   an_vs_xf_union.png         the union A_N vs x_F
#   an_vs_xf_<JP level>.png    A_N per JP level and allJP
#   an_pi0_<iso|noniso|incl>_<sig|sb>.root, .txt  AnAsymmetry -J -I: FMS-style
#                              pi0 candidates in the signal / sideband window
#   an_pi0_<iso|noniso|incl>.root, .txt, _mass.png  pi0_bgsub.C: their
#                              sideband-subtracted A_N vs pi0 x_F, mass fits
#   an_vs_xf_pi0_<class>[_fms].png  those A_N (FMS isolated / non-isolated overlaid)
#   pi0_peak_runs.png, .txt    the candidates' pi0 peak and width run by run
#                              (scale factor check; PeakShift / Wide / LowStat flags)
#   event_cuts_<var>.png, _all.png, .txt  each event-level cut variable (vertex z,
#                              vertex source, east BBC tiles, TOF multiplicity),
#                              N-1, with the cut (solid) and variations (dotted)
#   pi0_peak_runs_bad_runs.txt its flagged runs (D. Kapukchyan's format), masked in
#                              every A_N above with AnAsymmetry -Q (with his list)
#   pi0_peak_runs_scale_{global,per_period,per_run}.txt  pi0 scale tables for a rehist
#   an_vs_xf.png               final A_N vs x_F, statistical errors
#   an_projected.txt           the table with errors projected to the full
#                              good-run sample (project_errors.py)
#   an_vs_xf_projected.png     the points with projected errors, and
#   an_vs_xf_projected_zero.png the projected errors alone, at A_N = 0
#   an_vs_xf_<EM level>.png    A_N per EM level
#   an_vs_xf_np<N>.png         A_N for jets with exactly N = 1, 2 points,
#                              np12, 1 or 2, and np3plus, 3 or more (all A_N
#                              plots on the same axes: AN_YLO, AN_YHI, AN_XR)
#   an_vs_xf_allJP_np*.png,   the same photon classes for the JP levels alone and
#   an_vs_xf_union_np*.png     for the EM + JP union
#   an_means.png               <|cos phi|> and <P> per x_F bin and photon class
#   plots/                     yields and sanity plots, per-run trends (AnPlots)
#   resources.txt              wall time, memory, scratch per stage and job,
#                              and the full-pass extrapolation
#                              (estimate_resources.py on the timing tables)
#   condor_jobs.txt            condor's own record per job: wall clock,
#                              memory, disk, exit code (condor_history)
#
# Usage: ./smoke_report.sh <outdir> <good_runs.txt> [report name = report]
#   <good_runs.txt>: make_smoke_list.py's <prefix>_good_runs.txt, the
#   full sample the errors and resources are projected to
#   [report name]: subdirectory of <outdir>, e.g. report_2026-11-10, so
#   reports on a pass still running (whatever jobs have finished; yield files
#   written in the last 10 minutes are left out, in case condor is still
#   copying them) do not overwrite each other
# Compiles in the container; run it from the host (condor_history).
# SMOKE_PLOTS_ONLY=1: only the error projection and the A_N plots, from the
# AnAsymmetry outputs already in the report folder. For the plots alone (no
# projection, no pass folder needed): replot_report.sh <report dir>.
# AN_XF_EDGES=<e0,e1,...>: variable x_F bin edges for every A_N (AnAsymmetry -x).
# SMOKE_BIND=<dir>: also bind <dir> into the container (an outdir on a
# condor job's scratch disk; see report_job.sh).

set -e
OUTDIR="$(cd "$1" && pwd)"
GOOD="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
[ -d "${OUTDIR}/hists" ] && [ -f "$GOOD" ] || { sed -n '2,27p' "$0"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# a file of this script's directory, else of the one above (repo layout:
# asymmetry/sivers/ with the pass-level files in asymmetry/; a frozen report
# code directory has them all side by side)
src() { if [ -e "${HERE}/$1" ]; then echo "${HERE}/$1"; else echo "${HERE}/../$1"; fi; }
REP="${OUTDIR}/${3:-report}"
mkdir -p "${REP}/build"
cp "$(src AnAsymmetry.cpp)" "$(src AnPlots.cpp)" "$(src AnConfig.h)" "$(src AnRunTools.h)" "$(src AnMerge.cpp)" "$(src pi0_bgsub.C)" "$(src pi0_peak_runs.C)" "$(src plot_event_cuts.C)" "${REP}/build/"
find -L "${OUTDIR}/hists" -maxdepth 1 -name 'an_hists_*.root' -mmin +10 | sort > "${REP}/yields.list"
# jobs submitted so far, over all batches of the pass
NJOBS=$(cat "${OUTDIR}"/run_params_*.txt 2>/dev/null | awk -F= '$1 == "nProcesses" && $2 > 0 { n += $2 } END { if (n) print n }')
echo "$(wc -l < "${REP}/yields.list") yield files${NJOBS:+ of ${NJOBS} jobs} ($(date))" | tee "${REP}/completed.txt"

container() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 ${SMOKE_BIND:+-B ${SMOKE_BIND}} \
        /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd ${REP}; $1"
}

if [ -z "${SMOKE_PLOTS_ONLY}" ]; then
container "g++ -m64 -O2 build/AnAsymmetry.cpp -o build/AnAsymmetry \`root-config --cflags --libs\` -std=c++11 && g++ -m64 -O2 build/AnPlots.cpp -o build/AnPlots \`root-config --cflags --libs\` -std=c++11 && g++ -m64 -O2 build/AnMerge.cpp -o build/AnMerge \`root-config --cflags --libs\` -std=c++11"
# Merged flow (Sean 2026-11-16, submit_an_merge.sh): when <outdir>/merge holds
# the partial merges, the yields are not read again. AnMerge -a sums the
# partials, the pi0 peak runs on the per-run spectra, AnMerge -x subtracts the
# masked runs, and every A_N below reads the one masked file (AnAsymmetry -M).
# All partials must be there (SMOKE_MERGE_PARTIAL=1 accepts fewer);
# SMOKE_NO_MERGE=1 forces the per-file flow.
MERGED=""
if [ -z "${SMOKE_NO_MERGE}" ] && ls "${OUTDIR}"/merge/partial_*.root >/dev/null 2>&1; then
    NPART=$(ls "${OUTDIR}"/merge/partial_*.root | wc -l)
    NCHUNK=$(ls "${OUTDIR}"/merge/chunks/chunk_*.list 2>/dev/null | wc -l)
    if [ "${NPART}" != "${NCHUNK}" ] && [ -z "${SMOKE_MERGE_PARTIAL}" ]; then
        echo "ERROR: ${NPART} of ${NCHUNK} partial merges in ${OUTDIR}/merge (wait, or SMOKE_MERGE_PARTIAL=1)"; exit 1
    fi
    MERGED=1
    ls "${OUTDIR}"/merge/partial_*.root > "${REP}/partials.list"
    echo "merged flow: ${NPART} of ${NCHUNK} partial merges"
    container "./build/AnMerge -a partials.list merged.root" > "${REP}/anmerge.log" 2>&1 || { tail "${REP}/anmerge.log"; exit 1; }
    cat "${REP}/anmerge.log"
fi
PIN="yields.list"; [ -n "${MERGED}" ] && PIN="merged.root"
# run-by-run pi0 peak first (Sean 2026-11-16): checks the energy scale factor and
# flags runs; its bad-run list masks runs in every A_N below (-Q), together with
# D. Kapukchyan's list, until his covers the whole run. SMOKE_NO_PI0_MASK=1 skips the mask.
QOPT=""
rm -f "${REP}/pi0_peak_runs_bad_runs.txt"
container "root -b -q -l 'build/pi0_peak_runs.C+(\"${PIN}\",\"pi0_peak_runs\")'" > "${REP}/pi0_peak_runs.log" 2>&1
grep "^pi0_peak_runs:" "${REP}/pi0_peak_runs.log" || true
# event-level cut variables with the cut (solid) and its variations (dotted), N-1
# (Sean 2026-11-17; yield files from 2026-11-17 on, else it only warns)
container "root -b -q -l 'build/plot_event_cuts.C+(\"${PIN}\",\"event_cuts\")'" > "${REP}/event_cuts.log" 2>&1 || true
grep "^plot_event_cuts:" "${REP}/event_cuts.log" || true
if [ -s "${REP}/pi0_peak_runs_bad_runs.txt" ] && [ -z "${SMOKE_NO_PI0_MASK}" ]; then
    QOPT="-Q pi0_peak_runs_bad_runs.txt"
    echo "pi0-peak bad runs masked: $(grep -c '^[0-9]' "${REP}/pi0_peak_runs_bad_runs.txt")"
fi
# the A_N input: the masked merged file, or the per-file yields with the run selection
if [ -n "${MERGED}" ]; then
    container "./build/AnMerge -x merged.root ${OUTDIR}/hists $(src run22_polarization.txt) masked.root ${QOPT}" > "${REP}/anmerge_mask.log" 2>&1 || { tail "${REP}/anmerge_mask.log"; exit 1; }
    head -1 "${REP}/anmerge_mask.log"
    ANIN="-M masked.root"
else
    ANIN="${QOPT} -r yields.list $(src run22_polarization.txt)"
fi
# AN_XF_EDGES="e0,e1,...": variable x_F bins for every A_N below (AnAsymmetry -x,
# edges on the stored 0.025 grid); unset: the uniform 0.05 bins
if [ -n "${AN_XF_EDGES}" ]; then
    ANIN="-x ${AN_XF_EDGES} ${ANIN}"
    echo "x_F bin edges: ${AN_XF_EDGES}"
fi
container "./build/AnAsymmetry ${ANIN} an.root" > "${REP}/an_asymmetry.log"
grep -A20 "^all, all photon classes" "${REP}/an_asymmetry.log" || true
container "./build/AnAsymmetry -J ${ANIN} an_union.root" > "${REP}/an_asymmetry_union.log"
# FMS-style pi0 isolation (Sean 2026-11-16; AnConfig.h): isolated (z_em > 0.98),
# non-isolated (< 0.9) and all pi0 candidates, signal and sideband windows,
# then the sideband-subtracted A_N (pi0_bgsub.C); yield files from 2026-11-16 on
container "./build/AnAsymmetry -J -I incl:sig ${ANIN} an_pi0_incl_sig.root" > "${REP}/an_asymmetry_pi0_incl_sig.log" 2>&1 || true
if ! grep -q "has no h_xFpi_" "${REP}/an_asymmetry_pi0_incl_sig.log"; then
    for c in iso noniso incl; do
        for w in sig sb; do
            [ "${c}:${w}" = "incl:sig" ] && continue
            container "./build/AnAsymmetry -J -I ${c}:${w} ${ANIN} an_pi0_${c}_${w}.root" > "${REP}/an_asymmetry_pi0_${c}_${w}.log"
        done
        container "root -b -q -l 'build/pi0_bgsub.C+(\"an_pi0_${c}_sig.root\",\"an_pi0_${c}_sb.root\",\"an_pi0_${c}\",\"all\",\"${c}\")'" > "${REP}/pi0_bgsub_${c}.log" 2>&1
    done
else
    echo "no h_xFpi_* in the yield files (AnHistograms before 2026-11-16): no pi0 isolation A_N"
fi
else
    # SMOKE_PLOTS_ONLY=1: an*.root and an.txt are already in the
    # report folder (e.g. from condor jobs); redo the projection and the plots
    :
fi

# the full sample = all good-run events, disk and tape (run22_run_events_all.txt,
# from the tape catalog; Sean 2026-11-11), when that table is there
ALLEV="$(src run22_run_events_all.txt)"
"${HERE}/project_errors.py" "${OUTDIR}" "${REP}/an.txt" "${GOOD}" $([ -f "${ALLEV}" ] && echo "${ALLEV}") | tee "${REP}/projection.log"
# every A_N plot (list, axes, text): replot_report.sh, which also redraws them
# later from the report folder alone
"${HERE}/replot_report.sh" "${REP}"
if [ -n "${SMOKE_PLOTS_ONLY}" ]; then echo "report plots redone in ${REP}"; exit 0; fi
if [ -f "${REP}/masked.root" ] && [ -z "${SMOKE_NO_MERGE}" ]; then
    # merged flow: the summed pages show the masked sample (the runs in A_N); the per-run rows every run
    container "./build/AnPlots -m masked.root plots" > "${REP}/anplots.log" 2>&1 || { tail "${REP}/anplots.log"; exit 1; }
else
    container "./build/AnPlots yields.list plots" > "${REP}/anplots.log" 2>&1 || { tail "${REP}/anplots.log"; exit 1; }
fi

"$(src estimate_resources.py)" "${OUTDIR}/timing" --runs "${GOOD}" > "${REP}/resources.txt"
cat "${REP}/resources.txt"

REQ=$(grep '^request_id=' "${OUTDIR}"/run_params_*.txt | head -1 | cut -d= -f2)
if [ -n "$REQ" ]; then
    {
        echo "# condor_history for request ${REQ}: job wall_s cpu_s memory_MB (condor MemoryUsage) disk_kB (job sandbox only, not scratch) exit_code"
        condor_history seanp -constraint "regexp(\"${REQ}\", Arguments)" \
            -af:h ProcId RemoteWallClockTime RemoteUserCpu MemoryUsage DiskUsage ExitCode 2>/dev/null | sort -n
    } > "${REP}/condor_jobs.txt"
    awk '$1 ~ /^[0-9]+$/ && $2 > 0 { n++; w += $2; if ($2 > wm) wm = $2; if ($4 > mm) mm = $4 }
         END { if (n) printf "condor: %d jobs, mean wall %.0f min, max %.0f min, max memory %d MB\n", n, w / n / 60, wm / 60, mm }' "${REP}/condor_jobs.txt"
fi
echo "report in ${REP}"
