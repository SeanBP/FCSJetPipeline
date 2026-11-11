#!/bin/bash
# Offline A_N steps on a finished an_job pass (step 6 and the checks), into
# <outdir>/report/:
#   an.root, an.txt            AnAsymmetry -r (bad runs, per-run polarization;
#                              "all" = EM levels, "allJP" = JP levels)
#   an_union.root, .txt        the same with -J: "all" = EM + JP union, each
#                              jet once
#   an_vs_xf_union.png         the union A_N vs x_F
#   an_vs_xf_<JP level>.png    A_N per JP level and allJP
#   mgg_peak.png               two-photon jets' pair mass, the fitted pi0 peak
#   an_mass_in/out.root, .txt  AnAsymmetry -J -m in/out:<peak>:<sigma>: np2 rows =
#                              the two-photon jets inside / outside peak +- 2 sigma
#   an_vs_xf_mass.png          their A_N vs x_F, inside vs outside
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
# SMOKE_BIND=<dir>: also bind <dir> into the container (an outdir on a
# condor job's scratch disk; see report_job.sh).

set -e
OUTDIR="$(cd "$1" && pwd)"
GOOD="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
[ -d "${OUTDIR}/hists" ] && [ -f "$GOOD" ] || { sed -n '2,24p' "$0"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# a file of this script's directory, else of the one above (repo layout:
# asymmetry/sivers/ with the pass-level files in asymmetry/; a frozen report
# code directory has them all side by side)
src() { if [ -e "${HERE}/$1" ]; then echo "${HERE}/$1"; else echo "${HERE}/../$1"; fi; }
REP="${OUTDIR}/${3:-report}"
mkdir -p "${REP}/build"
cp "$(src AnAsymmetry.cpp)" "$(src AnPlots.cpp)" "$(src AnConfig.h)" "$(src mgg_peak.C)" "${REP}/build/"
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
container "g++ -m64 -O2 build/AnAsymmetry.cpp -o build/AnAsymmetry \`root-config --cflags --libs\` -std=c++11 && g++ -m64 -O2 build/AnPlots.cpp -o build/AnPlots \`root-config --cflags --libs\` -std=c++11"
container "./build/AnAsymmetry -r yields.list $(src run22_polarization.txt) an.root" > "${REP}/an_asymmetry.log"
grep -A20 "^all, all photon classes" "${REP}/an_asymmetry.log" || true
container "./build/AnAsymmetry -J -r yields.list $(src run22_polarization.txt) an_union.root" > "${REP}/an_asymmetry_union.log"
# two-photon jets in vs out of the pi0 mass window (fitted peak +- 2 sigma)
container "root -b -q -l 'build/mgg_peak.C(\"yields.list\",\"mgg_peak.png\")'" > "${REP}/mgg_peak.log" 2>&1
PEAK=$(sed -n 's/^peak=//p' "${REP}/mgg_peak.log" | tail -1)
SIGMA=$(sed -n 's/^sigma=//p' "${REP}/mgg_peak.log" | tail -1)
echo "pi0 peak: ${PEAK} GeV, sigma ${SIGMA} GeV"
if [ -n "$PEAK" ] && [ "$PEAK" != "0" ]; then
    container "./build/AnAsymmetry -J -m in:${PEAK}:${SIGMA} -r yields.list $(src run22_polarization.txt) an_mass_in.root" > "${REP}/an_asymmetry_mass_in.log"
    container "./build/AnAsymmetry -J -m out:${PEAK}:${SIGMA} -r yields.list $(src run22_polarization.txt) an_mass_out.root" > "${REP}/an_asymmetry_mass_out.log"
fi
else
    # SMOKE_PLOTS_ONLY=1: an*.root, an.txt and mgg_peak.log are already in the
    # report folder (e.g. from condor jobs); redo the projection and the plots
    PEAK=$(sed -n 's/^peak=//p' "${REP}/mgg_peak.log" | tail -1)
fi

# the full sample = all good-run events, disk and tape (run22_run_events_all.txt,
# from the tape catalog; Sean 2026-11-11), when that table is there
ALLEV="$(src run22_run_events_all.txt)"
"${HERE}/project_errors.py" "${OUTDIR}" "${REP}/an.txt" "${GOOD}" $([ -f "${ALLEV}" ] && echo "${ALLEV}") | tee "${REP}/projection.log"
# every A_N plot (list, axes, text): replot_report.sh, which also redraws them
# later from the report folder alone
"${HERE}/replot_report.sh" "${REP}"
if [ -n "${SMOKE_PLOTS_ONLY}" ]; then echo "report plots redone in ${REP}"; exit 0; fi
container "./build/AnPlots yields.list plots" > "${REP}/anplots.log" 2>&1 || { tail "${REP}/anplots.log"; exit 1; }

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
