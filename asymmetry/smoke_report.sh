#!/bin/bash
# Offline A_N steps on a finished an_job pass (step 6 and the checks), into
# <outdir>/report/:
#   an.root, an.txt            AnAsymmetry -r (bad runs, per-run polarization;
#                              "all" = EM levels, "allJP" = JP levels)
#   an_union.root, .txt        the same with -J: "all" = EM + JP union, each
#                              jet once
#   an_vs_xf_union.png         the union A_N vs x_F
#   an_vs_xf_<JP level>.png    A_N per JP level and allJP
#   an_vs_xf.png               final A_N vs x_F, statistical errors
#   an_projected.txt           the table with errors projected to the full
#                              good-run sample (project_errors.py)
#   an_vs_xf_projected.png     the points with projected errors, and
#   an_vs_xf_projected_zero.png the projected errors alone, at A_N = 0
#   an_vs_xf_<EM level>.png    A_N per EM level
#   an_vs_xf_np<N>.png         A_N for jets with exactly N = 1, 2 points, and
#                              np3plus, 3 or more
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

set -e
OUTDIR="$(cd "$1" && pwd)"
GOOD="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
[ -d "${OUTDIR}/hists" ] && [ -f "$GOOD" ] || { sed -n '2,24p' "$0"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REP="${OUTDIR}/${3:-report}"
mkdir -p "${REP}/build"
cp "${HERE}"/{AnAsymmetry.cpp,AnPlots.cpp,AnConfig.h,plot_an.C} "${REP}/build/"
find "${OUTDIR}/hists" -maxdepth 1 -name 'an_hists_*.root' -mmin +10 | sort > "${REP}/yields.list"
# jobs submitted so far, over all batches of the pass
NJOBS=$(cat "${OUTDIR}"/run_params_*.txt 2>/dev/null | awk -F= '$1 == "nProcesses" && $2 > 0 { n += $2 } END { if (n) print n }')
echo "$(wc -l < "${REP}/yields.list") yield files${NJOBS:+ of ${NJOBS} jobs} ($(date))" | tee "${REP}/completed.txt"

container() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
        /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd ${REP}; $1"
}

container "g++ -m64 -O2 build/AnAsymmetry.cpp -o build/AnAsymmetry \`root-config --cflags --libs\` -std=c++11 && g++ -m64 -O2 build/AnPlots.cpp -o build/AnPlots \`root-config --cflags --libs\` -std=c++11"
container "./build/AnAsymmetry -r yields.list ${HERE}/run22_polarization.txt an.root" > "${REP}/an_asymmetry.log"
grep -A20 "^all, all photon classes" "${REP}/an_asymmetry.log" || true
container "./build/AnAsymmetry -J -r yields.list ${HERE}/run22_polarization.txt an_union.root" > "${REP}/an_asymmetry_union.log"

"${HERE}/project_errors.py" "${OUTDIR}" "${REP}/an.txt" "${GOOD}" | tee "${REP}/projection.log"
SCALE=$(grep '^scale=' "${REP}/projection.log" | cut -d= -f2)
# what the plots show: "<test runs> of <good runs> good runs, <events> events = <percent> of the Run 22 good-run events"
SAMPLE=$(awk '/^test:/ { gsub(/[;,]/, ""); printf "%s of %s good runs (up to 40 files each), %.2g events = %.1f%% of the Run 22 good-run events", $2, $9, $6, 100 * $15 }' "${REP}/projection.log")
cat > "${REP}/plot_an.csh" <<PLOTS
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf.png","all","all",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_projected.png","all","all",${SCALE},0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_projected_zero.png","all","all",${SCALE},1,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_EM0.png","EM0","all",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_EM1.png","EM1","all",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_EM2.png","EM2","all",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_np1.png","all","np1",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_np2.png","all","np2",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_np3plus.png","all","np3plus",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an_union.root","an_vs_xf_union.png","all","all",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_allJP.png","allJP","all",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_JP0.png","JP0","all",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_JP1.png","JP1","all",1,0,10,"${SAMPLE}")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_JP2.png","JP2","all",1,0,10,"${SAMPLE}")' || exit 1
PLOTS
container "tcsh plot_an.csh" > "${REP}/plot_an.log" 2>&1 || { tail "${REP}/plot_an.log"; exit 1; }
container "./build/AnPlots yields.list plots" > "${REP}/anplots.log" 2>&1 || { tail "${REP}/anplots.log"; exit 1; }

"${HERE}/estimate_resources.py" "${OUTDIR}/timing" --runs "${GOOD}" > "${REP}/resources.txt"
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
