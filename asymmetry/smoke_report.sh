#!/bin/bash
# Offline A_N steps on a finished an_job pass (step 6 and the checks), into
# <outdir>/report/:
#   an.root, an.txt            AnAsymmetry -r (bad runs, per-run polarization)
#   an_vs_xf.png               final A_N vs x_F, statistical errors
#   an_projected.txt           the table with errors projected to the full
#                              good-run sample (project_errors.py)
#   an_vs_xf_projected.png     the points with projected errors, and
#   an_vs_xf_projected_zero.png the projected errors alone, at A_N = 0
#   an_vs_xf_<EM level>.png    A_N per EM level
#   plots/                     yields and sanity plots, per-run trends (AnPlots)
#   resources.txt              wall time, memory, scratch per stage and job,
#                              and the full-pass extrapolation
#                              (estimate_resources.py on the timing tables)
#   condor_jobs.txt            condor's own record per job: wall clock,
#                              memory, disk, exit code (condor_history)
#
# Usage: ./smoke_report.sh <outdir> <good_runs.txt>
#   <good_runs.txt>: make_smoke_list.py's <prefix>_good_runs.txt, the
#   full sample the errors and resources are projected to
# Compiles in the container; run it from the host (condor_history).

set -e
OUTDIR="$(cd "$1" && pwd)"
GOOD="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
[ -d "${OUTDIR}/hists" ] && [ -f "$GOOD" ] || { sed -n '2,24p' "$0"; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REP="${OUTDIR}/report"
mkdir -p "${REP}/build"
cp "${HERE}"/{AnAsymmetry.cpp,AnPlots.cpp,AnConfig.h,plot_an.C} "${REP}/build/"
ls "${OUTDIR}"/hists/an_hists_*.root > "${REP}/yields.list"
echo "$(wc -l < "${REP}/yields.list") yield files"

container() {
    singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
        /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c \
        "source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; module purge >& /dev/null; module load root-5.34.38 >& /dev/null; cd ${REP}; $1"
}

container "g++ -m64 -O2 build/AnAsymmetry.cpp -o build/AnAsymmetry \`root-config --cflags --libs\` -std=c++11 && g++ -m64 -O2 build/AnPlots.cpp -o build/AnPlots \`root-config --cflags --libs\` -std=c++11"
container "./build/AnAsymmetry -r yields.list ${HERE}/run22_polarization.txt an.root" > "${REP}/an_asymmetry.log"
grep -A20 "^all, all photon classes" "${REP}/an_asymmetry.log" || true

"${HERE}/project_errors.py" "${OUTDIR}" "${REP}/an.txt" "${GOOD}" | tee "${REP}/projection.log"
SCALE=$(grep '^scale=' "${REP}/projection.log" | cut -d= -f2)

cat > "${REP}/plot_an.csh" <<PLOTS
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf.png")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_projected.png","all","all",${SCALE},0)' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_projected_zero.png","all","all",${SCALE},1)' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_EM0.png","EM0","all")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_EM1.png","EM1","all")' || exit 1
root -b -q -l 'build/plot_an.C("an.root","an_vs_xf_EM2.png","EM2","all")' || exit 1
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
    awk 'NR > 2 && $2 > 0 { n++; w += $2; if ($2 > wm) wm = $2; if ($4 > mm) mm = $4 }
         END { if (n) printf "condor: %d jobs, mean wall %.0f min, max %.0f min, max memory %d MB\n", n, w / n / 60, wm / 60, mm }' "${REP}/condor_jobs.txt"
fi
echo "report in ${REP}"
