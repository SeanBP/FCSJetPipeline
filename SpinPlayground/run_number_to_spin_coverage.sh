#!/bin/bash
# Pure run-number -> spin-DB-coverage tool. No ROOT/StRoot/singularity
# needed -- just direct MySQL queries against the real DB servers (both
# reachable directly from this host).
#
# Usage: ./run_number_to_spin_coverage.sh <run_number>
#
# Given a real Run 22 run number, checks whether the STAR spin/polarization
# pattern DB (Calibrations_rhic.spinV124 / .spinStar, on
# dbx.star.bnl.gov:3316 -- same schema StSpinDbMaker reads via
# St_db_Maker's "Calibrations/rhic" path) has an entry actually covering
# that run's real start time (RunLog.runDescriptor on db04.star.bnl.gov,
# same mechanism GeometryPlayground's RunNumberToDateTime uses).
#
# IMPORTANT: unlike the FCS gain/geometry tables, spinV124/spinStar entries
# are narrow, fill-scoped windows (a few hours each, matching one RHIC
# fill) with real gaps between them -- NOT "most recent entry persists
# until superseded". A simple date-based DB lookup (St_db_Maker style) can
# easily land in a gap between fills and wrongly look empty even when the
# DB is well populated nearby, so this checks the exact run start time
# directly against each entry's [beginTime, endTime) window instead of
# sampling arbitrary dates. Investigated Aug 2026: spinV124 (254 entries)
# and spinStar (3362 entries) cover 2021-12-14 through 2022-04-18 with no
# gap over ~2.3 days; there is zero coverage of either table after
# 2022-04-18, for the rest of Run 22.
#
# spinBXmask is also checked but is a known red herring: its most recent
# real entry is from 2013-06-02 with an open-ended endTime (2037-01-01),
# so it "covers" any Run 22 timestamp trivially -- that's stale
# pre-Run-22 data resolving forward, not real coverage, and it doesn't
# carry polarization direction anyway (just a bunch-fill mask).
#
# spinStar has a similar trap: one entry (dataID 22639, entryTime
# 2025-01-17 -- inserted years after Run 22, comment "Run23091017, oleg")
# is backdated to beginTime=2022-04-01 with endTime=2037-01-01, so it
# "covers" every date from 2022-04-01 onward, masking the real gap. Any
# entry with an unusually long validity window (this script flags >7
# days) is reported as a suspect fallback/patch, not real per-fill data --
# spinV124 has no such entries (checked directly), so it's the more
# reliable of the two to trust at face value. The final ALL-TABLES verdict
# below mirrors StSpinDbMaker::isValid() (mNFound==3, i.e. ALL of
# spinV124/spinStar/spinBXmask must resolve) but additionally requires the
# resolving entry not be a flagged long-duration one, since isValid()
# itself has no way to know a "covering" entry is actually a stale
# placeholder.

set -e

RUN="$1"
if [ -z "$RUN" ] || ! [[ "$RUN" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 <run_number>"
    exit 1
fi

# ---- Step 1: resolve run number -> real start time via RunLog ----
YEAR_CODE=$((RUN/1000000))
FORMULA_PORT=$((3400 + (YEAR_CODE-1) - 1))
CANDIDATE_PORTS="3421 $FORMULA_PORT $((FORMULA_PORT-1)) $((FORMULA_PORT+1))"

STARTTIME=""
FOUND_PORT=""
SEEN_PORTS=""
for PORT in $CANDIDATE_PORTS; do
    case " $SEEN_PORTS " in *" $PORT "*) continue ;; esac
    SEEN_PORTS="$SEEN_PORTS $PORT"
    RESULT=$(mysql -h db04.star.bnl.gov --port="$PORT" --connect-timeout=10 -N -s \
        -e "SELECT startRunTime FROM RunLog.runDescriptor WHERE runNumber=$RUN LIMIT 1" 2>/dev/null || true)
    if [[ "$RESULT" =~ ^[0-9]+$ ]]; then
        STARTTIME="$RESULT"
        FOUND_PORT="$PORT"
        break
    fi
done

if [ -z "$STARTTIME" ]; then
    echo "Error: could not find run $RUN in RunLog on any candidate port ($CANDIDATE_PORTS)"
    exit 1
fi

RUN_DATETIME=$(date -u -d @"$STARTTIME" +"%Y-%m-%d %H:%M:%S")
echo "Run $RUN -> RunLog (port=$FOUND_PORT) start time = $RUN_DATETIME UTC"
echo ""

# Long-validity-window entries (>7 days) are flagged as suspect
# fallback/patch data rather than real per-fill measurements -- both known
# instances found so far (spinBXmask's 2013 default, spinStar's 2025-
# inserted "Run23091017, oleg" patch) are exactly this shape: an
# open-ended endTime far in the future, silently "covering" everything
# after their beginTime regardless of whether real data actually exists.
LONG_DURATION_DAYS=7

# ---- Step 2: check spinV124/spinStar/spinBXmask coverage at that exact timestamp ----
declare -A REAL_COVERAGE
for TABLE in spinV124 spinStar spinBXmask; do
    ROW=$(mysql -h dbx.star.bnl.gov --port=3316 -N -s -e "
        SELECT dataID, beginTime, endTime,
               TIMESTAMPDIFF(DAY, beginTime, endTime) AS dur_days,
               entryTime
        FROM $TABLE
        WHERE beginTime <= '$RUN_DATETIME' AND endTime > '$RUN_DATETIME'
          AND flavor='ofl' AND deactive=0
        ORDER BY beginTime DESC LIMIT 1" Calibrations_rhic 2>&1)

    if [ -n "$ROW" ]; then
        DUR_DAYS=$(echo "$ROW" | awk -F'\t' '{print $4}')
        if [ "$DUR_DAYS" -gt "$LONG_DURATION_DAYS" ]; then
            echo "$TABLE: COVERED, but by a SUSPECT long-duration entry ($DUR_DAYS days) -- likely a"
            echo "  fallback/patch, not real per-fill data. Treat as NOT reliably covered."
            echo "  dataID / beginTime / endTime / dur_days / entryTime: $ROW"
            REAL_COVERAGE[$TABLE]=0
        else
            echo "$TABLE: COVERED (real, fill-scoped entry, $DUR_DAYS days)"
            echo "  dataID / beginTime / endTime / dur_days / entryTime: $ROW"
            REAL_COVERAGE[$TABLE]=1
        fi
    else
        echo "$TABLE: NOT COVERED (no entry spans this run's start time)"
        PREV=$(mysql -h dbx.star.bnl.gov --port=3316 -N -s -e "
            SELECT beginTime,endTime FROM $TABLE WHERE beginTime <= '$RUN_DATETIME'
              AND flavor='ofl' AND deactive=0
            ORDER BY beginTime DESC LIMIT 1" Calibrations_rhic 2>&1)
        NEXT=$(mysql -h dbx.star.bnl.gov --port=3316 -N -s -e "
            SELECT beginTime,endTime FROM $TABLE WHERE beginTime > '$RUN_DATETIME'
              AND flavor='ofl' AND deactive=0
            ORDER BY beginTime ASC LIMIT 1" Calibrations_rhic 2>&1)
        echo "  nearest prior entry (beginTime, endTime): ${PREV:-none}"
        echo "  nearest next entry  (beginTime, endTime): ${NEXT:-none}"
        REAL_COVERAGE[$TABLE]=0
    fi
    echo ""
done

# ---- Step 3: two distinct verdicts ----
# spinBXmask has ZERO genuine Run-22 entries (checked directly, Aug 2026)
# -- it always resolves via the 2013 stale default for every single Run 22
# run, so requiring it to be "real" would make any combined verdict
# uniformly negative and useless as a signal. Two separate answers instead:
V124_RESOLVES=$(mysql -h dbx.star.bnl.gov --port=3316 -N -s -e "SELECT COUNT(*) FROM spinV124 WHERE beginTime<='$RUN_DATETIME' AND endTime>'$RUN_DATETIME' AND flavor='ofl' AND deactive=0" Calibrations_rhic 2>&1)
STAR_RESOLVES=$(mysql -h dbx.star.bnl.gov --port=3316 -N -s -e "SELECT COUNT(*) FROM spinStar WHERE beginTime<='$RUN_DATETIME' AND endTime>'$RUN_DATETIME' AND flavor='ofl' AND deactive=0" Calibrations_rhic 2>&1)
BXMASK_RESOLVES=$(mysql -h dbx.star.bnl.gov --port=3316 -N -s -e "SELECT COUNT(*) FROM spinBXmask WHERE beginTime<='$RUN_DATETIME' AND endTime>'$RUN_DATETIME' AND flavor='ofl' AND deactive=0" Calibrations_rhic 2>&1)

echo "VERDICT 1 -- what StSpinDbMaker::isValid() would literally report"
echo "  (mNFound==3, i.e. all 3 tables resolve to SOME entry -- it has no"
echo "   concept of a stale/fallback entry, so a suspect long-duration hit"
echo "   still counts):"
if [ "$V124_RESOLVES" -gt 0 ] && [ "$STAR_RESOLVES" -gt 0 ] && [ "$BXMASK_RESOLVES" -gt 0 ]; then
    echo "    -> isValid() = TRUE"
else
    echo "    -> isValid() = FALSE"
fi
echo ""
echo "VERDICT 2 -- genuine, non-fallback polarization PATTERN data available"
echo "  (spinV124 + spinStar only, both real/short-duration entries --"
echo "   spinBXmask excluded, since it never has genuine Run-22 data and"
echo "   doesn't carry polarization direction anyway; this is the check"
echo "   that actually matters for per-event spin info):"
if [ "${REAL_COVERAGE[spinV124]:-0}" = "1" ] && [ "${REAL_COVERAGE[spinStar]:-0}" = "1" ]; then
    echo "    -> YES"
else
    echo "    -> NO"
fi
