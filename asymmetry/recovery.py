#!/usr/bin/env python3
"""Recovery of an A_N pass's missing files (Sean 2026-11-16).

The distributed-disk nodes that hold the MuDsts drop out now and then and
usually come back within days. During a pass, an_job flags a node after
-L failed fetches (<outdir>/nodes/flagged/<node>) and skips its remaining
files. Every file's outcome is in the job's timing table (done, fetchfail,
nodeskip, mudstfail, jetfail). This tool collects what is missing and runs it
again once the nodes serve, as a recovery pass in <outdir>/recovery_<k>. The
yields go into the same sums (linked into hists/, merged by
submit_an_merge.sh), so the final results include them.

  recovery.py status <outdir>
      Every file of the pass and of its recovery passes: done, or why not:
        fetchfail  the copy failed twice
        nodeskip   its node was flagged
        mudstfail  runMudst.C wrote no SimpleTree (corrupt file?)
        jetfail    RecoJets / UnderlyingEvent failed
        lost       the file was done but its job's yield file is missing (the
                   job died later or was removed)
        notrun     no timing table: the job never ran (held, removed, still
                   queued)
      Writes <outdir>/recovery/status.txt (one line per file: status, node,
      list line) and prints the counts per status and per node. A file is
      done if any pass did it.
  recovery.py probe <outdir> [timeout s = 120]
      For every node with missing files, copies one of them (xrdcp, timeout):
      up or down. Writes <outdir>/recovery/probe.txt.
  recovery.py submit <outdir> [-n files per job = 40] [--mudst] [--down] [-N]
      A recovery pass over the missing files whose node probed up (and files
      with no node): submit_an_job.sh with the pass's own calibration,
      triggers, jet keeping and pi0 scale table, output
      <outdir>/recovery_<k>. --mudst also retries mudstfail / jetfail files
      (left out by default: they usually fail again); --down includes nodes
      that probed down; -N writes the job without submitting.
  recovery.py link <outdir> [-t <hists dir>]
      Links the finished recovery passes' yield (and jet) files into
      <outdir>/hists (or -t <dir>, for a combined pass directory) as
      an_hists_<run>_rec<k>_<job>.root. Then run submit_an_merge.sh again: it
      merges only the new files.
"""
import os
import re
import subprocess
import sys
import glob
import datetime
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL_STAGES = ('fetchfail', 'nodeskip', 'mudstfail', 'jetfail')
RUN_RE = re.compile(r'_(\d{8})_(?:raw_)?\d{7}\.MuDst')


def die(msg):
    sys.exit('recovery.py: ' + msg)


def read_params(d):
    ps = sorted(glob.glob(os.path.join(d, 'run_params_*.txt')))
    if len(ps) != 1:
        die('%s: expected one run_params_*.txt, found %d' % (d, len(ps)))
    p = {}
    for line in open(ps[0]):
        if '=' in line and not line.startswith('#'):
            k, v = line.rstrip('\n').split('=', 1)
            p[k] = v
    return p


def pass_files(d):
    """[(job index, tag, [(list line, local file name, node)])] of one submit_an_job.sh outdir"""
    p = read_params(d)
    label = p['label']
    sl = os.path.join(d, 'filelist_%s.sorted.list' % label)
    rg = os.path.join(d, 'filelist_%s.ranges' % label)
    nd = os.path.join(d, 'filelist_%s.sorted.nodes' % label)
    lines = [l.rstrip('\n') for l in open(sl)]
    nodes = [l.strip() for l in open(nd)] if os.path.exists(nd) else ['unknown'] * len(lines)
    jobs = []
    for j, r in enumerate(open(rg)):
        start, count = (int(x) for x in r.split())
        fl = []
        for k in range(start - 1, start - 1 + count):
            url = lines[k].rsplit('/', 1)[0]
            fl.append((lines[k], os.path.basename(url), nodes[k] if k < len(nodes) else 'unknown'))
        m = RUN_RE.search(fl[0][0])
        jobs.append((j, '%s_%d' % (m.group(1), j), fl))
    return jobs


def pass_status(d):
    """{list line: (status, node)} for one pass directory"""
    out = {}
    for j, tag, fl in pass_files(d):
        tim = os.path.join(d, 'timing', 'timing_%s.txt' % tag)
        hist = os.path.join(d, 'hists', 'an_hists_%s.root' % tag)
        if not os.path.exists(tim):
            for line, f, node in fl:
                out[line] = ('notrun', node)
            continue
        last = {}
        for row in open(tim):
            w = row.split()
            if len(w) >= 2 and w[0] in FAIL_STAGES + ('done',):
                last[w[1]] = w[0]
        have = os.path.exists(hist)
        for line, f, node in fl:
            st = last.get(f, 'lost' if have else 'notrun')
            if st == 'done' and not have:
                st = 'lost'
            out[line] = (st, node)
    return out


def recovery_dirs(d):
    ds = [x for x in glob.glob(os.path.join(d, 'recovery_*')) if re.search(r'recovery_\d+$', x)]
    return sorted(ds, key=lambda x: int(x.rsplit('_', 1)[1]))


def status(d, write=True):
    st = pass_status(d)
    for r in recovery_dirs(d):
        if not glob.glob(os.path.join(r, 'run_params_*.txt')):
            continue
        for line, v in pass_status(r).items():
            if line in st and st[line][0] != 'done':
                st[line] = v
    if write:
        os.makedirs(os.path.join(d, 'recovery'), exist_ok=True)
        with open(os.path.join(d, 'recovery', 'status.txt'), 'w') as f:
            f.write('# %s: status node list-line (recovery.py status)\n' % datetime.datetime.now().strftime('%F %T'))
            for line, (s, node) in st.items():
                f.write('%s %s %s\n' % (s, node, line))
    return st


def events(line):
    try:
        return int(line.rsplit('/', 1)[1])
    except ValueError:
        return 0


def cmd_status(d):
    st = status(d)
    n = Counter(s for s, _ in st.values())
    ev = Counter()
    for line, (s, _) in st.items():
        ev[s] += events(line)
    tot = sum(ev.values())
    print('%s: %d files, %d recovery passes' % (d, len(st), len(recovery_dirs(d))))
    for s in ('done', 'fetchfail', 'nodeskip', 'mudstfail', 'jetfail', 'lost', 'notrun'):
        if n[s]:
            print('  %-10s %7d files  %6.2f%% of events' % (s, n[s], 100.0 * ev[s] / max(1, tot)))
    per = defaultdict(Counter)
    for s, node in st.values():
        if s != 'done':
            per[node][s] += 1
    if per:
        print('missing files per node:')
        for node, c in sorted(per.items(), key=lambda x: -sum(x[1].values())):
            print('  %-24s %6d  %s' % (node, sum(c.values()), ', '.join('%s %d' % kv for kv in sorted(c.items()))))
    fl = sorted(glob.glob(os.path.join(d, 'nodes', 'flagged', '*')))
    if fl:
        print('flagged nodes:')
        for f in fl:
            print('  ' + open(f).read().strip())
    print('written: %s' % os.path.join(d, 'recovery', 'status.txt'))


def star(cmd, timeout=None):
    full = ('singularity exec -e -B /direct -B /star -B /afs -B /gpfs /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif '
            "tcsh -c 'source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null; %s'" % cmd)
    try:
        return subprocess.run(full, shell=True, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).returncode
    except subprocess.TimeoutExpired:
        return 124


def cmd_probe(d, timeout=120):
    st = status(d)
    pick = {}
    for line, (s, node) in st.items():
        if s in ('fetchfail', 'nodeskip', 'notrun', 'lost') and node != 'unknown' and node not in pick:
            pick[node] = line.rsplit('/', 1)[0]
    if not pick:
        print('no nodes with missing files')
        return
    tmp = os.path.join(os.environ.get('SCRATCH', '/tmp'), 'recovery_probe_%d' % os.getpid())
    os.makedirs(tmp, exist_ok=True)
    out = []
    for node, url in sorted(pick.items()):
        dst = os.path.join(tmp, 'probe.root')
        t0 = datetime.datetime.now()
        rc = star('timeout %d xrdcp -s -f %s %s' % (timeout, url, dst), timeout=timeout + 60)
        dt = (datetime.datetime.now() - t0).total_seconds()
        ok = rc == 0 and os.path.exists(dst) and os.path.getsize(dst) > 0
        if os.path.exists(dst):
            os.remove(dst)
        out.append('%s %s %.0f %s' % (node, 'up' if ok else 'down', dt, url))
        print('  %-24s %-4s (%.0f s)' % (node, 'up' if ok else 'down', dt))
    os.rmdir(tmp)
    os.makedirs(os.path.join(d, 'recovery'), exist_ok=True)
    with open(os.path.join(d, 'recovery', 'probe.txt'), 'w') as f:
        f.write('# %s: node up|down seconds probed-file (recovery.py probe)\n' % datetime.datetime.now().strftime('%F %T'))
        f.write('\n'.join(out) + '\n')


def jobs_in_condor(d):
    """jobs of a pass directory (its request_id) still idle, running or held"""
    req = read_params(d).get('request_id', '')
    if not req:
        return 0
    c = 'regexp("%s", Args) || regexp("%s", Arguments) || regexp("%s", Cmd)' % (req, req, req)
    try:
        out = subprocess.run(['condor_q', os.environ.get('USER', ''), '-constraint', c, '-af', 'ProcId'],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True, timeout=300).stdout
    except Exception:
        return -1
    return len([l for l in out.splitlines() if l.strip()])


def cmd_submit(d, args):
    nper, mudst, down, dry = '40', False, False, False
    i = 0
    while i < len(args):
        a = args[i]
        if a == '-n':
            nper = args[i + 1]; i += 1
        elif a == '--mudst':
            mudst = True
        elif a == '--down':
            down = True
        elif a == '-N':
            dry = True
        else:
            die('unknown option ' + a)
        i += 1
    # a file whose job is still queued would be done twice: wait for the pass
    for pd in [d] + recovery_dirs(d):
        nq = jobs_in_condor(pd)
        if nq != 0:
            die('%s still has %s jobs in condor; recover once they are done (or condor_rm them)' % (pd, nq if nq > 0 else 'unknown'))
    st = status(d)
    pf = os.path.join(d, 'recovery', 'probe.txt')
    up = None
    if not down:
        if not os.path.exists(pf):
            die('no %s: run "recovery.py probe %s" first (or --down to take every node)' % (pf, d))
        up = {l.split()[0] for l in open(pf) if not l.startswith('#') and l.split()[1] == 'up'}
    want = ('fetchfail', 'nodeskip', 'notrun', 'lost') + (('mudstfail', 'jetfail') if mudst else ())
    sel = [line for line, (s, node) in st.items() if s in want and (up is None or node == 'unknown' or node in up)]
    if not sel:
        print('nothing to recover (%d files missing; see recovery.py status)' % sum(1 for s, _ in st.values() if s != 'done'))
        return
    k = len(recovery_dirs(d)) + 1
    rd = os.path.join(d, 'recovery_%d' % k)
    os.makedirs(os.path.join(d, 'recovery'), exist_ok=True)
    lst = os.path.join(d, 'recovery', 'recovery_%d.list' % k)
    with open(lst, 'w') as f:
        f.write('\n'.join(sel) + '\n')
    p = read_params(d)
    cmd = [os.path.join(HERE, 'submit_an_job.sh'), '-f', lst, '-o', rd, '-n', nper, '-l', '%s_rec%d' % (p['label'], k),
           '-c', p['calib_json'], '-t', p['trigger_filter'], '-J', p.get('keep_jets_every', '0'), '-k', p.get('keep_logs', '0'),
           '-L', p.get('node_fail_limit', '3')]
    sc = os.path.join(d, 'pi0_scale_%s.txt' % p['label'])
    if os.path.exists(sc):
        cmd += ['-s', sc]
    if dry:
        cmd.append('-N')
    print('recovery pass %d: %d files (%s) -> %s' % (k, len(sel), ', '.join('%s %d' % kv for kv in Counter(st[l][0] for l in sel).items()), rd))
    print(' '.join(cmd))
    rc = subprocess.call(cmd)
    if rc != 0:
        die('submit_an_job.sh failed (%d)' % rc)


def cmd_link(d, args):
    target = os.path.join(d, 'hists')
    if len(args) >= 2 and args[0] == '-t':
        target = os.path.abspath(args[1])
    n = 0
    for r in recovery_dirs(d):
        k = r.rsplit('_', 1)[1]
        for sub, prefix in (('hists', 'an_hists_'), ('jets', 'jets_')):
            for f in sorted(glob.glob(os.path.join(r, sub, prefix + '*.root'))):
                m = re.match(re.escape(prefix) + r'(\d{8})_(\d+)\.root$', os.path.basename(f))
                if not m:
                    continue
                tdir = target if sub == 'hists' else os.path.join(os.path.dirname(target), 'jets')
                if sub == 'jets' and not os.path.isdir(tdir):
                    continue
                dst = os.path.join(tdir, '%s%s_rec%s_%s.root' % (prefix, m.group(1), k, m.group(2)))
                if not os.path.lexists(dst):
                    os.symlink(os.path.abspath(f), dst)
                    n += 1
    print('linked %d recovery files into %s (and its jets/); run submit_an_merge.sh %s to merge the new ones' % (n, target, os.path.dirname(target)))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    cmd, d = sys.argv[1], os.path.abspath(sys.argv[2])
    if cmd == 'status':
        cmd_status(d)
    elif cmd == 'probe':
        cmd_probe(d, int(sys.argv[3]) if len(sys.argv) > 3 else 120)
    elif cmd == 'submit':
        cmd_submit(d, sys.argv[3:])
    elif cmd == 'link':
        cmd_link(d, sys.argv[3:])
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main()
