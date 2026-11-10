#!/bin/tcsh
# $1 = run, $2 = output dir
set RUN = $1
set OUT = $2
set argv = ()
source /star/nfs4/AFS/star/group/star_cshrc.csh >& /dev/null
stardev >& /dev/null

set n = 0
foreach line (`cat list_${RUN}.txt`)
  set url = `echo $line | rev | cut -d/ -f2- | rev`
  set f = `basename $url`
  timeout 300 xrdcp --retry 3 -s $url $f
  if ( $status != 0 ) then
    echo "fetchfail $f"
    rm -f $f
    @ n++
    continue
  endif
  root4star -b -q runMudst.C\(\"$f\",-1,-2,\".\",1,0,\"fcsEM0,fcsEM1,fcsEM2,fcsEM0_tpc,fcsEM1_tpc,fcsEM2_tpc,fcsJP2,fcsJPA1,fcsJPBC1,fcsJPDE1,fcsJPA0,fcsJPBC0,fcsJPDE0\"\) >& mudst_$n.log
  set st = `ls -t SimpleTree_mudst_*.root | head -1`
  root4star -b -q slim.C\(\"$st\",\"${OUT}/vtx_${RUN}_$n.root\"\) >& slim_$n.log
  echo "file ${n}, $f -> $st"; tail -n 3 slim_$n.log
  rm -f $f $st
  @ n++
end
