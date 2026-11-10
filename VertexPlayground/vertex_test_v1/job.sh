#!/bin/bash
# Vertex test job: $1 = run. MuDst -> SimpleTree (production trigger filter) -> vertex branches only.
T=/star/data01/pwg/seanp/vertex_test_v1; RUN=$1
W=${_CONDOR_SCRATCH_DIR:-$PWD}/w; mkdir -p $W; cp -r $T/code/. $W/; cp $T/list_$RUN.txt $W/
cd $W
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B "$W" --pwd "$W" /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh inner.csh $RUN $T/out
ls -la $T/out/vtx_${RUN}_*.root
