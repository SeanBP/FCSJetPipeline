#!/bin/sh
# condor wrapper: ./run_job.sh <ckin3> <ckin4> <events> <seed> <mpi> <out.root>
exec singularity exec -e -B /direct -B /star -B /afs -B /gpfs /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif tcsh -c "source /star/nfs4/AFS/star/group/templates/cshrc >& /dev/null; stardev >& /dev/null; setenv LHAPDF_DATA_PATH \$PWD/lhapdf; ./recoil_truth $1 $2 $3 $4 $5 $6"
