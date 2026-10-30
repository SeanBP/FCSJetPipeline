source /star/nfs4/AFS/star/group/templates/cshrc >& /dev/null
stardev >& /dev/null
setenv LHAPDF_DATA_PATH /star/u/seanp/FCSJetPipeline/sim_to_jet/lhapdf_data
if ( ! -e pythia6428_nopdf.o ) gfortran -O2 -fno-automatic -c pythia6428_nopdf.f || exit 1
g++ -O2 recoil_truth.cpp pythia6428_nopdf.o -o recoil_truth `root-config --cflags --libs` -I$OPTSTAR/include -L$OPTSTAR/lib -lfastjet -lLHAPDF -lgfortran -Wl,-rpath,$OPTSTAR/lib || exit 1
if ( "$1" != "" ) then
  time ./recoil_truth $argv[1-6]
endif
