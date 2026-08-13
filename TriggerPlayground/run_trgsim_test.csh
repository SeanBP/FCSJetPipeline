#!/bin/csh
cd /star/u/seanp/FCSJetPipeline/TriggerPlayground
source /star/nfs4/AFS/star/group/star_cshrc.csh
stardev
root4star -b -q TrgSimTest.C
