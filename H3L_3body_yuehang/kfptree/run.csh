starver SL19e
mv StRoot StRoot_b
mv StRoot_v2 StRoot
cons
root4star -q -l -b readPicoDst.C | tee >& shift.log
ls -lhrt

mv StRoot StRoot_v2
mv StRoot_b StRoot
source setDEV2.csh
cons
root4star -b -q -l StRoot/macro/analysispicodst_ldxi_pol.C | tee a.log
ls -lhrt
