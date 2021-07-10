FILELIST="onetest.list"
starver SL20d
mv StRoot StRoot_b
mv StRoot_v2 StRoot
cons
root4star -q -l -b StRoot/macro/readPicoDst.C\(\"$FILELIST\",\"test\"\) 
ls -lhrt

mv StRoot StRoot_v2
mv StRoot_b StRoot
source setDEV2.csh
cons
root4star -b -q -l analysispicodst_h3l3b_tree_ME.C\(\"$FILELIST\",\"null\",true\) 

