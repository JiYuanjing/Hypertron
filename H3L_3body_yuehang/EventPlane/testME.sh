#! /bin/csh
set FILELIST="onetest.list"
starver SL20d
# cd EventPlane
cons
root -q -l -b readPicoDst.C\(\"../$FILELIST\",\"test\"\) 
#
# cd ../
# ls -lhrt
#
# source setDEV2.csh
# cons
# root4star -b -q -l analysispicodst_h3l3b_tree_ME.C\(\"$FILELIST\",\"kfout.root\",true\) 
