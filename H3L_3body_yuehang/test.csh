#! /bin/csh
starver .DEV 
source setDEV2.csh
cons
touch test.log
root4star -b -q -l analysispicodst_h3l3b_tree.C\(\"test_mc.list\",\"null\",false\) >> test.log
# ls -lhrt
