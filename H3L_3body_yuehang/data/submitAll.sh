#!/bin/bash
date
mkdir log list outHist script
rm -rf list/*.list_*
root -l -q splitlist.C\($1,\"$2\" \)
NUM=`ls list/$2* | wc -l`
echo $NUM
dir=$PWG
ifile=0
while [[ $ifile -le $NUM ]];
do
  ###############################################
  cp run.con run$2_${ifile}.con
  echo "Initialdir      = $dir" >> run$2_${ifile}.con
  echo "Executable      = runana.sh" >> run$2_${ifile}.con
  echo "Arguments      = list/$2_${ifile} outHist/$2_${ifile}.root"  >> run$2_${ifile}.con
  echo "Output         = log/job$2_${ifile}.out" >> run$2_${ifile}.con
  echo "Error          = log/job$2_${ifile}.err" >> run$2_${ifile}.con
  echo "Log            = log/job$2_${ifile}.log" >> run$2_${ifile}.con
  echo "Queue"         >> run$2_${ifile}.con

  condor_submit run$2_${ifile}.con

  mv run$2_${ifile}.con script/
  ##################################################
  let "ifile+=1";
done
