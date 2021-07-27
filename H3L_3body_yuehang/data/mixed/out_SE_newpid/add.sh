#! /bin/bash

#############################
#add all of the root files under the directory#
############################

tmp="tmplist.txt"
out="../fout_H3L_data_SE_Jul25.root"
N=`ls *.root | wc -l`
rm -rf ${tmp}
ls *.root > ${tmp}
if [ -f "${out}" ]; then
    echo "${out} already exit"
    exit 0
fi
#  root -l -q  myhadd.C\(`${out}`, `${tmp}`, `${N}`\) 

list=${out}" "
for i in `cat ${tmp}`; do
  list=${list}" "${i}
  # echo ${list}
done
echo ${list}
hadd ${list}
