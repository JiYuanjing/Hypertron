#!/bin/bash
date

if [ $# -lt 2 ]; then
echo "not enough arguements"
echo "./submitAll.sh <outDir> <xml>"
exit
fi

dir=$(pwd)
parentdir="$(dirname "$dir")"
tmp=${dir}/Tmp

mkdir -p ${tmp}
ver=$1
xml=$2
echo $dir
echo "do clean"
rm -rf $tmp/log_$1
rm -rf $tmp/submit/$1
rm -rf ${dir}/out/out_$1 
echo "rm -rf ${dir}/out/out_$1"

if [ ! -d ${dir}/out/out_$1 ]; then
mkdir -p ${dir}/out/out_$1
fi

if [ ! -d $tmp/log_$1 ]; then
mkdir -p $tmp/log_$1
fi

if [ ! -d $tmp/submit/$1 ]; then
mkdir -p $tmp/submit/$1
fi

cp $dir/${xml} $tmp/submit/$1/.
cd $tmp/submit/$1/

FILE=production_3p85GeV_fixedTarget_2018.txt

star-submit-template -template ${xml} -entities path=${parentdir},Tmp=$tmp,ver=$1,listOfFiles=${FILE}
