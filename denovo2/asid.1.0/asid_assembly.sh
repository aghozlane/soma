#!/bin/bash
# SOLiD(TM) Reads De Novo Assembly Project
# Copyright (c) 2010 Life Technologies
# Iterating through all gaps and calling assembly for each of them
tmp_dir="$1";
hsize=$2;
min_cont=$3;
exp_cov=$4;
cov_cutoff=$5;
numcores=$6
suffix=$7
FILES="${tmp_dir}/reads_*.de";
exeVELVET="$denovo2/velvet_0.7.55";
idx=0
idx2=0
numc=`expr ${numcores} - 1`;
rm -f ${tmp_dir}/multicore.flag;
for f in $FILES
do
  echo "Assembling $f file...";
  $denovo2/asid.1.0/asid_assembly_one.sh "${tmp_dir}/${idx}" $f ${hsize} ${min_cont} ${exp_cov} ${cov_cutoff} ${exeVELVET} ${suffix} "${tmp_dir}/multicore.flag" &
idx=`expr ${idx} + 1`;
idx2=`expr ${idx} % ${numcores}`;
if [ ${idx2} = 0 ] 
then 
xv=0;
while [ $xv -le ${numc} ]
do
  sleep 1;
  xv=`grep -c x ${tmp_dir}/multicore.flag`;
done
rm -f ${tmp_dir}/multicore.flag;
fi
done
echo ${idx2};
if [ ${idx2} -gt 0 ] 
then 
xv=0;
while [ $xv -lt ${idx2} ]
do
  sleep 1;
  xv=`grep -c x ${tmp_dir}/multicore.flag`;
done
rm -f ${tmp_dir}/multicore.flag;
fi
