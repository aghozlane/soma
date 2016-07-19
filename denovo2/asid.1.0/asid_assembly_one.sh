#!/bin/bash
# SOLiD(TM) Reads De Novo Assembly Project
# Copyright (c) 2010 Life Technologies
# Single Velvet call to assemble reads from one gap
tmp_dir=$1;
file=$2;
hsize=$3;
min_cont=$4;
exp_cov=$5;
cov_cutoff=$6;
exeVELVET="$7";
suffix="$8";
multicore_flag=$9;
  ${exeVELVET}/velveth_de ${tmp_dir} ${hsize} -fasta -short $file;
  ${exeVELVET}/velvetg_de ${tmp_dir} -min_contig_lgth ${min_cont} -exp_cov ${exp_cov} -cov_cutoff ${cov_cutoff} -read_trkg yes;
  mv -f ${tmp_dir}/LastGraph ${file}.${suffix}.graph;
  rm -r -f ${tmp_dir};
  echo "x" >> ${multicore_flag};



