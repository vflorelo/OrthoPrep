#!/bin/bash
cur_dir=$(pwd)
for dmnd in diamond diamond_def diamond_hard diamond_low diamond_med
do
  for inf in $(seq 1.00 0.01 1.10) $(seq 1.20 0.1 1.50)
  do
    sco_dir="sco_${dmnd}_inf_${inf}"
    cd ${sco_dir}
    sco_list=$(ls)
    sco_counts=$(echo "${sco_list}" | wc -l)
    echo "${sco_list}" | awk '{print "./get_og_len_diffs.sh "$1}' | parallel -j 24 | sort -V --parallel=24 | uniq > ${cur_dir}/${dmnd}.${inf}.sco_len_diffs.tsv
    echo "${sco_list}" | awk '{print "./get_og_mad.sh "$1}'       | parallel -j 24 | sort -V --parallel=24 | uniq > ${cur_dir}/${dmnd}.${inf}.sco_mad.tsv
    cd ${cur_dir}
  done
done
