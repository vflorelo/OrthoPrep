#!/bin/bash
work_dir=$1
species_list_file="${work_dir}/SpeciesIDs.txt"
species_num=$(cat ${species_list_file} | dos2unix | cut -d\: -f1)
for num in ${species_num}
do
    if [ -f "${work_dir}/Species${num}.fa" ]
    then
        cat "${work_dir}/Species${num}.fa" | dos2unix | perl -pe 'if(/\>/){s/$/\t/};s/\n//g;s/\>/\n/g' | tail -n+2 | sort -V | uniq | awk 'BEGIN{FS="\t"}{print $1 FS length($2)}' > "${work_dir}/Species${num}.sizes.tsv"
    fi
done