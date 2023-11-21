#!/bin/bash
work_dir=$1
cur_dir=$(pwd)
bckp_dir="${cur_dir}/bckp"
species_list_file="${work_dir}/SpeciesIDs.txt"
species_num=$(cat ${species_list_file} | dos2unix | cut -d\: -f1)
for num in ${species_num}
do
    if [ -f "${work_dir}/Species${num}.fa" ]
    then
        cp "${work_dir}/Species${num}.fa" "${bckp_dir}"
    fi
done