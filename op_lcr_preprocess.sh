#!/bin/bash
work_dir=$1
cur_dir=$2
threads=$3
no_masking=$4
species_list_file="${work_dir}/SpeciesIDs.txt"
species_num=$(cat ${species_list_file} | dos2unix | cut -d\: -f1)
uuid=$(uuidgen | cut -d\- -f5)
tmp_dir="${cur_dir}/tmp_${uuid}"
mkdir "${tmp_dir}"
cd ${tmp_dir}
for num in ${species_num}
do
    if [ -f "${work_dir}/Species${num}.fa" ]
    then
        sizes_tsv_file="${work_dir}/Species${num}.sizes.tsv"
        cp "${work_dir}/Species${num}.fa" "${tmp_dir}"
        run_flps_parallel.sh "Species${num}.fa" "${threads}" "Species${num}.flps.tsv" "Species${num}.flps.err"
        awk 'BEGIN{FS="\t";OFS="\t"}{lcr_len=($6-($5-1))}{if(lcr_len>=10){print $1,$5-1,$6,$9,$8,"+"}}' Species${num}.flps.tsv | sortBed -i - > Species${num}.flps.bed
        seg Species${num}.fa -l | grep \> | perl -pe 's/>//;s/ .*complexity=/\t/;s/ .*//;s/\(/\t/;s/\)//' | awk 'BEGIN{FS="\t";OFS="\t"}{gsub("-","\t",$2)}1' | awk 'BEGIN{FS="\t";OFS="\t"}{lcr_len=$3-($2-1)}{if(lcr_len>=10){print $1,$2-1,$3,"seg",$4,"+"}}' | sortBed -i - > Species${num}.seg.bed
        intersectBed -a Species${num}.flps.bed -b Species${num}.seg.bed | mergeBed > Species${num}.lcr.bed
        op_get_eff_len.sh "${sizes_tsv_file}" "Species${num}.lcr.bed" "${threads}"
        if [ "${no_masking}" == "FALSE" ]
        then
            maskFastaFromBed -fi Species${num}.fa -bed Species${num}.lcr.bed -fo ${work_dir}/Species${num}.fa -mc X
        fi
    fi
done
#rm -rf ${tmp_dir}