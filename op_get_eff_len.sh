#!/bin/bash
sizes_tsv_file=$1
lcr_bed_file=$2
threads=$3
uuid=$(uuidgen | cut -d- -f5)
get_eff_len() {
    prot_id=$1
    tsv_file=$2
    bed_file=$3
    lcr_len=$(grep -w ^"${prot_id}" "${bed_file}"  | awk 'BEGIN{FS="\t"}{lcr_len+=($3-$2)}END{print lcr_len}')
    prot_len=$(grep -w ^"${prot_id}" "${tsv_file}" | cut -f2)
    echo -e "${prot_id}\t${prot_len}\t${lcr_len}"  | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2-$3}'
    }
export -f get_eff_len
prot_list=$(cut -f1 ${sizes_tsv_file})
parallel -j ${threads} get_eff_len ::: "${prot_list[@]}" ::: ${sizes_tsv_file} ::: ${lcr_bed_file} | sort -V | uniq > ${uuid}.tsv
mv "${uuid}.tsv" "${sizes_tsv_file}"