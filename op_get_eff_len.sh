#!/bin/bash
prep_dir=$1
threads=$2
masking=$3
eff_len=$4
run_flps() {
    species=$1
    seq_id=$2
    seq_tsv_file=$3
    comp_file=$4
    fLPS2 -t 1e-6 -m 5 -M 25 <(grep -w "${seq_id}" "${seq_tsv_file}" | perl -pe 's/\t/\n/') -z thorough -c ${comp_file} 2>> ${species}.flps.err | grep -wv WHOLE
    }
get_eff_len() {
    seq_id=$1
    len_tsv_file=$2
    lcr_bed_file=$3
    lcr_len=$(grep -w ^"${seq_id}" "${lcr_bed_file}"  | awk 'BEGIN{FS="\t"}{lcr_len+=($3-$2)}END{print lcr_len}')
    eff_len=$(awk -v seq_id="${seq_id}" -v lcr_len="${lcr_len}" 'BEGIN{FS="\t"}{if($2==seq_id){print $3-lcr_len}}' "${len_tsv_file}" )
    echo -e "${seq_id}\t${eff_len}"
    }
export -f run_flps
export -f get_eff_len
species_list_file="${prep_dir}/SpeciesIDs.txt"
sequence_full_len_file="${prep_dir}/Sequence_full_len.tsv"
sequence_eff_len_file="${prep_dir}/Sequence_eff_len.tsv"
lcr_dir="${prep_dir}/LCR"
cat ${prep_dir}/Species*.fa \
    | dos2unix \
    | perl -pe 'if(/\>/){s/$/\t/;s/_/\t/};s/\n//g;s/\>/\n/g' \
    | tail -n+2 \
    | sort -V \
    | uniq \
    | awk 'BEGIN{FS="\t"}{print $1  FS $1"_"$2 FS length($3)}' > ${sequence_full_len_file}
if [ "${eff_len}" == "FALSE" ] || [ "${eff_len}" == "false" ]
then
    cat ${sequence_full_len_file} > "${sequence_eff_len_file}"
    exit 0
elif [ "${eff_len}" == "TRUE" ] || [ "${eff_len}" == "true" ]
then
    species_num=$(cat ${species_list_file} | dos2unix | cut -d\: -f1)
    cd ${lcr_dir}
    for num in ${species_num}
    do
        len_tsv_file="Species${num}.len.tsv"
        cp "${prep_dir}/Species${num}.fa" "${lcr_dir}"
        grep -w ^${num} ${sequence_full_len_file} > ${len_tsv_file}
        CompositionMaker Species${num}.fa
        perl -pe 'if(/\>/){s/$/\t/};s/\n//g;s/\>/\n\>/g' Species${num}.fa | tail -n+2 > Species${num}.seq.tsv
        seq_id_list=$(cut -f2 "${len_tsv_file}" | sort -V | uniq)
        seq_id_count=$(echo "${seq_id_list}" | wc -l)
        seq_id_iters=$(echo ${seq_id_count} | awk '{modulo=$1%1000}{if(modulo==0){print $1/1000}else{print (($1-modulo)/1000)+1 }}')
        echo > Species${num}.flps.tsv
        echo > Species${num}.lcr.tsv
        for iter in $(seq 1 ${seq_id_iters})
        do
            start_line=$(echo "${iter}" | awk '{print (($1-1)*1000)+1}')
            seq_id_sub_list=$(echo "${seq_id_list}" | tail -n+${start_line} | head -n1000)
            parallel -j ${threads} run_flps ::: ${num} ::: "${seq_id_sub_list[@]}" ::: Species${num}.seq.tsv ::: Species${num}.fa.COMPOSITION >> Species${num}.flps.tsv
        done
        echo "$(sort -V Species${num}.flps.tsv | uniq | grep -v ^$)" > Species${num}.flps.tsv
        awk 'BEGIN{FS="\t";OFS="\t"}{lcr_len=($6-($5-1))}{if(lcr_len>=10){print $1,$5-1,$6,$9,$8,"+"}}' Species${num}.flps.tsv | sortBed -i - > Species${num}.flps.bed
        seg Species${num}.fa -l \
            | grep \> \
            | perl -pe 's/>//;s/ .*complexity=/\t/;s/ .*//;s/\(/\t/;s/\)//' \
            | awk 'BEGIN{FS="\t";OFS="\t"}{gsub("-","\t",$2)}1' \
            | awk 'BEGIN{FS="\t";OFS="\t"}{lcr_len=$3-($2-1)}{if(lcr_len>=10){print $1,$2-1,$3,"seg",$4,"+"}}' \
            | sortBed -i - > Species${num}.seg.bed
        intersectBed -a Species${num}.flps.bed -b Species${num}.seg.bed | sortBed -i - | mergeBed > Species${num}.lcr.bed
        for iter in $(seq 1 ${seq_id_iters})
        do
            start_line=$(echo "${iter}" | awk '{print (($1-1)*1000)+1}')
            seq_id_sub_list=$(echo "${seq_id_list}" | tail -n+${start_line} | head -n1000)
            parallel -j ${threads} get_eff_len ::: "${seq_id_sub_list[@]}" ::: "${len_tsv_file}" ::: "Species${num}.lcr.bed" >> Species${num}.lcr.tsv
        done
        echo "$(sort -V Species${num}.lcr.tsv | uniq | grep -v ^$)" > Species${num}.lcr.tsv
        if [ "${masking}" == "TRUE" ] || [ "${masking}" == "true" ]
        then
            maskFastaFromBed -fi Species${num}.fa -bed Species${num}.lcr.bed -fo ${prep_dir}/Species${num}.fa -mc X
        fi
        rm Species${num}.seq.tsv Species${num}.len.tsv
        perl -pi -e "s/^/${num}\t/" Species${num}.lcr.tsv
    done
    cat Species*.lcr.tsv > "${sequence_eff_len_file}"
fi