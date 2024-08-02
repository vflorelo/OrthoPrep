#!/bin/bash
get_lengths() {
    fasta_dir=$1
    fasta_file=$2
    base_name=$(echo ${fasta_file} | rev | cut -d\. -f1 --complement | rev)
    lengths=$(fasta2tsv.sh ${fasta_dir}/${fasta_file} | awk 'BEGIN{FS="\t"}{print length($2)}' | sort -n | grep -v ^$ | grep . | perl -pe 's/\n/\,/g' | perl -pe 's/\,$//')
    echo -e "${base_name}\t${lengths}"
    }
export -f get_lengths
work_dir=$(pwd)
fasta_dir=$(realpath $1)
threads=$2
fasta_list=$(ls ${fasta_dir} | sort -V | uniq)
fasta_count=$(echo "${fasta_list}" | wc -l)
fasta_iters=$(echo ${fasta_count} | awk '{modulo=$1%1000}{if(modulo==0){print $1/1000}else{print (($1-modulo)/1000)+1 }}')
for iter in $(seq 1 ${fasta_iters})
do
    start_line=$(echo "${iter}" | awk '{print (($1-1)*1000)+1}')
    fasta_sub_list=$(echo "${fasta_list}" | tail -n+${start_line} | head -n1000)
    parallel -j ${threads} get_lengths ::: ${fasta_dir} ::: "${fasta_sub_list[@]}" >> ${cur_dir}/og_lengths.tmp.tsv
done
echo -e "og_id\tlengths\n$(sort -V ${cur_dir}/og_lengths.tmp.tsv | uniq | grep -v ^$ | grep .)" > ${work_dir}/og_lengths.tsv
rm ${cur_dir}/og_lengths.tmp.tsv
get_og_stats.py ${work_dir}/og_lengths.tsv ${work_dir}/og_stats.tsv