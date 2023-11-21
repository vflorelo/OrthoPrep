#!/bin/bash
fasta_file=$1
num_cores=$2
out_file=$3
tsv_file="${fasta_file}.tsv"
comp_file="${fasta_file}.COMPOSITION"
id_list=$(grep \> "${fasta_file}" | perl -pe 's/\ .*//;s/\>//' | sort -V | uniq)
CompositionMaker ${fasta_file}
fasta2tsv.sh ${fasta_file} > ${tsv_file}
echo "${id_list}" | awk -v tsv_file="${tsv_file}" -v comp_file="${comp_file}" '{print "run_flps.sh" , tsv_file , $1 , comp_file}' | parallel -j ${num_cores} 