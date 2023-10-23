#!/bin/bash
fasta_file=$1
size_list=$(perl -pe 'if(/\>/){s/$/\t/};s/\n//;s/\>/\n/g' "${fasta_file}" | tail -n+2 | awk 'BEGIN{FS="\t"}{print length($2)}')
min=$(echo "${size_list}" | sort -n  | head -n1)
max=$(echo "${size_list}" | sort -nr | head -n1)
med=$(echo "${size_list}" | sort -n | awk '{count[NR]=$1}END{if(NR%2){print count[(NR+1)/2]}else{print(count[(NR/2)]+count[(NR/2)+1])/2}}')
mad=$(echo "${size_list}" | sort -n | awk -v med="${med}" '{print sqrt(($1-med)^2)}' | sort -n | awk '{count[NR]=$1}END{if(NR%2){print count[(NR+1)/2]}else{print(count[(NR/2)]+count[(NR/2)+1])/2}}')
ran=$(echo $min $max | awk '{print $2-$1}')
echo -e "${fasta_file}\t${min}\t${max}\t${ran}\t${med}\t${mad}"
