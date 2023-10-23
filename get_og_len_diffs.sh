#!/bin/bash
fasta_file=$1
size_datablock=$(perl -pe 'if(/\>/){s/$/\t/};s/\n//;s/\>/\n/g' "${fasta_file}" | tail -n+2 | awk 'BEGIN{FS="\t"}{print $1 FS length($2)}')
id_list=$(echo  "${size_datablock}" | cut -f1 | sort -V | uniq | head -n-1)
sub_list=$(echo "${id_list}")
for query in ${id_list}
do
    sub_list=$(echo "${sub_list}" | grep -wv "${query}")
    for subject in ${sub_list}
    do
        qlen=$(echo "${size_datablock}" | grep -w ^${query}   | cut -f2)
        slen=$(echo "${size_datablock}" | grep -w ^${subject} | cut -f2)
        if [ "${qlen}" -ge "${slen}" ]
        then
            len_diff=$(echo -e "${qlen}\t${slen}" | awk '{print $1-$2}')
            short_diff_frac=$(echo -e "${len_diff}\t${slen}" | awk '{print $1/$2}')
            long_diff_frac=$(echo  -e "${len_diff}\t${qlen}" | awk '{print $1/$2}')
        else
            len_diff=$(echo -e "${slen}\t${qlen}" | awk '{print $1-$2}')
            short_diff_frac=$(echo -e "${len_diff}\t${qlen}" | awk '{print $1/$2}')
            long_diff_frac=$(echo  -e "${len_diff}\t${slen}" | awk '{print $1/$2}')
        fi
        echo -e "${query}\t${subject}\t${len_diff}\t${short_diff_frac}\t${long_diff_frac}"
    done
done