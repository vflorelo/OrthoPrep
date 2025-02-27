#!/bin/bash
function get_lengths() {
    fasta_dir=$1
    fasta_file=$2
    base_name=$(echo ${fasta_file} | rev | cut -d\. -f1 --complement | rev)
    lengths=$(fasta2tsv.sh ${fasta_dir}/${fasta_file} | awk 'BEGIN{FS="\t"}{print length($2)}' | sort -n | grep -v ^$ | grep . | perl -pe 's/\n/\,/g' | perl -pe 's/\,$//')
    echo -e "${base_name}\t${lengths}"
    }
export -f get_lengths
function usage(){
    echo "Utility to calculate the protein median length, median absolute deviation and size of a collection of orthogroups"
    echo
    echo "Use this script to get a quick summary of the orthogroups built before and after running OrthoPrep"
    echo
    echo "Options:"
    echo "  --fasta_dir  -> Directory containing protein sequences in fasta format (mandatory)"
    echo "  --threads    -> Number of CPU threads to use"
    echo
    echo "Example:"
    echo "  get_og_stats.sh --fasta_dir /path/to/my/sequences/folder --threads 16"
    echo
	}
export -f usage
work_dir=$(pwd)
while [ "$1" != "" ]
do
    case $1 in
        --fasta_dir    )
            shift
            fasta_dir=$(realpath $1)
            ;;
        --threads )
            shift
            threads=$1
            ;;
		--help         )
            usage
            exit 0
            ;;
	esac
	shift
done
if [ ! -d "${fasta_dir}" ] || [ -z "${fasta_dir}" ]
then
    echo "Missing fasta directory. Exiting"
    exit 0
fi
num_proc=$(nproc)
if [ -z "${threads}" ]
then
    threads=${num_proc}
else
    threads_test=$(echo -e "${threads}\t${num_proc}" | awk '{if((int($1)==$1) && ($1<=$2)){print $1}}')
    if [ -z "${threads_test}" ]
    then
        echo "Invalid threads value. Exiting"
        exit
    fi
fi
fasta_list=$(ls "${fasta_dir}" | sort -V | uniq)
fasta_count=$(echo "${fasta_list}" | wc -l)
fasta_iters=$(echo ${fasta_count} | awk '{modulo=$1%1000}{if(modulo==0){print $1/1000}else{print (($1-modulo)/1000)+1 }}')
for iter in $(seq 1 ${fasta_iters})
do
    start_line=$(echo "${iter}" | awk '{print (($1-1)*1000)+1}')
    fasta_sub_list=$(echo "${fasta_list}" | tail -n+${start_line} | head -n1000)
    parallel -j ${threads} get_lengths ::: ${fasta_dir} ::: "${fasta_sub_list[@]}" >> ${work_dir}/og_lengths.tmp.tsv
done
echo -e "og_id\tlengths\n$(sort -V ${work_dir}/og_lengths.tmp.tsv | uniq | grep -v ^$ | grep .)" > ${work_dir}/og_lengths.tsv
rm ${work_dir}/og_lengths.tmp.tsv
get_og_stats.py ${work_dir}/og_lengths.tsv ${work_dir}/og_stats.tsv