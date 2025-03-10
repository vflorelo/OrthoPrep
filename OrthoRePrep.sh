#!/bin/bash
usage="$(dirname $0)/reusage.sh"
source ${usage}
run_mode="normal"
while [ "$1" != "" ]
do
    case $1 in
        --prep_dir     )
            shift
            prep_dir=$(realpath $1)
            ;;
        --search_method )
            shift
            search_method=$1
            ;;
		--short_frac    )
            shift
            short_frac=$1
            ;;
		--long_frac     )
            shift
            long_frac=$1
            ;;
		--threads       )
            shift
            threads=$1
            ;;
        --control_file  )
            shift
            control_file=$1
            ;;
		--help          )
            usage
            exit 0
            ;;
	esac
	shift
done
if [ ! -d "${prep_dir}" ] || [ -z "${prep_dir}" ]
then
    echo "Missing OrthoPrep directory. Exiting"
    exit 0
fi
if [ -z "${search_method}" ]
then
    echo "Missing search method. Exiting"
    exit 0
fi
if [ -z "${control_file}" ]
then
    echo "No control file specified, proceeding in normal mode"
    run_mode="normal"
    if [ -z "${short_frac}" ] || [ -z "${long_frac}" ]
    then
        echo "Missing fraction values. Exiting"
        exit 0
    else
        frac_test=$(echo -e "${short_frac}\t${long_frac}" | awk 'BEGIN{FS="\t"}{if($1>=0 && $1<=1 && $2>=0 && $2<=1){print "pass"}else{print "fail"}}')
        if [ "${frac_test}" == "fail" ]
        then
            echo "Invalid fraction values. Exiting"
            exit 0
        fi
    fi
elif [ -f "${control_file}" ]
then
    if [ ! -z "${short_frac}" ] || [ ! -z "${long_frac}" ]
    then
        echo "Incompatible options. Exiting"
        exit 0
    else
        echo "Using ${control_file} for specific comparisons"
        run_mode="custom"
        frac_list=$(cut -f3,4 "${control_file}")
        frac_count=$(echo "${frac_list}" | awk 'BEGIN{FS="\t"}{if($1>=0 && $1<=1 && $2>=0 && $2<=1){print $0 FS NR}}' | wc -l)
        comp_num=$(cat "${control_file}" | wc -l)
        if [ ! "${comp_num}" -eq "${frac_count}" ]
        then
            echo "Invalid fraction values. Exiting"
            exit 0
        fi
    fi
elif [ ! -f "${control_file}" ]
then
    echo "Control file ${control_file} missing. Exiting"
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

check_deps=$(check_deps.sh)
dep_test=$(echo "${check_deps}" | grep -wc 0)
if [ "${dep_test}" -gt 0 ]
then
    echo "${check_deps}" | grep -w 0 | cut -f2
    echo "Exiting"
    exit 0
else
    orthofinder_cmd=$(echo "${check_deps}" | grep orthofinder | cut -f3)
fi
config_file=$(find $(dirname $(which ${orthofinder_cmd} )) -type f -name config.json | grep config.json$ )
makedb_command=$(cat ${config_file} | jq .${search_method}.db_cmd | sed -e 's/"//g')
if [ -z "${makedb_command}" ] || [ "${makedb_command}" == "null" ]
then
    echo "Invalid database command, Exiting"
    exit
fi
search_command=$(cat ${config_file} | jq .${search_method}.search_cmd | sed -e 's/"//g')
if [ -z "${search_command}" ] || [ "${search_command}" == "null" ]
then
    echo "Invalid search command, Exiting" 
    exit
fi
cur_date=$(date +%h%d)
cur_dir=$(pwd)
uuid=$(uuidgen | cut -d- -f5 )
prep_date=$(date +%Y-%m-%d)
tmp_dir="${prep_dir}/tmp"
log_dir="${prep_dir}/logs"
bckp_dir="${prep_dir}/resume_bckp"
log_file="OrthoPrep-${prep_date}.log"
log_file="${log_dir}/${log_file}"
mkdir -p "${prep_dir}" "${tmp_dir}" "${log_dir}" "${bckp_dir}"
if [ ! -f "${prep_dir}/SpeciesIDs.txt" ] || [ ! -d "${prep_dir}/masked_seqs" ] || [ ! -f "${prep_dir}/Sequence_eff_len.tsv" ]
then
    echo "Missing required files. Exiting" | tee -a ${log_file}
    exit
else
    dos2unix ${prep_dir}/SpeciesIDs.txt
    file_list=$(cut -d\: -f1 ${prep_dir}/SpeciesIDs.txt | grep -v "#" | sed -e "s|^|${prep_dir}/masked_seqs/Species|;s|$|.fa|")
    file_test=$(file $(echo "${file_list}") | grep -wc "No such file or directory" )
    if [ "${file_test}" -gt 0 ]
    then
        echo "Missing required files. Exiting" | tee -a ${log_file}
        out_test=$(file $(echo "${file_list}"))
        echo "${out_test}"
        exit
    fi
fi
echo "Running resume_OrthoPrep with the following options" | tee -a ${log_file}
echo -e "OrthoPrep directory\t->\t${prep_dir}"      | tee -a ${log_file}
echo -e "BLAST search method\t->\t${search_method}" | tee -a ${log_file}
if [ "${run_mode}" == "normal" ]
then
    echo -e "Shortest seq. fraction\t->\t${short_frac}" | tee -a ${log_file}
    echo -e "Longest seq. fraction\t->\t${long_frac}"   | tee -a ${log_file}
elif [ "${run_mode}" == "custom" ]
then
    echo -e "Custom control file\t->\t${control_file}"  | tee -a ${log_file}
fi
echo -e "Number of threads\t->\t${threads}"         | tee -a ${log_file}
echo "#############################################################"
species_list=$(cut -d\: -f1 ${prep_dir}/SpeciesIDs.txt | grep -v "#")
db_cmd_list=""
search_cmd_list=""
for query in ${species_list}
do
    if [ -f "${prep_dir}/${search_method}DBSpecies${query}.dmnd" ]
    then
        mv "${prep_dir}/${search_method}DBSpecies${query}.dmnd" "${bckp_dir}"
    fi
    cur_makedb_cmd=$(echo "${makedb_command}" | sed -e "s|INPUT|${prep_dir}/masked_seqs/Species${query}.fa|;s|OUTPUT|${prep_dir}/${search_method}DBSpecies${query}.dmnd|")
    db_cmd_list=$(echo -e "${db_cmd_list}\n${cur_makedb_cmd}")
    for subject in ${species_list}
    do
        cur_search_cmd=$(echo "${search_command}" | sed -e "s|DATABASE|${prep_dir}/${search_method}DBSpecies${query}.dmnd|;s|INPUT|${prep_dir}/masked_seqs/Species${subject}.fa|;s|OUTPUT|${prep_dir}/Blast${subject}_${query}.txt|")
        search_cmd_list=$(echo -e "${search_cmd_list}\n${cur_search_cmd}")
    done
done
echo "Step 1. Rebuilding diamond databases" | tee -a ${log_file}
db_cmd_list=$(echo "${db_cmd_list}" | sort -V | uniq | grep -v ^$ | grep .)
echo "${db_cmd_list}" > "${log_dir}/db_cmd_list.txt"
echo "${db_cmd_list}" | parallel -j ${threads} > ${log_dir}/makedb_commands.log 2> ${log_dir}/makedb_commands.err
if [ $? -gt 0 ]
then
    echo "Error rebuilding diamond databases. Restoring files" | tee -a ${log_file}
    rm ${prep_dir}/*.dmnd
    mv ${bckp_dir}/*.dmnd ${prep_dir}
    exit
fi
search_cmd_list=$(echo "${search_cmd_list}" | sort -V | uniq | grep -v ^$ | grep .)
echo "${search_cmd_list}" > "${log_dir}/search_cmd_list.txt"
num_commands=$(echo "${search_cmd_list}" | wc -l)
echo "Step 2. Running ${num_commands} diamond commands in parallel" | tee -a ${log_file}
mv "${prep_dir}"/Blast*.txt.gz "${bckp_dir}"
echo "${search_cmd_list}" | parallel -j ${threads} > ${log_dir}/search_commands.log 2> ${log_dir}/search_commands.err
if [ $? -gt 0 ]
then
    echo "Error running diamond searches. Restoring files" | tee -a ${log_file}
    rm ${prep_dir}/Blast*.gz
    mv ${bckp_dir}/Blast*.gz ${prep_dir}
    exit
fi
echo "Step 3. Filtering BLAST results based on size differences" | tee -a ${log_file}
eff_len_file="Sequence_eff_len.tsv"
blast_filter_cmd_list=""
for query in ${species_list}
do
    for subject in ${species_list}
    do
        if [ "${run_mode}" == "custom" ]
        then
            q_fasta=$(awk -v query="${query}"     'BEGIN{FS=": "}{if($1==query)  {print $2}}' ${prep_dir}/SpeciesIDs.txt)
            s_fasta=$(awk -v subject="${subject}" 'BEGIN{FS=": "}{if($1==subject){print $2}}' ${prep_dir}/SpeciesIDs.txt)
            short_frac=$(awk -v q="${q_fasta}" -v s="${s_fasta}" '{if((($1==q) && ($2==s)) || (($2==q)&&($1==s))){print $3}}' ${control_file})
            long_frac=$(awk  -v q="${q_fasta}" -v s="${s_fasta}" '{if((($1==q) && ($2==s)) || (($2==q)&&($1==s))){print $4}}' ${control_file})
        fi
        cur_bf_cmd="op_blast_filter.py ${prep_dir} ${query} ${subject} ${eff_len_file} ${short_frac} ${long_frac} > ${log_dir}/blast_filter_${query}_${subject}.log 2> ${log_dir}/blast_filter_${query}_${subject}.err"
        blast_filter_cmd_list=$(echo -e "${blast_filter_cmd_list}\n${cur_bf_cmd}")
    done
done
blast_filter_cmd_list=$(echo "${blast_filter_cmd_list}" | sort -V | uniq | grep -v ^$ | grep .)
echo "${blast_filter_cmd_list}" | parallel -j ${threads} > ${log_dir}/blast_filter.log 2> ${log_dir}/blast_filter.err
if [ $? -eq 0 ]
then
    mv ${prep_dir}/Blast*.txt.gz ${bckp_dir}
    mv ${tmp_dir}/Blast*.txt.gz ${prep_dir}
    echo ""  | tee -a ${log_file}
    echo "        Finished filtering BLAST results" | tee -a ${log_file}
else
    echo "Error filtering BLAST results"
    echo "Aborting" | tee -a ${log_file}
    exit 0
    echo "Check ${log_dir}/blast_filter.log and ${log_dir}/blast_filter.err"
fi
echo "        Removing temporary files" | tee -a ${log_file}
rm -rf ${tmp_dir}
echo | tee -a ${log_file}
echo "Now you can resume OrthoFinder by running:" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "orthofinder.py -b ${prep_dir} [other OrthoFinder options]" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "OrthoPrep v0.0.1" | tee -a ${log_file}