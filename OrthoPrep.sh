#!/bin/bash
usage="$(dirname $0)/usage.sh"
source ${usage}
no_masking="FALSE"
no_eff_len="FALSE"
run_mode="normal"
while [ "$1" != "" ]
do
    case $1 in
        --fasta_dir     )
            shift
            fasta_dir=$(realpath $1)
            ;;
        --search_method )
            shift
            search_method=$1
            ;;
        --no_masking    )
            no_masking="TRUE"
            ;;
        --no_eff_len    )
            no_eff_len="TRUE"
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
if [ ! -d "${fasta_dir}" ] || [ -z "${fasta_dir}" ]
then
    echo "Missing fasta directory. Exiting"
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
        fasta_file_list=$(cut -f1,2 "${control_file}" | perl -pe 's/\t/\n/' | sort -V | uniq | grep -v ^$)
        for fasta_file in ${fasta_file_list}
        do
            if [ ! -s "${fasta_dir}/${fasta_file}" ]
            then
                echo "Missing file ${fasta_file}. Exiting"
                exit 0
            fi
        done
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
fi
cur_date=$(date +%h%d)
cur_dir=$(pwd)
uuid=$(uuidgen | cut -d- -f5 )
prep_date=$(date +%y-%m-%d)
of_dir="${fasta_dir}/OrthoFinder/Results_${cur_date}/WorkingDirectory"
of_tmp_dir="${fasta_dir}/${uuid}"
prep_dir="${cur_dir}/OrthoPrep-${prep_date}"
bckp_dir="${prep_dir}/bckp"
lcr_dir="${prep_dir}/LCR"
tmp_dir="${prep_dir}/tmp"
log_dir="${prep_dir}/logs"
log_file="OrthoPrep-${prep_date}.log"
log_file="${log_dir}/${log_file}"
mkdir -p "${prep_dir}" "${tmp_dir}" "${bckp_dir}" "${lcr_dir}" "${log_dir}"
if [ -d "${fasta_dir}/OrthoFinder" ]
then
    echo "The folder ${fasta_dir}/OrthoFinder will be temporarily renamed to ${of_tmp_dir}"
    mv ${fasta_dir}/OrthoFinder ${of_tmp_dir}
    of_restore="TRUE"
else
    of_restore="FALSE"
fi
echo "Step 1. Preparing fasta files, and diamond commands" | tee -a ${log_file}
command_list=$(orthofinder.py -S ${search_method} -op -f ${fasta_dir} | grep -w ^diamond | grep blastp)
command_list=$(echo "${command_list}" | sed -e "s|${of_dir}|${prep_dir}|g")
num_commands=$(echo "${command_list}" | wc -l)
echo "${command_list}" > ${log_dir}/diamond.log 2> ${log_dir}/diamond.err
cp  "${of_dir}"/Species*.fa \
    "${of_dir}"/SequenceIDs.txt \
    "${of_dir}"/SpeciesIDs.txt \
    "${of_dir}"/${search_method}DBSpecies*.dmnd \
    "${bckp_dir}"
cp  "${of_dir}"/Species*.fa \
    "${of_dir}"/SequenceIDs.txt \
    "${of_dir}"/SpeciesIDs.txt \
    "${prep_dir}"
if [ $? -eq 0 ]
then
    rm -rf ${fasta_dir}/OrthoFinder
    if [ "${of_restore}" == "TRUE" ]
    then
        mv ${of_tmp_dir} ${fasta_dir}/OrthoFinder
    fi
fi
species_count=$(cat "${prep_dir}/SpeciesIDs.txt" | wc -l)
fasta_count=$(ls "${prep_dir}" | grep -c ".fa"$ )
if [ "${species_count}" -eq "${fasta_count}" ]
then
    echo "Step 2. Preprocessing sequence files" | tee -a ${log_file}
    op_get_eff_len.sh ${prep_dir} ${threads} ${no_masking} ${no_eff_len}
    if [ ! $? -eq 0 ]
    then
        echo "Error produced extracting LCRs" | tee -a ${log_file}
        echo "Aborting" | tee -a ${log_file}
        exit 0
    fi
    if [ "${no_masking}" == "FALSE" ]
    then
        for base_name in $(ls ${prep_dir} | grep fa$ | cut -d\. -f1)
        do
            diamond \
                makedb \
                --threads ${threads} \
                --db ${prep_dir}/${search_method}DB${base_name}.dmnd \
                --in ${prep_dir}/${base_name}.fa >> ${log_dir}/diamond.log 2>> ${log_dir}/diamond.err
        done
    elif [ "${no_masking}" == "TRUE" ]
    then
        cp ${bckp_dir}/*.dmnd ${prep_dir}
    fi
fi

echo "Step 3. Running ${num_commands} diamond commands in parallel" | tee -a ${log_file}
echo "${command_list}" | parallel -j ${threads}
echo "Step 4. Filtering BLAST results based on size differences" | tee -a ${log_file}
species_num=$(cat ${prep_dir}/SpeciesIDs.txt | dos2unix | cut -d\: -f1)
eff_len_file="Sequence_eff_len.tsv"
for query in ${species_num}
do
    for subject in ${species_num}
    do
        if [ "${run_mode}" == "custom" ]
        then
            q_fasta=$(awk -v query="${query}"     'BEGIN{FS=": "}{if($1==query)  {print $2}}' ${prep_dir}/SpeciesIDs.txt)
            s_fasta=$(awk -v subject="${subject}" 'BEGIN{FS=": "}{if($1==subject){print $2}}' ${prep_dir}/SpeciesIDs.txt)
            short_frac=$(awk -v q="${q_fasta}" -v s="${s_fasta}" '{if((($1==q) && ($2==s)) || (($2==q)&&($1==s))){print $3}}' ${control_file})
            long_frac=$(awk  -v q="${q_fasta}" -v s="${s_fasta}" '{if((($1==q) && ($2==s)) || (($2==q)&&($1==s))){print $4}}' ${control_file})
        fi
        echo "op_blast_filter.py ${prep_dir} ${query} ${subject} ${eff_len_file} ${short_frac} ${long_frac}" >> ${log_dir}/blast_filter.log 2>> ${log_dir}/blast_filter.err
        op_blast_filter.py ${prep_dir} ${query} ${subject} ${eff_len_file} ${short_frac} ${long_frac} >> ${log_dir}/blast_filter.log 2>> ${log_dir}/blast_filter.err
        if [ $? -eq 0 ]
        then
            mv ${prep_dir}/Blast${query}_${subject}.txt.gz ${bckp_dir}
            mv ${tmp_dir}/Blast${query}_${subject}.txt.gz ${prep_dir}
        else
            echo "Error filtering Blast${query}_${subject}.txt.gz."
            echo "Check ${log_dir}/blast_filter.log and ${log_dir}/blast_filter.err"
        fi
    done
done
if [ ! $? -eq 0 ]
then
    echo "Something went wrong while filtering BLAST results" | tee -a ${log_file}
    echo "Aborting" | tee -a ${log_file}
    exit 0
fi
echo "        Removing temporary files" | tee -a ${log_file}
rm -rf ${tmp_dir}
echo ""  | tee -a ${log_file}
echo "Finished filtering BLAST results" | tee -a ${log_file}
echo "Files in the ${prep_dir} directory were obtained by applying strict filters to the BLAST results" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "Now you can resume OrthoFinder by running:" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "orthofinder.py -b ${prep_dir} [other OrthoFinder options]" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "OrthoPrep v0.0.1" | tee -a ${log_file}