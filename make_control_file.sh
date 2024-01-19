#!/bin/bash
function usage(){
    echo "OrthoPrep takes a folder containing fasta files with protein sequences to be analysed by OrthoFinder 2+"
    echo
    echo "Use this script to prepare a control file to apply different filters to specific comparisons"
    echo
    echo "Options:"
    echo "  --fasta_dir    -> Directory containing protein sequences in fasta format (mandatory)"
    echo "  --same_species -> Short fraction value for same-species comparisons [0.30]"
    echo "  --same_genus   -> Short fraction value for within genus comparisons [0.35]"
    echo "  --diff_genus   -> Short fraction value for inter-genera comparisons [0.40]"
    echo "  --decrement    -> Amount to be subtracted to the short fraction to get long fraction values [0.05]"
    echo
    echo "Example:"
    echo "  make_control_file.sh --fasta_dir /path/to/my/sequences/folder --same_species 0.30 --same_genus 0.35 --diff_genus 0.4 --decrement 0.05"
    echo
    echo "Notes:"
    echo "  - For the script to work, fasta files should have the following structure:"
    echo "    Genus_species_strain.fasta"
	}
export -f usage
same_species="0.30"
same_genus="0.35"
diff_genus="0.40"
decrement="0.05"
while [ "$1" != "" ]
do
    case $1 in
        --fasta_dir    )
            shift
            fasta_dir=$(realpath $1)
            ;;
        --same_species )
            shift
            same_species=$1
            ;;
        --same_genus   )
            shift
            same_genus=$1
            ;;
        --diff_genus   )
            shift
            diff_genus=$1
            ;;
        --decrement    )
            shift
            decrement=$1
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
q_list=$(ls ${fasta_dir} | sort -V | uniq)
start_line=0
for query in ${q_list}
do
    q_genus=$(echo   "${query}" | cut -d_ -f1)
    q_species=$(echo "${query}" | cut -d_ -f2)
    start_line=$(echo ${start_line} | awk '{print $1+1}')
    s_list=$(echo "${q_list}" | tail -n+${start_line})
    for subject in ${s_list}
    do
        s_genus=$(echo   "${subject}" | cut -d_ -f1)
        s_species=$(echo "${subject}" | cut -d_ -f2)
        if [ "${query}" == "${subject}" ]
        then
            sf_value=${same_species}
        elif [ "${q_genus}" == "${s_genus}" ] && [ "${q_species}" == "${s_species}" ]
        then
            sf_value=${same_species}
        elif [ "${q_genus}" == "${s_genus}" ] && [ "${q_species}" != "${s_species}" ]
        then
            sf_value=${same_genus}
        elif [ "${q_genus}" != "${s_genus}" ]
        then
            sf_value=${diff_genus}
        fi
        lf_value=$(echo "${sf_value}" | awk -v decrement="${decrement}" '{print $1-decrement}')
        echo -e "${query}\t${subject}\t${sf_value}\t${lf_value}"
    done
done