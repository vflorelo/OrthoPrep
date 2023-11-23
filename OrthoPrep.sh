#!/bin/bash
function usage(){
    echo "OrthoPrep takes a folder containing fasta files with protein sequences to be analysed by OrthoFinder 2+"
    echo
    echo "The first step uses OrthoFinder to prepare the files in the expected format"
    echo "In the second step, the effective lengths of proteins are calculated"
    echo "The effective lengths are determined by detecting low complexity regions (LCR)"
    echo "Optional step, the sequences are then masked so that BLAST searches do not include LCRs"
    echo "The third step uses gnu-parallel to run the diamond commands prepared by OrthoFinder"
    echo "The fourth step filters the BLAST results based on the differences in effective length of each protein pair"
    echo "Optional step, if the sequences were masked, they are returned to the unmasked state"
    echo
    echo "The size difference is divided by a fraction of the length of each sequence and then compared against a manually set threshold"
    echo " - Pairs of sequences with a size difference lower than the set threshold are kept as comparable sequences"
    echo " - Pairs of sequences with a size difference greater than the set threshold are excluded"
    echo
    echo "Options:"
    echo "  --fasta_dir     -> Directory containing protein sequences in fasta format (mandatory)"
    echo "  --search_method -> Sequence similarity method to use in orthofinder (mandatory)"
    echo "  --no_eff_len    -> Turns off calculation of LCR based effective length, use at your own risk"
    echo "  --no_masking    -> Turns off masking of LCRs, use at your own risk"
    echo "  --short_frac    -> Maximum size difference fraction to be accepted (using smallest protein in pair as reference)"
    echo "  --long_frac     -> Maximum size difference fraction to be accepted (using longest protein in pair as reference)"
    echo "  --threads       -> Number of CPU threads to use"
    echo "  --control_file  -> Tab separated file containing short_frac and long_frac values for specific comparisons"
    echo "                     The following format should be used"
    echo "                     Species1.fasta Species2.fasta 0.30 0.25"
    echo "                     Species1.fasta Species3.fasta 0.25 0.25"
    echo "                     Species2.fasta Species3.fasta 0.15 0.25"
    echo
    echo "Notes:"
    echo "  --control_file and --long_frac/--short_frac options are mutually exclusive"
    echo "  --no_eff_len implies --no_masking"
    echo
    echo "Examples:"
    echo "  OrthoPrep.sh --fasta_dir /path/to/my/sequences/folder --short_frac 0.25 --long_frac 0.2 --threads 16"
    echo "  OrthoPrep.sh --fasta_dir /path/to/my/sequences/folder --control_file control.tsv --threads 16 "
    echo "  OrthoPrep.sh --fasta_dir /path/to/my/sequences/folder --control_file control.tsv --no_masking --threads 16"
    echo "  OrthoPrep.sh --fasta_dir /path/to/my/sequences/folder --control_file control.tsv --no_eff_len --threads 16"
	}
function check_files(){
    folder_name=$1
    file_name=$2
    file_path="${folder_name}/${file_name}"
    if [ -f "${file_path}" ]
    then
        echo "pass"
    else
        echo "fail"
    fi
    }
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
		* )
            usage
            exit 0
            ;;
	esac
	shift
done
########################################################
if [ ! -d "${fasta_dir}" ] || [ -z "${fasta_dir}" ]    #
then                                                   #
    echo "Missing fasta directory. Exiting"            #
    exit 0                                             #
fi                                                     #
########################################################

#####################################################################
if [ -z "${control_file}" ]                                         #
then                                                                #
    echo "No control file specified, proceeding in normal mode"     #
elif [ -f "${control_file}" ]                                       #
then                                                                #
    if [ ! -z "${short_frac}" ] || [ ! -z "${long_frac}" ]          #
    then                                                            #
        echo "Incompatible options. Exiting"                        #
        exit 0                                                      #
    else                                                            #
        echo "Using ${control_file} for specific comparisons"       #
        run_mode="custom"                                           #
    fi                                                              #
elif [ ! -f "${control_file}" ]                                     #
then                                                                #
    echo "Control file ${control_file} missing. Exiting"            #
    exit 0                                                          #
fi                                                                  #
#####################################################################


#######################################################################################
short_frac_test=$(echo "${short_frac}" | grep [0-9] | awk '{if($1>=0){print $1}}')    #
if [ -z "${short_frac_test}" ] && [ "${run_mode}" == "normal" ]                       #
then                                                                                  #
    echo "Invalid short fraction value. Exiting"                                      #
    exit 0                                                                            #
fi                                                                                    #
#######################################################################################

##############################################################################################
long_frac_test=$(echo "${long_frac}" | grep [0-9] | awk '{if($1>=0 && $1<=1){print $1}}')    #
if [ -z "${short_frac_test}" ] && [ "${run_mode}" == "normal" ]                              #
then                                                                                         #
    echo "Invalid long fraction value. Exiting"                                              #
    exit 0                                                                                   #
fi                                                                                           #
##############################################################################################

#############################################################################################################
if [ -z "${threads}" ]                                                                                      #
then                                                                                                        #
    threads=${num_proc}                                                                                     #
else                                                                                                        #
    num_proc=$(nproc)                                                                                       #
    threads_test=$(echo -e "${threads}\t${num_proc}" | awk '{if((int($1)==$1) && ($1<=$2)){print $1}}')     #
    if [ -z "${threads_test}" ]                                                                             #
    then                                                                                                    #
        echo "Invalid threads value. Exiting"                                                               #
        exit                                                                                                #
    fi                                                                                                      #
fi                                                                                                          #
#############################################################################################################

#####################################################
of_dep_test=$(which orthofinder.py 2> /dev/null)    #
if [ $? -gt 0 ] || [ -z "${of_dep_test}" ]          #
then                                                #
    echo "OrthoFinder not installed. Exiting"       #
    exit 0                                          #
fi                                                  #
#####################################################

#################################################
dmn_dep_test=$(which diamond 2> /dev/null)      #
if [ $? -gt 0 ] || [ -z "${dmn_dep_test}" ]     #
then                                            #
    echo "Diamond not installed. Exiting"       #
    exit 0                                      #
fi                                              #
#################################################

####################################################
par_dep_test=$(which parallel 2> /dev/null)        #
if [ $? -gt 0 ] || [ -z "${par_dep_test}" ]        #
then                                               #
    echo "GNU Parallel not installed. Exiting"     #
    exit 0                                         #
fi                                                 #
####################################################

#################################################################
cur_date=$(date +%h%d)                                          #
cur_dir=$(pwd)                                                  #
prep_date=$(date +%y-%m-%d)                                     #
log_file="OrthoPrep-${prep_date}.log"                           #
res_dir="Results_${cur_date}"                                   #
work_dir="${fasta_dir}/OrthoFinder/${res_dir}/WorkingDirectory" #
prep_dir="${cur_dir}/OrthoPrep-${prep_date}"                    #
tmp_dir="${prep_dir}/tmp"                                       #
bckp_dir="${prep_dir}/bckp"                                     #
#################################################################

############################################################################################################
if [ "${run_mode}" == "custom" ]                                                                           #
then                                                                                                       #
    fasta_file_list=$(cut -f1,2 "${control_file}" | perl -pe 's/\t/\n/' | sort -V | uniq | grep -v ^$)     #
    for fasta_file in ${fasta_file_list}                                                                   #
    do                                                                                                     #
        file_present=$(check_files ${fasta_dir} ${fasta_file})                                            #
        if [ "${file_present}" == "fail" ]                                                                 #
        then                                                                                               #
            echo "Missing file ${fasta_file}. Exiting"                                                     #
            exit 0                                                                                         #
        fi                                                                                                 #
    done                                                                                                   #
fi                                                                                                         #
############################################################################################################

##############################################################################################
if [ -d "${fasta_dir}/OrthoFinder" ]                                                         #
then                                                                                         #
    echo "We're about to remove ${fasta_dir}/OrthoFinder that may contain previous runs"     #
    echo "Proceed? [N]/Y"                                                                    #
    read rm_check                                                                            #
    if [ "${rm_check}" == "Y" ] || [ "${rm_check}" == "y" ]                                  #
    then                                                                                     #
        rm -rf ${fasta_dir}/OrthoFinder                                                      #
    elif [ "${rm_check}" == "N" ] || [ "${rm_check}" == "n" ] || [ -z "${rm_check}" ]        #
    then                                                                                     #
        echo "Directory not removed, aborting"                                               #
        exit 0                                                                               #
    else                                                                                     #
        echo "Invalid option, aborting"                                                      #
        exit 0                                                                               #
    fi                                                                                       #
fi                                                                                           #
##############################################################################################


#########################################################################################################
echo "Step 1. Preparing fasta files, and diamond commands" | tee -a ${log_file}                         #
mkdir -p "${prep_dir}"                                                                                  #
command_list=$(orthofinder.py -S ${search_method} -op -f ${fasta_dir} | grep -w ^diamond | grep blastp) #
command_list=$(echo "${command_list}" | sed -e "s|${work_dir}|${prep_dir}|g")                           #
cp  "${work_dir}"/Species*.fa \
    "${work_dir}"/SequenceIDs.txt \
    "${work_dir}"/SpeciesIDs.txt \
    "${work_dir}"/${search_method}DBSpecies*.dmnd \
    "${prep_dir}"
#######################################################################################################


#############################################################################
species_count=$(cat "${prep_dir}/SpeciesIDs.txt" | wc -l)                   #
fasta_count=$(ls "${prep_dir}" | grep -c ".fa"$ )                           #
if [ "${species_count}" -eq "${fasta_count}" ]                              #
then                                                                        #
    echo "Step 2. Preprocessing sequence files" | tee -a ${log_file}        #
    if [ "${no_masking}" == "TRUE" ]                                        #
    then                                                                    #
        mkdir ${bckp_dir}                                                   #
        cp ${prep_dir}/Species*.fa ${bckp_dir}                              #
        cp ${prep_dir}/${search_method}DBSpecies*.dmnd ${bckp_dir}          #
    fi                                                                      #
    op_get_eff_len.sh ${prep_dir} ${threads} ${no_masking} ${no_eff_len}    #
    if [ ! $? -eq 0 ]                                                       #
    then                                                                    #
        echo "Error produced extracting LCRs" | tee -a ${log_file}          #
        echo "Aborting" | tee -a ${log_file}                                #
        exit 0                                                              #
    fi                                                                      #
    if [ "${no_masking}" == "FALSE" ]                                       #
    then                                                                    #
        rm ${prep_dir}/${search_method}DBSpecies*.dmnd                      #
        for base_name in $(ls ${prep_dir} | grep fa$ | cut -d\. -f1)        #
        do
            diamond \
                makedb \
                --threads ${threads} \
                --db ${prep_dir}/${search_method}DB${base_name}.dmnd \
                --in ${prep_dir}/${search_method}DB${base_name}.fa          #
        done                                                                #
    fi                                                                      #
fi                                                                          #
#############################################################################

##########################################################################################
num_commands=$(echo "${command_list}" | wc -l)                                           #
echo "Step 3. Running ${num_commands} diamond commands in parallel" | tee -a ${log_file} #
echo "${command_list}" > ${prep_dir}/diamond_command_list.txt                            #
echo "${command_list}" | parallel -j ${threads}                                          #
##########################################################################################

##########################################################################################
echo "Step 4. Filtering BLAST results based on size differences" | tee -a ${log_file}
species_num=$(cat ${prep_dir}/SpeciesIDs.txt | dos2unix | cut -d\: -f1)
eff_len_file="Sequence_eff_len.tsv"
mkdir "${tmp_dir}"
for query in ${species_num}
do
    for subject in ${species_num}
    do
        if [ "${run_mode}" == "normal" ]
        then
            op_blast_filter.py ${prep_dir} ${short_frac} ${long_frac}
        elif [ "${run_mode}" == "custom" ]
        then
            q_fasta=$(awk -v query="${query}"     'BEGIN{FS=": "}{if($1==query)  {print $2}}' ${prep_dir}/SpeciesIDs.txt)
            s_fasta=$(awk -v subject="${subject}" 'BEGIN{FS=": "}{if($1==subject){print $2}}' ${prep_dir}/SpeciesIDs.txt)
            short_frac=$(awk -v query="${q_fasta}" -v subject="${s_fasta}" '{if(($1==query) && ($2==subject)){print $3}}' ${control_file})
            long_frac=$(awk  -v query="${q_fasta}" -v subject="${s_fasta}" '{if(($1==query) && ($2==subject)){print $4}}' ${control_file})
        fi
        op_blast_filter.py ${prep_dir} ${query} ${subject} ${eff_len_file} ${short_frac} ${long_frac}
        mv ${prep_dir}/Blast${query}_${subject}.txt.gz ${bckp_dir}
        mv ${tmp_dir}/Blast${query}_${subject}.txt.gz ${prep_dir}
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