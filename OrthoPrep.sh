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
    echo "  --fasta_dir    -> Directory containing protein sequences in fasta format (mandatory)"
    echo "  --no_eff_len   -> Turns off calculation of LCR based effective length, use at your own risk"
    echo "  --no_masking   -> Turns off masking of LCRs, use at your own risk"
    echo "  --short_frac   -> Maximum size difference fraction to be accepted (using smallest protein in pair as reference)"
    echo "  --long_frac    -> Maximum size difference fraction to be accepted (using longest protein in pair as reference)"
    echo "  --threads      -> Number of CPU threads to use"
    echo "  --control_file -> Tab separated file containing short_frac and long_frac values for specific comparisons"
    echo "                    The following format should be used"
    echo "                    Species1.fasta Species2.fasta 0.30 0.25"
    echo "                    Species1.fasta Species3.fasta 0.25 0.25"
    echo "                    Species2.fasta Species3.fasta 0.15 0.25"
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
        --fasta_dir    )
            shift
            fasta_dir=$(realpath $1)
            ;;
        --no_masking   )
            no_masking="TRUE"
            ;;
        --no_eff_len   )
            no_eff_len="TRUE"
            ;;
		--short_frac   )
            shift
            short_frac=$1
            ;;
		--long_frac    )
            shift
            long_frac=$1
            ;;
		--threads      )
            shift
            threads=$1
            ;;
        --control_file )
            shift
            control_file=$1
            ;;
		--help )
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

#####################################################################
cur_date=$(date +%h%d)                                              #
cur_dir=$(pwd)                                                      #
prep_date=$(date +%y-%m-%d)                                         #
log_file="OrthoPrep-${prep_date}.log"                               #
res_dir="Results_${cur_date}"                                       #
work_dir="${fasta_dir}/OrthoFinder/${res_dir}/WorkingDirectory"     #
prep_dir="${cur_dir}/OrthoPrep-${prep_date}"                        #
#####################################################################

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


#######################################################################################################
echo "Step 1. Preparing fasta files, and diamond commands" | tee -a ${log_file}                       #
mkdir -p "${prep_dir}"                                                                                #
#command_list=$(orthofinder.py -S diamond_ulow -op -f ${fasta_dir} | grep -w ^diamond | grep blastp)  # <- 1e-3
command_list=$(orthofinder.py -S diamond_vlow -op -f ${fasta_dir} | grep -w ^diamond | grep blastp)   # <- 1e-6
#command_list=$(orthofinder.py -S diamond_low  -op -f ${fasta_dir} | grep -w ^diamond | grep blastp)  # <- 1e-12
#command_list=$(orthofinder.py -S diamond_med  -op -f ${fasta_dir} | grep -w ^diamond | grep blastp)  # <- 1e-15
#command_list=$(orthofinder.py -S diamond_hard -op -f ${fasta_dir} | grep -w ^diamond | grep blastp)  # <- 1e-18
#command_list=$(orthofinder.py -S diamond_def  -op -f ${fasta_dir} | grep -w ^diamond | grep blastp)  # <- 1e-3 no masking
command_list=$(echo "${command_list}" | sed -e "s|${work_dir}|${prep_dir}|g")                         #
cp "${work_dir}"/Species*.fa "${work_dir}"/SequenceIDs.txt "${work_dir}"/SpeciesIDs.txt "${prep_dir}" #
#######################################################################################################


########################################################################################
species_count=$(cat "${prep_dir}/SpeciesIDs.txt" | wc -l)                              #
fasta_count=$(ls "${prep_dir}" | grep -c ".fa"$ )                                      #
if [ "${species_count}" -eq "${fasta_count}" ]                                         #
then                                                                                   #
    echo "Step 2. Preprocessing sequence files" | tee -a ${log_file}                   #
    echo "        Calculating sequence lengths" | tee -a ${log_file}                   #
    cat ${prep_dir}/Species*.fa \
        | dos2unix \
        | perl -pe 'if(/\>/){s/$/\t/};s/\n//g;s/\>/\n/g' \
        | tail -n+2 \
        | sort -V \
        | uniq \
        | awk 'BEGIN{FS="\t"}{print $1 "@" $1 FS length($2)}' \
        | perl -pe 's/_.*\@/\t/' | tee -a ${prep_dir}/Sequence_len.tsv                        #
    if [ ! $? -eq 0 ]                                                                  #
    then                                                                               #
        echo "Error produced getting sequence lengths" | tee -a ${log_file}            #
        echo "Aborting" | tee -a ${log_file}                                           #
        exit 0                                                                         #
    fi                                                                                 #
    if [ "${no_masking}" == "FALSE" ]                                                  #
    then                                                                               #
        mkdir ${prep_dir}/bckp                                                         #
        cp ${prep_dir}/Species*.fa ${prep_dir}/bckp                                    #
    fi                                                                                 #
    if [ "${no_eff_len}" == "FALSE" ]                                                  #
    then                                                                               #
        echo "        Calculating effective lengths" | tee -a ${log_file}              #
        op_lcr_preprocess.sh "${prep_dir}" "${threads}" "${no_masking}"                #
        if [ ! $? -eq 0 ]                                                              #
        then                                                                           #
            echo "Error produced extracting LCRs" | tee -a ${log_file}                 #
            echo "Aborting" | tee -a ${log_file}                                       #
            exit 0                                                                     #
        fi                                                                             #
    fi                                                                                 #
fi                                                                                     #
########################################################################################

##########################################################################################
num_commands=$(echo "${command_list}" | wc -l)                                           #
echo "Step 3. Running ${num_commands} diamond commands in parallel" | tee -a ${log_file} #
echo "${command_list}" | parallel -j ${threads}                                          #
##########################################################################################

##########################################################################################
echo "Step 4. Filtering BLAST results based on size differences" | tee -a ${log_file}
echo "        Copying required files to destination directories" | tee -a ${log_file}
if [ "${no_masking}" == "TRUE" ]
then
    seq_origin_dir="${prep_dir}"
elif [ "${no_masking}" == "FALSE" ]
then
    seq_origin_dir="${cur_dir}/bckp"
fi
mkdir ${prep_dir}/WorkingDirectory
cp ${prep_dir}/SequenceIDs.txt ${prep_dir}/SpeciesIDs.txt ${seq_origin_dir}/Species*.fa ${prep_dir}/WorkingDirectory
echo "        Filtering pairs" | tee -a ${log_file}
species_num=$(cat ${prep_dir}/SpeciesIDs.txt | dos2unix | cut -d\: -f1)
sequence_len_file="Sequence_len.tsv"
for query in ${species_num}
do
    for subject in ${species_num}
    do
        if [ "${run_mode}" == "normal" ]
        then
            op_blast_filter.py ${prep_dir} ${short_frac} ${long_frac}
        elif [ "${run_mode}" == "custom" ]
        then
            q_fasta=$(awk -v query="${query}"     'BEGIN{FS=": "}{if($1==query)  {print $2}}' ${prep_dir}/SequenceIDs.txt)
            s_fasta=$(awk -v subject="${subject}" 'BEGIN{FS=": "}{if($1==subject){print $2}}' ${prep_dir}/SequenceIDs.txt)
            short_frac=$(awk -v query="${q_fasta}" -v subject="${s_fasta}" '{if(($1==query) && ($2==subject)){print $3}}' ${control_file})
            long_frac=$(awk  -v query="${q_fasta}" -v subject="${s_fasta}" '{if(($1==query) && ($2==subject)){print $4}}' ${control_file})
        fi
        op_blast_filter.py ${prep_dir} ${query} ${subject} ${sequence_len_file} ${short_frac} ${long_frac}
    done
fi
if [ ! $? -eq 0 ]
then
    echo "Something went wrong while filtering BLAST results" | tee -a ${log_file}
    echo "Aborting" | tee -a ${log_file}
    exit 0
fi
echo "        Removing temporary files" | tee -a ${log_file}
#rm ${work_dir}/Species*.sizes.tsv
#rm -rf ${cur_dir}/bckp
echo ""  | tee -a ${log_file}
echo "Finished filtering BLAST results" | tee -a ${log_file}
echo "Files in the ${prep_dir}/WorkingDirectory directory were obtained by applying strict filters to the BLAST results" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "Files in the ${prep_dir}/WorkingDirectory directory were obtained by applying more relaxed filters to the BLAST results" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "Now you can resume OrthoFinder by running:" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "orthofinder.py -b ${prep_dir}/WorkingDirectory [other OrthoFinder options]" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "OrthoPrep v0.0.1"
