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
            fasta_dir=$1
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
log_date=$(date +%y-%m-%d)                                          #
log_file="OrthoPrep-${log_date}.log"                                #
res_dir="Results_${cur_date}"                                       #
work_dir="${fasta_dir}/OrthoFinder/${res_dir}/WorkingDirectory"     #
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


######################################################################################################
echo "Step 1. Preparing fasta files, and diamond commands" | tee -a ${log_file}                      #
#command_list=$(orthofinder.py -S diamond_ulow -op -f ${fasta_dir} | grep -w ^diamond | grep blastp) # <- 1e-3
command_list=$(orthofinder.py -S diamond_vlow -op -f ${fasta_dir} | grep -w ^diamond | grep blastp)  # <- 1e-6
#command_list=$(orthofinder.py -S diamond_low  -op -f ${fasta_dir} | grep -w ^diamond | grep blastp) # <- 1e-12
#command_list=$(orthofinder.py -S diamond_med  -op -f ${fasta_dir} | grep -w ^diamond | grep blastp) # <- 1e-15
#command_list=$(orthofinder.py -S diamond_hard -op -f ${fasta_dir} | grep -w ^diamond | grep blastp) # <- 1e-18
#command_list=$(orthofinder.py -S diamond_def  -op -f ${fasta_dir} | grep -w ^diamond | grep blastp) # <- 1e-3 no masking
cd "${work_dir}"                                                                                     #
work_dir=$(pwd)                                                                                      #
cd "${cur_dir}"                                                                                      #
######################################################################################################

########################################################################################
echo "Step 2. Getting sequence lengths and effective lengths" | tee -a ${log_file}     #
op_get_lengths.sh "${work_dir}"                                                        #
if [ ! $? -eq 0 ]                                                                      #
then                                                                                   #
    echo "Something went wrong while extracting sequence lengths" | tee -a ${log_file} #
    echo "Aborting" | tee -a ${log_file}                                               #
    exit 0                                                                             #
fi                                                                                     #
if [ "${no_masking}" == "FALSE" ]                                                      #
then                                                                                   #
    mkdir ${cur_dir}/bckp                                                              #
    op_bckp_seqs.sh "${work_dir}"                                                      #
    if [ ! $? -eq 0 ]                                                                  #
    then                                                                               #
        echo "Something went wrong backing up sequences" | tee -a ${log_file}          #
        echo "Aborting" | tee -a ${log_file}                                           #
        exit 0                                                                         #
    fi                                                                                 #
fi                                                                                     #
########################################################################################
if [ "${no_eff_len}" == "FALSE" ]                                                      #
then                                                                                   #
    op_lcr_preprocess.sh "${work_dir}" "${threads}" "${no_masking}"                    #
    if [ ! $? -eq 0 ]                                                                  #
    then                                                                               #
        echo "Something went wrong extracting LCRs from proteins" | tee -a ${log_file} #
        echo "Aborting" | tee -a ${log_file}                                           #
        exit 0                                                                         #
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
    seq_origin_dir="${work_dir}"
elif [ "${no_masking}" == "FALSE" ]
then
    seq_origin_dir="${cur_dir}/bckp"
fi
mkdir ${work_dir}/hq ${work_dir}/mq
cp ${work_dir}/SequenceIDs.txt ${work_dir}/SpeciesIDs.txt ${seq_origin_dir}/Species*.fa ${work_dir}/hq
cp ${work_dir}/SequenceIDs.txt ${work_dir}/SpeciesIDs.txt ${seq_origin_dir}/Species*.fa ${work_dir}/mq
echo "        Filtering pairs" | tee -a ${log_file}
if [ "${run_mode}" == "normal" ]
then
    op_blast_filter.py ${work_dir} ${short_frac} ${long_frac}
elif [ "${run_mode}" == "custom" ]
then
    num_comparisons=$(cat ${control_file} | wc -l)
    for i in $(seq 1 ${num_comparisons})
    do
        query=$(tail   -n+${i} ${control_file} | head -n1 | cut -f1)
        subject=$(tail -n+${i} ${control_file} | head -n1 | cut -f2)
        cur_sf=$(tail  -n+${i} ${control_file} | head -n1 | cut -f3)
        cur_lf=$(tail  -n+${i} ${control_file} | head -n1 | cut -f4)
        query_id=$(grep   -w "${query}"   ${work_dir}/SpeciesIDs.txt | cut -d\: -f1)
        subject_id=$(grep -w "${subject}" ${work_dir}/SpeciesIDs.txt | cut -d\: -f1)
        blast_file="Blast${query_id}_${subject_id}.txt.gz"
        query_size_file="Species${query_id}.sizes.tsv"
        subject_size_file="Species${subject_id}.sizes.tsv"
        op_c_blast_filter.py ${work_dir} ${blast_file} ${query_size_file} ${subject_size_file} ${cur_sf} ${cur_lf}
    done
fi
if [ ! $? -eq 0 ]
then
    echo "Something went wrong while filtering BLAST results" | tee -a ${log_file}
    echo "Aborting" | tee -a ${log_file}
    exit 0
fi
echo "        Removing temporary files" | tee -a ${log_file}
rm ${work_dir}/Species*.sizes.tsv
rm -rf ${cur_dir}/bckp
rm -rf ${cur_dir}/
echo ""  | tee -a ${log_file}
echo "Finished filtering BLAST results" | tee -a ${log_file}
echo "Files in the ${work_dir}/hq directory were obtained by applying strict filters to the BLAST results" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "Files in the ${work_dir}/hq directory were obtained by applying more relaxed filters to the BLAST results" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "Now you can resume OrthoFinder by running:" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "orthofinder.py -b ${work_dir}/hq [other OrthoFinder options]" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "or" | tee -a ${log_file}
echo | tee -a ${log_file}
echo "orthofinder.py -b ${work_dir}/mq [other OrthoFinder options]" | tee -a ${log_file}
echo  | tee -a ${log_file}
echo "OrthoPrep.sh v0.0.1"
