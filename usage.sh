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
    echo "  --blast_filter  -> Filters BLAST results based on total- or effective-length [TRUE]/FALSE"
    echo "  --eff_len       -> Calculates LCR based effective length [TRUE]/FALSE"
    echo "  --masking       -> Hardmasking of LCRs for BLAST searches [TRUE]/FALSE"
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
    echo
    echo "Examples:"
    echo "  OrthoPrep.sh --fasta_dir /path/to/my/sequences/folder --short_frac 0.25 --long_frac 0.2 --threads 16"
    echo "  OrthoPrep.sh --fasta_dir /path/to/my/sequences/folder --control_file control.tsv --threads 16 "
    echo "  OrthoPrep.sh --fasta_dir /path/to/my/sequences/folder --control_file control.tsv --masking FALSE --threads 16"
    echo "  OrthoPrep.sh --fasta_dir /path/to/my/sequences/folder --control_file control.tsv --eff_len FALSE --threads 16"
	}
export -f usage