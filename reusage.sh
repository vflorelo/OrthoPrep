function usage(){
    echo "resumeOrthoPrep takes a folder containing sequences previously processed by OrthoPrep"
    echo "to be re-run with OrthoFinder (or similar) with different options for BLAST comparisons"
    echo
    echo "Diamond databases are re-build based on masked sequences"
    echo "BLAST results are then updated and filtered with the pre-calculated effective lengths"
    echo
    echo "The size difference is divided by a fraction of the length of each sequence and then compared against a manually set threshold"
    echo " - Pairs of sequences with a size difference lower than the set threshold are kept as comparable sequences"
    echo " - Pairs of sequences with a size difference greater than the set threshold are excluded"
    echo
    echo "Options:"
    echo "  --prep_dir      -> Directory containing OrthoPrep results (mandatory)"
    echo "  --search_method -> Sequence similarity method to use in orthofinder (mandatory)"
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
    echo "  resume_OrthoPrep.sh --prep_dir /path/to/previous/orthoprep/run --short_frac 0.25 --long_frac 0.2 --threads 16 --search_method diamond"
    echo "  resume_OrthoPrep.sh --prep_dir /path/to/previous/orthoprep/run --control_file control.tsv --threads 16 --search_method diamond"
	}
export -f usage