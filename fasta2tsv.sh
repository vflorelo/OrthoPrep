#!/bin/bash
fasta_file=$1
format_test=$(grep -c \> "${fasta_file}")
if [ -z "${format_test}" ] || [ "${format_test}" -eq 0 ]
then
	echo "Not a fasta file"
	exit 1
fi
perl -pe 'if(/\>/){s/$/\t/};s/\n//g;s/\>/\n\>/g' "${fasta_file}" | tail -n+2
