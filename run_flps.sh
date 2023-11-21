#!/bin/bash
tsv_file=$1
seq_id=$2
comp_file=$3
fLPS2 -t 1e-6 -m 5 -M 25 <(grep -w "${seq_id}" "${tsv_file}" | perl -pe 's/\t/\n/') -z thorough -c ${comp_file} | grep -wv WHOLE
