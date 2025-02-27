#!/bin/bash
orthofinder_test=$(which orthofinder.py 2> /dev/null)
if [ $? -gt 0 ] || [ -z "${orthofinder_test}" ]
then
    orthofinder_test=$(which orthofinder 2> /dev/null)
    if [ $? -gt 0 ] || [ -z "${orthofinder_test}" ]
    then
        orthofinder_msg="OrthoFinder not installed"
        orthofinder_test="0"
        orthofinder_cmd="orthofinder"
    else
        orthofinder_msg="OrthoFinder OK"
        orthofinder_test="1"
        orthofinder_cmd="orthofinder"
    fi
else
    orthofinder_msg="OrthoFinder OK"
    orthofinder_test="1"
    orthofinder_cmd="orthofinder.py"
fi
diamond_test=$(which diamond 2> /dev/null)
if [ $? -gt 0 ] || [ -z "${diamond_test}" ]
then
    diamond_msg="Diamond not installed"
    diamond_test="0"
else
    diamond_msg="Diamond OK"
    diamond_test="1"
fi
parallel_test=$(which parallel 2> /dev/null)
if [ $? -gt 0 ] || [ -z "${parallel_test}" ]
then
    parallel_msg="GNU Parallel not installed"
    parallel_test="0"
else
    parallel_msg="GNU Parallel OK"
    parallel_test="1"
fi
flps_test=$(which fLPS2 2> /dev/null)
if [ $? -gt 0 ] || [ -z "${flps_test}" ]
then
    flps_msg="fLPS2 not installed"
    flps_test="0"
else
    flps_msg="fLPS2 OK"
    flps_test="1"
fi
seg_test=$(which seg 2> /dev/null)
if [ $? -gt 0 ] || [ -z "${seg_test}" ]
then
    seg_msg="seg not installed"
    seg_test="0"
else
    seg_msg="seg OK"
    seg_test="1"
fi
bed_test=$(which bedtools 2> /dev/null)
if [ $? -gt 0 ] || [ -z "${bed_test}" ]
then
    bed_msg="bedTools not installed"
    bed_test="0"
else
    bed_msg="bedTools OK"
    bed_test="1"
fi
echo -e "${orthofinder_test}\t${orthofinder_msg}\t${orthofinder_cmd}\n${diamond_test}\t${diamond_msg}\n${parallel_test}\t${parallel_msg}\n${flps_test}\t${flps_msg}\n${seg_test}\t${seg_msg}\n${bed_test}\t${bed_msg}"