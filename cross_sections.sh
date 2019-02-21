#!/bin/bash

#Set up parameters.
log2f=0.59 

type1="cbp20"
type2="cbp80"
type3="se1"

file1="deseq2_WT_cbp20_005.txt"
file2="deseq2_WT_cbp80_005.txt"
file3="deseq2_WT_se1_005.txt"

# 3 outputs' cross section.
log2_output3="cross_sectionslog2_${log2f}_${type1}_${type2}_${type3}.txt"
cat $file1 $file2 $file3 | sed 1,1d | awk -v var=$log2f '{if (($3>=var)||($3<=-var)) print $1}' | sort | uniq -c | awk '{if ($1>2) print $2}' > $log2_output3
echo $log2f
wc -l $log2_output3

# 2 outputs' cross section.
log2_output12="cross_sectionslog2_${log2f}_${type1}_${type2}.txt"
cat ${file1} ${file2} | sed 1,1d | awk -v var=$log2f '{if (($3>=var)||($3<=-var)) print $1}' | sort | uniq -c | awk '{if ($1>1) print $2}'  > $log2_output12
wc -l $log2_output12 

log2_output13="cross_sectionslog2_${log2f}_${type1}_${type3}.txt"
cat ${file1} ${file3} | sed 1,1d | awk -v var=$log2f '{if (($3>=var)||($3<=-var)) print $1}' | sort | uniq -c | awk '{if ($1>1) print $2}'  > $log2_output13
wc -l $log2_output13

log2_output23="cross_sectionslog2_${log2f}_${type2}_${type3}.txt"
cat ${file2} ${file3} | sed 1,1d | awk -v var=$log2f '{if (($3>=var)||($3<=-var)) print $1}' | sort | uniq -c | awk '{if ($1>1) print $2}'  > $log2_output23
wc -l $log2_output23
