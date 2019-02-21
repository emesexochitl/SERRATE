#!/bin/bash

for stress in dehydt salt
do

for num in 1 6 12
do

echo $stress $num

cat deseq1_soy_${stress}_${num}h.csv | sed -e 2d | awk -F"," '{if (($6>=2)&&($8<0.05)) print $1}' | tr '[:lower:]' '[:upper:]' > deseq1_soy_${stress}_${num}h_clean.txt

done

done
