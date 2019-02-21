#!/bin/bash

while read l
do

#hit=`grep "Mapped" tophat_out_${l}/align_summary.txt`

#echo $l $hit

#echo $l

#echo tophat_out_${l}/${l}_accepted_hits.bam | transpose -t

mv tophat_out_${l}/accepted_hits.bam tophat_out_${l}/${l}_accepted_hits.bam 

done < SRR_Acc_List.txt
