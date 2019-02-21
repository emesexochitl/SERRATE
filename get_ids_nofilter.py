#!/usr/bin/python

import datetime

ori="arabidopsis_thaliana_iso1_intronnum.txt"
ori_list=[]
with open(ori) as inputfile:
    for line in inputfile:
       ori_list.append(line.strip().split('\t'))

ori1=[item[0] for item in ori_list] # gene ID.
ori2=[item[1] for item in ori_list] # Intron number.

deseq= "../diff_expression/deseq2_WT_se1_005.txt"
deseq_list=[]
with open(deseq) as inputfile:
    for line in inputfile:
        deseq_list.append(line.strip().split(' '))

del deseq_list[:1]
print deseq_list[:5]

deseq1=[item[0] for item in deseq_list] #gene_ID
deseq2=[float(item[2]) for item in deseq_list] #log2 FC
deseq3=[item[6] for item in deseq_list] # padj

output="deseq2_WT_se1_intronnum_005_log2_all_nofiler.txt"
output1="deseq2_WT_se1_intronnum_005_log2_1up_nofiler.txt"
output2="deseq2_WT_se1_intronnum_005_log2_1down_nofiler.txt"

#fcnum = float(1)
fcnum= 0
#print ori1[0], download1[0]
myfile= open(output, "w")
myfile1= open(output1, "w")
myfile2= open(output2, "w")

for i in xrange(0, len(ori_list)):
    for j in xrange(0, len(deseq_list)):

        if ori1[i] == deseq1[j] and deseq2[j] >= fcnum:
            myfile1.write("%s_%s\t%s\t%s\t%s\n" % ("Arabidopsis", "se1_total_up", deseq1[j], ori2[i], deseq2[j]))
        if ori1[i] == deseq1[j] and deseq2[j] <= (fcnum)*-1:
            myfile2.write("%s_%s\t%s\t%s\t%s\n" % ("Arabidopsis", "se1_total_down", deseq1[j], ori2[i], deseq2[j]))

        if ori1[i] == deseq1[j]:
            myfile.write("%s_%s\t%s\t%s\t%s\n" % ("Arabidopsis", "se1_total_all", deseq1[j], ori2[i], deseq2[j]))


        else:
            continue

myfile.close()
myfile1.close()
myfile2.close()

