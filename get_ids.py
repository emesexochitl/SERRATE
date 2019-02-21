#!/usr/bin/python

import datetime

ori="sra_rice.txt"
ori_list=[]
with open(ori) as inputfile:
    for line in inputfile:
       ori_list.append(line.strip().split('\t'))

ori1=[item[0] for item in ori_list]

download="downloaded.txt"
download_list=[]
with open(download) as inputfile:
    for line in inputfile:
       download_list.append(line.strip().split('\t'))

download1=[item[0] for item in download_list]

output="sra_rice_new3.txt"
print ori1[0], download1[0]

myfile= open(output, "w")

print len(list(set(ori1) - set(download1)))

for i in list(set(ori1) - set(download1)):
   myfile.write("%s\n" % (i))

myfile.close()


