#!/usr/bin/python
import sys
import os
import commands
import time
import re
import csv
import string


input_path = sys.argv[1]   # /hi-c-human/10_histone_modification
output_path = sys.argv[2]   # ./
bam_file = sys.argv[3]   #the .bam file
bed_file = sys.argv[4]  #the BS bed file
readcount_file = sys.argv[5]  # ../10_histone_modification/histone_rep_count.txt
output_file = sys.argv[6]  # ./histone_rep_BS_RPM.txt
pwm_length = sys.argv[7]
flank_length = sys.argv[8]
IsBpWise= sys.argv[9]

input_path = os.path.normpath(input_path)
input_path = input_path + '/'
output_path = os.path.normpath(output_path)
output_path = output_path + '/'


f = open(readcount_file,"r+")
line = f.readline()
f.close()

mappedReads = line.rstrip()

if int(IsBpWise) == 0:   
    scalingFactor = float(1000000)/(int(mappedReads)*(int(pwm_length)+2*int(flank_length)))
    cmd = 'coverageBed -abam '+ bam_file + ' -b ' + bed_file + ' -d | ' + 'awk \'{a[$1"\t"$2"\t"$3"\t"]+=$5}END{for(i in a){print i,a[i]*' + str(scalingFactor) + '}}\' - >> ' + output_file
    os.system(cmd)
elif int(IsBpWise) == 1:
    scalingFactor = float(1000000)/(int(mappedReads))
    cmd = 'coverageBed -abam '+ bam_file + ' -b ' + bed_file + ' -d | ' + 'awk \'{print $1"\t"$2"\t"$3"\t"$4"\t"$5*' + str(scalingFactor) + '}\' - ' + '>> ' + output_file
    os.system(cmd)
else:
    sys.exit("Wrong IsBpWise option!")
