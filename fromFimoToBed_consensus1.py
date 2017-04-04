import sys
import os
import commands
import time
import re
import csv
import string


input_path = sys.argv[1] # ./fimo_out
output_path = sys.argv[2] # ./$tf

input_path = os.path.normpath(input_path)
input_path = input_path + '/'
output_path = os.path.normpath(output_path)
output_path = output_path + '/'

input_file = input_path + 'unique_fimo_out.txt'
output_file1 = output_path + 'BSs_regions_for_getfasta.bed'
f = open(input_file, "r+")
lines = f.readlines()
f.close()


out = open(output_file1, "w")
out2 = open(output_path + 'BSs_motif_score.bed',"w")

for line in lines:
    items = re.split('\t', line)
    items2 = re.split(':',items[1])
    out.write(items2[0])
    out.write('\t')
    items3 = re.split('-', items2[1])
    out.write(str(int(items3[0])+int(items[2])-3))
    out.write('\t')
    out.write(str(int(items3[0])+int(items[3])+2))
    out2.write(items2[0])
    out2.write('\t')
    out2.write(str(int(items3[0])+int(items[2])))
    out2.write('\t')
    out2.write(str(int(items3[0])+int(items[3])))
    out2.write('\t')
    out2.write(items[5]+'\n')
    out.write('\t')
    out.write("forward" + '\t' + str(1)+ '\t')
    out.write(items[4]+'\n')
out.close()
out2.close()

output_file2 = output_path + 'BSs_seq_with_2_flanks.fa'
cmd = 'bedtools getfasta -fi /home/cmb-04/rr/bxin/hi-c-human/hg19.fa -bed ' + output_file1 + ' -s -fo ' + output_file2
os.system(cmd)

out = open(output_path + 'BSs_regions_from_fasta.bed',"w")

with open(output_file2, "r+") as f:
	for line in f:
		items = line.rstrip()
		sub = items[1:-3]
		subs = re.split(':', sub)
		subss = re.split('-', subs[1])
		out.write(subs[0]+'\t'+str(int(subss[0])+3)+'\t'+str(int(subss[1])-2)+'\t')
		out.write(f.next().upper())
out.close()
f.close()

