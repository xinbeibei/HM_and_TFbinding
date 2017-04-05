#!/usr/bin/python
import sys
import os
import commands
import time
import csv

if len(sys.argv) != 6:
    print 'Incorrect number of arguments:', len(sys.argv)
    print 'Usage:', sys.argv[0], '<tf_name> <input_path> <path_to_r_script> <output_path> <feature_list_location>'
    sys.exit(1)

tf = sys.argv[1]
input_path = sys.argv[2]
r_script = sys.argv[3]
output_path = sys.argv[4]
feature_list_location = sys.argv[5]

# create needed folders if not existed |||
if not os.path.exists(output_path):
    os.makedirs(output_path)
# |||

input_path = os.path.normpath(input_path)
input_path = input_path + '/'

output_path = os.path.normpath(output_path)
output_path = output_path + '/'

# get feature naming from file ++
with open(feature_list_location) as feature_list:
    cmb = list(csv.reader(feature_list, delimiter='\t', quoting=csv.QUOTE_NONE))
header = [row[0] for row in cmb]
header.insert(0, '#seqs')
header.insert(0, 'Dataset')
#header = ['Dataset', '#sequence', '1mer', '1mer+shape', '1/2mer', '1/2mer+shape', '1/3mer', '1/3mer+shape', '1/2/3mer'] # this header is for the matrix r2 in summarize.R
# ++
f = open(output_path+'summary_auroc.txt', 'w') # fixed filename "summary.txt"
f.write(' '.join(header))
f.write('\n')
f.close()

# construct command
#cmd = 'R --no-restore --no-save --args ' + input_path + identifier + ' ' + output_path + ' ' + feature_list_location + ' < ' + r_script + ' 2>&1 >/dev/null'
#os.system(cmd)
cmd = 'Rscript ' + r_script + ' ' + input_path + tf + ' ' + output_path + ' ' + feature_list_location
os.system(cmd)