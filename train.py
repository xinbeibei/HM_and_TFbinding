#!/usr/bin/python
import sys
import os
import commands
import time
import csv

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

if len(sys.argv) != 8:
    print 'Incorrect number of arguments:', len(sys.argv)
    print 'Usage:', sys.argv[0], '<tf_name> <encoded_path> <index_path> <output_path> <k> <Rscript_location> <feature_list_location>'
    sys.exit(1)

tf = sys.argv[1]
encoded_path = sys.argv[2]
index_path = sys.argv[3]
output_path = sys.argv[4]
k = int(sys.argv[5])
rscript = sys.argv[6]
feature_list_location = sys.argv[7]

# create needed folders if not existed |||
if not os.path.exists(output_path):
    os.makedirs(output_path)

output_path = os.path.normpath(output_path)
output_path = output_path + '/'

index_path = os.path.normpath(index_path)
index_path = index_path + '/'

encoded_path = os.path.normpath(encoded_path)
encoded_path = encoded_path + '/'

## get feature combinations from file ++
with open(feature_list_location) as feature_list:
    cmb = list(csv.reader(feature_list, delimiter='\t', quoting=csv.QUOTE_NONE))
combinations = [row[1] for row in cmb]

for cb in combinations:
    for subi in range(1, k+1):
        # generate .qsub file ...
        encoded_filename = encoded_path + tf + '.' + cb
        train_index_filename = index_path + tf + '.fold' + str('{0:02}'.format(subi)) + '.train'
        log_filename = output_path + tf + '.' + cb + '.fold' + str('{0:02}'.format(subi)) + '.train.log'
        out_direction = output_path + tf + '.' + cb + '.fold' + str('{0:02}'.format(subi)) + '.train'

        cmd = 'R --no-restore --no-save --args ' + encoded_filename + ' ' + train_index_filename + ' ' + out_direction  + ' ' + str(k) + ' < ' + rscript + '\n'
        os.system(cmd)