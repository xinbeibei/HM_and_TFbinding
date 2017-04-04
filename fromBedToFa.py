#!/usr/bin/python
import sys
import os
import commands
import time
import re
import csv
import string


input_file = sys.argv[1] # ./hi-c-human/$tf/.bed
output_file = sys.argv[2] # ./$tf/.fa

f = open(input_file, "r+")
lines = f.readlines()
f.close()

out = open(output_file, "w")
for line in lines:
    items = re.split('\t', line)
    out.write('>' + items[0] + ':' + str(items[1]) + '-' + str(items[2]) + '\n')
    out.write(str(items[3]))
out.close()


