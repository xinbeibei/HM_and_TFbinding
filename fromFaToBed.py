#!/usr/bin/python
import sys
import os
import commands
import time
import re
import csv
import string


input_file = sys.argv[1] # ./hi-c-human/$tf/.fa
output_file = sys.argv[2] # ./$tf/.bed
#flankinglen = sys.argv[3]  # can be 0 if no flanking region is added

out = open(output_file, "w")
with open(input_file, "r+") as f:
	for line in f:
		items = line.rstrip()
		sub = items[1:]
		subs = re.split(':', sub)
		subss = re.split('-', subs[1])
		out.write(subs[0]+'\t'+str(int(subss[0]))+'\t'+str(int(subss[1]))+'\n')
		f.next()
out.close()
f.close()

