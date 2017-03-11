# Load ravi as first arg, key as second arg, and generate a new ravi file with fragments in first column.

import sys

flip = {}

a = open(sys.argv[1] + '.flip','w')

with open(sys.argv[2]) as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		fragment = values[0]
		num = values[1]
		flip[num] = fragment

with open(sys.argv[1]) as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		num = values[0]
		nca = values[1].split('.')[0]
		score = values[2]
		a.write(flip[num] + '\t' + nca + '\t' + score + '\n')

a.close()
