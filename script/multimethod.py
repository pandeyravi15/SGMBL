import math
import sys

# First arg is key
# Second arg is blast
# Third arg is flip

NC = {}

with open('taxonomy.txt') as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		nca = values[0]
		genus = values[1].split('_')[0]
		NC[nca] = genus

Frag2Genus = {}

with open(sys.argv[1]) as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		Frag2Genus[values[0]] = [values[1]]

blast = {}

with open(sys.argv[2]) as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		fragment = values[0]
		#nca = values[1].split('|')[3].split('.')[0]
		nca = values[1]
		evalue = float(values[10])
		if evalue == 0.0:
			evalue = 2.2e-308
		blast[fragment + ':' + nca] = math.log(evalue,10)
ravi = {}

with open(sys.argv[3]) as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		fragment = values[0]
		nca = values[1]
		score = float(values[2])
		ravi[fragment + ':' + nca] = score
		
results = {}

for entry in ravi:
	fragment = entry.split(':')[0]
	nca = entry.split(':')[1]
	if entry in blast:
		finalscore = ((ravi[entry]) + (1.2*(4-blast[entry])))
		#print(finalscore)
	else:
		finalscore = ravi[entry]
	if fragment not in results:
		results[fragment] = [finalscore, nca]
	else:
		if finalscore > results[fragment][0]:
			results[fragment] = [finalscore, nca]
		else:
			pass

a = open(sys.argv[3] + '.results','w')

for rez in results:
	a.write(rez + '\t' + results[rez][1] + '\n')

a.close()

blast = {}

with open(sys.argv[2]) as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		fragment = values[0]
		#nca = values[1].split('|')[3].split('.')[0]
		nca = values[1]
		if fragment not in blast:
			blast[fragment] = nca
		else:
			pass

results = {}

for entry in ravi:
	fragment = entry.split(':')[0]
	nca = entry.split(':')[1]
	if fragment in blast:
		results[fragment] = [10000, blast[fragment]]
	if fragment not in results:
		results[fragment] = [finalscore, nca]
	else:
		if finalscore > results[fragment][0]:
			results[fragment] = [finalscore, nca]
		else:
			pass
b = open(sys.argv[3] + '.tophit_results','w')

for rez in results:
	b.write(rez + '\t' + results[rez][1] + '\n')

b.close()


