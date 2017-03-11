import sys

# Load in the files from arguments.

blastfile = sys.argv[1]
ravifile = sys.argv[2]
keyfile = sys.argv[3]

# Generate the result file.

resultfile = sys.argv[2] + '.results'

# Load in BLAST results, and choose the lowest evalue result(s)
# Store the evalue into the zeroth position, and then all matching accessions in the following.
# Pop off the evalue at the end.

blast = {}

with open(blastfile) as infile:
	for lines in infile:
		lines = lines.rstrip()
		fragment = lines.split('\t')[0]
		evalue = float(lines.split('\t')[10])
		dirtync = lines.split('\t')[1]
		#cleanernc = dirtync.split('|')[3]
		cleanernc = dirtync
		#nc = cleanernc.split('.')[0]
		nc = cleanernc
		if fragment in blast:
			if blast[fragment][0] > evalue:
				blast[fragment] = []
				blast[fragment].append(evalue)
				blast[fragment].append(nc)
			elif blast[fragment][0] == evalue:
				blast[fragment].append(nc)
			elif blast[fragment][0] < evalue:
					pass
		else:
			blast[fragment] = []
			blast[fragment].append(evalue)
			blast[fragment].append(nc)

for bs in blast:
	blast[bs].pop(0)

# Load in the key file, which allows the reading of Ravi's output.
# This will convert Ravi's sequential numbering with the actual fragment ID.
# Key reads number to fragment ID. (Key[1] = r2937)

key = {}

with open(keyfile) as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		key[values[1]] = values[0]
	
# Using the key file, find all fragments that are not included in the BLAST results.  These will be handled by Ravi alone.
# For all others, we will need to compare their associated scores from Ravi output to find the best match.

noblast = []
yesblast = []

for keys in key:
	if key[keys] not in blast:
		noblast.append(key[keys])
	else:
		yesblast.append(key[keys])


# Open up Ravi's output, and parse.
# If the fragment had blast matches, keep a dictionary entry for that fragment and all accessions in the blast output with scores.
# Do this by concat the frag and accession as key, with score as value.

# If the fragment had no blast output, simply choose the best value for that fragment (highest) and add to a solo dictionary.

bestofblast = {}
ravisolo = {}

with open(ravifile) as infile:
	for lines in infile:
		lines = lines.rstrip()
		values = lines.split('\t')
		frag = key[values[0]]
		acc = values[1].split('.')[0]
		score = float(values[2])
		if (frag in yesblast) and (acc in blast[frag]):
			tempkey = frag + acc
			bestofblast[tempkey] = score
		elif (frag in noblast):
			if frag in ravisolo:
				if ravisolo[frag][0] < score:
					ravisolo[frag] = []
					ravisolo[frag].append(score)
					ravisolo[frag].append(acc)
				elif ravisolo[frag][0] == score:
					#print('Match Encountered.')
					ravisolo[frag].append(acc)
				else:
					pass
			else:
				ravisolo[frag] = []
				ravisolo[frag].append(score)
				ravisolo[frag].append(acc)


# Go through blast matches, and pick the highest scoring value and write this to a results output file.
# If there is only one unique accession in the blast output (using sets), simply choose this and attribute BLAST.
# Otherwise, go through the BLAST accessions and choose the highest score.

a = open(resultfile,'w')
a.write('Fragment\tMatching_Accession\tMethod\n')

for frags in blast:
	stb = -100000000000000000000000000000000000000000.0
	if len(list(set(blast[frags]))) == 1:
		a.write(frags + '\t' + blast[frags][0] + '\tBLAST\n')
	else:
		for acc in blast[frags]:
			tempkey = frags + acc
			if bestofblast[tempkey] > stb:
				stb = bestofblast[tempkey]
				bestacc = acc
			else:
				pass
		a.write(frags + '\t' + bestacc + '\t' + 'BLAST+Ravi\n')


# Go through the BLASTless fragments, and choose the best scoring accession using the ravisolo dictionary.

for frags in ravisolo:
	a.write(frags + '\t' + ','.join(ravisolo[frags][1:]) + '\tRavi\n')


