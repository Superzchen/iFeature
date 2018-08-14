#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import readFasta
import saveCode
import checkFasta

USAGE = """
USAGE:
	python CKSAAGP.py input.fasta <k_space> <output>

	input.fasta:      the input protein sequence file in fasta format.
	k_space:          the gap of two amino acids, integer, defaule: 5
	output:           the encoding file, default: 'encodings.tsv'
"""

def generateGroupPairs(groupKey):
	gPair = {}
	for key1 in groupKey:
		for key2 in groupKey:
			gPair[key1+'.'+key2] = 0
	return gPair

def CKSAAGP(fastas, gap = 5, **kw):
	if gap < 0:
		print('Error: the gap should be equal or greater than zero' + '\n\n')
		return 0

	if checkFasta.minSequenceLength(fastas) < gap+2:
		print('Error: all the sequence length should be greater than the (gap value) + 2 = ' + str(gap+2) + '\n\n')
		return 0

	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}

	AA = 'ARNDCQEGHILKMFPSTWYV'

	groupKey = group.keys()

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	gPairIndex = []
	for key1 in groupKey:
		for key2 in groupKey:
			gPairIndex.append(key1+'.'+key2)

	encodings = []
	header = ['#']
	for g in range(gap + 1):
		for p in gPairIndex:
			header.append(p+'.gap'+str(g))
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		for g in range(gap + 1):
			gPair = generateGroupPairs(groupKey)
			sum = 0
			for p1 in range(len(sequence)):
				p2 = p1 + g + 1
				if p2 < len(sequence) and sequence[p1] in AA and sequence[p2] in AA:
					gPair[index[sequence[p1]]+'.'+index[sequence[p2]]] = gPair[index[sequence[p1]]+'.'+index[sequence[p2]]] + 1
					sum = sum + 1

			if sum == 0:
				for gp in gPairIndex:
					code.append(0)
			else:
				for gp in gPairIndex:
					code.append(gPair[gp] / sum)

		encodings.append(code)

	return encodings

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print(USAGE)
		sys.exit(1)
	fastas = readFasta.readFasta(sys.argv[1])
	gap = int(sys.argv[2]) if len(sys.argv) >= 3 else 5
	output = sys.argv[3] if len(sys.argv) >= 4 else 'encoding.tsv'
	encodings = CKSAAGP(fastas, gap)
	saveCode.savetsv(encodings, output)
