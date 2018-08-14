#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import readFasta
import saveCode
import checkFasta

USAGE = """
USAGE:
	python CKSAAP.py input.fasta <k_space> <output>

	input.fasta:      the input protein sequence file in fasta format.
	k_space:          the gap of two amino acids, integer, defaule: 5
	output:           the encoding file, default: 'encodings.tsv'
"""

def CKSAAP(fastas, gap=5, **kw):
	if gap < 0:
		print('Error: the gap should be equal or greater than zero' + '\n\n')
		return 0

	if checkFasta.minSequenceLength(fastas) < gap+2:
		print('Error: all the sequence length should be larger than the (gap value) + 2 = ' + str(gap+2) + '\n\n')
		return 0

	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
	encodings = []
	aaPairs = []
	for aa1 in AA:
		for aa2 in AA:
			aaPairs.append(aa1 + aa2)
	header = ['#']
	for g in range(gap+1):
		for aa in aaPairs:
			header.append(aa + '.gap' + str(g))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], i[1]
		code = [name]
		for g in range(gap+1):
			myDict = {}
			for pair in aaPairs:
				myDict[pair] = 0
			sum = 0
			for index1 in range(len(sequence)):
				index2 = index1 + g + 1
				if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[index2] in AA:
					myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
					sum = sum + 1
			for pair in aaPairs:
				code.append(myDict[pair] / sum)
		encodings.append(code)
	return encodings

if __name__ == '__main__':
	myAAorder = {
		'alphabetically': 'ACDEFGHIKLMNPQRSTVWY',
		'polarity': 'DENKRQHSGTAPYVMCWIFL',
		'sideChainVolume': 'GASDPCTNEVHQILMKRFYW',
	}
	kw = {'order': 'ACDEFGHIKLMNPQRSTVWY'}

	if len(sys.argv) == 1:
		print(USAGE)
		sys.exit(1)
	fastas = readFasta.readFasta(sys.argv[1])
	gap = int(sys.argv[2]) if len(sys.argv) >= 3 else 5
	output = sys.argv[3] if len(sys.argv) >= 4 else 'encoding.tsv'

	if len(sys.argv) >= 5:
		if sys.argv[4] in myAAorder:
			kw['order'] = myAAorder[sys.argv[4]]
		else:
			tmpOrder = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', sys.argv[4])
			kw['order'] = tmpOrder if len(tmpOrder) == 20 else 'ACDEFGHIKLMNPQRSTVWY'
	encodings = CKSAAP(fastas, gap, **kw)
	saveCode.savetsv(encodings, output)