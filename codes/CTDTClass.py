#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import readFasta
import saveCode

USAGE = """
USAGE:
	python CTDTClass.py input.fasta output amino_acids_group_1 amino_acids_group_2 ... amino_acids_group_N

	input.fasta:                 the input protein sequence file in fasta format.	
	output:                      the encoding file.
	amino_acids_group_x          the amino acids groups.

EXAMPLE:
	python CTDTClass.py example/test-protein.txt CTDTClass.tsv RKEDQN GASTPHY CLVIMFW
"""

def CTDTClass(fastas, groups):
	encodings = []
	header = ['#']
	myOrder = []
	for g in range(len(groups) - 1):
		for g1 in range(g+1, len(groups)):
			header.append('Tr.%d%d%d%d' %(g+1, g1+1, g1+1, g+1))
			myOrder.append('%d%d%d%d' % (g + 1, g1 + 1, g1 + 1, g + 1))
	encodings.append(header)

	myDict = {}
	for g in range(len(groups)):
		for aa in groups[g]:
			myDict[aa] = g + 1

	myValues = list(set(myDict.values()))
	myKey = {}
	for v in range(len(myValues)):
		for v1 in range(v, len(myValues)):
			myKey['%d%d' % (myValues[v], myValues[v1])] = '%d%d%d%d' % (
			myValues[v], myValues[v1], myValues[v1], myValues[v])
			myKey['%d%d' % (myValues[v1], myValues[v])] = '%d%d%d%d' % (
			myValues[v], myValues[v1], myValues[v1], myValues[v])

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]

		myCount = {}
		for pair in aaPair:
			tmp = myKey['%d%d' %(myDict[pair[0]], myDict[pair[1]])]
			if tmp in myCount:
				myCount[tmp] = myCount[tmp] + 1
			else:
				myCount[tmp] = 1

		for key in myOrder:
			if key in myCount:
				code.append(myCount[key] / len(aaPair))
			else:
				code.append(0)
		encodings.append(code)
	return encodings

if __name__ == '__main__':
	if len(sys.argv) < 5:
		print(USAGE)
		sys.exit(1)

	groups = sys.argv[3:]
	myStr = ''.join(groups)
	myStr = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', myStr)
	if len(myStr) != 20 or len(set(myStr)) != 20:
		print('\nERROR: The amino acid must be no-repeat in each groups and the sum is 20!\n')
	fastas = readFasta.readFasta(sys.argv[1])
	encodings = CTDTClass(fastas, groups)
	saveCode.savetsv(encodings, sys.argv[2])