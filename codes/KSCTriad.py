#!/usr/bin/env python
#_*_coding:utf-8_*_

USAGE = """
USAGE:
	python KSCTriad.py input.fasta <K> output
	
	input.fasta:  the input protein sequence file in fasta format.
	K:            the max space number, integer, defaule: 5
	output:       the encoding file, default: 'encodings.tsv'
"""
import re, sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import readFasta
import saveCode

def CalculateKSCTriad(sequence, gap, features, AADict):
	res = []
	for g in range(gap+1):
		myDict = {}
		for f in features:
			myDict[f] = 0

		for i in range(len(sequence)):
			if i+g+1 < len(sequence) and i+2*g+2<len(sequence):
				fea = AADict[sequence[i]] + '.' + AADict[sequence[i+g+1]]+'.'+AADict[sequence[i+2*g+2]]
				myDict[fea] = myDict[fea] + 1

		maxValue, minValue = max(myDict.values()), min(myDict.values())
		for f in features:
			res.append((myDict[f] - minValue) / maxValue)

	return res

def KSCTriad(fastas, gap = 0, **kw):
	AAGroup = {
		'g1': 'AGV',
		'g2': 'ILFP',
		'g3': 'YMTS',
		'g4': 'HNQW',
		'g5': 'RK',
		'g6': 'DE',
		'g7': 'C'
	}

	myGroups = sorted(AAGroup.keys())

	AADict = {}
	for g in myGroups:
		for aa in AAGroup[g]:
			AADict[aa] = g

	features = [f1 + '.'+ f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]

	encodings = []
	header = ['#']
	for g in range(gap+1):
		for f in features:
			header.append(f+'.gap'+str(g))
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		if len(sequence) < 2*gap + 3:
			print('Error: for "KSCTriad" encoding, the input fasta sequences should be greater than (2*gap+3). \n\n')
			return 0
		code = code + CalculateKSCTriad(sequence, gap, features, AADict)
		encodings.append(code)

	return encodings

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print(USAGE)
		sys.exit(1)
	fastas = readFasta.readFasta(sys.argv[1])
	k = int(sys.argv[2]) if len(sys.argv) >= 3 else 5
	output = sys.argv[3] if len(sys.argv) >= 4 else 'encoding.tsv'
	encodings = KSCTriad(fastas, k)
	saveCode.savetsv(encodings, output)
