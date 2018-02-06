#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import checkFasta

def DisorderB(fastas, **kw):
	if checkFasta.checkFasta(fastas) == False:
		print('Error: for "DisorderB" encoding, the input fasta sequences should be with equal length. \n\n')
		return 0

	disDir = kw['path']
	if disDir == None:
		print('Error: please specify the directory of predicted protein disorder files by "--path"')
		return 0

	encodings = []
	header = ['#']
	for p in range(1, 2*len(fastas[0][1])+1):
		header.append('disorderB.F' + str(p))

	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], i[1]
		code = [name]
		if os.path.exists(disDir + '/' + name + '.dis') == False:
			print('Error: the predicted disorder information file (.dis) for protein ' + name + ' does not exist.')
			return 0
		with open(disDir + '/' + name + '.dis') as f:
			records = f.readlines()
		tag = 0
		for i in range(len(records)):
			if re.search('^-------', records[i]):
				tag = i
				break
		records = records[tag+1:-1]

		proteinSeq = ''
		disValue = []
		myDict = {'D':[0, 1], 'O':[1, 0]}
		for line in records:
			array = line.rstrip().split() if line.rstrip() != '' else None
			key = array[3] if array[3] == 'D' else 'O'
			proteinSeq = proteinSeq + array[1]
			disValue.append(key)

		pos = proteinSeq.find(sequence)
		if pos == -1:
			print('Warning: could not find the peptide in proteins.\n\n')
		else:
			for p in range(pos, pos+len(sequence)):
				code = code + myDict[disValue[p]]
		encodings.append(code)

	return encodings