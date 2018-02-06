#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re

def calculateDicorderContent(pos, endPos, disValue):
	newValues = disValue[pos: endPos]
	return [newValues.count('D')/(endPos - pos), newValues.count('O')/(endPos-pos)]

def DisorderC(fastas, **kw):
	encodings = []
	header = ['#', 'disorder-content', 'order-content']
	encodings.append(header)

	disDir = kw['path']
	if disDir == None:
		print('Error: please specify the directory of predicted protein disorder files by "--path"')
		return 0

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
		for line in records:
			array = line.rstrip().split() if line.rstrip() != '' else None
			key = array[3] if array[3] == 'D' else 'O'
			proteinSeq = proteinSeq + array[1]
			disValue.append(key)

		pos = proteinSeq.find(sequence)
		if pos == -1:
			print('Warning: could not find the peptide in proteins.\n\n')
		else:
			code = code + calculateDicorderContent(pos, pos+len(sequence), disValue)
		encodings.append(code)

	return encodings