#!/usr/bin/env python
#_*_coding:utf-8_*_

import os

def calculateSSE(pos, end, SSE):
	newValues = SSE[pos:end]
	return [newValues.count('H')/(end-pos), newValues.count('E')/(end-pos), newValues.count('C')/(end-pos)]

def SSEC(fastas, **kw):
	ssDir = kw['path']
	if ssDir == None:
		print('Error: please specify the directory of predicted protein disorder files by "--path"')
		return 0

	encodings = []
	header = ['#', 'H', 'E', 'C']
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], i[1]
		code = [name]

		if os.path.exists(ssDir + '/' + name + '.ss2') == True:
			with open(ssDir + '/' + name + '.ss2') as f:
				records = f.readlines()[2:]
		elif os.path.exists(ssDir + '/' + name + '.spXout') == True:
			with open(ssDir + '/' + name + '.spXout') as f:
				records = f.readlines()[1:]
		else:
			print('Error: the predicted secondary structure (.ss2 or .spXout) for protein ' + name + ' does not exist.')
			return 0

		proteinSeq = ''
		SSE = []
		for line in records:
			array = line.strip().split() if line.rstrip() != '' else None
			proteinSeq = proteinSeq + array[1]
			SSE.append(array[2])

		pos = proteinSeq.find(sequence)
		if pos == -1:
			print('Warning: could not find the peptide in proteins.\n\n')
		else:
			code = code + calculateSSE(pos, pos+len(sequence), SSE)
		encodings.append(code)

	return encodings