#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import checkFasta


def ASA(fastas, **kw):
	if checkFasta.checkFasta(fastas) == False:
		print('Error: for "ASA" encoding, the input fasta sequences should be with equal length. \n\n')
		return 0

	encodings = []
	header = ['#']
	for p in range(1, len(fastas[0][1])+1):
		header.append('ASA.F' + str(p))
	encodings.append(header)

	disDir = kw['path']
	if disDir == None:
		print('Error: please specify the directory of predicted protein ASA file by "--path"')
		return 0
	for i in fastas:
		name, sequence = i[0], i[1]
		code = [name]
		if os.path.exists(disDir + '/' + name + '.dis') == False:
			print('Error: the predicted ASA information file (.spXout) for protein ' + name + ' does not exist.')
			return 0

		with open(disDir + '/' + name + '.spXout') as f:
			records = f.readlines()[1:]

		proteinSeq = ''
		asaValue = []
		for line in records:
			array = line.strip().split() if line.strip() != '' else None
			proteinSeq = proteinSeq + array[1]
			asaValue.append(array[10])
		pos = proteinSeq.find(sequence)
		if pos == -1:
			print('Warning: could not find the peptide in proteins.\n\n')
		else:
			for p in range(pos, pos+len(sequence)):
				code.append(asaValue[p])
		encodings.append(code)

	return encodings