#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, platform, os, re
import numpy as np
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import checkFasta
import readFasta
import saveCode

USAGE = """
USAGE:
	python SOCNumber.py input.fasta <nlag> <output>

	input.fasta:      the input protein sequence file in fasta format.
	nlag:             the nlag value, integer, defaule: 30
	output:           the encoding file, default: 'encodings.tsv'
"""

def SOCNumber(fastas, nlag=30, **kw):
	if checkFasta.minSequenceLengthWithNormalAA(fastas) < nlag + 1:
		print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
		return 0

	dataFile = re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + r'\data\Schneider-Wrede.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/Schneider-Wrede.txt'
	dataFile1 = re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + r'\data\Grantham.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/Grantham.txt'
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	AA1 = 'ARNDCQEGHILKMFPSTWYV'

	DictAA = {}
	for i in range(len(AA)):
		DictAA[AA[i]] = i

	DictAA1 = {}
	for i in range(len(AA1)):
		DictAA1[AA1[i]] = i

	with open(dataFile) as f:
		records = f.readlines()[1:]
	AADistance = []
	for i in records:
		array = i.rstrip().split()[1:] if i.rstrip() != '' else None
		AADistance.append(array)
	AADistance = np.array(
		[float(AADistance[i][j]) for i in range(len(AADistance)) for j in range(len(AADistance[i]))]).reshape((20, 20))

	with open(dataFile1) as f:
		records = f.readlines()[1:]
	AADistance1 = []
	for i in records:
		array = i.rstrip().split()[1:] if i.rstrip() != '' else None
		AADistance1.append(array)
	AADistance1 = np.array(
		[float(AADistance1[i][j]) for i in range(len(AADistance1)) for j in range(len(AADistance1[i]))]).reshape(
		(20, 20))

	encodings = []
	header = ['#']
	for n in range(1, nlag + 1):
		header.append('Schneider.lag' + str(n))
	for n in range(1, nlag + 1):
		header.append('gGrantham.lag' + str(n))
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		for n in range(1, nlag + 1):
			code.append(sum(
				[AADistance[DictAA[sequence[j]]][DictAA[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]) / (
						len(sequence) - n))

		for n in range(1, nlag + 1):
			code.append(sum([AADistance1[DictAA1[sequence[j]]][DictAA1[sequence[j + n]]] ** 2 for j in
							 range(len(sequence) - n)]) / (len(sequence) - n))
		encodings.append(code)
	return encodings

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print(USAGE)
		sys.exit(1)
	fastas = readFasta.readFasta(sys.argv[1])
	nlag = int(sys.argv[2]) if len(sys.argv) >= 3 else 30
	output = sys.argv[3] if len(sys.argv) >= 4 else 'encoding.tsv'
	encodings = SOCNumber(fastas, nlag)
	saveCode.savetsv(encodings, output)