#!/usr/bin/env python
#_*_coding:utf-8_*_

import re, sys, os
from collections import Counter
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import readFasta
import saveCode
import checkFasta

USAGE = """
USAGE:
	python EAAC.py input.fasta <sliding_window> <output>

	input.fasta:      the input protein sequence file in fasta format.
	sliding_window:   the sliding window, integer, defaule: 5
	output:           the encoding file, default: 'encodings.tsv'
	order:            the out order, select from ['alphabetically', 'polarity', 'sideChainVolume' or userDefined] 
"""

def EAAC(fastas, window=5, **kw):
	if checkFasta.checkFasta(fastas) == False:
		print('Error: for "EAAC" encoding, the input fasta sequences should be with equal length. \n\n')
		return 0

	if window < 1:
		print('Error: the sliding window should be greater than zero' + '\n\n')
		return 0

	if checkFasta.minSequenceLength(fastas) < window:
		print('Error: all the sequence length should be larger than the sliding window :' + str(window) + '\n\n')
		return 0

	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
	#AA = 'ARNDCQEGHILKMFPSTWYV'
	encodings = []
	header = ['#']
	for w in range(1, len(fastas[0][1]) - window + 2):
		for aa in AA:
			header.append('SW.'+str(w)+'.'+aa)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], i[1]
		code = [name]
		for j in range(len(sequence)):
			if j < len(sequence) and j + window <= len(sequence):
				count = Counter(re.sub('-', '', sequence[j:j+window]))
				for key in count:
					count[key] = count[key] / len(re.sub('-', '', sequence[j:j+window]))
				for aa in AA:
					code.append(count[aa])
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
	sw = int(sys.argv[2]) if len(sys.argv) >= 3 else 5
	output = sys.argv[3] if len(sys.argv) >= 4 else 'encoding.tsv'

	if len(sys.argv) >= 5:
		if sys.argv[4] in myAAorder:
			kw['order'] = myAAorder[sys.argv[4]]
		else:
			tmpOrder = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', sys.argv[4])
			kw['order'] = tmpOrder if len(tmpOrder) == 20 else 'ACDEFGHIKLMNPQRSTVWY'

	encodings = EAAC(fastas, sw, **kw)
	saveCode.savetsv(encodings, output)

