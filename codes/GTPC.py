#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def GTPC(fastas, **kw):
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}

	groupKey = group.keys()
	baseNum = len(groupKey)
	triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	encodings = []
	header = ['#'] + triple
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])

		code = [name]
		myDict = {}
		for t in triple:
			myDict[t] = 0

		sum = 0
		for j in range(len(sequence) - 3 + 1):
			myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1
			sum = sum +1

		if sum == 0:
			for t in triple:
				code.append(0)
		else:
			for t in triple:
				code.append(myDict[t]/sum)
		encodings.append(code)

	return encodings