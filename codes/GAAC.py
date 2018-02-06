#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
from collections import Counter

def GAAC(fastas, **kw):
	group = {
		'alphatic': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharge': 'KRH',
		'negativecharge': 'DE',
		'uncharge': 'STCPNQ'
	}

	groupKey = group.keys()

	encodings = []
	header = ['#']
	for key in groupKey:
		header.append(key)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		count = Counter(sequence)
		myDict = {}
		for key in groupKey:
			for aa in group[key]:
				myDict[key] = myDict.get(key, 0) + count[aa]

		for key in groupKey:
			code.append(myDict[key]/len(sequence))
		encodings.append(code)

	return encodings