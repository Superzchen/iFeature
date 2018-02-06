#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def TPC(fastas, **kw):
	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
	encodings = []
	triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
	header = ['#'] + triPeptides
	encodings.append(header)

	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		tmpCode = [0] * 8000
		for j in range(len(sequence) - 3 + 1):
			tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j+1]]*20 + AADict[sequence[j+2]]] = tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j+1]]*20 + AADict[sequence[j+2]]] +1
		if sum(tmpCode) != 0:
			tmpCode = [i/sum(tmpCode) for i in tmpCode]
		code = code + tmpCode
		encodings.append(code)
	return encodings
