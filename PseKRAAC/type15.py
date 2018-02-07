#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

AAGroup = {
	2:['MFILVAW', 'CYQHPGTSNRKDE'],
	3:['MFILVAW', 'CYQHPGTSNRK', 'DE'],
	4:['MFILV', 'ACW', 'YQHPGTSNRK', 'DE'],
	5:['MFILV', 'ACW', 'YQHPGTSN', 'RK', 'DE'],
	6:['MFILV', 'A', 'C', 'WYQHPGTSN', 'RK', 'DE'],
	7:['MFILV', 'A', 'C', 'WYQHP', 'GTSN', 'RK', 'DE'],
	8:['MFILV', 'A', 'C', 'WYQHP', 'G', 'TSN', 'RK', 'DE'],
	9:['MF', 'ILV', 'A', 'C', 'WYQHP', 'G', 'TSN', 'RK', 'DE'],
	10:['MF', 'ILV', 'A', 'C', 'WYQHP', 'G', 'TSN', 'RK', 'D', 'E'],
	11:['MF', 'IL', 'V', 'A', 'C', 'WYQHP', 'G', 'TSN', 'RK', 'D', 'E'],
	12:['MF', 'IL', 'V', 'A', 'C', 'WYQHP', 'G', 'TS', 'N', 'RK', 'D', 'E'],
	13:['MF', 'IL', 'V', 'A', 'C', 'WYQHP', 'G', 'T', 'S', 'N', 'RK', 'D', 'E'],
	14:['MF', 'I', 'L', 'V', 'A', 'C', 'WYQHP', 'G', 'T', 'S', 'N', 'RK', 'D', 'E'],
	15:['MF', 'IL', 'V', 'A', 'C', 'WYQ', 'H', 'P', 'G', 'T', 'S', 'N', 'RK', 'D', 'E'],
	16:['MF', 'I', 'L', 'V', 'A', 'C', 'WYQ', 'H', 'P', 'G', 'T', 'S', 'N', 'RK', 'D', 'E'],
	20:['M', 'F', 'I', 'L', 'V', 'A', 'C', 'W', 'Y', 'Q', 'H', 'P', 'G', 'T', 'S', 'N', 'R', 'K', 'D', 'E'],
}
'''
for key in AAGroup:
	aa = set(list(''.join(AAGroup[key])))
	if len(aa) != 20:
		print(key, 1)
	if key != len(AAGroup[key]):
		print(key, 2)
	print(''.join(AAGroup[key]))
'''
def gapModel(fastas, myDict, gDict, gNames, ktuple, glValue):
	encodings = []
	header = ['#']

	if ktuple == 1:
		header = header + [g + '_gap' + str(glValue) for g in gNames]
		encodings.append(header)
		for i in fastas:
			name, sequence = i[0], re.sub('-', '', i[1])
			code = [name]
			numDict = {}
			for j in range(0, len(sequence), glValue+1):
				numDict[gDict[myDict[sequence[j]]]] = numDict.get(gDict[myDict[sequence[j]]], 0) + 1

			for g in gNames:
				code.append(numDict.get(g, 0))
			encodings.append(code)

	if ktuple == 2:
		header = header + [g1 + '_' + g2 + '_gap' + str(glValue) for g1 in gNames for g2 in gNames]
		encodings.append(header)
		for i in fastas:
			name, sequence = i[0], re.sub('-', '', i[1])
			code = [name]
			numDict = {}
			for j in range(0, len(sequence), glValue + 1):
				if j+1 < len(sequence):
					numDict[gDict[myDict[sequence[j]]]+'_'+gDict[myDict[sequence[j+1]]]] = numDict.get(gDict[myDict[sequence[j]]]+'_'+gDict[myDict[sequence[j+1]]], 0) + 1

			for g in [g1+'_'+g2 for g1 in gNames for g2 in gNames]:
				code.append(numDict.get(g, 0))
			encodings.append(code)

	if ktuple == 3:
		header = header + [g1 + '_' + g2 + '_' + g3 + '_gap' + str(glValue) for g1 in gNames for g2 in gNames for g3 in gNames]
		encodings.append(header)
		for i in fastas:
			name, sequence = i[0], re.sub('-', '', i[1])
			code = [name]
			numDict = {}
			for j in range(0, len(sequence), glValue + 1):
				if j + 1 < len(sequence) and j + 2 < len(sequence):
					numDict[gDict[myDict[sequence[j]]] + '_' + gDict[myDict[sequence[j + 1]]] + '_' + gDict[myDict[sequence[j + 2]]]] = numDict.get(gDict[myDict[sequence[j]]] + '_' + gDict[myDict[sequence[j + 1]]] + '_' + gDict[myDict[sequence[j + 2]]], 0) + 1

			for g in [g1 + '_' + g2 + '_' +g3 for g1 in gNames for g2 in gNames for g3 in gNames]:
				code.append(numDict.get(g, 0))
			encodings.append(code)

	return encodings

def lambdaModel(fastas, myDict, gDict, gNames, ktuple, glValue):
	if glValue == 0:
		print('Warning: the lambda value should not be zero in "lambda correlation" model')
		return 0

	encodings = []
	header = ['#']

	if ktuple == 1:
		header = header + [g + '_LC' + str(glValue) for g in gNames]
		encodings.append(header)
		for i in fastas:
			name, sequence = i[0], re.sub('-', '', i[1])
			code = [name]
			numDict = {}
			for j in range(0, len(sequence)):
				numDict[gDict[myDict[sequence[j]]]] = numDict.get(gDict[myDict[sequence[j]]], 0) + 1

			for g in gNames:
				code.append(numDict.get(g, 0))
			encodings.append(code)

	if ktuple == 2:
		header = header + [g1 + '_' + g2 + '_LC' + str(glValue) for g1 in gNames for g2 in gNames]
		encodings.append(header)
		for i in fastas:
			name, sequence = i[0], re.sub('-', '', i[1])
			code = [name]
			numDict = {}
			for j in range(0, len(sequence)):
				if j + glValue < len(sequence):
					numDict[gDict[myDict[sequence[j]]] + '_' + gDict[myDict[sequence[j + glValue]]]] = numDict.get(
						gDict[myDict[sequence[j]]] + '_' + gDict[myDict[sequence[j + glValue]]], 0) + 1

			for g in [g1 + '_' + g2 for g1 in gNames for g2 in gNames]:
				code.append(numDict.get(g, 0))
			encodings.append(code)

	if ktuple == 3:
		header = header + [g1 + '_' + g2 + '_' + g3 + '_LC' + str(glValue) for g1 in gNames for g2 in gNames for g3 in
						   gNames]
		encodings.append(header)
		for i in fastas:
			name, sequence = i[0], re.sub('-', '', i[1])
			code = [name]
			numDict = {}
			for j in range(0, len(sequence)):
				if j + glValue < len(sequence) and j + 2*glValue < len(sequence):
					numDict[gDict[myDict[sequence[j]]] + '_' + gDict[myDict[sequence[j + glValue]]] + '_' + gDict[
						myDict[sequence[j + 2*glValue]]]] = numDict.get(
						gDict[myDict[sequence[j]]] + '_' + gDict[myDict[sequence[j + glValue]]] + '_' + gDict[
							myDict[sequence[j + 2*glValue]]], 0) + 1

			for g in [g1 + '_' + g2 + '_' + g3 for g1 in gNames for g2 in gNames for g3 in gNames]:
				code.append(numDict.get(g, 0))
			encodings.append(code)

	return encodings

def type1(fastas, subtype, raactype, ktuple, glValue):
	if raactype not in AAGroup:
		print('Error: the "--raactype" value is not correct.')
		return 0

	# index each amino acids to their group
	myDict = {}
	for i in range(len(AAGroup[raactype])):
		for aa in AAGroup[raactype][i]:
			myDict[aa] = i

	gDict = {}
	gNames = []
	for i in range(len(AAGroup[raactype])):
		gDict[i] = 'T15.G.'+str(i+1)
		gNames.append('T15.G.'+str(i+1))

	encodings = []
	if subtype == 'g-gap':
		encodings = gapModel(fastas, myDict, gDict, gNames, ktuple, glValue)
	else:
		encodings = lambdaModel(fastas, myDict, gDict, gNames, ktuple, glValue)

	return encodings