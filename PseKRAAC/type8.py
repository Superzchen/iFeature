#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

AAGroup = {
	2:['ADEGKNPQRST', 'CFHILMVWY'],
	3:['ADEGNPST', 'CHKQRW', 'FILMVY'],
	4:['AGNPST', 'CHWY', 'DEKQR', 'FILMV'],
	5:['AGPST', 'CFWY', 'DEN', 'HKQR', 'ILMV'],
	6:['APST', 'CW', 'DEGN', 'FHY', 'ILMV', 'KQR'],
	7:['AGST', 'CW', 'DEN', 'FY', 'HP', 'ILMV', 'KQR'],
	8:['AST', 'CG', 'DEN', 'FY', 'HP', 'ILV', 'KQR', 'MW'],
	9:['AST', 'CW', 'DE', 'FY', 'GN', 'HQ', 'ILV', 'KR', 'MP'],
	10:['AST', 'CW', 'DE', 'FY', 'GN', 'HQ', 'IV', 'KR', 'LM', 'P'],
	11:['AST', 'C', 'DE', 'FY', 'GN', 'HQ', 'IV', 'KR', 'LM', 'P', 'W'],
	12:['AST', 'C', 'DE', 'FY', 'G', 'HQ', 'IV', 'KR', 'LM', 'N', 'P', 'W'],
	13:['AST', 'C', 'DE', 'FY', 'G', 'H', 'IV', 'KR', 'LM', 'N', 'P', 'Q', 'W'],
	14:['AST', 'C', 'DE', 'FL', 'G', 'H', 'IV', 'KR', 'M', 'N', 'P', 'Q', 'W', 'Y'],
	15:['AST', 'C', 'DE', 'F', 'G', 'H', 'IV', 'KR', 'L', 'M', 'N', 'P', 'Q', 'W', 'Y'],
	16:['AT', 'C', 'DE', 'F', 'G', 'H', 'IV', 'KR', 'L', 'M', 'N', 'P', 'Q', 'S', 'W', 'Y'],
	17:['AT', 'C', 'DE', 'F', 'G', 'H', 'IV', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'W', 'Y'],
	18:['A', 'C', 'DE', 'F', 'G', 'H', 'IV', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'Y'],
	19:['A', 'C', 'D', 'E', 'F', 'G', 'H', 'IV', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'Y'],
	20:['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'V', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'Y'],
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
		gDict[i] = 'T8.G.'+str(i+1)
		gNames.append('T8.G.'+str(i+1))

	encodings = []
	if subtype == 'g-gap':
		encodings = gapModel(fastas, myDict, gDict, gNames, ktuple, glValue)
	else:
		encodings = lambdaModel(fastas, myDict, gDict, gNames, ktuple, glValue)

	return encodings