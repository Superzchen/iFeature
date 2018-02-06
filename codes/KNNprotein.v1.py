#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
import math
import numpy as np
pPath = re.sub(r'codes$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from codes import readFasta
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

def CalculateSimilarity(sequence1, sequence2):
	alignment = pairwise2.align.globalds(sequence1, sequence2, blosum62, -10, -0.5)[0]
	sum = 0
	for i in range(len(alignment[0])):
		if alignment[0][i] == alignment[1][i]:
			sum = sum + 1
	return 2 * sum / (len(sequence1) + len(sequence2))

def CalculateContent(mySimilarity, j, myLabelSets):
	content = []
	myDict = {}
	for i in myLabelSets:
		myDict[i] = 0
	for i in range(j):
		myDict[mySimilarity[i][0]] = myDict[mySimilarity[i][0]] + 1
	for i in myLabelSets:
		content.append(myDict[myLabelSets[i]] / j)
	return content

def KNNprotein(fastas, **kw):
	trainFile = kw['train']
	labelFile = kw['label']

	if os.path.exists(labelFile) == False:
		print('Error: the label file does not exist.')
		sys.exit(1)

	if trainFile == None or labelFile == None:
		print('Error: please specify the directory of train file ["--train"] and the label file ["--label"]')
		sys.exit(1)
	trainData = readFasta.readFasta(trainFile)
	with open(labelFile) as f:
		records = f.readlines()
	myLabel = {}
	for i in records:
		array = i.rstrip().split() if i.strip() != '' else None
		myLabel[array[0]] = int(array[1])
	myLabelSets = list(set(myLabel.values()))

	if len(trainData) != len(myLabel):
		print('ERROR: the inconsistent sample number in train and label file.')
		sys.exit(1)

	kValues = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15,
			   0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30]
	kNum = []
	for i in kValues:
		kNum.append(math.ceil(len(trainData) * i))

	encodings = []
	header = ['#']
	for k in kValues:
		for l in myLabelSets:
			header.append('Top' + str(k) + '.label' + str(l))
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', i[1])
		code = [name]
		mySimilarity = []
		for j in range(len(trainData)):
			if name != trainData[j][0]:
				mySimilarity.append([myLabel[trainData[j][0]], CalculateSimilarity(re.sub('[^ARNDCQEGHILKMFPSTWYV]', '', trainData[j][1]), sequence)])
		mySimilarity = np.array(mySimilarity)
		mySimilarity = mySimilarity[np.lexsort(-mySimilarity.T)]
		for j in kNum:
			code = code + CalculateContent(mySimilarity, j, myLabelSets)
		encodings.append(code)
	return encodings

