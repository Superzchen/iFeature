#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
import pandas as pd
import math

binBox = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

def calProb(array):
	myProb = {}
	myClass = set(array)
	for i in myClass:
		myProb[i] = array.count(i) / len(array)
	return myProb

def jointProb(newArray, labels):
	myJointProb = {}
	for i in range(len(labels)):
		myJointProb[str(newArray[i]) + '-' + str(labels[i])] = myJointProb.get(str(newArray[i]) + '-' + str(labels[i]), 0) + 1

	for key in myJointProb:
		myJointProb[key] = myJointProb[key] / len(labels)
	return myJointProb

def MIC(encodings, labelFile):
	features = encodings[0][1:]
	encodings = np.array(encodings)[1:]
	data = encodings[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)

	e = ''
	if shape[0] < 5 or shape[1] < 2:
		return 0, e

	with open(labelFile) as f:
		records = f.readlines()
	myDict = {}
	try:
		for i in records:
			array = i.rstrip().split() if i.strip() != '' else None
			myDict[array[0]] = int(array[1])
	except IndexError as e:
		print(e)
		return 0, e

	labels = []
	for i in encodings:
		labels.append(myDict.get(i[0], 0))

	dataShape = data.shape

	if dataShape[1] != len(features):
		print('Error: inconsistent data shape with feature number.')
		return 0, 'Error: inconsistent data shape with feature number.'

	if dataShape[0] != len(labels):
		print('Error: inconsistent data shape with sample number.')
		return 0, 'Error: inconsistent data shape with sample number.'

	sampleNumber = len(data)
	labelClass = set(labels)
	probY = calProb(labels)
	myFea = {}
	for i in range(len(features)):
		array = data[:, i]
		newArray = list(pd.cut(array, len(binBox), labels= binBox))
		binBoxClass = set(newArray)
		probX = calProb(newArray)
		probXY = jointProb(newArray, labels)
		mic = 0
		for x in probX.keys():
			for y in probY.keys():
				if str(x) + '-' + str(y) in probXY:
					mic = mic + probXY[str(x) + '-' + str(y)] * math.log(probXY[str(x) + '-' + str(y)]/(probX[x] * probY[y]), 2)
		myFea[features[i]] = mic

	res = []
	res.append(['feature', 'MIC-value'])
	for key in sorted(myFea.items(), key=lambda item:item[1], reverse=True):
		res.append([key[0], '{0:.3f}'.format(myFea[key[0]])])
	return res, e
