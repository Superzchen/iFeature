#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
from math import sqrt

def multipl(a,b):
	sumofab=0.0
	for i in range(len(a)):
		temp=a[i]*b[i]
		sumofab+=temp
	return sumofab

def corrcoef(x,y):
	n=len(x)
	sum1=sum(x)
	sum2=sum(y)
	sumofxy=multipl(x,y)
	sumofx2 = sum([pow(i,2) for i in x])
	sumofy2 = sum([pow(j,2) for j in y])
	num=sumofxy-(float(sum1)*float(sum2)/n)
	den=sqrt((sumofx2-float(sum1**2)/n)*(sumofy2-float(sum2**2)/n))
	if den != 0:
		return num/den
	else:
		return 0

def pearsonr(encodings, labelFile):
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
	print(labels)

	dataShape = data.shape

	if dataShape[1] != len(features):
		print('Error: inconsistent data shape with feature number.')
		return 0, 'Error: inconsistent data shape with feature number.'
	if dataShape[0] != len(labels):
		print('Error: inconsistent data shape with sample number.')
		return 0, 'Error: inconsistent data shape with sample number.'

	myFea = {}
	for i in range(len(features)):
		array = list(data[:, i])
		try:
			myFea[features[i]] = corrcoef(array, labels)
		except (ValueError, RuntimeWarning) as e:
			return 0, e

	res = []
	res.append(['feature', 'pcc'])
	for key in sorted(myFea.items(), key=lambda item: item[1], reverse=True):
		res.append([key[0], '{0:.3f}'.format(myFea[key[0]])])
	return res, e