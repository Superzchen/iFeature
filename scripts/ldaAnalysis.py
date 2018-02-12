#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
import pandas as pd
pPath = re.sub(r'scripts$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from clusters import readCode
import argparse
import numpy as np
from sklearn.decomposition import LatentDirichletAllocation
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def ldaAnalysis(encodings, labelFile, n_components = 2):
	with open(labelFile) as f:
		records = f.readlines()[1:]
	myDict = {}
	try:
		for i in records:
			array = i.rstrip().split() if i.strip() != '' else None
			myDict[array[0]] = int(array[1])
	except IndexError as e:
		return 0, e

	labels = []
	for i in encodings:
		labels.append(myDict.get(i[0], 0))

	encodings = np.array(encodings)[1:]
	data = encodings[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)
	newData = LatentDirichletAllocation(n_components=n_components, max_iter=5,
                                learning_method='batch',
                                learning_offset=50.,
                                random_state=0).fit_transform(data, labels)
	lda = []
	for i in range(len(data)):
		lda.append([encodings[i][0]] + list(newData[i]))
	return lda

def saveLDA(pca, file='LDA.txt'):
	with open(file, 'w') as f:
		f.write('protein\t')
		for i in range(1, len(pca[0])):
			f.write('pc.' + str(i) + '\t')
		f.write('\n')
		for i in pca:
			for j in i:
				f.write(str(j) + '\t')
			f.write('\n')
	return None

def ldaPlot(lda, label, file='lda.png'):
	lda = np.array(lda)[:]
	data = lda[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)

	df = pd.DataFrame({'X': data[:, 0], 'Y': data[:, 1], 'L': label})
	mySet = set(label)
	for l in mySet:
		newData = df.loc[df.loc[:, "L"] == l, :]
		plt.scatter(np.array(newData.X), np.array(newData.Y), 20, label="%s" % l)
		plt.legend(loc='best')
	plt.savefig(file)
	return None

if __name__ == '__main__':
	parser = argparse.ArgumentParser(usage="it's usage tip.", description="PCA and plot")
	parser.add_argument("--file", required=True, help="data file for LDA")
	parser.add_argument("--label", required=True, help="sample label file")
	parser.add_argument("--out", help="output file of LDA analysis.")
	parser.add_argument("--ncomponents", help="number of n components, default 2")
	args = parser.parse_args()
	n_components = int(args.ncomponents) if args.ncomponents != None else 2
	ldaFile = args.out if args.out != None else 'lda.txt'

	data = readCode.readCode(args.file)
	label = []
	if args.label != None:
		with open(args.label) as f:
			records = f.readlines()
		for line in records:
			array = line.strip().split() if line.strip() != '' else None
			label.append(array[1])

	lda = ldaAnalysis(data, args.label, n_components)
	saveLDA(lda, ldaFile)
	ldaPlot(lda, label, '%s.png' %ldaFile)



