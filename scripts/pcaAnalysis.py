#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
import pandas as pd
pPath = re.sub(r'scripts$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from clusters import readCode
import argparse
import numpy as np
from sklearn.decomposition import PCA
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def pacAnalysis(encodings, n_components = 2):
	encodings = np.array(encodings)[1:]
	data = encodings[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)
	newData = PCA(n_components = n_components).fit_transform(data)
	pca = []
	for i in range(len(data)):
		pca.append([encodings[i][0]] + list(newData[i]))
	return pca

def savePCA(pca, file='PCA.txt'):
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

def pcaPlot(pca, label, file='pca.png'):
	pca = np.array(pca)[:]
	data = pca[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)

	if len(label) == 0:
		plt.scatter(data[:, 0], data[:, 1], 20, c='r')
	else:
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
	parser.add_argument("--file", required=True, help="data file for PCA")
	parser.add_argument("--out", help="output file of PCA analysis.")
	parser.add_argument("--ncomponents", help="number of n components, default 2")
	parser.add_argument("--label", help="optimal, the label file for the input file")
	args = parser.parse_args()
	output = args.out if args.out != None else 'PCA.txt'
	n_components = int(args.ncomponents) if args.ncomponents != None else 2
	data = readCode.readCode(args.file)
	label = []
	if args.label != None:
		with open(args.label) as f:
			records = f.readlines()
		for line in records:
			array = line.strip().split() if line.strip() != '' else None
			label.append(array[1])

	pca = pacAnalysis(data, n_components)
	savePCA(pca, output)
	pcaPlot(pca, label, '%s.png' %output)



