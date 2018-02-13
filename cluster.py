#!/usr/bin/env python
#_*_coding:utf-8_*_

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import pylab
from clusters import *

if __name__ == '__main__':
	parser = argparse.ArgumentParser(usage="it's usage tip.",
									 description="cluster for the generated numerical represention")
	parser.add_argument("--file", required=True, help="input encoding file")
	parser.add_argument("--type", required=True,
						choices=['kmeans', 'hcluster', 'apc', 'meanshift', 'dbscan'], help="select cluster method")
	parser.add_argument("--sof", choices = ['sample', 'feature'], help="cluster for sample or feature, default: sample")
	parser.add_argument("--nclusters", help="specify the cluster number for kmeans cluster method. default: 3")
	parser.add_argument("--out", help="output file")
	args = parser.parse_args()

	kw = {'nclusters': args.nclusters if args.nclusters != None else 3, 'sof':args.sof if args.sof !=None else 'sample'}

	output = args.out if args.out != None else '%s_cluster.txt' %args.type
	encodings = readCode.readCode(args.file)
	myFun = args.type + '.' + args.type + '(encodings, **kw)'
	print('Cluster type: ' + args.type)
	myCluster, e = eval(myFun)
	saveCluster.saveCluster(myCluster, e, output)
	## t-sne plot

	if myCluster != 0:
		if kw['sof'] == 'sample':
			data = np.array(encodings)[1:, 1:].astype(float)
		else:
			data = np.array(encodings).T[1:, 1:].astype(float)
		labels = np.array(myCluster)[0:, 1:].reshape(-1, )
		e = ''
		try:
			Y = tsne.tsne(data, 2, 50, 20.0)
		except RuntimeWarning as e:
			Y = pca.pca(data, n_components = 2)

		df = pd.DataFrame({'X': Y[:, 0], 'Y': Y[:, 1], 'L': labels})

		mySet = set(labels)
		if len(mySet) > 5:
			pylab.scatter(Y[:, 0], Y[:, 1], 20, labels)
		else:
			for l in mySet:
				newData = df.loc[df.loc[:, "L"] == l, :]
				pylab.scatter(np.array(newData.X), np.array(newData.Y), 20, label="Cluster_%s" % l)
		pylab.legend(loc='best')
		#pylab.show()
		pylab.savefig('%s.png' %output)
		pylab.close()