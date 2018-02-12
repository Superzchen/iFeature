#!/usr/bin/env python
#_*_coding:utf-8_*_

import scipy.cluster.hierarchy as sch
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt

def hcluster(encodings, image='hcluster.png', **kw):
	if kw['sof'] == 'sample':
		encodings = np.array(encodings)[1:]
	else:
		encodings = np.array(encodings).T[1:]

	if len(encodings) < 3:
		return 0, 'sample number should be greater than 3.'

	#encodings = np.array(encodings)[1:]
	data = encodings[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)

	names = []
	for i in encodings:
		names.append(i[0])

	e = ''
	try:
		disMat = sch.distance.pdist(data, 'euclidean')
		Z = sch.linkage(disMat, method='average')
	except ValueError as e:
		print(e)
		return 0, e
	P = sch.dendrogram(Z, labels=names, leaf_rotation=270, leaf_font_size=8)
	plt.savefig(image)
	plt.close()
	cluster = sch.fcluster(Z, 1, 'inconsistent')
	res = []
	for i in range(len(data)):
		res.append([encodings[i][0], cluster[i]])
	return res, e