#!/usr/bin/env python
#_*_coding:utf-8_*_

from sklearn.cluster import AffinityPropagation
import numpy as np

def apc(encodings, **kw):
	if kw['sof'] == 'sample':
		encodings = np.array(encodings)[1:]
	else:
		encodings = np.array(encodings).T[1:]
	#encodings = np.array(encodings)[1:]

	if len(encodings) < 3:
		return 0, 'sample number should be greater than 3.'

	data = encodings[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)

	e = ''
	try:
		af = AffinityPropagation().fit(data)
	except ValueError as e:
		print(e)
		return 0, e

	cluster_centers_indices = af.cluster_centers_indices_
	labels = af.labels_
	#n_clusters_ = len(cluster_centers_indices)

	res = []
	for i in range(len(data)):
		res.append([encodings[i][0],labels[i]])
	return res, 0