#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
from sklearn.cluster import MeanShift, estimate_bandwidth

def meanshift(encodings, **kw):
	if kw['sof'] == 'sample':
		encodings = np.array(encodings)[1:]
	else:
		encodings = np.array(encodings).T[1:]

	if len(encodings) < 3:
		return 0, 'sample number should be greater than 3.'

	data = encodings[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)

	e = ''
	try:
		bandwidth = estimate_bandwidth(data)
		ms = MeanShift(bandwidth=bandwidth, bin_seeding=True).fit(data)
	except ValueError as e:
		print(e)
		return 0, e

	labels = ms.labels_
	cluster_centers = ms.cluster_centers_

	labels_unique = np.unique(labels)
	n_clusters_ = len(labels_unique)

	res = []
	for i in range(len(data)):
		res.append([encodings[i][0], labels[i]])
	return res, e