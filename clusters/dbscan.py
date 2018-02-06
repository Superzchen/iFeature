#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

def dbscan(encodings, **kw):
	if kw['sof'] == 'sample':
		encodings = np.array(encodings)[1:]
	else:
		encodings = np.array(encodings).T[1:]

	if len(encodings) < 5:
		return 0, 'sample number should be greater than 5.'

	data = encodings[:, 1:]
	shape = data.shape
	data = np.reshape(data, shape[0] * shape[1])
	data = np.reshape([float(i) for i in data], shape)

	e = ''
	dataNew = StandardScaler().fit_transform(data)
	db = DBSCAN().fit(dataNew)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_

	res = []
	for i in range(len(data)):
		res.append([encodings[i][0], labels[i]])
	return res, e