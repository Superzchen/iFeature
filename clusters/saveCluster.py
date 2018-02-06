#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
import pandas as pd

def saveCluster(cluster, e, file='Clusters.txt'):
	with open(file, 'w') as f:
		if cluster == 0:
			f.write(str(e))
		else:
			myCluster = np.array(cluster)
			df = pd.DataFrame({'name': myCluster[:, 0], 'cluster': myCluster[:, 1]})
			mySet = set(np.array(df.cluster).tolist())
			f.write('# The sample/feature can be clustered into %d clusters:\n' %len(mySet))
			for l in sorted(mySet):
				newData = np.array(df.loc[df.loc[:, "cluster"] == l, :].name)
				f.write('Cluster_%s:\t' % l)
				for i in newData:
					f.write(i + '\t')
				f.write('\n')
			f.write('\n==============================================================\n')

			f.write('Protein/Feature\tcluster\n')
			for i in cluster:
				f.write(i[0] + '\t' + str(i[1]) + '\n')
	return None