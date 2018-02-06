#!/usr/bin/env python
#_*_coding:utf-8_*_

def saveFeature(feature, e, method, file='featureRank.txt'):
	with open(file, 'w') as f:
		if feature == 0:
			f.write(str(e))
		else:
			f.write('# Feature selection method: %s\n' %method)
			f.write('# The features were ranked according to their importance, the topper the more important the feature is\n');
			f.write('=======================\n')
			for i in feature:
				f.write(i[0] + '\t' + str(i[1]) + '\n')
	return None