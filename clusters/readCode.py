#!/usr/bin/env python
#_*_coding:utf-8_*_

import os, sys
import numpy as np

def readCode(file):
	encodings = []
	if not os.path.exists(file):
		print('Error: file does not exist.')
		sys.exit(1)
	with open(file) as f:
		records = f.readlines()
	for i in records:
		array = i.rstrip().split() if i.strip() != '' else None
		encodings.append(array)
	return np.array(encodings)