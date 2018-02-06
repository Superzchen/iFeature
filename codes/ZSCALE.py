#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import checkFasta

def ZSCALE(fastas, **kw):
	if checkFasta.checkFasta(fastas) == False:
		print('Error: for "ZSCALE" encoding, the input fasta sequences should be with equal length. \n\n')
		return 0

	zscale = {
		'A': [0.24,  -2.32,  0.60, -0.14,  1.30], # A
		'C': [0.84,  -1.67,  3.71,  0.18, -2.65], # C
		'D': [3.98,   0.93,  1.93, -2.46,  0.75], # D
		'E': [3.11,   0.26, -0.11, -0.34, -0.25], # E
		'F': [-4.22,  1.94,  1.06,  0.54, -0.62], # F
		'G': [2.05,  -4.06,  0.36, -0.82, -0.38], # G
		'H': [2.47,   1.95,  0.26,  3.90,  0.09], # H
		'I': [-3.89, -1.73, -1.71, -0.84,  0.26], # I
		'K': [2.29,   0.89, -2.49,  1.49,  0.31], # K
		'L': [-4.28, -1.30, -1.49, -0.72,  0.84], # L
		'M': [-2.85, -0.22,  0.47,  1.94, -0.98], # M
		'N': [3.05,   1.62,  1.04, -1.15,  1.61], # N
		'P': [-1.66,  0.27,  1.84,  0.70,  2.00], # P
		'Q': [1.75,   0.50, -1.44, -1.34,  0.66], # Q
		'R': [3.52,   2.50, -3.50,  1.99, -0.17], # R
		'S': [2.39,  -1.07,  1.15, -1.39,  0.67], # S
		'T': [0.75,  -2.18, -1.12, -1.46, -0.40], # T
		'V': [-2.59, -2.64, -1.54, -0.85, -0.02], # V
		'W': [-4.36,  3.94,  0.59,  3.44, -1.59], # W
		'Y': [-2.54,  2.44,  0.43,  0.04, -1.47], # Y
		'-': [0.00,   0.00,  0.00,  0.00,  0.00], # -
	}
	encodings = []
	header = ['#']
	for p in range(1, len(fastas[0][1])+1):
		for z in ('1', '2', '3', '4', '5'):
			header.append('Pos'+str(p) + '.ZSCALE' + z)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], i[1]
		code = [name]
		for aa in sequence:
			code = code + zscale[aa]
		encodings.append(code)
	return encodings