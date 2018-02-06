#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import readFasta
import saveCode

USAGE = """
USAGE:
	python CTDCClass.py input.fasta output amino_acids_group_1 amino_acids_group_2 ... amino_acids_group_N
	
	input.fasta:                 the input protein sequence file in fasta format.	
	output:                      the encoding file.
	amino_acids_group_x          the amino acids groups.

EXAMPLE:
	python CTDCClass.py example/test-protein.txt CTDCClass.tsv RKEDQN GASTPHY CLVIMFW
"""

def Count(seq1, seq2):
	sum = 0
	for aa in seq1:
		sum = sum + seq2.count(aa)
	return sum

def CTDCClass(fastas, groups):
	encodings = []
	header = ['#']
	for g in range(len(groups)):
		header.append('g.'+str(g+1))
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		for group in groups:
			code.append(Count(group, sequence) / len(sequence))
		encodings.append(code)
	return encodings


if __name__ == '__main__':
	if len(sys.argv) < 5:
		print(USAGE)
		sys.exit(1)

	groups = sys.argv[3:]
	myStr = ''.join(groups)
	myStr = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', myStr)
	if len(myStr) != 20 or len(set(myStr)) != 20:
		print('\nERROR: The amino acid must be no-repeat in each groups and the sum is 20!\n')
	fastas = readFasta.readFasta(sys.argv[1])
	encodings = CTDCClass(fastas, groups)
	saveCode.savetsv(encodings, sys.argv[2])