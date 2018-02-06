#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
pPath = re.sub(r'scripts$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from codes import readFasta
import argparse


def generateDisorder(fastas, outDir, vsl2):
	"""
	Generate predicted protein disorder information by bysing VSL2 software.

	Parameters
	----------
	file : file
		the file, which include the protein sequences in fasta format.
	
	vsl2: the VSL2 software path
		VSL2 software path, which could be download at: http://www.dabi.temple.edu/disprot/download/VSL2.tar.gz

	Returns
	-------
	a string: 
		A directory name, which include the predicted protein disorder information.
	"""

	if os.path.exists(outDir) == False:
		os.mkdir(outDir)

	for i in fastas:
		name, sequence = re.sub('\|', '', i[0]), i[1]
		with open(name + '.txt', 'w') as f:
			f.write(sequence + '\n')
		myCmd = 'java -jar ' + vsl2  + ' -s:' + name + '.txt >' + outDir + '/' + name + '.dis'
		if os.path.exists(outDir + '/' + name + '.dis') == False:
			os.system(myCmd)
	return outDir

if __name__ == '__main__':
	parser = argparse.ArgumentParser(usage="it's usage tip.", description="generate protein secondary structure profile")
	parser.add_argument("--file", required=True, help="protein sequence file in fasta format")
	parser.add_argument("--vsl2", help="the path of vsl2 program")
	args = parser.parse_args()

	vsl2 = args.vsl2 if args.vsl2 != None else 'VSL2.jar'
	fastas = readFasta.readFasta(args.file)
	outputDir = generateDisorder(fastas, 'out', vsl2)
	print('The predicted disorder are stored in: ' + outputDir)