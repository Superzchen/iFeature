#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
import shutil
pPath = re.sub(r'scripts$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from codes import readFasta
import argparse


def generateSecondaryStructure(fastas, outDir, psipred):
	"""
	Generate predicted protein secondary structure by using the 'PSIPRED Version 4.01' program.

	Parameters
	----------
	file : file
		the file, which include the protein sequences in fasta format.

	psipred: string
		the path of PSIPRED Version 4.01 program, the PSIPRED Version 4.01 program could be download at: http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/.

	Returns
	-------
	a string: 
		A directory name, which include the predicted protein secondary structure information.
	"""

	if os.path.exists(outDir) == False:
		os.mkdir(outDir)

	for i in fastas:
		name, sequence = re.sub('\|', '', i[0]), i[1]
		with open(name + '.txt', 'w') as f:
			f.write('>'+name+'\n'+sequence + '\n')
		myCmd = psipred + ' ' + name + '.txt'
		if os.path.exists(outDir + '/' + name + '.ss2') == False:
			os.system(myCmd)
			os.remove(name + '.txt')
			os.remove(name + '.ss')
			os.remove(name + '.horiz')
			shutil.move(name+'.ss2', outDir)
	return outDir

if __name__ == '__main__':
	parser = argparse.ArgumentParser(usage="it's usage tip.", description="generate protein secondary structure profile")
	parser.add_argument("--file", required=True, help="protein sequence file in fasta format")
	parser.add_argument("--psipred", help="the path of psipred program")
	args = parser.parse_args()

	psipred = args.psipred if args.psipred != None else 'runpsipred'
	fastas = readFasta.readFasta(args.file)
	outputDir = generateSecondaryStructure(fastas, 'out', psipred)
	print('The predicted secodnary structure are stored in: ' + outputDir)