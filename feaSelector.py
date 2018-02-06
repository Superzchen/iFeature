#!/usr/bin/env python
#_*_coding:utf-8_*_

import argparse
from featureSelection import *
import sys, os, re
pPath = re.sub(r'scripts$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from clusters import readCode

if __name__ == '__main__':
	parser = argparse.ArgumentParser(usage="it's usage tip.",
									 description="feature selection")
	parser.add_argument("--file", required=True, help="input encoding file")
	parser.add_argument("--label", required=True, help="sample label")
	parser.add_argument("--type", required=True,
						choices=['CHI2', 'IG', 'MIC', 'pearsonr'], help="select cluster method")
	parser.add_argument("--out", help="output file")
	args = parser.parse_args()

	encodings = readCode.readCode(args.file)
	myFun = args.type + '.' + args.type + '(encodings, args.label)'
	print('Feature selection method: ' + args.type)
	selectedFeatures, e = eval(myFun)
	output = args.out if args.out != None else 'featureRank.txt'
	saveFeature.saveFeature(selectedFeatures, e, args.type, output)


