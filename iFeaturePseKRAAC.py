#!/usr/bin/env python
#_*_coding:utf-8_*_

import argparse, sys
from PseKRAAC import *
from codes import readFasta, saveCode

USAGE = """The 'raactype' value for each subtype descriptor could be chosen from:
	type1    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type2	 [2, 3, 4, 5, 6,    8,                        15,                 20]
	type3A   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type3B   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type4	 [         5,       8, 9,     11,     13,                         20]
	type5    [   3, 4,          8,    10,                 15,                 20]
	type6A   [      4, 5,                                                     20]
	type6B   [         5,                                                       ]
	type6C   [         5,                                                       ]
	type7    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type8    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type9    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type10   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type11   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type12   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,     20]
	type13   [      4,                        12,                 17,         20]
	type14   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
	type15   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20]
	type16   [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20]	
"""

USAGEHASH = {
	'type1' :[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type2' :[2, 3, 4, 5, 6,    8,                        15,                 20],
	'type3A':[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type3B':[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type4'	:[         5,       8, 9,     11,     13,                         20],
	'type5' :[   3, 4,          8,    10,                 15,                 20],
	'type6A':[      4, 5,                                                     20],
	'type6B':[         5,                                                       ],
	'type6C':[         5,                                                       ],
	'type7' :[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type8' :[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type9' :[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type10':[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type11':[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type12':[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,     20],
	'type13':[      4,                        12,                 17,         20],
	'type14':[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
	'type15':[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20],
	'type16':[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,             20],
}



if __name__ == '__main__':
	parser = argparse.ArgumentParser(usage="it's usage tip.",
									 description="Generating PseKRAAC descriptors:")
	parser.add_argument("--file", help="input fasta file")
	parser.add_argument("--type",
						choices=['type1', 'type2', 'type3A', 'type3B', 'type4', 'type5', 'type6A', 'type6B', 'type6C',
								 'type7', 'type8', 'type9', 'type10', 'type11', 'type12', 'type13', 'type14', 'type15',
								 'type16'], help="the descriptor type")
	parser.add_argument("--subtype", choices=['g-gap', 'lambda-correlation'], help="the subtype of the descriptor type, default is 'g-gap'")
	parser.add_argument("--ktuple", choices=[1, 2, 3], help="k-tuple peptide, default is 2", type=int)
	parser.add_argument("--gap_lambda", choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], help="the gap value or lambda value for the 'g-gap' model or 'lambda-correlation' model", type=int)
	parser.add_argument("--raactype", help="the reduced amino acids cluster type", type=int)
	parser.add_argument("--show", help="show detatiled available '--raactype' value for each type",
						action="store_true")
	parser.add_argument("--out", dest='outFile', help="the generated descriptor file")
	args = parser.parse_args()

	if args.show:
		print(USAGE)
		sys.exit(1)

	if args.file == None:
		print('The following arguments are required: --file\n')
		sys.exit(1)
	if args.type == None:
		print('The following arguments are required: --type\n')
		sys.exit(1)
	if args.raactype == None:
		print('The following arguments are required: --raactype\n')
		sys.exit(1)
	if args.gap_lambda == None:
		print('The following arguments are required: --gap_lambda\n')
		sys.exit(1)

	subtype = args.subtype if args.subtype != None else 'g-gap'
	ktuple = args.ktuple if args.ktuple != None else 2


	if args.raactype not in USAGEHASH[args.type]:
		print('The "raactype" value error. For detailed parameter, please see:\n\n')
		print(USAGE)
		sys.exit(1)


	fastas = readFasta.readFasta(args.file)
	myFun = args.type + '.type1' + '(fastas, subtype, args.raactype, ktuple, args.gap_lambda)'
	print(myFun)
	print('Descriptor type: ' + args.type)
	print('Subtype model: ' + subtype)
	print('reduced amino acids cluster type: ' + str(args.raactype))
	print('k_tuple peptide number: ' + str(ktuple))
	print('gap or lambda value: ' + str(args.gap_lambda) + '\n\n')
	encodings = eval(myFun)
	outFile = args.outFile if args.outFile != None else 'encoding.tsv'
	saveCode.savetsv(encodings, outFile)