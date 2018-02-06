#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def CTDT(fastas, **kw):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}

	groups = [group1, group2, group3]
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	encodings = []
	header = ['#']
	for p in property:
		for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
			header.append(p + '.' + tr)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
		for p in property:
			c1221, c1331, c2332 = 0, 0, 0
			for pair in aaPair:
				if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
					c1221 = c1221 + 1
					continue
				if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
					c1331 = c1331 + 1
					continue
				if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
					c2332 = c2332 + 1
			code = code + [c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)]
		encodings.append(code)
	return encodings