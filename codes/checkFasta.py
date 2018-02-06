#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def checkFasta(fastas):
	status = True
	lenList = set()
	for i in fastas:
		lenList.add(len(i[1]))
	if len(lenList) == 1:
		return True
	else:
		return False

def minSequenceLength(fastas):
	minLen = 10000
	for i in fastas:
		if minLen > len(i[1]):
			minLen = len(i[1])
	return minLen

def minSequenceLengthWithNormalAA(fastas):
	minLen = 10000
	for i in fastas:
		if minLen > len(re.sub('-', '', i[1])):
			minLen = len(re.sub('-', '', i[1]))
	return minLen