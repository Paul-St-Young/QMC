#!/usr/bin/env python
"""
GAMESS output analysis script
Author: Paul Young
Last Modified: July 9 2014
"""

import sys

class GAMESS:
	def __init__(self,gmsout):
		self.gmsout = gmsout
			
	def grab(self,entry):
	# prints out the first entry it finds in the output
		for line in open(self.gmsout,'r'):
			if line.find(entry)!=-1:
				print entry,line.split()[-1]
				
	def count_basis(self):
	# count # basis determinants, 
		in_basis = False
		MAX_N_BASIS = 1000000
		n_basis = 0
		for line in open(self.gmsout,'r'):
			if line.find("ALPHA")!=-1 and line.find("BETA")!=-1:
				in_basis = True
			if line.find("DENSITY MATRIX OVER ACTIVE MO")!=-1:
				in_basis = False
				break
			if in_basis and n_basis<MAX_N_BASIS and line.find("|")!=-1:
				n_basis+=1
		return n_basis/2-2


beh=GAMESS(sys.argv[1])
print "============================================="
print " Analyzing GAMESS output " + sys.argv[1]
beh.grab("TOTAL ENERGY =")
print "TOTAL NUMBER OF DET =",beh.count_basis()
