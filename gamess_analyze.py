#!/usr/bin/env python
"""
GAMESS output analysis script
Author: Paul Young
Last Modified: July 28 2014
"""

ENTRIES=["TOTAL ENERGY =", "FINAL MCSCF ENERGY IS","COMPUTING THE HAMILTONIAN FOR"]

import sys

class GAMESS:
	def __init__(self,gmsout):
		self.gmsout = gmsout
			
	def grab(self,entry,line):
	# prints out the first entry it finds in the output
		if line.find(entry)!=-1:
			#print entry,line.split()[-1]
			print line.lstrip()[:-1]
				
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
		
	def count_CSF(self):
	# count # CSF 
		in_basis = False
		n_basis = 0
		for line in open(self.gmsout,'r'):
			if line.find("CSF")!=-1 and line.find("OCCUPANCY")!=-1:
				in_basis = True
			if in_basis:
				if len(line.lstrip(" ").split(" "))!=1:
					n_basis+=1
			if line.find("END OF CI-MATRIX DIAGONALIZATION")!=-1:
				in_basis=False
				break
		n_basis-=3
		return n_basis
		
	def main_loop(self):
		prt=0
		for line in open(self.gmsout,'r'):
			for entry in ENTRIES:
				self.grab(entry,line)
			# check for errors
			if line.find("WARNING")!=-1:
				prt=3
				if line.find("NORMALIZATION")!=-1:
					prt=0
				if line.find("MCSCF POPULATIONS")!=-1:
					prt=0
			if prt>0:
				print line.rstrip("\n")
				prt-=1

beh=GAMESS(sys.argv[1])
print "============================================="
print " Analyzing GAMESS output " + sys.argv[1]
beh.main_loop()
print "TOTAL NUMBER OF DET =",beh.count_basis()
print "TOTAL NUMBER OF CSF =",beh.count_CSF()

# check errors in output
	# CI eigenvector not found?
	# Occupation of orbital too small?
