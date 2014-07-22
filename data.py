#!/usr/bin/env python
# Author: Paul Young
# Last Modified: 07/12/2014
# Objectified autorun script for GAMESS->QMCPACK calculation on the molecule specified

template_dir = "QMC_Templates/" # remember the trailing /

"""
Prerequisites:
====
	1. GAMESS is compiled and "rungms" script location is added to PATH
	2. QMCPACK is compiled and binary location is added to PATH
	3. rungms outputs the .dat file to cwd
	4. $(molecule).inp is present in cwd
	5. Directory containing QMC templates is specified
	
To do:
----
	1. change all os.system to subprocess.call
"""

import os
import subprocess

class MrData:
	def __init__(self,thisMolecule):
		# Tell Mr. Data the name of the molecule you would like to analyze
		# for example, if the GAMESS input is LiH.inp, then "molecule=LiH"
		self.molecule=thisMolecule
		
		# some file names
		self.inp = self.molecule + ".inp"
		self.out = self.molecule + ".out"
		self.reinp = self.molecule+"_rerun.inp"
		self.reout = self.molecule+"_rerun.out"
		self.temp = "gamess2qmcpack.tmp"
		self.gamess_err = "gamess.err"
		self.gamess_dat = self.molecule+".dat "+self.molecule+"_rerun.dat"

		self.wfs = self.molecule+"_wfs.xml"
		self.opt_wfs=self.molecule+"_opt_wfs.xml"
		self.ptcl = self.molecule+"_ptcl.xml"
		self.convert_out = self.molecule+"_convert.out"

		self.template_main = template_dir+"main.xml"
	
	# ======================= rungms ======================= #
	def rungms(self,i):
		# REM use ISPHER=1 at $CONTROL in inp
		gamess_err_file = open(self.gamess_err,'w') # log error from GAMESS runs
		if (i==1 or i==0):
			# run GAMESS for the first time
			# ====
			print("============================================")
			print("Executing 1st GAMESS caulation with "+self.inp)
			print("Check progress in "+self.out)
			with open(self.out,'w') as outfile:
				subprocess.call(["rungms",self.inp],stdout=outfile,stderr=gamess_err_file)
		if (i==2 or i==0):
			# creat reinp with MOREAD as GUESS
			# ====
				# find NORB
				# ----
			print("Looking for OPTIMIZED MCSCF orbits")
			gamess_out=open(self.molecule+".dat",'r') # make sure GAMESS outputs to cwd
			os.system("cp "+self.inp +" "+ self.reinp)
			reinpfile=open(self.reinp,'a')
			in_orbit=False
			in_vec=False
			count_NORB=0
			cur=0
			prev=-1
			for line in gamess_out:
				if line.find("OPTIMIZED MCSCF")!=-1:
					in_orbit=True
				if in_orbit and line.find("$VEC")!=-1:
					in_vec=True
				if in_orbit and in_vec:
					cur = line[:2]
					if cur!=prev:
						count_NORB+=1
					prev=cur
					reinpfile.write(line)
					if line.find("$END")==-1:
						NORB=line[:2]
				if line.find("$END")!=-1:
					in_vec=False
					in_orbit=False		
			count_NORB-=2 # first transition is cur=$ prev=-1, last transition is cur=$,prev=last
			if (int(NORB)!=int(count_NORB)%100):
                        	print "Mismatch! NORB="+str(NORB)+" <> count_NORB="+str(count_NORB)
				exit(0)
			else:
				NORB=str(count_NORB)
			reinpfile.close()
				# creat reinp
				# ----
			print("Found "+NORB+" orbits, augmenting GUESS to 'GUESS=MOREAD NORB="+NORB+"' in "+self.reinp)
			os.system("sed 's/GUESS=HCORE/GUESS=MOREAD NORB="+
				NORB+"/g' <"+self.reinp+" >"+self.temp)
			os.system("rm "+self.reinp+"; mv "+self.temp+" "+self.reinp)
	
			# rerun GAMESS with MO as GUESS
			# ====
			print("============================================")
			print("Executing 2nd GAMESS caulation with "+self.reinp)
			print("Check progress in "+self.reout)
	
			with open(self.reout,'w') as outfile:
				subprocess.call(["rungms",self.reinp],stdout=outfile,stderr=gamess_err_file)
				
		gamess_err_file.close()

		print("GAMESS calculations completed")
		print("=========================================")
		
		# blindly reorganize folder for now, should check for completion first.
		os.system("mkdir gamess;mv *.inp *.out *.dat "+self.gamess_err+" gamess")
		os.system("cp gamess/"+self.reout+" .")
	# ======================= convert4qmc ======================= #
	def convert4qmc(self):
		# Convert GAMESS output
		print("Converting GAMESS output to QMCPACK input")
		with open(self.convert_out,"w") as outfile:
			subprocess.call(["convert4qmc","-gamessAscii",self.reout,"-ci",self.reout,"-threshold","0.0001","-add3BodyJ"],stdout=outfile)
	
		# folder management
		# ====
		os.system("mv sample.Gaussian-G2.xml "+self.wfs+"; mv sample.Gaussian-G2.ptcl.xml "+self.ptcl)
		os.system("mkdir convert; mv *.xml *.combo.dat "+self.reout+" "+self.convert_out+" convert")

	# ======================= cuspCorrection ======================= #

	# --------------- addOption2QMCline --------------- #
	def addOption2QMCline(self,linehead,option):
		outfile=open(self.temp,'w')
		with open(self.wfs,'r') as infile:
			for line in infile:
				if line.find(linehead)!=-1:
					outfile.write(line[:-2]+option)
				else:
					outfile.write(line)
		outfile.close()
		subprocess.call(["rm",self.wfs])
		subprocess.call(["mv",self.temp,self.wfs])
	
	# --------------- mvQMCblock --------------- #
	def mvQMCblock(self,blockname,filename):
		block_file = open(filename,'w')
		in_block=False
		outfile=open(self.temp,'w')
		with open(self.wfs,'r') as infile:
			for line in infile:
				if in_block or line.find(blockname)!=-1:
					block_file.write(line)
				else:
					outfile.write(line)
				if line.find(blockname)!=-1:
					in_block=not in_block
		block_file.close()
		outfile.close()
		subprocess.call(["rm",self.wfs])
		subprocess.call(["mv",self.temp,self.wfs])
	
	# --------------- fileline2file --------------- #
	def fileline2file(self,templatefile,linehead,linefile):
		tmp_file=open(self.temp,'w')
		with open(templatefile,'r') as infile:
			for line in infile:
				if line.find(linehead)!=-1:
					with open(linefile,'r') as replacement:
						tmp_file.write(replacement.read())
				else:
					tmp_file.write(line)
		tmp_file.close()
		subprocess.call(["rm",templatefile])
		subprocess.call(["mv",self.temp,templatefile])

	def add_QMCblock_before(self,linehead,block_file):
		tmp_file=open(self.temp,'w')
		with open(self.wfs,'r') as infile:
			for line in infile:
				if line.find(linehead)!=-1:
					with open(block_file,'r') as block:
						tmp_file.write(block.read())
				tmp_file.write(line)
		tmp_file.close()
		subprocess.call(["rm",self.wfs])
		subprocess.call(["mv",self.temp,self.wfs])
	
	def cuspCorrection(self):
		# add cuspCorrection and remove jastrow in the wavefunction
		os.system("cp convert/"+self.wfs+" convert/"+self.ptcl+" .")
		self.addOption2QMCline("<determinantset"," cuspCorrection=\"yes\">\n")
		self.mvQMCblock("jastrow","jastrow.xml")
		qmcblock = "cuspCorrection"
		os.system("sed 's/molecule/"+self.molecule+"/g' "+self.template_main+" > "+qmcblock+".xml")
		self.fileline2file(qmcblock+".xml","qmcblock",template_dir+qmcblock+".xml")
	
		print("============================================")
		print("Executing Cusp Correction")
		print("Check progress in "+qmcblock+".out")
		os.system("export OMP_NUM_THREADS=1")
		with open(qmcblock+".out",'w') as outfile:
			subprocess.call(["qmcapp",qmcblock+".xml"],stdout=outfile)
		
		# folder management
		os.system("mkdir cusp;mv newOrbs* *.xml wftest.000 eloc.dat *.out cusp")
	
	# ======================= optJastrow ======================= #
	def optJastrow(self):
		# modify the wavefunction
		os.system("cp cusp/"+self.wfs+" cusp/"+self.ptcl+" cusp/*.cuspInfo.xml cusp/jastrow.xml .")
		self.addOption2QMCline("updet"," cuspInfo=\"spo-up.cuspInfo.xml\">\n")
		self.addOption2QMCline("downdet"," cuspInfo=\"spo-dn.cuspInfo.xml\">\n")
		self.add_QMCblock_before("</wavefunction>","jastrow.xml")
	
		qmcblock="optJastrow"
		os.system("sed 's/molecule/"+self.molecule+"/g' "+self.template_main+" > "+qmcblock+".xml")
		self.fileline2file(qmcblock+".xml","qmcblock",template_dir+qmcblock+".xml")
	
		print("============================================")
		print("Optimizing Jastrows")
		print("Check progress in "+qmcblock+".out")
		
		
		#with open(qmcblock+".out",'w') as outfile:
			#subprocess.call(["qmcapp",qmcblock+".xml"],stdout=outfile)
			
		# folder management
	
	
	# ======================= rundmc ======================= #
	def rundmc(self):
		# write molecule-specified data to qmc template
		qmcblock="dmc"
		os.system("sed 's/molecule/"+self.molecule+"/g' "+self.template_main+" > "+qmcblock+".xml")
		self.fileline2file(qmcblock+".xml","qmcblock",template_dir+qmcblock+".xml")
		os.system("sed -i 's/"+self.wfs+"/"+self.opt_wfs+"/' "+qmcblock+".xml")
	
		# run dmc
		#print("============================================")
		#print("Running DMC")
		#print("Check progress in "+qmcblock+".out")
		#with open(qmcblock+".out",'w') as outfile:
		#	subprocess.call(["qmcapp",qmcblock+".xml"],stdout=outfile)
		
		# folder management
		#new_folder="qmc/"+qmcblock
		#os.system("mkdir -p "+new_folder+";mv newOrbs* J1* J2* "+new_folder+";mv "+qmcblock+".xml "+new_folder+"; mv "+qmcblock+".out "+new_folder)
	

# ======================= main ======================= #
import argparse

def main():
	parser = argparse.ArgumentParser(description='LtCdr. Data at your service. I can analyze the properties of unknown materials if sufficient information regarding its constituents is available.')
	parser.add_argument("molecule", help="molecule to analyze. (name of the rungms input file should be 'molecule.inp')")
	parser.add_argument("-g", "--gms", type=int, choices=[0,1,2], help="run gamess calculation (you may specify if it's the 1st or 2nd calculation, if 0 is given I'll run both)" )
	parser.add_argument("-c", "--convert", action="store_true", help="convert gamess output to qmc input" )
	parser.add_argument("-cusp", "--cuspCorrection", action="store_true", help="perform cusp correction" )
	parser.add_argument("-j", "--optJastrow", action="store_true", help="perform jastrow optimization" )
	parser.add_argument("-d", "--runDMC", action="store_true", help="run Diffusion Monte Carlo" )
	args = parser.parse_args()
	args = parser.parse_args()
	
	Data=MrData(args.molecule)
	
	if args.gms==1:
		Data.rungms(1)
	elif args.gms==2:
		Data.rungms(2)
	elif args.gms==0:
		Data.rungms(0)
	else:
		pass
	
	if args.convert:
		Data.convert4qmc()
	if args.cuspCorrection:
		Data.cuspCorrection()
	if args.optJastrow:
		Data.optJastrow()
	if args.runDMC:
		Data.rundmc()
main()


