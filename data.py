#!/usr/bin/env python
# Author: Paul Young
# Last Modified: 08/08/2014
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
	2. Read # of printed determinants correctly
"""

import os
import subprocess

class MrData:
	def __init__(self,thisMolecule):
		# Tell Mr. Data the name of the molecule you would like to analyze
		# for example, if the GAMESS input is LiH.inp, then "molecule=LiH"
		self.known_molecules=[["H","HYDROGEN"],["He","HELIUM"],["Li","LITHIUM"],["Be","BERYLLIUM"],["B","BORON"],["C","CARBON"],["N","NITROGEN"],["O","OXYGEN"],["F","FLUORINE"],["Ne","NEON"]] # this list is ONLY used in setupGMS function
		self.best_NACT={'H':1,'He':5,'Li':14,'Be':14,'B':13,'C':10,'N':16,'O':14,'F':14,'Ne':9} # empirical, criteria: SOCI has ~1000 CSF with CI coefficient > 0.0001
		self.MULT={'H':2,'He':1,'Li':2,'Be':1,'B':2,'C':3,'N':4,'O':3,'F':2,'Ne':1} # Hund's rule
		self.molecule=thisMolecule
		self.basis_name=""
		self.basis_set=""
		
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
		self.basis_file="BSE_cc-pV5Z-H"
		self.z=-1
		
	def setupGMS(self,GMS_Template):
		f=open(self.basis_file,'r')
		
		# split the file into blocks
		block=[]
		newblock=""
		for line in f:
			newblock+=line
			if line.lstrip(" ")=="\n":
				block.append(newblock)
				newblock=""
		self.basis_name=block[0].split("\n")[0].split(" ")[2]
		
		# identify the beginning of basis data
		for i in range(len(block)):
			if block[i].split("\n")[0]=="$DATA":
				start_id=i
				break

		# extract basis data
		basis={}
		for i in range(start_id,len(block)):
			if i==start_id:
				basis[block[i].split("\n")[1]]="\n".join( block[i].split("\n")[2:-1] )
			else:
				basis[block[i].split("\n")[0]]="\n".join( block[i].split("\n")[1:-1] )
		
		# save relevant basis sets
		found=False
		z=0
		for mol in self.known_molecules:
			z+=1
			if self.molecule in mol:
				self.basis_set=basis[mol[-1]]
				found=True
				self.z=z
				break
		if not found:
			print "One of the constituents in this material does not match any known terrestrial element"
			exit(0)
		
		# write GMS input
		with open(self.inp,'w') as inp:
			inp.write(" ")
			for line in open(GMS_Template,'r'):
				line=line.replace("molecule",self.molecule)
				line=line.replace("basis_name",self.basis_name)
				line=line.replace("basis_set",self.basis_set)
				line=line.replace("zcharge",str(self.z)+".0")
				line=line.replace("NELS=3","NELS="+str(self.z))
				line=line.replace("NACT=12","NACT="+str(self.best_NACT[self.molecule]))
				line=line.replace("MULT=2","MULT="+str(self.MULT[self.molecule]))
				inp.write(line)
				
		subprocess.call(["mkdir",self.molecule])
		subprocess.call(["mv",self.inp,self.molecule])
	# end def setupGMS
				
	def countNORB(self,gamess_dat):
		in_orbit=False
		in_vec=False
		count_NORB=0
		cur=0
		prev=-1
		for line in open(gamess_dat,'r'):
			if line.find("OPTIMIZED MCSCF")!=-1:
				in_orbit=True
			if in_orbit and line.find("$VEC")!=-1:
				in_vec=True
			if in_orbit and in_vec:
				cur = line[:2]
				if cur!=prev:
					count_NORB+=1
				prev=cur
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
			return count_NORB
			
	# ======================= rungms ======================= #
	def rungms(self,i):
		# REM use ISPHER=1 at $CONTROL in inp
		gamess_err_file = open(self.gamess_err,'w') # log error from GAMESS runs
		if (i==1 or i==0):
			# run GAMESS for the first time
			# ====
			print("============================================")
			print("Executing 1st GAMESS calcuation with "+self.inp)
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
				NORB+" PRTMO=.T./g' <"+self.reinp+" >"+self.temp)
			os.system("rm "+self.reinp+"; mv "+self.temp+" "+self.reinp)
	
			# rerun GAMESS with MO as GUESS
			# ====
			print("Executing 2nd GAMESS caulation with "+self.reinp)
			print("Check progress in "+self.reout)
	
			with open(self.reout,'w') as outfile:
				subprocess.call(["rungms",self.reinp],stdout=outfile,stderr=gamess_err_file)
		# end if i==2
				
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
			#subprocess.call(["convert4qmc","-gamessAscii",self.reout,"-ci",self.reout,"-readInitialGuess",str(self.countNORB("gamess/"+self.molecule+".dat")),"-threshold","0.0001","-add3BodyJ"],stdout=outfile)
			subprocess.call(["convert4qmc","-gamessAscii",self.reout,"-ci",self.reout,"-threshold","0.0001","-add3BodyJ"],stdout=outfile)
	
		# folder management
		# ====
		os.system("mv sample.Gaussian-G2.xml "+self.wfs+"; mv sample.Gaussian-G2.ptcl.xml "+self.ptcl)
		os.system("mkdir convert; mv *.xml *.combo.dat "+self.reout+" "+self.convert_out+" convert")

	# ======================= cuspCorrection ======================= #

	# --------------- addOption2QMCline --------------- #
	def addOption2QMCline(self,file2edit,linehead,option):
		outfile=open(self.temp,'w')
		with open(file2edit,'r') as infile:
			for line in infile:
				if line.find(linehead)!=-1:
					outfile.write(line[:-2]+option)
				else:
					outfile.write(line)
		outfile.close()
		subprocess.call(["rm",file2edit])
		subprocess.call(["mv",self.temp,file2edit])
	
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

	def add_QMCblock_before(self,file2edit,linehead,block_file):
		tmp_file=open(self.temp,'w')
		with open(file2edit,'r') as infile:
			for line in infile:
				if line.find(linehead)!=-1:
					with open(block_file,'r') as block:
						tmp_file.write(block.read())
				tmp_file.write(line)
		tmp_file.close()
		subprocess.call(["rm",file2edit])
		subprocess.call(["mv",self.temp,file2edit])
		
	def fill_QMCblock(self,qmcblock):
		os.system("sed 's/molecule/"+self.molecule+"/g' "+self.template_main+" > "+qmcblock+".xml")
		os.system("sed -i \"s/today/$(date)/\" "+qmcblock+".xml")
		host_computer=subprocess.check_output(["uname","-a"]).split(" ")[1]
		os.system("sed -i \"s/host_computer/"+host_computer+"/\" "+qmcblock+".xml")
		self.fileline2file(qmcblock+".xml","qmcblock",template_dir+qmcblock+".xml")
	
	def cuspCorrection(self):
		# add cuspCorrection and remove jastrow in the wavefunction
		os.system("cp convert/"+self.wfs+" convert/"+self.ptcl+" .")
		self.addOption2QMCline(self.wfs,"<determinantset"," cuspCorrection=\"yes\">\n")
		self.mvQMCblock("jastrow","jastrow.xml")
		qmcblock = "cuspCorrection"
		
		self.fill_QMCblock(qmcblock)
		print("Executing Cusp Correction")
		print("Check progress in "+qmcblock+".out")
		os.system("export OMP_NUM_THREADS=1")
		with open(qmcblock+".out",'w') as outfile:
			subprocess.call(["qmcapp",qmcblock+".xml"],stdout=outfile)
			
		# modify the wavefunction
		self.addOption2QMCline(self.wfs,"name=\"spo-up\""," cuspInfo=\"spo-up.cuspInfo.xml\">\n")
		self.addOption2QMCline(self.wfs,"name=\"spo-dn\""," cuspInfo=\"spo-dn.cuspInfo.xml\">\n")
		self.add_QMCblock_before(self.wfs,"</wavefunction>","jastrow.xml")
		
		# folder management
		os.system("mkdir cusp;mv newOrbs* *.xml wftest.000 eloc.dat *.out cusp")
	
	def vmc(self,tag):
		vmc_dir="vmc"+tag
		os.system("mkdir "+vmc_dir)
		print("Setting up VMC run in "+vmc_dir)
		
		# copy necessary files
		#os.system("cp cusp/"+self.wfs+" cusp/"+self.ptcl+" cusp/*.cuspInfo.xml "+vmc_dir)
		os.system("cp cusp/"+self.wfs+" cusp/"+self.ptcl+" cusp/spo*"+" .")
		self.mvQMCblock("jastrow","jastrow.xml")
		os.system("mv "+self.wfs+" "+self.ptcl+" spo*"+" jastrow.xml "+vmc_dir)
		
		qmcblock="vmc"
		self.fill_QMCblock(qmcblock)
		os.system("mv "+qmcblock+".xml "+vmc_dir)
		print("============================================")
	
	# ======================= optJastrow ======================= #
	def optJastrow(self,tag):
		opt_dir="opt"+tag
		os.system("mkdir "+opt_dir)
		print("Setting up Jastrow Optimization in "+opt_dir)
		
		# copy necessary files
		os.system("cp cusp/"+self.wfs+" cusp/"+self.ptcl+" cusp/*.cuspInfo.xml "+opt_dir)
	
		qmcblock="optJastrow"
		self.fill_QMCblock(qmcblock)
		os.system("mv "+qmcblock+".xml "+opt_dir)
		print("============================================")
		
		
		#with open(qmcblock+".out",'w') as outfile:
			#subprocess.call(["qmcapp",qmcblock+".xml"],stdout=outfile)
	
	# ======================= rundmc ======================= #
	def rundmc(self,tag):
		dmc_dir="dmc"+tag
		os.system("mkdir "+dmc_dir)
		
		# write molecule-specified data to qmc template
		print("Setting up DMC run")
		qmcblock="dmc"
		self.fill_QMCblock(qmcblock)
		os.system("sed -i 's/"+self.wfs+"/"+self.opt_wfs+"/' "+qmcblock+".xml")
		os.system("mv "+qmcblock+".xml "+dmc_dir)
		
		# !!!!! Dangerous hard code !!!!!
		# by default use opt0/$molecule.s005.opt.xml
		os.system("cd opt0;cp spo-* "+self.molecule+".s005.opt.xml "+self.molecule+"_ptcl.xml ../"+dmc_dir)
		os.system("cd "+dmc_dir+";mv "+self.molecule+".s005.opt.xml "+self.molecule+"_opt_wfs.xml")

# ======================= main ======================= #
import argparse

def main():
	parser = argparse.ArgumentParser(description='LtCdr. Data at your service. I can analyze the properties of unknown materials if sufficient information regarding its constituents is available.')
	parser.add_argument("molecule", help="molecule to analyze. (name of the rungms input file should be 'molecule.inp')")
	parser.add_argument("-g", "--gms", type=int, choices=[0,1,2,3], help="run gamess calculation 1=CAS, 2=CAS_rerun, 3=SOCI, 0=all" )
	parser.add_argument("-c", "--convert", action="store_true", help="convert gamess output to qmc input" )
	parser.add_argument("-cusp", "--cuspCorrection", action="store_true", help="perform cusp correction" )
	parser.add_argument("-j", "--optJastrow", type=str, help="\'-j tag\' perform jastrow optimization in folder opt$tag" )
	parser.add_argument("-v", "--runVMC", type=str, help="\'-v tag\' perform VMC run in folder vmc$tag" )
	parser.add_argument("-d", "--runDMC", type=str, help="Setup Diffusion Monte Carlo" )
	parser.add_argument("-a", "--all", type=str, help="run everything with a given GAMESS input Template. For example, 'data.py -a GMS_Template.inp H' will run everything for Hydrogen atom" )
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
		Data.optJastrow(args.optJastrow)
	if args.runVMC:
		Data.vmc(args.runVMC)
	if args.runDMC:
		Data.rundmc(args.runDMC)
	if args.all:
		Data.setupGMS(args.all)
		#os.chdir(Data.molecule)
		#Data.rungms(0)
main()


