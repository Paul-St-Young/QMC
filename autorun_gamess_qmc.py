#!/usr/bin/env python
# Author: Paul Young
# Last Modified: 07/02/2014
# Autorun script for GAMESS->QMCPACK calculation on the molecule specified

import sys
molecule = sys.argv[1] # just the name of the molecule
# so the GAMESS input should read "molecule+'.inp'"
template_dir = "/home/yyang173/Templates/QMC_Templates/" # remember the trailing /

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

# some file names
inp = molecule + ".inp"
out = molecule + ".out"
reinp = molecule+"_rerun.inp"
reout = molecule+"_rerun.out"
temp = "gamess2qmcpack.tmp"
gamess_err = "gamess.err"
gamess_dat = molecule+".dat "+molecule+"_rerun.dat"


wfs = molecule+"_wfs.xml"
opt_wfs=molecule+"_opt_wfs.xml"
ptcl = molecule+"_ptcl.xml"
convert_out = molecule+"_convert.out"

template_main = template_dir+"main.xml"

import os
import subprocess
# ======================= rungms ======================= #
def rungms():
	# REM use ISPHER=1 at $CONTROL in inp
	gamess_err_file = open(gamess_err,'w') # log error from GAMESS runs
	
	# run GAMESS
	# ====
	print("============================================")
	print("Executing 1st GAMESS caulation with "+inp)
	print("Check progress in "+out)
	with open(out,'w') as outfile:
		subprocess.call(["rungms",inp],stdout=outfile,stderr=gamess_err_file)

	# creat reinp with MOREAD as GUESS
	# ====
		# find NORB
		# ----
	print("Looking for OPTIMIZED MCSCF orbits")
	gamess_out=open(molecule+".dat",'r') # make sure GAMESS outputs to cwd
	os.system("cp "+inp +" "+ reinp)
	reinpfile=open(reinp,'a')
	in_orbit=False
	in_vec=False
	for line in gamess_out:
		if line.find("OPTIMIZED MCSCF")!=-1:
			in_orbit=True
		if in_orbit and line.find("$VEC")!=-1:
			in_vec=True
		if in_orbit and in_vec:
			reinpfile.write(line)
			if line.find("$END")==-1:
				NORB=line[:2]
		if line.find("$END")!=-1:
			in_vec=False
			in_orbit=False		
	reinpfile.close()
		# creat reinp
		# ----
	print("Found "+NORB+" orbits, augmenting GUESS to 'GUESS=MOREAD NORB="+NORB+"' in "+reinp)
	os.system("sed 's/GUESS=HCORE/GUESS=MOREAD NORB="+
		NORB+"/g' <"+reinp+" >"+temp)
	os.system("rm "+reinp+"; mv "+temp+" "+reinp)
	
	# rerun GAMESS with MO as GUESS
	# ====
	print("============================================")
	print("Executing 2nd GAMESS caulation with "+reinp)
	print("Check progress in "+reout)
	
	with open(reout,'w') as outfile:
		subprocess.call(["rungms",reinp],stdout=outfile,stderr=gamess_err_file)
		
	gamess_err_file.close()
	
	print("GAMESS calculations completed")
	print("=========================================")

# ======================= convert4qmc ======================= #
def convert4qmc():

	# Convert GAMESS output
	print("Converting GAMESS output to QMCPACK input")
	with open(convert_out,"w") as outfile:
		subprocess.call(["convert4qmc","-gamessAscii",reout,"-ci",reout,"-threshold","0.0001"],stdout=outfile)
	
	# folder management
	# ====
	os.system("mkdir gamess;mv *.inp *.out "+gamess_err+" *.dat gamess; mv sample.Gaussian-G2.xml "+wfs+"; mv sample.Gaussian-G2.ptcl.xml "+ptcl)

# ======================= cuspCorrection ======================= #

	# --------------- addOption2QMCline --------------- #
def addOption2QMCline(linehead,option):
	outfile=open(temp,'w')
	with open(wfs,'r') as infile:
		for line in infile:
			if line.find(linehead)!=-1:
				outfile.write(line[:-2]+option)
			else:
				outfile.write(line)
	outfile.close()
	subprocess.call(["rm",wfs])
	subprocess.call(["mv",temp,wfs])
	
	# --------------- mvQMCblock --------------- #
def mvQMCblock(blockname,filename):
	block_file = open(filename,'w')
	in_block=False
	outfile=open(temp,'w')
	with open(wfs,'r') as infile:
		for line in infile:
			if in_block or line.find(blockname)!=-1:
				block_file.write(line)
			else:
				outfile.write(line)
			if line.find(blockname)!=-1:
				in_block=not in_block
	block_file.close()
	outfile.close()
	subprocess.call(["rm",wfs])
	subprocess.call(["mv",temp,wfs])
	
	# --------------- fileline2file --------------- #
def fileline2file(templatefile,linehead,linefile):
	tmp_file=open(temp,'w')
	with open(templatefile,'r') as infile:
		for line in infile:
			if line.find(linehead)!=-1:
				with open(linefile,'r') as replacement:
					tmp_file.write(replacement.read())
			else:
				tmp_file.write(line)
	tmp_file.close()
	subprocess.call(["rm",templatefile])
	subprocess.call(["mv",temp,templatefile])

def add_QMCblock_before(linehead,block_file):
	tmp_file=open(temp,'w')
	with open(wfs,'r') as infile:
		for line in infile:
			if line.find(linehead)!=-1:
				with open(block_file,'r') as block:
					tmp_file.write(block.read())
			tmp_file.write(line)
	tmp_file.close()
	subprocess.call(["rm",wfs])
	subprocess.call(["mv",temp,wfs])
	
def cuspCorrection():
	# add cuspCorrection and remove jastrow in the wavefunction
	addOption2QMCline("<determinantset"," cuspCorrection=\"yes\">\n")
	mvQMCblock("jastrow","jastrow.xml")
	qmcblock = "cuspCorrection"
	os.system("sed 's/molecule/"+molecule+"/g' "+template_main+" > "+qmcblock+".xml")
	fileline2file(qmcblock+".xml","qmcblock",template_dir+qmcblock+".xml")
	
	print("============================================")
	print("Executing Cusp Correction")
	print("Check progress in "+qmcblock+".out")
	os.system("export OMP_NUM_THREADS=1")
	with open(qmcblock+".out",'w') as outfile:
		subprocess.call(["qmcapp",qmcblock+".xml"],stdout=outfile)
		
	# folder management
	os.system("mkdir -p qmc/cusp;mv newOrbs* *.qmc.xml wftest.000 eloc.dat qmc/cusp;mv *.cont.xml "+qmcblock+".xml qmc/cusp; mv "+qmcblock+".out qmc/cusp")
	
# ======================= optJastrow ======================= #
def optJastrow():
	addOption2QMCline("updet"," cuspInfo=\"spo-up.cuspInfo.xml\">\n")
	addOption2QMCline("downdet"," cuspInfo=\"spo-dn.cuspInfo.xml\">\n")
	add_QMCblock_before("</wavefunction>","jastrow.xml")
	
	qmcblock="optJastrow"
	os.system("sed 's/molecule/"+molecule+"/g' "+template_main+" > "+qmcblock+".xml")
	fileline2file(qmcblock+".xml","qmcblock",template_dir+qmcblock+".xml")
	
	print("============================================")
	print("Optimizing Jastrows")
	print("Check progress in "+qmcblock+".out")
	with open(qmcblock+".out",'w') as outfile:
		subprocess.call(["qmcapp",qmcblock+".xml"],stdout=outfile)
	
	qmcblock="optJastrow"
	# folder management
	new_folder="qmc/"+qmcblock
	
	os.system("cp "+molecule+".s005.opt.xml "+opt_wfs)
	
	os.system("mkdir -p "+new_folder+";mv newOrbs* jastrow.xml J1* J2* "+new_folder+";mv "+molecule+".s* "+qmcblock+".xml "+new_folder+"; mv "+qmcblock+".out "+new_folder)
	
	
# ======================= rundmc ======================= #
def rundmc():
	# write molecule-specified data to qmc template
	qmcblock="dmc"
	os.system("sed 's/molecule/"+molecule+"/g' "+template_main+" > "+qmcblock+".xml")
	fileline2file(qmcblock+".xml","qmcblock",template_dir+qmcblock+".xml")
	os.system("sed -i 's/"+wfs+"/"+opt_wfs+"/' "+qmcblock+".xml")
	
	# run dmc
	print("============================================")
	print("Running DMC")
	print("Check progress in "+qmcblock+".out")
	with open(qmcblock+".out",'w') as outfile:
		subprocess.call(["qmcapp",qmcblock+".xml"],stdout=outfile)
		
	# folder management
	new_folder="qmc/"+qmcblock
	os.system("mkdir -p "+new_folder+";mv newOrbs* J1* J2* "+new_folder+";mv "+qmcblock+".xml "+new_folder+"; mv "+qmcblock+".out "+new_folder)
	

# ======================= main ======================= #
rungms()
convert4qmc()
cuspCorrection()
optJastrow()
#rundmc()


