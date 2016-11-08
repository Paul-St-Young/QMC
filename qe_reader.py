import numpy as np
from mmap import mmap

def retrieve_occupations(nscf_outfile, max_nbnd_lines=10):
    """ read the eigenvalues and occupations of DFT orbitals at every available kpoint in an non-scf output produced by pwscf """ 

    fhandle = open(nscf_outfile,'r+')
    mm = mmap(fhandle.fileno(),0)

    # read number of k points
    nk_prefix = "number of k points="
    idx = mm.find(nk_prefix)
    mm.seek(idx)

    nk_line = mm.readline()
    nk = int( nk_line.strip(nk_prefix).split()[0] )

    # skip to the end of band structure calculation
    idx = mm.find('End of band structure calculation')
    mm.seek(idx)

    # read the eigenvalues and occupations at each kpoint
    kpt_prefix = "k ="
    data = []
    for ik in range(nk):
        idx = mm.find(kpt_prefix)
        mm.seek(idx)
        kpt_line = mm.readline()
        kpt = map(float,kpt_line.strip(kpt_prefix).split()[:3])

        mm.readline() # skip empty line
        eval_arr = np.array([])
        for iline in range(max_nbnd_lines):
            tokens = mm.readline().split()
            if len(tokens)==0:
                break
            # end if
            eval_arr = np.append(eval_arr, map(float,tokens))
        # end for iline
        
        idx = mm.find('occupation numbers')
        mm.seek(idx)
        mm.readline() # skip current line
        occ_arr = np.array([])
        for iline in range(4):
            tokens = mm.readline().split()
            if len(tokens)==0:
                break
            # end if
            occ_arr = np.append(occ_arr, map(float,tokens))
        # end for iline

        entry = {'ik':ik,'kpt':list(kpt),'eval':list(eval_arr),'occ':list(occ_arr)}
        data.append(entry)
    # end for
    mm.close()
    fhandle.close()
    return data
# end def

import os
import h5py
def retrieve_psig(h5_file,only_occupied=False,occupations=None):
    """ return a list dictionaries of DFT orbital coefficients in PW basis by reading an hdf5 file written by pw2qmcpack. If only_occupied=True and a database of occupied orbitals are given, then only read orbitals that are occupied. """
    if only_occupied and (occupations is None):
        raise NotImplementedError("no occupation database is given")
    # end if

    ha = 27.21138602 # ev from 2014 CODATA

    orbitals = []

    h5handle = h5py.File(h5_file)
    electron = h5handle['electrons']

    kpt_labels = []
    for key in electron.keys():
        if key.startswith('kpoint'):
            kpt_labels.append(key)
        # end if
    # end for key

    nk = electron['number_of_kpoints'].value
    assert nk==len(kpt_labels)

    for label in kpt_labels:

        # get kpoint index
        kpt_idx = int( label.split('_')[-1] )

        # get plane wave wave numbers
        if kpt_idx == 0:
            mypath = os.path.join(label,'gvectors')
            gvecs = electron[mypath].value
        # end if

        # verify eigenstates at this kpoint
        kpt_ptr = electron[os.path.join(label,'spin_0')]
        nbnd = kpt_ptr['number_of_states'].value

        evals = kpt_ptr['eigenvalues'].value

        # compare to nscf output (eigenvalues and occupation)
        if occupations is not None:
            mydf = occupations[occupations['ik']==kpt_idx]
            myval= mydf['eval'].values[0]
            myocc= mydf['occ'].values[0]
            assert nbnd == len(myval), "expect %d bands, nscf has %d bands" % (nbnd,len(myval))
            assert np.allclose(evals*ha,myval,atol=1e-4), str(evals*ha-myval)
        # end if

        for iband in range(nbnd):
            if only_occupied and (np.isclose(myocc[iband],0.0)):
                continue
            # end if
            psig  = kpt_ptr['state_%d/psi_g'%iband].value
            entry = {'ik':kpt_idx,'iband':iband,'eval':evals[iband],'psig':psig}
            orbitals.append(entry)
        # end for iband

    # end for label

    h5handle.close()
    return orbitals
# end def retrieve_psig

def retrieve_system(h5_file):
    h5handle = h5py.File(h5_file)
    lattice  = h5handle['supercell/primitive_vectors'].value
    #elem     = h5handle['atoms/species_ids'].value # why so complicated?
    pos      = h5handle['atoms/positions'].value
    gvecs    = h5handle['electrons/kpoint_0/gvectors']
    h5handle.close()
    return {'axes':lattice,'pos':pos,'gvecs':gvecs}
# end def
