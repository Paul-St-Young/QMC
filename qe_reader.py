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
    gvecs    = h5handle['electrons/kpoint_0/gvectors'].value
    h5handle.close()
    return {'axes':lattice,'pos':pos,'gvecs':gvecs}
# end def

import subprocess as sp
def find_pwscf_io(path,infile_subfix='-scf.in',outfile_subfix='.out',use_last=False):
    # assuming there is only 1 pair of pw.x input and output in path
    #  return the names of the input and output files
    out = sp.check_output(['ls',path])
    
    infile  = ''
    outfile = ''
    
    found_in  = False
    found_out = False
    for fname in out.split('\n')[:-1]:
        if fname.endswith(infile_subfix):
            if found_in and not use_last:
                raise NotImplementedError('multiple inputs found in %s'%path)
            # end if
            infile   = fname
            found_in = True
        elif fname.endswith(outfile_subfix):
            if found_out and not use_last:
                raise NotImplementedError('multiple outputs found in %s'%path)
            # end if
            outfile   = fname
            found_out = True
        # end if
    # end for fname
    
    if not found_in:
        raise IOError('infile not found in %s'%path)
    elif not found_out:
        raise IOError('outfile not found in %s'%path)
    # end if
    
    return infile,outfile
# end def find_pwscf_io

def all_lines_with_tag(mm,tag,nline_max):
    """ return a list of memory indices pointing to the start of tag """
    mm.seek(0) # rewind file
    all_idx = []
    for iline in range(nline_max):
        idx = mm.find(tag)
        if idx == -1:
            break
        # end if
        mm.seek(idx)
        all_idx.append(idx)
        mm.readline()
    # end for iline
    
    # guard
    if iline >= nline_max-1:
        raise NotImplementedError('may need to increase nline_max')
    # end if
    return all_idx
# end def all_lines_with_tag

import struct
def available_structures(pw_out,nstruct_max=2000,natom_max=1000,ndim=3
        ,variable_cell=False):
    """ find all available structures in a pwscf output """
    fhandle = open(pw_out,'r+')
    mm = mmap(fhandle.fileno(),0)

    idx = mm.find('lattice parameter')
    mm.seek(idx)
    lat_line = mm.readline()
    alat = float( lat_line.split()[-2] )
    
    # locate all axes
    axes_tag = 'CELL_PARAMETERS ('.encode()
    axes_starts = all_lines_with_tag(mm,axes_tag,nstruct_max)
    naxes = len(axes_starts)
    if (naxes != 0) and (not variable_cell):
        raise NotImplementedError('CELL_PARAMETERS found, are you sure this is not a variable cell run?')
    # end if
    
    # locate all atomic cd positions
    pos_tag    = 'ATOMIC_POSITIONS'.encode()
    pos_starts = all_lines_with_tag(mm,pos_tag,nstruct_max)
    npos = len(pos_starts)
    
    if variable_cell and (npos != naxes):
        raise NotImplementedError('expect same number of cells as atomic positions in a variable cell calculation. got (naxes,npos)=(%d,%d)'%(naxes,npos))
    # end if
    
    # count number of atoms
    mm.seek(pos_starts[0])
    mm.readline() # skip tag line
    natom = 0
    for iatom in range(natom_max):
        line = mm.readline()
        tokens = line.split()
        if len(tokens) != 4:
            break
        # end if
        natom += 1
    # end for iatom
    
    # read initial crystal axes
    axes = np.zeros([ndim,ndim])
    if not variable_cell: 
        idx = all_lines_with_tag(mm,'crystal axes'.encode(),nstruct_max)[0]
        mm.seek(idx)
        tag_line = mm.readline()
        unit_text= tag_line.split()[-1].strip('()')
        for idim in range(ndim):
            line = mm.readline()
            axes[idim,:] = line.split()[3:3+ndim]
            if 'alat' in unit_text:
                axes[idim,:] *= alat
            else:
                raise NotImplementedError('crystal axes: what unit is %s?'%unit_text)
            # end if
        # end for 
    # end if
    
    bohr = 0.52917721067 # angstrom (CODATA 2014)
    angstrom = False # assume bohr
    nstructs = max(naxes,npos)
    all_axes = np.zeros([nstructs,ndim,ndim])
    all_pos  = np.zeros([nstructs,natom,ndim])
    for istruct in range(nstructs):
        
        if variable_cell: # read cell parameters
            cell_idx = axes_starts[istruct]
            mm.seek(cell_idx)
            mm.readline() # skip tag line
            
            axes_text = ''
            for idim in range(ndim):
                axes[idim,:] = mm.readline().split()
            # end for idim
        # end if variable_cell
        
        all_axes[istruct,:,:] = axes
        
        pos_idx = pos_starts[istruct]
        mm.seek(pos_idx)
        tag_line = mm.readline()
        unit_text= tag_line.split()[-1]
        if 'angstrom' in unit_text:
            angstrom = True
        elif 'crystal' in unit_text:
            raise NotImplementedError('crsytal units')
        else:
            raise NotImplementedError('what unit is this? %s' % unit_text)
        # end if
        
        for iatom in range(natom):
            line = mm.readline()
            name = line.split()[0]
            pos_text = line.strip(name)
            try:
                name,xpos,ypos,zpos = struct.unpack('4sx14sx14sx13s',pos_text)
                all_pos[istruct,iatom,:] = [xpos,ypos,zpos]
                if angstrom:
                    all_pos[istruct,iatom,:] /= bohr
            except:
                print 'failed to read (istruct,iatom)=(%d,%d)' % (istruct,iatom)
            # end try
        # end for iatom
        
    # end for istruct
        
    fhandle.close()
    return all_axes,all_pos
# end def available_structures

def md_traces(md_out,nstep=2000):
    """ extract scalar traces from pwscf md output md_out 
     look for tags defined in line_tag_map """
    fhandle = open(md_out,'r+')
    mm = mmap(fhandle.fileno(),0)

    line_tag_map = { # unique identifier of the line that contains each key
        'fermi energy':'the Fermi energy is',
        'total energy':'!',
        'kinetic energy':'kinetic energy',
        'temperature':'temperature',
        'econst':'Ekin + Etot'
    }
    val_idx_map  = {} # assume -2
    val_type_map = {} # assume float

    mm.seek(0)
    data = []
    for istep in range(nstep):
        if mm.tell() >= mm.size():
            break
        # end if
        found_stuff = False
        
        entry = {'istep':istep}
        for label in line_tag_map.keys():
            
            # locate line with value for label
            idx = mm.find(line_tag_map[label].encode())
            if idx == -1:
                continue
            # end if
            found_stuff = True
            mm.seek(idx)
            line = mm.readline()
            
            # locate value in line
            rval_idx = -2 # assume patten "label = value unit"
            if label in val_idx_map.keys():
                rval_idx = val_idx_map[label]
            # end if
            rval = line.split()[rval_idx]
            
            # convert value
            val_type = float
            if label in val_type_map.keys():
                val_type = val_type_map[key]
            # end if
            value = val_type(rval)
            
            entry[label] = value # !!!! assume float value
        # end for
        
        if found_stuff:
            data.append(entry)
        else:
            break
        # end if
    # end for istep
    if istep >= nstep-1:
        print "WARNING: %d structures found, nstep may need to be increased" %istep
    # end if
    fhandle.close()
    return data
# end def md_traces# end def md_traces