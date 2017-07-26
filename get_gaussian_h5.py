#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
from itertools import product

def naive_cubic_grid(alat,grid_shape,func):
  """ put 3D function 'func' on grid, one point at a time """
  orb_grid   = np.zeros(grid_shape)
  dx,dy,dz   = alat/grid_shape
  for ix in range(grid_shape[0]):
    for iy in range(grid_shape[1]):
      for iz in range(grid_shape[2]):
        myr = np.array([ix,iy,iz]) * np.array([dx,dy,dz])
        orb_grid[ix,iy,iz] = func(myr)
      # end for iz
    # end for iy
  # end for ix
  return orb_grid
# end def naive_cubic_grid

def img_gauss(rvec,center,axes,alpha):
  dr = rvec-center # naive displacement

  # minimum image
  dists = []
  images = product([-1,0,1],repeat=3) # check all neighboring images
  for image_shift in images:
    disp  = rvec - center + np.dot(image_shift,axes)
    dists.append( np.linalg.norm(disp) )
  # end for

  pdr = min(dists)
  r2 = np.dot(pdr,pdr)
  return np.exp(-alpha*r2)
# end def img_gauss

def spline_orbitals(aoR,mo_coeff,gs,kpt=[0,0,0],ikpt=0,ispin=0):
  # e.g. gs=[8 8 8]; fftgrid_shape = 2*gs+1
  from pyscf_orbital_routines import mo_orbitals

  ngrid,nao0 = aoR.shape
  nao,nmo    = mo_coeff.shape
  assert nao==nao0

  # energies of the orbitals are use to sort the orbitals only
  fake_mo_energy = np.arange(nmo)

  fft_normalization = 1.0
  gvecs,eig_df = mo_orbitals(fake_mo_energy,mo_coeff,aoR,gs,fft_normalization)

  # finish dataframe
  eig_df['ikpt']  = ikpt
  eig_df['ispin'] = ispin
  kpt_data = np.zeros([len(eig_df),len(kpt)])
  kpt_data[:] = kpt # copy kpt to every row of kpt_data
  eig_df['reduced_k'] = kpt_data.tolist()
  eig_df.set_index(['ikpt','ispin','istate'],inplace=True,drop=True)
  eig_df.sort_index()
  return gvecs,eig_df
# end def spline_orbitals

def get_gaussian_aoR(axes,pos,gs):
  alat = axes[0,0]
  assert np.allclose(axes,np.eye(3)*alat), 'assuming cubic cell'
  inv_axes = np.linalg.inv(axes)
  # assuming x,y,z = [0,alat)
  grid_shape = 2*gs+1

  natom = len(pos)
  aoR = np.zeros([np.prod(grid_shape),natom])
  for iatom in range(natom):
    center = pos[iatom]
    func   = lambda x:img_gauss(x,center,axes,expo)
    orb_grid = naive_cubic_grid(alat,grid_shape,func) # put orbital on grid
    aoR[:,iatom] = orb_grid.flatten()
  # end for iatom
  return aoR
# end def get_gaussian_aoR

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('inp_xml',type=str,help='QMCPACk input xml file to read axes,pos from')
  parser.add_argument('expo',type=float,help='gaussian exponent in a.u. e.g. 9.0')
  parser.add_argument('gs0',type=int,help='gs=(gs0,gs0,gs0)')
  args = parser.parse_args()

  from input_xml import InputXml
  inp = InputXml()
  inp.read(args.inp_xml)
  expo = args.expo
  gs = np.array([args.gs0]*3) # FFT grid shape is 2*gs+1

  wf_h5_fname = 'exp%3.2f_gs%d.h5' % (expo,gs[0])

  # get structure (i.e. axes,pos)
  axes = inp.lattice_vectors()
  pos  = inp.atomic_coords(pset_name='wf_centers')
  natom = len(pos)

  # put gaussian orbitals on real-space grid
  aoR = get_gaussian_aoR(axes,pos,gs)
  mo_coeff = np.eye(natom,natom)

  # FFT to put orbitals on PW basis
  gvecs,eig_df = spline_orbitals(aoR,mo_coeff,gs)

  # make wavefunction h5 file
  elem = ['H']*natom
  from pyscf_orbital_routines import generate_pwscf_h5, atom_text
  from pyscf.pbc.gto.cell import Cell
  cell = Cell()
  cell.build(a=axes,atom=atom_text(elem,pos),unit='B',gs=gs)
  generate_pwscf_h5(cell,gvecs,eig_df,h5_fname=wf_h5_fname)

# end __main__
