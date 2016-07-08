#!/usr/bin/env python

import numpy as np

# Fourier grid
# =====
def kint_to_kvec(kint,nk,alat):
    """ convert an integer kint, to kpoint given:
         1. real space box length: alat
         2. number of grid points: nk """

    # In discrete fourier transform, if implicit real space size is alat,
    #  then k-space grid point spacing is 2pi/alat. This is because finite
    #  real space is implicitly periodic.
    k_space_unit = 2*np.pi/alat
    # ! warning: real space grid spacing is not physical, the number of 
    #  real space grid point does NOT change k-space unit. Only the total
    #  size of real space change k-space unit/spacing.

    # even vs. odd entry mapping, consult numpy.fft manual
    """
    if nk%2==0:
        if kint <= nk/2:
            return kint * k_space_unit
        else:
            return (kint-nk) * k_space_unit
        # end if
    else:
        print "code in odd grid size you lazy!"
        raise IndexError
    # end if
    """
    # turns out even and odd treatments are the same
    if kint <= nk/2:
        return kint * k_space_unit
    else:
        return (kint-nk) * k_space_unit
    # end if

# end def

def get_3d_kgrid(fk,alat,nk,exclude_gamma=False):
    """ put a k-space function fk onto an FFT grid for inverse Fourier transform
     fr = np.fft.ifft(fk)*(alat/nk)**ndim ! Don't forget the normalization! """
    kgrid = np.zeros([nk,nk,nk],dtype=complex)
    kint_unit = lambda x:kint_to_kvec(x,nk,alat)
    for kx in range(nk):
        for ky in range(nk):
            for kz in range(nk):
                if exclude_gamma and np.all([kx,ky,kz]==[0,0,0]):
                    continue
                # end if
                gvec = np.array(map(kint_unit,[kx,ky,kz]))
                kgrid[kx,ky,kz] = fk(gvec)
            # end for
        # end for
    # end for
    return kgrid
# end def

# inverse Fourier transform
# =====
def inv_fft(kgrid,alat):
    # !!! assume isotropic uniform cubic grid
    ndim = len(kgrid.shape)
    nk   = kgrid.shape[0]
    dx   = float(alat)/nk
    half = float(alat)/2.
    if nk%2==0:
        x = np.arange(-half,half,dx)
    else:
        x = np.linspace(-half,half,nk)
    # end if
    fr = np.fft.ifftn(kgrid)/dx**ndim
    fr = np.fft.fftshift(fr)
    return x,fr
# end def

def fft(fr,alat):

    # !!! assume isotropic uniform cubic grid
    ndim = len(fr.shape)
    nk   = fr.shape[0]
    dx   = float(alat)/nk

    kgrid = np.fft.fftn( np.fft.ifftshift(fr)*dx**ndim )
    return kgrid

# end def
