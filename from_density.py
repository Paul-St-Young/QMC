#!/usr/bin/env python

import h5py
from lxml import etree
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# helper functions
# ====================

def grab_density_parameters(qmc_input_xml):
    root = etree.parse(qmc_input_xml)
    dens_est_list = root.xpath('//estimator[@type="density"]')
    if len(dens_est_list) != 1:
        print "expected 1 density estimator, found: ", len(dens_est_list)
        raise InputError()
    # end if
    dens_est = dens_est_list[0]
    nx,ny,nz = [int(1./float(x)) for x in dens_est.attrib["delta"].split()]

    x = np.linspace( float(dens_est.attrib["x_min"]),
                    float(dens_est.attrib["x_max"]),nx)
    y = np.linspace( float(dens_est.attrib["y_min"]),
                    float(dens_est.attrib["y_max"]),ny)
    z = np.linspace( float(dens_est.attrib["z_min"]),
                    float(dens_est.attrib["z_max"]),nz)
    
    return x,y,z
# end def grab_density_parameters

def grab_density(h5file):
    # get density
    f = h5py.File(h5file)
    for name,quantity in f.items():
        if name.startswith('density'):
            density = quantity.get("value")[:]
        # end if name.startswith
    # end for name,quantity
    f.close()

    return density
# end def grab_density

def a2w(alpha):
    return 1./np.sqrt(2.*alpha)
def w2a(sigma):
    return 1./(4.*sigma**2.)


# fitting procedure
# ====================

def gauss_1d(x,A,xo,alpha):
    return A*np.exp(-alpha*(x-xo)**2.)
# end def gauss_1d

def fit_1d_gaussian(x,y,z,density,ro,plot=True):
    
    # find density maximum
    max_density = max( density[density.nonzero()] )
    max_idx = np.where(density==max_density)
    max_idx = [item[0] for item in max_idx]
    max_idx[2] += 1

    max_xyz = x[max_idx[0]],y[max_idx[1]],z[max_idx[2]]
    
    # x slice
    xslice = density[:,max_idx[1],max_idx[2]]
    xcoeff,var = curve_fit(gauss_1d,x,xslice,p0=(max(xslice),ro[0],10))
    alphax = xcoeff[-1]
    
    
    # y slice
    yslice = density[max_idx[0],:,max_idx[2]]
    ycoeff,var = curve_fit(gauss_1d,y,yslice,p0=(max(yslice),ro[1],10))
    alphay = ycoeff[-1]
    
    
    # z slice
    zslice = density[max_idx[0],max_idx[1],:]
    zcoeff,var = curve_fit(gauss_1d,z,zslice,p0=(max(zslice),ro[2],10))
    alphaz = zcoeff[-1]

    if plot:
        fig = plt.figure()

        # x slice
        ax  = fig.add_subplot(131)
        ax.set_xlabel("x (bohr)",fontsize=14)
        ax.set_ylabel("density slice (arb. u.)",fontsize=14)

        ax.plot(x,xslice,".")
        ax.plot(x,gauss_1d(x,*xcoeff),lw=2,color="black")
        ax.set_xticks( np.linspace(min(x),max(x),4) )
        ax.get_xaxis().set_major_formatter( FormatStrFormatter("%1.2f") )

        # y slice
        ax  = fig.add_subplot(132)
        ax.set_xlabel("y (bohr)",fontsize=14)

        ax.plot(y,yslice,".")
        ax.plot(y,gauss_1d(y,*ycoeff),lw=2,color="black")
        ax.set_xticks( np.linspace(min(y),max(y),4) )
        ax.get_xaxis().set_major_formatter( FormatStrFormatter("%1.2f") )

        # z slice
        ax  = fig.add_subplot(133)
        ax.set_xlabel("z (bohr)",fontsize=14)

        ax.plot(z,zslice,".")
        ax.plot(z,gauss_1d(z,*zcoeff),lw=2,color="black")
        ax.set_xticks( np.linspace(min(z),max(z),4) )
        ax.get_xaxis().set_major_formatter( FormatStrFormatter("%1.2f") )

        fig.tight_layout()

    # end if plot

    return alphax,alphay,alphaz
# end def fit_1d_gaussian

def density_from_1d_fit(qmc_input_xml,qmc_stat_h5,nequil,ro,plot=True):
    
    # load data,extract equilibrated piece
    x,y,z = grab_density_parameters(qmc_input_xml)
    density = grabDensity(qmc_stat_h5)
    equil = density[nequil:].sum(axis=0)
    
    alphas = fit_1d_gaussian(x,y,z,equil,ro,plot)
    
    return np.array(alphas)
# end def
