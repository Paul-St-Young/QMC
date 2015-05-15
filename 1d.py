#!/usr/bin/env python

# finite difference for solving 1D Schrodinger

import numpy as np

def getGroundState(x,V):
    # using Hartree atomic units x MUST be in Bohr and V MUST be in Hatree
    dx = x[1]-x[0]
    # set up Hamiltonian
    m = 1822.88839
    t = 1./(2*m*dx**2)
    Hdiag = V+2*t*np.ones(len(x))
    Hoffd = -t*np.ones(len(x)-1)
    H     = np.diag(Hdiag) + np.diag(Hoffd,1) + np.diag(Hoffd,-1)
    
    # get ground state
    w,v = np.linalg.eig(H)
    minIdx = list(w).index(w.min())
    psi = v[:,minIdx]
    
    return w[minIdx], psi
# end def getGroundState

import matplotlib.pyplot as plt
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description='Solve for ground state proton density on given cold curve')
    parser.add_argument('coldFile',type=str
        ,help='File containing cold curve in x:y:yerr format')
    parser.add_argument('-n','--npoints',type=int,default=100
        ,help='Number of grid points for finite difference')
    parser.add_argument('-o','--order',type=int,default=5
        ,help='Order of polynomial in fit')
    parser.add_argument('-x1','--xmin',type=float,default=None
        ,help='xmin to use from fitted curve') 
    parser.add_argument('-x2','--xmax',type=float,default=None
        ,help='xmax to use from fitted curve')    
    args = parser.parse_args()
    
    # get cold curve
    X, Y, Ye = np.loadtxt(args.coldFile).T
    V = np.poly1d( np.polyfit(X,Y,args.order) )
    
    # get ground state proton density on fitted cold curve with fine grid
    xmin, xmax = args.xmin, args.xmax
    if args.xmin==None:
        xmin=X[0]
    # end if xmin
    if args.xmax==None:
        xmax=X[-1]
    # end if xmax
    x=np.arange(xmin,xmax,(xmax-xmin)/args.npoints)
    E,psi = getGroundState(x,V(x))
    density = psi*psi
    
    # calculate Re, Ro
    print "Re =", round( x[list(V(x)).index(V(x).min())] ,4), " by fit min"
    print "Ro =", round( sum(density*x)/sum(density) ,4), " by integrate"
    
    fig, ax1 = plt.subplots()

    
    # plot density
    ax1.plot(x,density,label="1D Grid")
    ax1.set_ylabel(r'$\psi_o^*\psi_o$',fontsize=16)
    ax1.set_xlabel(r'$r_{CH} (Bohr)$',fontsize=16)
    ax1.legend(loc=0)
    
    # plot curve
    ax2 = ax1.twinx()
    ax2.plot(X,Y,'o',markersize=2)
    if len(Ye)==len(X):
        ax2.errorbar(X,Y,yerr=Ye,fmt=None,color='b')
    # end if Ye
    
    # plot fit
    ax2.plot(x,V(x),'-r')
    ax2.set_ylabel('V (Hatree)',fontsize=14)
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    # end for
    
    plt.title("HCN+ Proton Density in 1D", fontsize=16)
    plt.show()
    
# end __main__



