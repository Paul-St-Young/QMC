#!/usr/bin/env python

# finite difference for solving 1D Schrodinger

import numpy as np

def bestFit(xmin=0.7,xmax=1.6,dx=1e-2):
    # best fit V(x) = a+b*(x-c)*(x-c)+e*(x-c)**3+f*(x-c)**4
    a    = -92.8291 
    b    = 0.0205672
    c    = 1.82335  
    e    = 0.770105 
    f    = 0.760898

    # set up potential
    x = np.arange(xmin,xmax,dx)
    V = a+b*(x-c)*(x-c)+e*(x-c)**3+f*(x-c)**4
    
    return x,V
# end def bestFit

def integral(x,dx,density):
    A = 0.0
    for i in range(len(x)):
        A += density[i]*dx
    # end for i
    return A
# end def integral

def normalize(x,dx,density):
    A = 0.0
    for i in range(len(x)):
        A += density[i]*dx
    # end for i
    for i in range(len(x)):
        density[i] /= A
    # end for i
# end def normalize

def getGroundState(x,V):
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

def loadCurve(curvefile,mol="HCN+",dmcline=2):
    datline = 0
    curveX = []
    curveY = []
    curveYe = []
    f = open(curvefile,'r')
    for line in f:
        tokens = line.split()
        if len(tokens)==1:
            curveX.append( float(tokens[0]) )
        elif len(tokens)>1 and tokens[0]==mol:
            datline += 1
            if datline==dmcline:
                curveY.append( float(tokens[3]) )
                curveYe.append( float(tokens[5]) )
                datline=0
        # end if
    # end for 
    f.close()
    return curveX, curveY, curveYe
# end def

import matplotlib.pyplot as plt
import mix_fix
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='plot HCN+ curves')
    parser.add_argument('VMC',type=str,help='h5 file containing VMC gofr')
    parser.add_argument('DMC',type=str,help='h5 file containing VMC gofr')
    args = parser.parse_args()

    xmin = 0.7
    xmax = 1.6

    r, fixed_gr = mix_fix.perform_quick_fix(args.VMC,args.DMC)
    ang_per_bohr = 0.529177249
    r *= ang_per_bohr
    dr = r[1]-r[0]
    normalize(r,dr,fixed_gr)
    start,end = mix_fix.get_start_and_end(fixed_gr,1e-2)
    coeff = mix_fix.fit_gaussian(r,fixed_gr,start,end)

    x,V = bestFit()
    E,psi = getGroundState(x,V)
    curveX, curveY, curveYe = loadCurve("HCN+-cas.dat")

    # get expectation
    xavg = 0.0
    norm = 0.0
    dx = x[1]-x[0]
    for i in range(len(x)):
        xavg += psi[i]*psi[i]*x[i]*dx
        norm += psi[i]*psi[i]*dx
    # end for i

    # plot rCH density
    fig, ax1 = plt.subplots()
    density = psi*psi
    normalize(x,dx,density)
    ax1.plot(x,density,label="1D Grid")
    ax1.plot(r, mix_fix.gaussian(r,*coeff), '--',label="non-adiabatic fit")
    #ax1.plot(r,fixed_gr,'--')
    ax1.set_ylabel(r'$\psi_o^*\psi_o$',fontsize=20)
    ax1.set_xlabel(r'$r_{CH} (\AA)$',fontsize=20)
    ax1.legend(loc=0)

    # plot curve and fit
    ax2 = ax1.twinx()
    ax2.plot(curveX,curveY,'o',markersize=2)
    ax2.errorbar(curveX,curveY,yerr=curveYe,fmt=None,color='b')
    ax2.plot(x,V,'-r')
    ax2.set_ylabel('V (Hatree)',fontsize=16)
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    # end for
    
    plt.title("HCN+ Proton Density in 1D", fontsize=20) #$\langle r_{CH}\rangle$="+str( round(xavg/norm,4) ) )
    plt.xlim(xmin,xmax)
    plt.show()
# end __main__



