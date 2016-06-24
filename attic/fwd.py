#!/usr/bin/env python
import numpy as np

def fwd(scalarfile,fwdStart,fwdNumber,equil):
    data = np.loadtxt(scalarfile).T
    X=[]
    V=[]
    Ve=[]
    for i in range(fwdNumber):
        X.append(i)
        curdat = data[fwdStart+i,equil:]
        V.append( curdat.mean() )
        Ve.append( error(curdat) )
    # end for i
    return X,V,Ve
# end def fwd

from col import error
import matplotlib.pyplot as plt
import argparse
from scipy import stats
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Plot time evolution of g(r)')
    parser.add_argument("scalarfile", type=str
        , help="file name in the format ${molecule}.s00x.scalar.dat")
    parser.add_argument('-e','--equil', type=int, default=0
        , help="number of equilibration blocks")
    parser.add_argument('-s','--start', type=int, default=10
        , help="first column that contains forwadwalking data")
    parser.add_argument('-n','--number', type=int, default=10
        , help="number of forwad walking data points")
    parser.add_argument('-b','--best', type=int, default=6
        , help="best data point")
    parser.add_argument('-p','--plot', action='store_true'
        , help="plot convergence")
    #parser.add_argument('-d','--dump', action='store_true'
        #, help="dump test data")
    args = parser.parse_args()
  
    # read data 
    X,V,Ve = fwd(args.scalarfile,args.start,args.number,args.equil)

    # process best value
    print V[args.best],"+-",Ve[args.best]

    # plot
    if args.plot:
        fig,ax=plt.subplots()
        ax.plot(V,'-x')
        ax.errorbar(X,V,yerr=Ve,fmt=None,color='b')
        plt.show() 
    # end if plot

# end __main__
