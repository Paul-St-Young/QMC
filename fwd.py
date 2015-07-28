#!/usr/bin/env python
import numpy as np

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
    parser.add_argument('-n','--number', type=int, default=20
        , help="number of forwad walking data")
    #parser.add_argument('-d','--dump', action='store_true'
        #, help="dump test data")
    args = parser.parse_args()
  
    # read data 
    data = np.loadtxt(args.scalarfile).T

    # process best value
    datStart = args.start
    datMax   = args.number
    print "after",args.equil,"blocks of equilibration, the data has length", len(data[datStart,args.equil:])
    dat = data[datStart+datMax-1,args.equil:]
    print dat.mean(),"+-",error(dat)

    """
    testdat = data[9,args.equil:]
    print testdat.mean(),"+-",stats.sem(testdat)
    if (args.dump):
        with open("dist.dat",'w') as f:
            for i in range(len(testdat)):
                f.write( str(testdat[i])+"\n" )
            # end for i
        # end with open
    # end if dump
    plt.plot(testdat,'-')
    plt.show()
    """

    # plot
    X=[]
    V=[]
    Ve=[]
    for i in range(datMax):
        X.append(i)
        curdat = data[datStart+i,args.equil:]
        V.append( curdat.mean() )
        Ve.append( error(curdat) )
    # end for i
    fig,ax=plt.subplots()
    ax.plot(V,'-x')
    ax.errorbar(X,V,yerr=Ve,fmt=None,color='b')
    plt.show() 

    """
    dat = data[27,args.equil:]
    plt.plot(dat)
    plt.show()
    """
# end __main__
