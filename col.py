#!/usr/bin/env python
import numpy as np
import h5py

def corr(myg):
    # autocorrelation
    g = np.array(myg)
    mu=g.mean()
    s=g.std()

    sumR=0.0
    for k in range(1,len(g)):
        R=0.0

        # calculate autocorrelation
        for t in range(len(g)-k):
            R+=(g[t]-mu)*(g[t+k]-mu)
        #end for t
        R/=(len(g)-k)*s**2.

        # accumulate until R<=0
        if R>0:
            sumR+=R
        else:
            break
        #end if
    #end for k
    
    return 1+2.*sumR
#end def corr

import matplotlib.pyplot as plt
import argparse
from scipy import stats
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Extract a column of scalar data')
    parser.add_argument("scalar", type=str, help="scalar file name")
    parser.add_argument('-e','--equil', type=int, default=0, help="number of equilibration steps")
    parser.add_argument('-c','--column', type=int, default=0, help="column to analyze")
    args = parser.parse_args()
    
    # list options
    with open(args.scalar) as f:
        headline = f.readline()
        tokens = headline.strip('#').split()
    # end with open
    if (args.column==0):
        print np.array( [range(len(tokens)),tokens] ).T
    # end if args.column==0
  
    # read data
    data = np.loadtxt(args.scalar).T
    data = data[args.column,args.equil:]
    correlation=corr(data)
    print tokens[args.column],data.mean(),"+-",data.std()/np.sqrt(len(data)/correlation),correlation
    # plot
    plt.plot(data)
    plt.show()

# end __main__
