#!/usr/bin/env python
import numpy as np
import h5py

def reblock(data,block):
    n=int(np.floor(float(len(data))/block))
    newdata=[]
    for i in range(n):
        newdata.append(data[i*block:(i+1)*block].mean())
    # end for i
    return np.array(newdata)
# end def reblock

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

def error(data):
    return data.std()/np.sqrt(len(data)/corr(data))
# end def error

import matplotlib.pyplot as plt
import argparse
from scipy import stats
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Extract a column of scalar data')
    parser.add_argument("scalar", type=str, help="scalar file name")
    parser.add_argument('-e','--equil', type=int, default=0, help="number of equilibration steps")
    parser.add_argument('-c','--column', type=int, default=0, help="column to analyze")
    parser.add_argument('-rb','--reblock', type=int, default=1, help="reblock data")
    parser.add_argument('-p','--plot', action='store_true', help="plot data")
    parser.add_argument('-d','--dump', action='store_true', help="dump data")
    args = parser.parse_args()
    
  
    # read data
    data = np.loadtxt(args.scalar).T

    if not isinstance(data[0], (float,int,long,complex)):
        # if no column specified, list options
        with open(args.scalar) as f:
            headline = f.readline()
            tokens = headline.strip('#').split()
        # end with open
        if (args.column==0):
            print np.array( [range(len(tokens)),tokens] ).T
        # end if args.column==0
        data = data[args.column,args.equil:]
    else:
        data = data[args.equil:]
    # if len(data)!=1

    data = reblock(data,args.reblock)
    correlation=corr(data)

    if not isinstance(data[0], (float,int,long,complex)):
        print(tokens[args.column]),
    # end if not isinstance
    print data.mean(),"+-",data.std()/np.sqrt(len(data)/correlation),correlation

    # plot
    if args.plot:
        plt.plot(data)
        plt.show()
    # end if

    # dump
    if args.dump:
        np.savetxt("data.dat",data)
    # end if dump

# end __main__
