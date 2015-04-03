#!/usr/bin/env python

import numpy as np

def rehist(r,gr,group=2):
    # coarsen gr by a factor of group, but improve accuracy on each point
    newr  = np.array([  r[i-group:i] for i in range(group,len(r),group)  ]).mean(axis=1)
    newgr = np.array([ gr[i-group:i] for i in range(group,len(gr),group) ]).mean(axis=1)
    return newr,newgr
# end def
def nonzero(r,vmc,dmc,thres):
    # return r,vmc,dmc only on the interval where both vmc and dmc are nonzero
    start   = 0
    end     = len(r)-1
    started = False
    for i in range(len(r)):
        if (not started) and (vmc[i]>thres or dmc[i]>thres):
            started = True
            start = i
        # end if not started
        if started and (vmc[i]<thres and dmc[i]<thres):
            end = i
            break
        #end if started
    # end for
    return r[start:end],vmc[start:end],dmc[start:end]
# end def nonzero

from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import argparse
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Use g(r) to fix wf')
    parser.add_argument('VMC',type=str,help='text file containing VMC gofr')
    parser.add_argument('DMC',type=str,help='text file containing VMC gofr')
    parser.add_argument("-t", "--threshold", type=float, default=1e-2
                        , help="threshold for small number, default 1e-2" )
    parser.add_argument("-cv", "--coarsenVMC", type=int, default=1
                        , help="by what factor to corasen VMC data" )
    parser.add_argument("-cd", "--coarsenDMC", type=int, default=1
                        , help="by what factor to corasen DMC data" )
    args = parser.parse_args()
    
    r,vmc= np.loadtxt(args.VMC).T
    r,dmc= np.loadtxt(args.DMC).T
    r,vmc,dmc=nonzero(r,vmc,dmc,args.threshold)
    
    # corasen if needed
    newrv,newvmc=rehist(r,vmc,args.coarsenVMC)
    newrd,newdmc=rehist(r,dmc,args.coarsenDMC)
    
    # spline fit to refine
    r=np.linspace(max(newrv[0],newrd[0]),min(newrv[-1],newrd[-1]),1000)
    vmc=interp1d(newrv,newvmc,kind='cubic')(r)
    dmc=interp1d(newrd,newdmc,kind='cubic')(r)
    fixed=[ dmc[i]/vmc[i]*dmc[i] for i in range(len(r)) ]
    print "Ro = ",sum(fixed*r)/sum(fixed)
    
    plt.plot(r,vmc,'--',label="VMC")
    plt.plot(r,dmc,'-.',label="DMC")
    plt.plot(r,fixed,label="Fix")
    plt.legend(loc=0)
    plt.show()

    
# end __main__
