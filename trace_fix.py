#!/usr/bin/env python
# Plot or dump proton density for HCN+

import h5py
import numpy as np

from gofr import grabGR,normalize2PDF
from fix  import nonzero

import argparse
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Plot time evolution of g(r)')
    parser.add_argument('VMC',type=str,help='h5 file containing VMC gofr')
    parser.add_argument('DMC',type=str,help='h5 file containing DMC gofr')
    parser.add_argument("-vd", "--VMCData", type=int, default=8000, 
        help="how much VMC data to keep")
    parser.add_argument("-w", "--windowSize", type=int, default=100, 
        help="size of time slice window in steps")
    parser.add_argument("-r","--roll", type=int, default=100, 
        help="step size for rolling average")
    parser.add_argument("-e","--equil", type=int, default=100, 
        help="number of equilibration blocks")
    parser.add_argument("-t","--threshold", type=float, default=1e-2, 
        help="small number cutoff used in fix")
    args = parser.parse_args()

    mol = args.VMC.split(".")[0] # name of molecule
    r,GR = grabGR(args.DMC,'0','2')
    r,GRV= grabGR(args.VMC,'0','2')
    grv=GRV[len(GRV)-args.VMCData:].mean(axis=0)
    normalize2PDF(r,grv)
    # ---- grv now contain probability distributions
    
    idx=[]
    peak_trace=[]
    for i in range(0,(len(GR)-args.windowSize)/args.roll):
        endIdx   = len(GR)-i*args.roll
        
        startIdx = endIdx-args.windowSize
        gr = GR[startIdx:endIdx].mean(axis=0)
        normalize2PDF(r,gr)
        # ---- gr now contain probability distributions
        rf,vmc,dmc=nonzero(r,grv,gr,args.threshold)
        fixed=[ dmc[i]/vmc[i]*dmc[i] for i in range(len(rf)) ]
        idx.append(endIdx)
        peak_trace.append(sum(fixed*rf)/sum(fixed))
        vm = sum(vmc*rf)/sum(vmc)
        dm = sum(dmc*rf)/sum(dmc)
        #peak_trace.append(dm*dm/vm)
    # end for i
    peak_trace.reverse()
    idx.reverse()
    with open('grt.dat','w') as f:
        for i in range(len(peak_trace)):
            f.write( str(idx[i])+" "+str(peak_trace[i])+"\n"  )
        # end for i
    # end with open
    print "Using block", args.equil, "and beyond",np.array(peak_trace).mean()

    
# end __main__
