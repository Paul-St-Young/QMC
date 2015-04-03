#!/usr/bin/env python
# Plot or dump proton density for HCN+

import h5py
import numpy as np

def grabGR(h5file,myi,myj):
    # get gofr list with label "gofr_ion0_myi_myj"
    r = []
    GR = []
    
    f = h5py.File(h5file)
    for name,quantity in f.items():
        if name.startswith('gofr'):
            g,p,i,j = name.split('_')
            if i==str(myi) and j==str(myj):
                dr = quantity.get("delta")[0]
                maxr = quantity.get("cutoff")[0]
                r = np.arange( int(maxr/dr) ) * dr
                for k in range(len(r)):
                    r[k] += dr/2.0
                # end for k
                GR = quantity.get("value")[:]
        # end if name.startswith
    # end for name,quantity
    f.close()
    
    return r,GR
# end def grabGR
def normalize2PDF(r,gr):
    A = 0.0
    dr=r[1]-r[0]
    for i in range(len(r)):
        gr[i] *= 4*np.pi*r[i]*r[i]*dr
        A += gr[i]*dr
    # end for
    for i in range(len(r)):
        gr[i] /= A
    # end for
# end def

import matplotlib.pyplot as plt
import pandas
import argparse
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Plot time evolution of g(r)')
    parser.add_argument("h5file", type=str, help="file name in the format ${molecule}.s00x.stat.h5")
    parser.add_argument("-w", "--windowSize", type=int, default=1000, help="size of time slice window in steps" )
    parser.add_argument("-d","--dump", action='store_true', help="dump last time")
    args = parser.parse_args()

    mol = args.h5file.split(".")[0] # name of molecule
    r,GR = grabGR(args.h5file,'0','2')
    
    nslices=len(GR)/args.windowSize
    slices = []
    sGR = []
    for i in range(nslices,0,-1):
        sGR.append( GR[(i-1)*args.windowSize:i*args.windowSize].mean(axis=0) )
        slices.append(i)
    # end for i
    sGR.reverse()
    slices.reverse()
    
    for i in range(len(sGR)):
        normalize2PDF(r,sGR[i])
    # end for i
    # ---- sGR now contain probability distributions
    if args.dump:
        with open('gr.dat','w') as f:
            for i in range(len(r)):
                f.write( str(r[i])+" "+str(sGR[-1][i])+"\n" )
            # end for i
        # end with open
    else:
        ar = np.array(sGR)
        df = pandas.DataFrame(ar.T,index=r,columns=slices)
        df.plot(title=mol+" run split (from end) into "+str(nslices)+" windows each of size "+str(args.windowSize))
        plt.show()
    # end if dump
    
# end __main__
