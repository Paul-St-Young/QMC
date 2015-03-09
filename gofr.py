#!/usr/bin/env python

import h5py
import pandas
import matplotlib.pyplot as plt
import numpy as np
import argparse

slices = [0,-1]
species_map={0:'C',1:'N',2:'H'}

def timeSlice(gofr, windowSize):
    nbins = int(gofr.get("cutoff")[0]/gofr.get("delta")[0]) # number of bins
    r=np.arange( nbins ) * gofr.get("delta")[0]
    GR=[]
    for nwindow in range( 0,len(gofr.get("value"))/windowSize ):
    
       GR.append( gofr.get("value")[nwindow*windowSize:(nwindow+1)*windowSize].mean(axis=0) )

    # end for nwindow
    return r,GR
# end def timeSlice

def normalize(r,dr,gr):
    for i in range(len(r)):
        gr[i] *= 4*np.pi*r[i]*r[i]*dr
    # end for
# end def normalize

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Plot time evolution of g(r)')
    parser.add_argument("h5file", type=str, help="file name in the format ${molecule}.s00x.stat.h5")
    parser.add_argument("-w", "--windowSize", type=int, default=50, help="size of time slice window in steps" )
    args = parser.parse_args()

    # get name of molecule
    mol = args.h5file.split(".")[0] # name of molecule
    #species_map={}
    #for i in range(len(mol.strip("+").strip("-"))):
        #atom = mol[i]
        #species_map[i]=atom
    # end i

    f = h5py.File(args.h5file)
    for name,quantity in f.items():
        if name.startswith('gofr'):
            g,p,i,j = name.split('_')
            i = int(i)
            j = int(j)
            dr = quantity.get("delta")[0]
            
            if (i!=j):
                windowSize = args.windowSize
                nwindow = len(quantity.get("value"))/windowSize
                r,GR = timeSlice(f.get(name), windowSize)
                
                ar = []
                for s in slices:
                    normalize(r,dr,GR[s])
                    ar.append(GR[s])
                # end for
                ar = np.array(ar)
                
                label = slices[:]
                for k in range(len(slices)):
                    label[k] = "window "+str(slices[k]+1)
                    print species_map[i]+"-"+species_map[j] + " window " + str(slices[k]+1) + " " + str( (ar[k]*r).mean()/ar[k].mean() )
                # end for k
                
                df = pandas.DataFrame(ar.T,index=r,columns=label)
                ax = df.plot(title=mol+" run split into "+str(nwindow)+" windows each of size "+str(windowSize))
                ax.set_xlabel(species_map[i]+"-"+species_map[j]+" (a.u.)")
                ax.set_ylabel("g(r)")
            # end if i!=j
        # end if name.startswith
    # end for name,quantity
    
    
    plt.show()
    
    
# end __main__
