#!/usr/bin/env python

import h5py
import pandas
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import argparse

slices = [0,1]
species_map={'0':'C','1':'N','2':'H'}
#species_map={'0':'Li','1':'H'}

def timeSlice(gofr, windowSize):
    nbins = int(gofr.get("cutoff")[0]/gofr.get("delta")[0]) # number of bins
    
    # get bin centers
    r=np.arange( nbins ) * gofr.get("delta")[0]
    for i in range(len(r)):
        r[i] += dr/2.0
    # end i
    
    # get all time slices of g(r)
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

def average(r,gr):
    return sum(r*gr)/sum(gr)
# end def average

def gaussian(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
# end def gaussian

def fit_gaussian(r,gr,start,end):
    p0 = [1., average(r,ar[k]) , 1.]
    coeff, var_matrix = curve_fit(gaussian, r[start:end], gr[start:end], p0=p0)
    return coeff
# end def fit_gaussian

def get_start_and_end(gr,thres):
    started = False
    start = 0
    end = len(gr)-1
    for i in range(len(gr)):
        if (not started) and gr[i]>thres:
            started = True
            start = i
        # end if not started
        if started and gr[i]<thres:
            end = i
            break
        #end if started
    # end for
    return start,end
# end def

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Plot time evolution of g(r)')
    parser.add_argument("h5file", type=str, help="file name in the format ${molecule}.s00x.stat.h5")
    parser.add_argument("-w", "--windowSize", type=int, default=50, help="size of time slice window in steps" )
    args = parser.parse_args()

    start = 0 
    end = 0

    # get name of molecule
    mol = args.h5file.split(".")[0] # name of molecule

    f = h5py.File(args.h5file)
    for name,quantity in f.items():
        if name.startswith('gofr'):
            g,p,i,j = name.split('_')
            dr = quantity.get("delta")[0]
            if (i!=j and i=='0' and j=='2'):
            #if (i!=j):
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
                    if slices[k]==-1:
                        slices[k]=nwindow-1
                    # end if
                    start,end = get_start_and_end(ar[k],1e-2)
                    coeff = fit_gaussian(r,ar[k],start,end)
                    A,mean,sigma = coeff
                    print species_map[i]+"-"+species_map[j] + " window " + str(slices[k]+1) + " mu=" + str( round(mean,3) ) + " var=" + str( round(sigma*sigma,3) )
                    #print sum(ar[k]*r)/sum(ar[k])
                    
                # end for k
                
                df = pandas.DataFrame(ar.T,index=r,columns=label)
                ax = df.plot(title=mol+" run split into "+str(nwindow)+" windows each of size "+str(windowSize))
                ax.set_xlabel(species_map[i]+"-"+species_map[j]+" (a.u.)")
                ax.set_ylabel("P(r)")

                #plt.plot(r[start:end], gaussian(r[start:end],*coeff), label='Fitted data')

            # end if i!=j
            
        # end if name.startswith
    # end for name,quantity
    
    plt.show() 
# end __main__
