#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit

def linear(x,a,b):
    return a+b*x
# end def linear

import matplotlib.pyplot as plt
import argparse
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument("scalar", type=str, help="scalar file name")
    parser.add_argument('-p','--plot', action='store_true', help="plot extrapolation")
    #parser.add_argument('-c','--column', type=int, default=0, help="column to analyze")
    args = parser.parse_args()
    data = np.loadtxt(args.scalar).T
    T=data[0];E=data[1];Ee=data[2]
    
    popt,pcov = curve_fit(linear,T,E,sigma=Ee,absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    Efit = popt[0]
    Efite= perr[0]
    print Efit,Efite
    
    if args.plot:
        x = np.linspace(0,T[0])
        
        fig, ax = plt.subplots()
        fig.tight_layout()
        plt.subplots_adjust(left=0.15)
        ax.set_xlim(-T[-1],T[0]*1.5)
        
        ax.plot(T,E,'x',ms=10,color='black')
        ax.errorbar(T,E,yerr=Ee,fmt=None,ecolor='black')
        
        ax.plot(x,linear(x,*popt),'--')
        ax.errorbar(0,Efit,fmt='x',color='r',ms=10,yerr=Efite,ecolor='r')
        
        ax.set_ylabel("E(DMC) (Ha)",fontsize=16)
        plt.show()
    # end if plot
    
# end __main__
