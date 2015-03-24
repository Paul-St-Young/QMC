#!/usr/bin/python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas
import argparse

from scipy.optimize import curve_fit

def gofrGrabber(h5file,myi,myj):
    # get gofr list with label "gofr_ion0_myi_myj"
    dr = 0.0
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
    
    return r,dr,GR
# end def gofrGrabber

def normalize2PDF(r,dr,gr):
    A = 0.0
    for i in range(len(r)):
        gr[i] *= 4*np.pi*r[i]*r[i]*dr
        A += gr[i]*dr
    # end for
    for i in range(len(r)):
        gr[i] /= A
    # end for
# end def 

def average(r,gr):
    return sum(r*gr)/sum(gr)
# end def average

def gaussian(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
# end def gaussian
def fit_gaussian(r,gr,start,end):
    p0 = [1., average(r,gr) , 1.]
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

def mix(vgr,dgr,eps):
    vstart, vend = get_start_and_end(vgr,eps)
    dstart, dend = get_start_and_end(dgr,eps)
    start = max(vstart,dstart)
    end = min(vend,dend)
    
    lnfix = list( np.zeros( min(vstart,start) ) )
    if start > vstart:
        lnfix += list( np.zeros(start-vstart) )
    # end if
    
    for i in range(start,end):
        lnfix.append(dgr[i]/vgr[i])
    # end i
    
    lnfix += list( np.zeros(len(vgr)-end) )
    
    return lnfix  
# end def mix

def quick_fix(vgr,dgr,eps):
    fixed_gr = []
    for i in range(len(vgr)):
        if vgr[i]>eps:
            fixed_gr.append(dgr[i]*dgr[i]/vgr[i])
        else:
            fixed_gr.append(dgr[i])
        # end if
    # end for
    return fixed_gr
# end def quick_fix

def cut(r,f,r1,r2): # cut out a piece of f and pad with 0
    cutgr = []
    for i in range(len(r)):
        if r[i]>r1 and r[i]<r2:
            cutgr.append(f[i])
        else:
            cutgr.append(0)
        # end if
    # end for
    return cutgr
# end def cut

def perform_quick_fix(vmcfile,dmcfile,vw=4000,dw=500):
    r,dr,VGR = gofrGrabber(vmcfile,'0','2')
    vgr = VGR[len(VGR)-vw:].mean(axis=0)
    normalize2PDF(r,dr,vgr)
    r,dr,DGR = gofrGrabber(dmcfile,'0','2')
    dgr = DGR[len(DGR)-dw:].mean(axis=0)
    normalize2PDF(r,dr,dgr)
    # ---- vgr,dgr now contain probability distributions

    mean_vgr = []
    for i in range(0,len(vgr)-5,6):
        mean_vgr.append( sum(vgr[i:i+6])/6. )
    # end for i
    vgr = mean_vgr[:]

    fixed_gr = quick_fix(vgr,dgr,1e-2)
    normalize2PDF(r,dr,fixed_gr)
    return r,fixed_gr
# end def perform_quick_fix

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Use g(r) to fix wf')
    parser.add_argument('VMC',type=str,help='h5 file containing VMC gofr')
    parser.add_argument('DMC',type=str,help='h5 file containing VMC gofr')
    parser.add_argument("-vw", "--VMCwindowSize", type=int, default=4000
                        , help="size of time slice window in steps" )
    parser.add_argument("-dw", "--DMCwindowSize", type=int, default=500
                        , help="size of time slice window in steps" )
    args = parser.parse_args()

    '''
    start,end = get_start_and_end(vgr,1e-3)
    coeff = fit_gaussian(r,vgr,start,end)
    print "VMC fit ", coeff[1]
    start,end = get_start_and_end(dgr,1e-2)
    coeff = fit_gaussian(r,dgr,start,end)
    print "DMC fit ", coeff[1]
    '''

    # quick fixed
    r,fixed_gr = perform_quick_fix(args.VMC,args.DMC,args.VMCwindowSize,args.DMCwindowSize)
    start,end = get_start_and_end(fixed_gr,1e-2)
    coeff = fit_gaussian(r,fixed_gr,start,end)
    print "fix fit ", coeff[1]
    print "fix int ", sum(fixed_gr*r)/sum(fixed_gr)
   
    ar = np.array( fixed_gr )
    df = pandas.DataFrame(ar.T,index=r)
    ax = df.plot(title=r"Probability Distribution of $r_{CH}$")
    ax.set_xlabel(r"$r_{CH}$ (bohr)",fontsize=16)
    ax.set_ylabel("P(r)",fontsize=16)
    plt.plot(r, gaussian(r,*coeff), label='Fitted data')
    
    plt.show()
# end __main__
