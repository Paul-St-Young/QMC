#!/usr/bin/env python

import h5py
import numpy as np

def grabDensity(h5file):
    # get density
    f = h5py.File(h5file)
    for name,quantity in f.items():
        if name.startswith('density'):
            density = quantity.get("value")[:]
        # end if name.startswith
    # end for name,quantity
    f.close()
    
    return density
# end def grabDensity

import xml.etree.ElementTree as ET
def grabLimits(inputfile):
    tree = ET.parse(inputfile)
    root = tree.getroot()
    for section in root:
        if section.tag=='hamiltonian':
            for ham in section:
                if ham.attrib['name']=='density':
                    xmin=ham.attrib['x_min']
                    xmax=ham.attrib['x_max']
                    ymin=ham.attrib['y_min']
                    ymax=ham.attrib['y_max']
                    zmin=ham.attrib['z_min']
                    zmax=ham.attrib['z_max']
                    delta=ham.attrib['delta']
                # end if ham==density
            # end for ham in section
        # end if section.tag
    # end for section in root
    return [xmin,xmax,ymin,ymax,zmin,zmax,delta]
# end def grabInput

import argparse
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Plot proton density')
    parser.add_argument('XML', type=str, default=None, help="input XML")
    parser.add_argument("DMC", type=str, help="h5 file with DMC density")
    parser.add_argument('-e','--equil', type=int, help="number of equilibration steps")
    args = parser.parse_args()

    # get density grid parameters
    limits = grabLimits(args.XML)
    xmin,xmax,ymin,ymax,zmin,zmax = map(float,limits[:-1])
    d1,d2,d3 = map(float,limits[-1].split())
    dx = (xmax-xmin)/int(1/d1)
    dy = (ymax-ymin)/int(1/d2)
    dz = (zmax-zmin)/int(1/d3)

    # get density on grid
    density = grabDensity(args.DMC)[args.equil:]
    avgdens = density.mean(axis=0)
    #avgdens = avgdens.transpose()

    print xmin,xmax,ymin,ymax,zmin,zmax 
    print (xmax-xmin)/dx,(ymax-ymin)/dy,(zmax-zmin)/dz
        
    print 'dumpping to file'
    f = open('density.dat','w')
    for i in range(len(avgdens)):
        x = i*dx+xmin+dx/2.0
        for j in range(len(avgdens[0])):
            y = j*dy%(ymax-ymin)+ymin+dy/2.0
            for k in range(len(avgdens[0][0])):
                z = k*dz%(zmax-zmin)+zmin+dz/2.0
                f.write( "%2.3f %2.3f %2.3f %1.5f\n" % (x,y,z,avgdens[i][j][k]) )
            # end for k
        # end for j
    # end for i
    f.close()

# end __main__
