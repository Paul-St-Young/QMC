#!/usr/bin/env python

import numpy as np
def combine(m1,e1,m2,e2):
    m=(e2*e2*m1+e1*e1*m2)/(e1*e1+e2*e2)
    e=e1*e2/np.sqrt(e1*e1+e2*e2)
    return (m,e)
# end def combine

import argparse
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='combine 2 data points')
    parser.add_argument("m1", type=float, help="mean1")
    parser.add_argument("e1", type=float, help="error1")
    parser.add_argument("m2", type=float, help="mean2")
    parser.add_argument("e2", type=float, help="error2")
    args = parser.parse_args()

    print combine(args.m1,args.e1,args.m2,args.e2)

# end __main__
