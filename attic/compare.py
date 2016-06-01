#!/usr/bin/env python

import math,numpy
import argparse
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='compare 2 data points')
    parser.add_argument("m1", type=float, help="mean1")
    parser.add_argument("e1", type=float, help="error1")
    parser.add_argument("m2", type=float, help="mean2")
    parser.add_argument("e2", type=float, help="error2")
    args = parser.parse_args()
    
    diff = numpy.abs(args.m1-args.m2)
    e = math.sqrt(args.e1*args.e1+args.e2*args.e2)
    x = diff/e
    
    print "data points are",x, "standard error away"
    print "1-erf(x/sqrt(2))=",1-math.erf(x/math.sqrt(2))
    
# end __main__
