#!/usr/bin/env python
from mmap import mmap
def get_madelung(qmcout,keyword='Madelung',val_type=float,val_loc=-1):
    """ find the value of the line 'keyword = value' in qmcout """
    
    # open file
    fhandle = open(qmcout,'r+')
    mm = mmap(fhandle.fileno(),0)
    
    # find first line with the keyword
    idx = mm.find(keyword)
    if idx == -1:
        raise IOError(keyword + ' not found')
    # end 
    
    # go to line and read
    mm.seek(idx)
    line = mm.readline().strip('\n')
    val_text = line.split('=')[val_loc]
    
    # try to obtain the value of keyword
    val = val_type(val_text) # rely on internal exception
    
    # close file
    fhandle.close()
    
    return val
# end def get_madelung
