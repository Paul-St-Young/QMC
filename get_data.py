#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import nexus_addon as na
import subprocess as sp

import xml.etree.ElementTree as ET
def twist_options(twist_input,orb_options = ['meshfactor','precision','twistnum']):
    tree = ET.parse( twist_input )
    dft_orb_builder_node = tree.findall(
        './/sposet_builder[@type="bspline"]')[0]
    
    myoptns = {}
    for optn in orb_options:
        myoptns[optn] = dft_orb_builder_node.attrib[optn]
    # end for
    return myoptns
# end def 

if __name__ == '__main__':
    # initialize analyzer
    from qmca import QBase
    options = {"equilibration":"auto"}
    QBase.options.transfer_from(options)

    paths = sp.check_output(['find','../default','-name','dmc_4']).split('\n')[:-1]

    for path in paths:
        # the purpose of this loop is to generate raw data base twists.json, one for each folder

        if os.path.isfile( os.path.join(path,'twists.json') ):
            continue
        # end if

	
	# locate all inputs in this folder
        cmd = 'ls %s | grep in.xml | grep twistnum'%path
        proc = sp.Popen(cmd,shell=True,stdout=sp.PIPE,stderr=sp.PIPE)
        out,err = proc.communicate()

        inputs = out.split('\n')[:-1]
	
	# make a database of all the scalar files
	data = []
	for qmc_input in inputs:
	    infile = os.path.join(path,qmc_input)
	    mydf = pd.DataFrame(na.scalars_from_input(infile))
	    toptn= twist_options(infile)
	    for optn in toptn.keys():
		mydf[optn] = toptn[optn]
	    # end for
	    data.append( mydf )
	# end for
	df = pd.concat(data).reset_index().drop('index',axis=1)
	
	# save raw data in local directory
	pd.concat([df,df['settings'].apply(pd.Series)],axis=1).to_json(
	    os.path.join(path,'twists.json'))

    # end for path

    for path in paths:
        # the purpose of this loop is to generate analyzed database 'scalars.json', one for each folder 

        if os.path.exists( os.path.join(path,'scalars.json') ):
            continue
        # end if

	# load local data
	df = pd.read_json(os.path.join(path,'twists.json'))
	
	# !!!! only analyze real twists
	#real_twists = [ 0,  8, 10, 32, 34, 40, 42,  2]
	#df = df[ df['twistnum'].apply(lambda x:x in real_twists) ]
	
	# exclude columns that don't need to be averaged, add more as needed
	special_colmns = ['iqmc','method','path','settings','vol_unit','volume']
	columns_to_average = df.drop(special_colmns,axis=1).columns

	mean_names = []
	error_names= []
	for col_name in columns_to_average:
	    if col_name.endswith('_mean'):
		mean_names.append(col_name)
	    elif col_name.endswith('_error'):
		error_names.append(col_name)
	    # end if
	# end for col_name

	col_names   = []
	for iname in range(len(mean_names)):
	    mname = mean_names[iname].replace('_mean','')
	    ename = error_names[iname].replace('_error','')
	    assert mname == ename
	    col_names.append(mname)
	# end for i
	
	# perform twist averaging
	new_means  = df.groupby('iqmc')[mean_names].apply(np.mean)
	ntwists = len(df[df['iqmc']==0]) # better way to determine ntwists?
	new_errors = df.groupby('iqmc')[error_names].apply(
	    lambda x:np.sqrt((x**2.).sum())/ntwists)
	
	# make averaged database
	dfev = pd.merge(new_means.reset_index(),new_errors.reset_index())
	extras = df[special_colmns].groupby('iqmc').apply(lambda x:x.iloc[0])
	newdf = pd.merge( extras.drop('iqmc',axis=1).reset_index(), dfev)
	
	newdf.to_json(os.path.join(path,'scalars.json'))
	
    # end for 

    # congregate data
    import dmc_databse_analyzer as dda
    data = []
    for path in paths:
        jfile = path + '/scalars.json'
        data.append( dda.process_dmc_data_frame( pd.read_json(jfile) ) )
    # end for path

    df = pd.concat(data).reset_index().drop('index',axis=1)
    df.to_json('scalars.json')

# end __main__
