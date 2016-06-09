from fileio import TextFile

# extend MyTextFile to recognize block pattern
class MyTextFile(TextFile):
    def block_between(self,header,trailer):
        
        header_pos = self.mm.find(header)
        trailer_pos = self.mm.find(trailer)
        
        block = self.mm[header_pos:trailer_pos]
        if header_pos >= trailer_pos:
            print "header {:d} appears after trailer {:d}".format(header_pos,trailer_pos)
        # end if
        
        trim_head_tail = "\n".join( block.split("\n")[1:-2] )
        return trim_head_tail.strip("\n")
    # end def 
# end class

def extract_scalar( qmc_run
    ,extract_names   = ["LocalEnergy","LocalEnergyVariance"], extract_all=False
    ,extract_attribs = ["mean","error"], warn=True ):

    """ Given a qmc run, extract desired scalar quantities. Return a dictionary of the extracted results. Available attributes can be found by printing qmc_run.scalars.keys(). """

    entry = {}

    if extract_all:
        for name in qmc_run.scalars.keys():
            if name not in entry.keys():
                scalar = qmc_run.scalars[name]
            # end if 
            for attrib in scalar.keys():
                if attrib in extract_attribs:
                    entry[name+"_"+attrib] = scalar[attrib]
                # end if
            # end for
        # end for
        for name in qmc_run.scalars_hdf.keys():
            if name not in entry.keys():
                scalar = qmc_run.scalars_hdf[name]
            # end if 
            for attrib in scalar.keys():
                if attrib in extract_attribs:
                    entry[name+"_"+attrib] = scalar[attrib]
                # end if
            # end for
        # end for name
    else:

        for name in extract_names:
            scalar = {}
            if name in qmc_run.scalars.keys():
                scalar = qmc_run.scalars[name]
            elif name in qmc_run.scalars_hdf.keys():
                scalar = qmc_run.scalars_hdf[name]
            else:
                if warn: print( "%s not found for %s!" % (name,qmc_run.info.file_prefix) )
            # end if
            for attrib in scalar.keys():
                if attrib in extract_attribs:
                    entry[name+"_"+attrib] = scalar[attrib]
                # end if
            # end for
        # end for

    # end if
    
    return entry
# end def

def scalar_rename(name):
    my_names = {
        "LocalEnergyVariance_error":"Ve",
        "LocalEnergyVariance_mean":"V",
        "LocalEnergy_mean":"E",
        "LocalEnergy_error":"Ee"
    }
    
    
    if name in my_names:
        new_name = my_names[name] 
    else:
        new_name = name
    # end if
    
    return new_name
# end def

import os
import pandas as pd
from nexus import QmcpackAnalyzer
def qmcpack_scalars(qmc_inputs,img_name = "image.p",save_image=True
    # default inputs to extract_scalar
    ,extract_names   = ["LocalEnergy","LocalEnergyVariance"],extract_all=False
    ,extract_attribs = ["mean","error"],warn=True):
    
    """ build on extract_scalar, pull out useful scalar outputs from QMCPACK run folder
    qmc_inputs: a list of locations of QMCPACK inputs
    return: a database containing basic scalar data
    
    save an image of the analyzer in the raw data directory to save time if a second pass is needed.
    useful = ["ElecElec","ElecIon","IonIon","KEcorr","Kinetic","LocalEnergy"
        ,"LocalEnergy_mpc_kc","LocalEnergy_sq","LocalPotential","MPC"
        ,"Pressure","AcceptRatio","BlockCPU","BlockWeight"]
    """
    
    data = []
    for qmc_input in qmc_inputs:

        # initialize an analyzer
        qa = QmcpackAnalyzer(qmc_input)
        
        # find the folder containing the input file
        path = "/".join( qmc_input.split("/")[:-1] )
        input_name = "/".join( qmc_input.split("/")[-1] )
        if os.path.isfile(path+"/"+img_name):
            qa.load(path+"/"+img_name)
        else:
            qa.analyze()
            if save_image:
                qa.save(path+"/"+img_name)
            # end if save_image
        # end if

        # loop through all runs in the <qmc="blah"/> sections
        for iqmc in qa.qmc.keys():
            entry = extract_scalar(qa.qmc[iqmc], extract_names=extract_names
                    ,extract_all=extract_all, warn=False)
            entry["iqmc"] = iqmc
            entry["path"] = path
            entry["input_name"] = input_name
            data.append(entry)
        # end for iqmc
    # end for qmc_input
    df = pd.DataFrame.from_dict(data)
    
    return df
# end def 
