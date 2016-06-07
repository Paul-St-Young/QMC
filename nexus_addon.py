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
    ,extract_names   = ["LocalEnergy","LocalEnergyVariance"]
    ,extract_attribs = ["mean","error"] ):

    """ Given a qmc run, extract desired scalar quantities. Return a dictionary of the extracted results. Available attributes can be found by printing qmc_run.scalars.keys(). """

    entry = {}

    for name in qmc_run.scalars.keys():
        if name in extract_names:
            scalar = qmc_run.scalars[name]
            for attrib in scalar.keys():
                if attrib in extract_attribs:
                    entry[name+"_"+attrib] = scalar[attrib]
                # end if
            # end for
        # end if
    # end for
    
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
    ,extract_names   = ["LocalEnergy","LocalEnergyVariance"]
    ,extract_attribs = ["mean","error"]):
    
    """ pull out useful scalar outputs from QMCPACK runs.
    qmc_inputs: a list of locations of QMCPACK inputs
    return: a database containing basic scalar data
    
    save an image of the analyzer in the raw data directory to save time if a second pass is needed.
    """
    
    data = []
    for qmc_input in qmc_inputs:

        # initialize an analyzer
        qa = QmcpackAnalyzer(qmc_input)
        
        # find the folder containing the input file
        path = "/".join( qmc_input.split("/")[:-1] )
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
            entry = extract_scalar(qa.qmc[iqmc])
            entry["iqmc"] = iqmc
            entry["path"] = path
            data.append(entry)
        # end for iqmc
    # end for qmc_input
    df = pd.DataFrame.from_dict(data)
    
    return df
# end def 
