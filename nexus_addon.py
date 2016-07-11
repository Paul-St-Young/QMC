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

    try:
        qmc_run.scalars
    except:
        return entry
    # end try

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
from nexus import QmcpackAnalyzer
def qmcpack_scalars(qmc_inputs,img_name = "image.p",save_image=True
    # default inputs to extract_scalar
    ,extract_names   = ["LocalEnergy","LocalEnergyVariance"],extract_all=False
    ,extract_attribs = ["mean","error"],warn=True):
    
    """ build on extract_scalar, pull out useful scalar outputs from a list of QMCPACK inputs
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

        # extract some input data
        qi = qa.info.input
        calculations = []
        for key in qi.simulation.keys():
            if key=="calculations":
                calculations = qi.simulation[key].to_dict()
            elif key=="qmc":
                calculations = qi.simulation[key].to_dict()
            # end if
        # end for
        
        # find the folder containing the input file
        path = "/".join( qmc_input.split("/")[:-1] ) + "/"
        input_name = qmc_input.split("/")[-1]
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
            if calculations != []:
                entry["settings"] = calculations[iqmc]
            # end if
            entry["path"] = path
            entry["input_name"] = input_name
            data.append(entry)
        # end for iqmc
    # end for qmc_input
    
    return data
# end def 

from nexus import GamessAnalyzer

def collect_gamess_from_inputs(inputs):
    data = []
    for gamess_input in inputs:
        
        # analyze run
        ga = GamessAnalyzer(gamess_input)
        ga.analyze()
        
        # initialize entry
        entry = {"path":gamess_input}
        
        # record inputs
        input_essentials = set(["system","contrl","basis","det","data"])
        for key in ga.info.input.keys():
            if key in input_essentials:
                entry.update( ga.info.input[key].to_dict() )
            # end if
        # end for
        
        # record energies
        entry.update( ga.energy.to_dict() )
        
        # save entry
        data.append(entry)
    # end for
    return data
# end def collect_gamess_from_inputs

from nexus import QmcpackInput
import numpy as np
def get_jastrow_list(qmcpack_input,one_body=True,two_body=True,just_wf=False):
    
    def two_body_jastrow(correlation):
        x = np.linspace(0,correlation.rcut,correlation.size)
        y = correlation.coefficients.coeff
        return x,y
    # end def
    
    qi = QmcpackInput(qmcpack_input)
    if just_wf:
        # if qmcpack_input is an output of optimization (prefix.opt.xml)
        jastrows = qi.qmcsystem.wavefunction.jastrows
    else:
        jastrows = qi.simulation.qmcsystem.wavefunction.jastrows
    # end if
    
    jastrow_list = []
    if one_body:
        j1_keys = jastrows.J1.keys()
        if "correlations" in j1_keys:
            for name in jastrows.J1.correlations.keys():
                qij1 = jastrows.J1.correlations[name]
                x,y  = two_body_jastrow( qij1 )
                entry = {"name":name,"x":x,"y":y
                         ,"type":"one-body"
                         ,"path":qmcpack_input}
                jastrow_list.append(entry)
            # end for name
        else:
            qij1 = jastrows.J1.correlation
            x,y  = two_body_jastrow( qij1 )
            entry = {"name":qij1.elementtype,"x":x,"y":y
                     ,"type":"one-body"
                     ,"path":qmcpack_input}
            jastrow_list.append(entry)
        # end if
    # end if

    if two_body:
        for name in jastrows.J2.correlations.keys():
            qij2 = jastrows.J2.correlations[name]
            x,y  = two_body_jastrow( qij2 )
            entry = {"name":name,"x":x,"y":y
                     ,"type":"two-body"
                     ,"path":qmcpack_input}
            jastrow_list.append(entry)
        # end for name
    # end if
    
    return jastrow_list
# end def 
