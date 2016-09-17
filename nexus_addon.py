from fileio import TextFile
import numpy as np
from mmap import mmap

def read_xsf_datagrid(filename):
    
    lat_skip = 4 # first 4 lines after DATAGRID_3D_ORBITAL are lattice descriptions
    end_skip = 4 # 4 extra lines after data grid
    
    myfile = open(filename,"r+")
    mm = mmap(myfile.fileno(),0)
    idx = mm.find("DATAGRID_3D_ORBITAL")
    mm.seek(idx)
    mm.readline()

    grid_size = map(int,mm.readline().strip().strip("\n").split())
    for i in range(lat_skip):
        mm.readline()
    # end for
    
    data_grid_text = mm[mm.tell():mm.size()]
    data_grid = [map(float,line.strip().split()) for line in data_grid_text.split("\n")[:-end_skip]]
    data_grid = np.array(data_grid).astype(float)
    data_grid = data_grid.reshape(tuple(grid_size))
    return data_grid
# end def

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
        try:
            qa = QmcpackAnalyzer(qmc_input)
        except:
            if warn:
                raise AttributeError("failed at %s"%qmc_input)
            else:
                continue
            # end if
        # end try

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

def to_dict(iterable):
    d = dict()
    for k,v in iterable.iteritems():
        try:
            len(v.keys())
            d[k] = to_dict(v)
        except:
            d[k] = v
        #end if
    #end for
    return d
#end def to_dict

from qmca import DatAnalyzer
from lxml import etree

# specify equilibration length before use
#from qmca import QBase
#options = {"equilibration":"auto"}
#QBase.options.transfer_from(options)

def scalars_from_input(qmcinput,extract = ["mean","error"],skip_failed=False):
    """ Look for calculations specified in qmcinput xml, analyze each scalar.dat file and return a dictionary of data. Need to specify equilibration length with:
    from qmca import QBase
    options = {"equilibration":"auto"}
    QBase.options.transfer_from(options) """
    
    # current directory
    subdir = "/".join( qmcinput.split("/")[:-1] ) +"/"
    
    # parse input
    xml = etree.parse(qmcinput)
    
    # used for QMCPACK naming
    proj_id   = xml.xpath("//project")[0].get("id") 
    num_start = int( xml.xpath("//project")[0].get("series") )
    
    data = []
    calcs = xml.xpath("//qmc")
    for iqmc in range(len(calcs)):

        entry = {}

        # get input parameters
        method = calcs[iqmc].attrib["method"]
        params = calcs[iqmc].xpath(".//parameter")
        param_dict = {}
        for param in params:

            value = param.text
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass
                # end try
            # end try

            param_dict[ param.get("name") ] = value 

        # end for param
        entry["settings"] = param_dict

        # get scalar values
        scalar_file = subdir + ".".join([proj_id
            ,"s"+str(num_start+iqmc).zfill(3),"scalar","dat"]).strip("/")
        entry["path"] = scalar_file
        entry["iqmc"] = iqmc
        entry["method"] = method

        try:
            qa = DatAnalyzer(scalar_file,0)
        except:
            if skip_failed:
                #data.append(entry)
                continue
            # end if
            print "failed to read ",scalar_file, " did you initialize QBase? Read docstring"
            print " if there are known failed simulations, set skip_failed=True. "
        # end try 
        
        scalar_attribs = to_dict(qa.stats)
        for key in scalar_attribs:
            for name in extract:
                entry[key+"_"+name] = scalar_attribs[key][name]
            # end for
        # end for

        data.append(entry)
    # end for i
    return data
# end def 

from nexus import PwscfAnalyzer
def qe_scalars(qeinputs,warn=True):
    """ given a list of quantum espresso inputs, return a dictionary of simulation data.
     use warn=False if you know which simulations did not finish. """
    
    data = []
    for pwinput in qeinputs:
        pa = PwscfAnalyzer(pwinput)
        pa.analyze()

        failed = False
        try:
            pa.E
        except AttributeError:
            failed = True
            if warn:
                print pwinput," failed/did not finish"
            # end if
        # end try

        if failed:
            entry = {
                "path":pa.path,
                "system":pa.input.system.to_dict(),
                "kgrid":pa.input.k_points.grid,
                "failed":True
            }
        else:
            entry = {
                "path":pa.path,
                "energy":pa.E,
                "pressure":pa.pressure,
                "volume":pa.volume,
                "system":pa.input.system.to_dict(),
                "forces":pa.forces,
                "stress":pa.stress,
                "kgrid":pa.input.k_points.grid,
                "walltime":pa.walltime,
                "failed":False
            }
        # end if failed
        data.append(entry)
    # end for
    return data

# end for
