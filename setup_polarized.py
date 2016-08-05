#!/usr/bin/env python

class codata2016:
    def __init__(self):
        self.mp = 1.672621898e-27 # kg
        self.me = 9.10938356e-31  # kg
    # end def
# end class

# HELPER FUNCTIONS
# ======================================
from lxml import etree
def xml_print(node):
    print etree.tostring(node,pretty_print=True)
# end def

from copy import deepcopy
import numpy as np

def matrix_to_text(matrix):
    text = "\n"+\
      "\n".join([" ".join( map(str,pos) ) for pos in matrix])\
    +"\n"
    return text
# end def matrix_to_text

# START EDIT QMCPACK INPUT XML 
# ======================================
def change_ion0_to_wf_ceneters(ion0):

    # change ion0 particle set to be artificial centers around which
    #  to expand the wave function
    # ======================================

    ion0.attrib['name'] = 'wf_centers'
    protons = ion0.find("group")
    protons.attrib.pop("mass")
    protons.attrib.pop("size")
    #protons.attrib["name"] = "centers"

    for child in protons:
        if child.attrib["name"] != "position" and child.attrib["name"] != "charge":
            protons.remove(child)
        elif child.attrib["name"] == "charge":
            child.text = child.text.replace("1","0")
        # end if
    # end for

    # record and return the wf_centers for later use 
    # ======================================
    centers_text = ion0.find("group/attrib").text
    centers = [map(float,center.split()) for center in centers_text.split("\n")[1:-1]]
    nions = len(centers)
    ion0.find("group").attrib["size"] = str(nions)

    return centers

# end def change_ion0_to_wf_ceneters

def edit_quantum_particleset(e_particleset,centers,rs):
    # centers: wf_centers
    nions = len(centers)

    # initialize electrons manually: remove random, add attrib positions
    # ======================================

    # do not initialize eletrons at random
    # --------------------------------------
    e_particleset.attrib.pop("random")

    def gaussian_move(centers,sigma):
        move = np.random.randn( *np.shape(centers) )*sigma
        return np.copy(centers + move)
    # end def

    # sprinkle particles around ion positions
    # --------------------------------------
    electron_pos = gaussian_move(centers,0.2*rs)
    proton_pos = gaussian_move(centers,0.01*rs)

    # should be 2 groups in electronic wf, one u, one d
    spin_groups = e_particleset.findall("group")

    for i in range(2): # loop through u,d
        
        # 0 for up spin, 1 for down spin
        e_section = spin_groups[i]
        
        # assuming adjacent protons have adjacent indexing, the following will
        #  place u on the even sublattice, d on the odd sublattice
        epos = electron_pos[i::2]
        epos_text = "\n"+"\n".join([" ".join( map(str,pos) ) for pos in epos])+"\n"
        
        new_section = etree.Element("attrib",
                {"name":"position","datatype":"posArray","condition":"0"})
        new_section.text = epos_text
        e_section.append( new_section )
        
    # end for i

    # add protons to the particleset 
    # ======================================
    #  steal from electron section
    unit = codata2016()
    p_section = deepcopy( e_section )
    p_section.attrib["name"] = "H"
    p_section.attrib['mass'] = str(unit.mp/unit.me)
    p_section.attrib['size'] = str(nions)
    p_section.xpath('.//parameter[@name="charge"]')[0].text =  '    1    '
    p_section.xpath('.//parameter[@name="mass"]')[0].text = '    %s    ' % str(unit.mp/unit.me)
    ppos_text = "\n"+"\n".join([" ".join( map(str,pos) ) for pos in proton_pos])+"\n"
    p_section.find("attrib").text = ppos_text

    e_particleset.append(p_section)
# end def

def edit_jastrows(wf):

    # 1.  grab Uep and remove ep jastrow
    j1  = wf.xpath('//jastrow[@name="J1"]')[0]
    Uep = j1.xpath(".//coefficients")[0].text
    wf.remove(j1)

    # 2. edit 2-body jastrow to add e-p correlation
    j2   = wf.xpath('//jastrow[@name="J2"]')[0]
    term = j2.find("correlation")

    etypes = {0:"u",1:"d"}
    for ie in range(2):
         
        eHterm = deepcopy(term)

        etype = etypes[ie] # 0:u, 1:d
        hname = "H"

        eHterm.attrib["speciesA"] = etype
        eHterm.attrib["speciesB"] = hname

        coeff = eHterm.find("coefficients")
        coeff.attrib["id"] = etype+hname 
        coeff.text = Uep # initialize e-p Jastrow to clamped values
        j2.append(eHterm)

    # end for ie

# end def edit_jastrows

def edit_hamiltonian(ham):
    # remove all interactions that requires the "ion0" particleset
    for interaction in ham: 
        use_ion0 = False
        for key in interaction.attrib.keys():
            if interaction.attrib[key]=="ion0":
                use_ion0 = True
            # end if
        # end for key
        if use_ion0:
            ham.remove(interaction)
        # end if
    # end for
# end def edit_hamiltonian

def edit_determinantset(wf,centers,ion_width):
    nions = len(centers)
    # build single-particle orbitals around "wf_centers" instead of "ion0"
    wf.find("sposet_builder").attrib['source'] = 'wf_centers'

    slater = wf.find("determinantset/slaterdeterminant")

    # add basis to first determinant
    basis = etree.Element("basisset",{"ref":"wf_centers"})

    atomic_basis = etree.Element("atomicBasisSet",{
        "name":"G2","angular":"cartesian","type":"GTO",
        "elementType":"H","normalized":"no"
    })
    atomic_basis.append(etree.Element("grid",{
        "type":"log","ri":"1.e-6","rf":"1.e2","npts":"1001"
    }))
    basis_group = etree.Element("basisGroup",{
        "rid":"R0","n":"1","l":"0","type":"Gaussian"
    })
    basis_group.append(etree.Element("radfunc",{
        "exponent":str(ion_width),"contraction":"1"}))

    atomic_basis.append(basis_group)
    basis.append(atomic_basis)

    # add a determinant for protons
    Hdet = etree.Element("determinant",
            {"id":"Hdet","group":"H","spin":"0","size":str(nions)
             ,"type":"mo","source":"wf_centers"})
    coeff = etree.Element("coefficient",
            {"id":"HdetC","type":"constArray","size":str(nions)})
    coeff.text = matrix_to_text(np.eye(nions))
    Hdet.append(basis)
    Hdet.append(coeff)
    slater.append(Hdet)
# end def edit_determinantset

#  Main Routine
# ======================================
def bo_to_nobo(bo_input_name,nobo_input_name,ion_width=10.0,rs=1.31):

    xml = etree.parse(bo_input_name)

    ion0 = xml.xpath('//particleset[@name="ion0"]')[0]
    centers = change_ion0_to_wf_ceneters(ion0)

    e_particleset = xml.xpath('//particleset[@name="e"]')[0]
    edit_quantum_particleset(e_particleset,centers,rs)

    wf = xml.xpath("//wavefunction")[0]
    edit_jastrows(wf)
    edit_determinantset(wf,centers,ion_width)

    ham = xml.xpath("//hamiltonian")[0]
    edit_hamiltonian(ham)

    xml.write(nobo_input_name,pretty_print=True)
# end def bo_to_nobo

if __name__ == "__main__":
    import sys
    prefix = sys.argv[1]
    bo_to_nobo(prefix+"-dmc.in.xml",prefix+"-nobo.in.xml")
# end __main__
