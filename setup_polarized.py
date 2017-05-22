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
    if 'mass' in protons.attrib:
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

def edit_quantum_particleset(e_particleset,centers,rs,hname='p'):
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
    p_section.attrib["name"] = hname
    p_section.attrib['mass'] = str(unit.mp/unit.me)
    p_section.attrib['size'] = str(nions)
    p_section.xpath('.//parameter[@name="charge"]')[0].text =  '    1    '
    p_section.xpath('.//parameter[@name="mass"]')[0].text = '    %s    ' % str(unit.mp/unit.me)
    ppos_text = "\n"+"\n".join([" ".join( map(str,pos) ) for pos in proton_pos])+"\n"
    p_section.find("attrib").text = ppos_text

    e_particleset.append(p_section)
# end def

def edit_jastrows(wf,hname='p'):

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

        eHterm.attrib["speciesA"] = etype
        eHterm.attrib["speciesB"] = hname

        coeff = eHterm.find("coefficients")
        coeff.attrib["id"] = etype+hname 
        coeff.text = Uep # initialize e-p Jastrow to clamped values
        j2.append(eHterm)

    # end for ie

# end def edit_jastrows

def edit_backflows(wf):

    bf_node = wf.find('.//backflow')
    if bf_node is None:
        return # nothing to do

    # 1.  grab eta_ei and remove e-I backflow
    bf1  = wf.find('.//transformation[@type="e-I"]')
    if bf1 is None:
        return # nothing to do
    eta1 = bf1.xpath('.//correlation')
    bf_node.remove(bf1)

    # 2. edit 2-body backflow to add e-p correlation
    bf2  = wf.find('.//transformation[@type="e-e"]')
    etypes = {0:"u",1:"d"}
    for eta in eta1:
        ion_name = eta.attrib.pop('elementType')
        ion_name = 'p' # !!!! override ion name
        for ie in etypes.keys():
            ename = etypes[ie]
            eta.set('speciesA',ename)
            eta.set('speciesB',ion_name)
            coeff = eta.find('.//coefficients')
            coeff.set('id',ename+ion_name+'B')
            bf2.append(deepcopy(eta))
        # end for ie
    # end for eta

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
    """ construct <sposet_builder type="mo"> for protons and add one determinant to slater """
    nions = len(centers)

    # get electronic sposet_builder
    ebuilder = wf.find('sposet_builder[@source="ion0"]')
    if ebuilder is None:
        raise RuntimeError('electronic sposet_builder with source="ion0" not found.')
    # end if 

    # build electronic single-particle orbitals around "wf_centers" instead of "ion0"
    ebuilder.set('source','wf_centers')

    # start <sposet_builder> for protons
    pbuilder = etree.Element('sposet_builder',attrib={
	'type':'mo',            # use MolecularOrbitalBuilder
	'source':'wf_centers',  # use static lattice sites for proton w.f.
	'transform':'yes',      # use numerical radial function and NGOBuilder
	'name':'proton_builder'
    }) # !!!! transformOpt flag forces cuspCorr="yes" and key="NMO"

    # construct <basisset>
    pbasis = etree.Element("basisset")
    pao_basis = etree.Element('atomicBasisSet',attrib={
	'elementType':'H',     # identifier for aoBuilder
	'angular':'cartesian', # Use Gamess-style order of angular functions
	'type':'GTO',          # use Gaussians in NGOBuilder::addRadialOrbital()
	'normalized':'yes'     # do not mess with my input coefficients
    })

    # build <grid>
    bgrid = etree.Element('grid',{
	'npts':'1001',
	'rf':'100',
	'ri':'1.e-6',
	'type':'log'
    })

    # build <basisGroup>
    bgroup = etree.Element('basisGroup',{
	'l':'0',
	'n':'1',
	'rid':'R0'
    })
    bgroup.append(etree.Element('radfunc',{'contraction':'1.0','exponent':'9.0'}))

    pao_basis.append(bgrid)
    pao_basis.append(bgroup)
    pbasis.append(pao_basis)
    # finished construct </basisset>

    # build <sposet>
    psposet = etree.Element('sposet',attrib={
	'name':'spo_p',
	'id':'proton_orbs',
	'size':str(nions) # SPOSetBase::put
	# no 'cuspInfo'
    })
    pbuilder.append(pbasis)  # build basis set
    pbuilder.append(etree.Comment("Identity coefficient matrix by default SPOSetBase::setIdentity"))
    pbuilder.append(psposet) # build single-particle orbitals
    # end </sposet_builder>

    # build <determinant>
    pdet = etree.Element('determinant',{
	'group':'p',
	'id':'pdet',
	'size':str(nions),
	'sposet':psposet.get('name')
    })

    slater = wf.find('determinantset/slaterdeterminant')
    ebuilder_idx = wf.index(ebuilder)
    wf.insert(ebuilder_idx+1,pbuilder) # insert after electronic sposet_builder
    slater.append(pdet)

# end def edit_determinantset

#  Main Routine
# ======================================
def bo_to_nobo(bo_input_name,nobo_input_name,ion_width=10.0,rs=1.31):

    parser = etree.XMLParser(remove_blank_text=True)
    xml = etree.parse(bo_input_name,parser)

    ion0 = xml.xpath('//particleset[@name="ion0"]')[0]
    centers = change_ion0_to_wf_ceneters(ion0)

    e_particleset = xml.xpath('//particleset[@name="e"]')[0]
    edit_quantum_particleset(e_particleset,centers,rs)

    wf = xml.xpath("//wavefunction")[0]
    edit_jastrows(wf)
    edit_backflows(wf)
    edit_determinantset(wf,centers,ion_width)

    ham = xml.xpath("//hamiltonian")[0]
    edit_hamiltonian(ham)

    xml.write(nobo_input_name,pretty_print=True)
# end def bo_to_nobo

if __name__ == "__main__":
    import sys
    #prefix = sys.argv[1]
    #bo_to_nobo(prefix+"-dmc.in.xml",prefix+"-nobo.in.xml")
    inp_xml = sys.argv[1]
    out_xml = 'nobo-'+inp_xml
    bo_to_nobo(inp_xml,out_xml)
# end __main__
