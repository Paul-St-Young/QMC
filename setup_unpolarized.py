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

def find_dimers(axes,pos,sep_min=1.0,sep_max=1.4):
  natom,ndim = pos.shape
  import pymatgen as mg
  elem = ['H'] * natom
  struct = mg.Structure(axes,elem,pos,coords_are_cartesian=True)
  dtable = struct.distance_matrix

  xidx,yidx = np.where( (sep_min<dtable) & (dtable<sep_max) )
  pairs = set()
  for myx,myy in zip(xidx,yidx):
    mypair1 = (myx,myy)
    mypair2 = (myy,myx)
    if (not mypair2 in pairs) and (not mypair2 in pairs):
      pairs.add(mypair1)
    # end if
  # end for
  return pairs
# end def find_dimers
def Hu_idx(axes,pos,sep_min=1.0,sep_max=1.4):
  ''' locate atoms to be set to spin up '''
  pairs = find_dimers(axes,pos,sep_min,sep_max)

  Hu_idx_set = set()
  for pair in pairs:
    Hu_idx_set.add(pair[0])
  # end for pair
  return Hu_idx_set
# end def Hu_idx

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

def edit_quantum_particleset(e_particleset,centers,rs,axes,ion_width):
    # centers: wf_centers
    nions = len(centers)
    # locate even and odd sub-lattices
    Hu_idx_set = Hu_idx(axes,np.array(centers))

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
    # average <R^2>
    sig2 = 1./(4.*ion_width) # not 3.* for each direction
    ion_sig = np.sqrt(sig2)
    proton_pos = gaussian_move(centers,ion_sig) # protons are slow, initialize well

    Hu_centers = []
    Hd_centers = []
    eu_centers = []
    ed_centers = []
    for i in range(nions):
        if i in Hu_idx_set:
            Hu_centers.append(proton_pos[i])
            eu_centers.append(electron_pos[i])
        else:
            Hd_centers.append(proton_pos[i])
            ed_centers.append(electron_pos[i])
        # end if
    # end for i

    # should be 2 groups in electronic wf, one u, one d
    spin_groups = e_particleset.findall("group")

    for i in range(2): # loop through u,d
        
        # 0 for up spin, 1 for down spin
        e_section = spin_groups[i]
        
        #  place u on Hu, d on Hd
        if i==0:
            epos = eu_centers
        else: # i==1
            epos = ed_centers
        # end if

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
    pu_section = deepcopy( e_section )
    pu_section.attrib["name"] = "Hu"
    pu_section.attrib['mass'] = str(unit.mp/unit.me)
    pu_section.attrib['size'] = str(nions/2)
    pu_section.xpath('.//parameter[@name="charge"]')[0].text =  '    1    '
    pu_section.xpath('.//parameter[@name="mass"]')[0].text = '    %s    ' % str(unit.mp/unit.me)
    ppos_text = "\n"+"\n".join([" ".join( map(str,pos) ) for pos in Hu_centers])+"\n"
    pu_section.find("attrib").text = ppos_text
    e_particleset.append(pu_section)

    p_section = deepcopy( e_section )
    p_section.attrib["name"] = "Hd"
    p_section.attrib['mass'] = str(unit.mp/unit.me)
    p_section.attrib['size'] = str(nions/2)
    p_section.xpath('.//parameter[@name="charge"]')[0].text =  '    1    '
    p_section.xpath('.//parameter[@name="mass"]')[0].text = '    %s    ' % str(unit.mp/unit.me)
    ppos_text = "\n"+"\n".join([" ".join( map(str,pos) ) for pos in Hd_centers])+"\n"
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
        for hname in ["Hu","Hd"]:
         
            eHterm = deepcopy(term)

            etype = etypes[ie] # 0:u, 1:d

            eHterm.attrib["speciesA"] = etype
            eHterm.attrib["speciesB"] = hname

            coeff = eHterm.find("coefficients")
            coeff.attrib["id"] = etype+hname 
            coeff.text = Uep # initialize e-p Jastrow to clamped values
            j2.append(eHterm)

        # end for hname
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

def edit_determinantset(wf,centers,ion_width,axes):
    nions = len(centers)
    # get electronic sposet_builder
    ebuilder = wf.find('.//sposet_builder[@source="ion0"]')
    if ebuilder is None:
        raise RuntimeError('electronic sposet_builder with source="ion0" not found.')
    # end if 

    # build single-particle orbitals around "wf_centers" instead of "ion0"
    wf.find("sposet_builder").attrib['source'] = 'wf_centers'

    # build electronic single-particle orbitals around "wf_centers" instead of "ion0"
    ebuilder.set('source','wf_centers')
    # use double precision
    ebuilder.set('precision','double')

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
    bgroup.append(etree.Element('radfunc',{'contraction':'1.0','exponent':str(ion_width)}))
    pao_basis.append(bgrid)
    pao_basis.append(bgroup)
    pbasis.append(pao_basis)
    # finished construct </basisset>
    pbuilder.append(pbasis)

    # build <sposet>

    # split into up and down protons
    Hu_idx_set = Hu_idx(axes,np.array(centers))
    identity = np.eye(nions)
    Hup_det  = []
    Hdown_det= []
    for i in range(nions):
        if i in Hu_idx_set:
            Hup_det.append(identity[i])
        else:
            Hdown_det.append(identity[i])
        # end if
    # end for i

    pup_sposet = etree.Element('sposet',attrib={
        'name':'spo_uH',
        'id':'proton_orbs_up',
        'size':str(nions/2)
    })
    coeff = etree.Element("coefficient",
            {"id":"HudetC","type":"constArray","size":str(nions/2)})
    coeff.text = matrix_to_text(Hup_det)
    pup_sposet.append(coeff)

    pdn_sposet = etree.Element('sposet',attrib={
        'name':'spo_dH',
        'id':'proton_orbs_down',
        'size':str(nions/2)
    })
    coeff1 = etree.Element("coefficient",
            {"id":"HddetC","type":"constArray","size":str(nions/2)})
    coeff1.text = matrix_to_text(Hdown_det)
    pdn_sposet.append(coeff1)
    # done construct </sposet>
    pbuilder.append(pup_sposet)
    pbuilder.append(pdn_sposet)
    # end </sposet_builder>

    # build <determinant>
    uHdet = etree.Element('determinant',{
        'group':'uH',
        'id':'uHdet',
        'size':str(nions/2),
        'sposet':pup_sposet.get('name'),
        'no_bftrans':'yes' # do not transform proton coordinates
    })
    dHdet = etree.Element('determinant',{
        'group':'dH',
        'id':'dHdet',
        'size':str(nions/2),
        'sposet':pdn_sposet.get('name'),
        'no_bftrans':'yes' # do not transform proton coordinates
    })
    # end </determinant>
    slater = wf.find('determinantset/slaterdeterminant')
    slater.append(uHdet)
    slater.append(dHdet)

    ebuilder_idx = wf.index(ebuilder)
    wf.insert(ebuilder_idx+1,pbuilder) # insert proton_builder after electronic sposet_builder
# end def edit_determinantset

#  Main Routine
# ======================================
def bo_to_nobo(bo_input_name,nobo_input_name,ion_width=9.0,rs=1.31):

    parser = etree.XMLParser(remove_blank_text=True)
    xml = etree.parse(bo_input_name,parser)
    from input_xml import InputXml
    inp = InputXml()
    inp.read(bo_input_name)
    axes = inp.lattice_vectors()
    del inp

    ion0 = xml.xpath('//particleset[@name="ion0"]')[0]
    centers = change_ion0_to_wf_ceneters(ion0)
    # !! minimum atom spacing, special to cubic cell
    alat = np.sort(np.unique(np.array(centers).flatten()))[1]

    e_particleset = xml.xpath('//particleset[@name="e"]')[0]
    edit_quantum_particleset(e_particleset,centers,rs,axes,ion_width)

    wf = xml.xpath("//wavefunction")[0]
    edit_jastrows(wf)
    edit_determinantset(wf,centers,ion_width,axes)

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
