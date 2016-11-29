import numpy as np
from mmap import mmap
# Aldous Huxley

class MoldenNormalMode:
    
    def __init__(self,filename,nmode=None,natom=None):
        self.fhandle   = open(filename,"r+")
        self.mm        = mmap(self.fhandle.fileno(),0)
        self.sections  = None
        self.nmode     = nmode
        self.natom     = natom
    # end def
    
    def find_sections(self,start_marker="[",nsection_max=20):
    
        self.mm.seek(0)
        sections = {}
        for isection in range(nsection_max):
            # each section is marked by a []
            idx = self.mm.find(start_marker)
            if idx == -1:
                break
            # end if

            self.mm.seek(idx)
            section_name = self.mm.readline().strip()[1:-1]
            sections[section_name] = idx
        # end for
        
        self.sections = sections
        return sections
    # end def find_sections
    
    def section_text(self,section_label,nline_max = 2048):
        if self.sections is None:
            raise NotImplementedError("call find_sections() first")
        # end if
        
        self.mm.seek( self.sections[section_label] )
        self.mm.readline() # skip section label line
        block = ''
        for iat in range(nline_max):
            line = self.mm.readline()
            if line.startswith("[") or self.mm.size() == self.mm.tell():
                break
            # end if
            block += line
        # end for
        return block
    # end def
    
    def find_modes(self, nmode_max = 256, start_marker = "vibration"):
        
        self.mm.seek(0)
        
        vib_start_idx = []
        for imode in range(nmode_max):
            idx = self.mm.find(start_marker)
            vib_start_idx.append(idx)
            if idx == -1:
                break
            # end if
            self.mm.seek(idx)
            self.mm.readline()
        # end for

        return vib_start_idx
    # end def
    
    def read_modes(self,ndim=3):
        
        sections = self.find_sections()
        atom_lines = self.section_text("FR-COORD").split("\n")[:-1] 
        natom = len(atom_lines)
        self.natom = natom
        
        vib_start_idx = self.find_modes()
        nmode = len(vib_start_idx)-1
        self.nmode = nmode

        normal_modes = np.zeros([nmode,natom,ndim])
        for imode in range(nmode):
            vibration_block = self.mm[vib_start_idx[imode]:vib_start_idx[imode+1]]
            disp_vecs_texts = vibration_block.split("\n")[1:]
            disp_vecs = [map(float,line.split()) for line in disp_vecs_texts]
            if len(disp_vecs) == natom +1:
                disp_vecs = disp_vecs[:-1]
            # end if
            normal_modes[imode,:,:] = disp_vecs
        # end for imode
        return normal_modes
    # end def

    def read_freqs(self):
        if self.nmode is None:
          self.read_modes()
        # end if

        freq_lines = self.section_text("FREQ").split("\n")
        if len(freq_lines) == self.nmode+1:
            freq_lines = freq_lines[:-1]
        # end if
        
        freqs = map(float,freq_lines)
        return freqs
    # end def
# end class MoldenNormalMode

from copy import deepcopy
def read_two_body_jastrows(jastrows):
    
    if (jastrows.attrib["type"] != "Two-Body"):
        raise TypeError("input is not a two-body Jastrow xml node")
    elif (jastrows.attrib["function"].lower() != "bspline"):
        raise NotImplementedError("can only handle bspline Jastrows for now")
    # end if
    
    data = []
    for corr in jastrows.xpath('./correlation'):
        coeff = corr.xpath('./coefficients')[0]
        entry = deepcopy( corr.attrib )
        entry.update(coeff.attrib)
        entry['coeff'] = np.array(coeff.text.split(),dtype=float)
        data.append(entry)
    # end for corr
    
    return data
# end def read_two_body_jastrows

class SearchableFile:

    def __init__(self,fname):
        self.fhandle = open(fname,'r+')
        self.mm = mmap(self.fhandle.fileno(),0)
    # end def __init__

    def __del__(self):
        self.fhandle.close()
    # end def

    def rewind(self):
        self.mm.seek(0)
    # end def

    def find(self,string):
        idx = self.mm.find(string.encode())
        return idx
    # end def

    def find_first(self,string):
        self.rewind()
        idx = self.mm.find(string.encode())
        return idx
    # end def

    def locate_block(self,header,trailer):
        begin_idx = self.mm.find(header.encode())
        end_idx   = self.mm.find(trailer.encode())
        return begin_idx,end_idx
    # end def

    def block_text(self,header,trailer,skip_header=True,skip_trailer=True):
        bidx,eidx = self.locate_block(header,trailer)
        if skip_header:
            self.mm.seek(bidx)
            self.mm.readline()
            bidx = self.mm.tell()
        # end if
        if not skip_trailer:
            self.mm.seek(eidx)
            self.mm.readline()
            eidx = self.mm.tell()
        # end if
        return self.mm[bidx:eidx]
    # end def block_text

# end class SearchableFile

class BlockInterpreter:

    def __init__(self):
        pass
    # end def

    def matrix(self,text,columns_begin=0,columns_end=-1):
        rows = text.split('\n')
        if len(rows[-1])==0:
            rows.pop()
        # end if
        if columns_end == -1:
            matrix_text = [row.split()[columns_begin:] for row in rows]
        else: 
            matrix_text = [row.split()[columns_begin:columns_end] for row in rows]
        # end if
        return np.array(matrix_text,dtype=float)
    # end def matrix

# end class BlockInterpreter
