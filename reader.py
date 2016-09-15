import numpy as np
from mmap import mmap
# Aldous Huxley

class MoldenNormalMode:
    
    def __init__(self,filename):
        self.fhandle   = open(filename,"r+")
        self.mm        = mmap(self.fhandle.fileno(),0)
        self.sections  = None
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
    
    def find_modes(self, nmode_max = 50, start_marker = "vibration"):
        
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
        
        vib_start_idx = self.find_modes()
        nmode = len(vib_start_idx)-1

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
# end class
