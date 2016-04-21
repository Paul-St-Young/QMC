def qmca_ev_parse(qmca_output):
    # parse output from command:
    #  qmca -q ev $scalar_file
    evline=qmca_output[1].split()
    e=float(evline[3])
    ee=float(evline[5])
    v=float(evline[6])
    ve=float(evline[8])
    return e,ee,v,ve
# end def qmca_ev_parse

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
