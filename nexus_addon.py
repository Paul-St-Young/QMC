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
