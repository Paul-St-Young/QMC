from lxml import etree
def write_j2section(entry):

    j2section = etree.Element('jastrow'
        ,{'type':'Two-Body','name':'J2','function':'bspline','print':'yes'})
    for dcorr in entry:

        corr_dict = {}
        for key in ['speciesA','speciesB','size','rcut']:
            corr_dict[key] = dcorr[key]
        # end for 
        corr = etree.Element('correlation',corr_dict)

        coeff = etree.Element('coefficients',{'id':dcorr['id'],'type':'Array'})
        coeff_text = "\n" + " ".join(map(str,dcorr['coeff'])) + "\n"
        coeff.text = coeff_text

        corr.append(coeff)
        j2section.append(corr)
    # end for dcorr
    
    return j2section
# end def write_j2section
