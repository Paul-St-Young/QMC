#!/usr/bin/env python
import os
import pandas as pd

import qmcpack_reader as qpr
import nexus_addon as na

if __name__ == '__main__':

    subdir = './'

    folders = ['../opt']

    for subdir in folders:
        fname  = 'cmca4-blyp-3600-ecut60-opt.in.xml'

        # get jastrows (jas.json)
        myinput = os.path.join(subdir,fname)
        qpr.extract_jastrows(myinput)
    # end for subdir

    for subdir in folders:
        # get scalars (opt_scalar.json)
        scalar_json = os.path.join(subdir,'opt_scalar.json')
        if os.path.isfile(scalar_json):
            continue
        # end if

        # initialize analyzer
        from qmca import QBase
        options = {"equilibration":"auto"}
        QBase.options.transfer_from(options)

        entry = na.scalars_from_input(myinput)
        pd.DataFrame(entry).to_json(scalar_json)
    # end for subdir

    data = []
    for subdir in folders:
        # collect best jastrow 
        best_jas = qpr.collect_best_jastrow_set(subdir)

        data.append(best_jas)
    # end for 
    print data

# end __main__
