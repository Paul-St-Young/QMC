import h5py
import numpy as np

def grab_stat_entries(stat_file_name,name):
    """ stat_file_name: .stat.h5 file name, name: ex. gofr """
    stat_file = h5py.File(stat_file_name)
    
    data = []
    # find each entry and extract data
    for estimator in stat_file.keys():
        if estimator.startswith(name):
            entry = {"name":estimator}
            for key,value in stat_file[estimator].iteritems():
                entry[key] = value[:]
                if entry[key].shape == (1,):
                    entry[key] = entry[key][0]
                # end if
            # end for
            data.append(entry)
        # end if
    # end for
    return data
# end def grab_stat_entries

def get_gofr_with_name(df,myname,nequil):
    mydf = df[ df.name.apply(lambda x:myname in x) ]
    rcut,dx,name,val,valsq = mydf.iloc[0]
    x = np.arange(0,rcut,dx)
    y = val[nequil:].mean(axis=0)
    return x,y
# end def

def max_overlap_ratio(y_vmc,y_dmc,eps):
    ratio = y_dmc/y_vmc
    small_idx = np.where(y_vmc<eps)
    ratio[small_idx] = 1
    return ratio
# end def
