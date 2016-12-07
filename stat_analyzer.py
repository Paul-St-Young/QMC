#!/usr/bin/env python
import h5py
import numpy as np

def grab_stat_entries(stat_file_name,name,exact_name=False):
    """ stat_file_name: .stat.h5 file name, name: ex. gofr """
    stat_file = h5py.File(stat_file_name)
    
    data = []
    # find each entry and extract data
    for estimator in stat_file.keys():
        if estimator.startswith(name):
            if exact_name and estimator != name:
                continue
            # end if
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

def dft_orbital(fhandle,ik,ibnd,nx=0):
    
    # get orbital in reciprocal space
    gvecs = fhandle["electrons/kpoint_0/gvectors"]
    orbg  = fhandle["electrons/kpoint_%d/spin_0/state_%d/psi_g" % (ik,ibnd)]
    
    if nx == 0:
        try:
            orbr = fhandle["electrons/kpoint_0/spin_0/state_0/psi_r"]
            nx   = orbr.shape[0]
        except:
            raise IndexError("psi_r does not exist, please provide FFT grid size nx")
        # end try
    # end if
    
    # populate FFT grid
    kgrid = np.zeros([nx,nx,nx],dtype=complex)
    for i in range(len(gvecs)):
        gvec = gvecs[i]
        kgrid[tuple(gvec)] = orbg[i,0] + 1j*orbg[i,1]
    # end for i
    
    # invfft to obtain orbital in real space
    orbr = np.fft.fftshift( np.fft.ifftn(kgrid) ) * nx**3.
    
    # get real space grid points
    # !!!! specific to cubic grids for now
    lattice = fhandle["supercell/primitive_vectors"].value
    alat    = lattice[0,0]
    x       = np.arange(-alat/2,alat/2,alat/nx)
    orbr = np.abs(orbr)
    
    return x,orbr
# end def

import pandas as pd
def raw_stats(stat_files,observable,exact_name=False):
    """ given a list of files, save all raw data related to observable """

    data = []
    twistnum = 0 # !!!! need better way to determine twistnum
    # does not matter if twistnum never used
    for fname in stat_files:
        mydf = pd.DataFrame( grab_stat_entries(fname,observable,exact_name) )
        mydf['twistnum'] = twistnum
        twistnum += 1
        data.append(mydf)
    # end for
    if len(data) == 0:
        raise NotImplementedError('no observable %s found in %s etc.' % (observable,stat_files[0]))
    # end if
    df = pd.concat(data)
    
    return df
# end def raw_stats

def avg_twists(raw_df,mat_entries=['value'],drop_columns=[],skip_array=False):
    """ raw_df must at least contain the 'name', 'value', and 'twistnum' columns
     twistnum column does not have to be meaningful, it is simply there to distinguish
     data from different twists, all twisted are weighted 1. """
    
    if len(mat_entries) == 0:
        raise NotImplementedError('must provide at least one column to average.')
    # end if
    
    mat_avg = lambda x:np.mean(x,axis=0)
    groups2avg = raw_df.groupby('name')[mat_entries[0]]

    avg_df = pd.DataFrame( groups2avg.apply(mat_avg) )
    for col_name in mat_entries[1:]:
        avg_df[col_name] = raw_df.groupby('name')[col_name].apply(mat_avg)
    # end for col_name
    
    default_drops = ['name','twistnum']
    
    # get others parameters, which are assumed to be same across twists
    params = raw_df.iloc[0].drop(default_drops+mat_entries+drop_columns).to_dict()
    
    for key,val in params.iteritems():
        if (type(val) is np.ndarray) and (not skip_array):
            raise NotImplementedError(key+' is array')
        # end if
        avg_df[key] = params[key]
    # end for key
    return avg_df.reset_index()
# end def avg_twists

def equil_spectrum(matrix_data,nequil,warn=True):
    """ assuming a matrix of (nrow,ncol) = (nstep,nval) of simulation data
     each row is a snapshot of a spectrum of data, return average spectrum
     with the first nequil steps thrown out """
    
    if type(matrix_data) is not np.ndarray:
        raise NotImplementedError('first argument must be a numpy array')
    # end if
    
    nstep,nval = matrix_data.shape
    
    if (nval < nequil):
        raise IndexError('number of data points %d less than equilibration length %d' % (nval,nequil))
    elif (nval-nequil<10) and warn:
        print "CAUTION: less than 10 blocks of data remaining"
    # end if
    
    avg_data = matrix_data[nequil:].mean(axis=0)
    return avg_data
# end def equil_spectrum

import os
def analyze_stat_h5s(stat_files,nequil,prefix,exact_name=False
    ,mat_entries = {
    'gofr':['value','value_squared'],
    'sk':['value','value_squared','kpoints']
}):
    """ analyze a list of stat.h5 files generated by QMCPACK, 
     the default is to analyze g(r) and s(k), assuming the names
     'gofr*', 'sk*'. Analysis steps:
      1. load raw data from each stat.h5 file
      2. perform twist average on matrix entries
      3. throw out equilibration
      4. save each analyzed observable in a local database named
       prefix_observable.json
    """
    
    for observable in mat_entries.keys():

        target_json = '%s_%s.json'%(prefix,observable)
        if os.path.isfile(target_json):
            continue
        # end if

        # load raw data
        df = raw_stats(stat_files,observable,exact_name=exact_name)

        # average over twists
        avg_df = avg_twists(df,mat_entries[observable])

        # throw out equilibration (!!!! assume value and value_squared)
        avg_df['value'] = avg_df['value'].apply(
            lambda x:equil_spectrum(x,nequil))
        avg_df['value_squared'] = avg_df['value_squared'].apply(
            lambda x:equil_spectrum(x,nequil))
        avg_df.to_json(target_json)
        
    # end for observable
    
# end def analyze_stat_h5s

import subprocess as sp
def find_stat_h5_by_series(path,iseries,igroup=0):
    proc = sp.Popen(['ls',path],stdout=sp.PIPE,stderr=sp.PIPE)
    out,err = proc.communicate()

    h5files = []
    for fname in out.strip('\n').split('\n'):
        if fname.endswith('h5'):
            ftokens = fname.split('.')
            # is this a stat.h5 file?
            if ftokens[-2] != 'stat':
                continue
            # end if
            # determine group and series indices
            mygroup  = -1 # !!!! never used, only for debug
            myseries = -1
            group_text = ftokens[1]
            series_text= ftokens[2]
            has_group = True
            if group_text.startswith('s'):
                mygroup = 0
                myseries= int(group_text[1:])
                has_group=False
            elif group_text.startswith('g'):
                mygroup = int(group_text[1:])
            else:
                raise NotImplementedError('unrecognized group label %s'%group_text)
            # end if
            if has_group:
                if not series_text.startswith('s'):
                    raise NotImplementedError('unrecognized series label %s'%series_text)
                # end if
                myseries = int(series_text[1:])
            # end if
            if (mygroup==-1) or (myseries==-1):
                raise RuntimeError('failed to determine group/series %s %s for %s'%(mygroup,myseries,fname))
            # end if

            if myseries != iseries:
                continue
            # end if

            h5files.append(fname)
        # end if
    # end for fname
    return h5files
# end def find_h5_by_group_and_series

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='analyze g(r) and S(k) and plot if asked. QMCPACK stat.h5 files must be available in given path.')
    parser.add_argument('path',type=str,help='path to QMCPACK run directory')
    parser.add_argument('-is','--iseries',type=int,default=2,help='series index 1 -> s001')
    parser.add_argument('-exact','--use_exact_name',action='store_true',help='only collect observable with the exact name given. Default is to use name.startswith(name_given).')
    parser.add_argument('-eq','--nequil',type=int,default=5,help='number of equilibration steps')
    parser.add_argument('-psk','--plot_sk',action='store_true',help='plot S(k)')
    parser.add_argument('-pgr','--plot_gr',action='store_true',help='plot g(r)')
    parser.add_argument('--gr_name',type=str,default='gofr_e_0_1',help='plot g(r)')
    parser.add_argument('-nogr','--no_gr',action='store_true',help='do not collect g(r)')
    parser.add_argument('-nosk','--no_sk',action='store_true',help='do not collect S(k)')
    args = parser.parse_args()
    
    qmc_path= args.path
    iseries = args.iseries
    nequil  = args.nequil
    plot_sk = args.plot_sk
    plot_gr = args.plot_gr
    gr_name = args.gr_name
    exact_name = args.use_exact_name

    stat_files = find_stat_h5_by_series(qmc_path,iseries)
    nfile = len(stat_files)
    print "found %d stat.h5 files in %s" % (nfile,qmc_path)
    if (nfile<=0):
        raise RuntimeError('no file to read')
    # end if
    prefix = stat_files[0].split('.')[0]+'.equil%d'%nequil

    mat_entries = {
        'gofr':['value','value_squared'],
        'sk':['value','value_squared','kpoints']
    }
    if args.no_sk:
        mat_entries.pop('sk')
    elif argsno_gr:
        mat_entries.pop('gofr')
    # end if

    stat_file_locs = [os.path.join(qmc_path,fname) for fname in stat_files]
    analyze_stat_h5s(stat_file_locs,nequil,prefix,exact_name=exact_name
        ,mat_entries = mat_entries)
        #,mat_entries = {'sk':['kpoints','value','value_squared']})

    if plot_sk:
        import matplotlib.pyplot as plt

        # analyze S(k)
        sk_df = pd.read_json(prefix+'_sk.json')
        kvecs = sk_df['kpoints'][0]
        sk_arr= np.array(sk_df['value'][0])

        kmags = np.linalg.norm(kvecs,axis=1)
        # average S(k) if multiple k-vectors have the same magnitude
        unique_kmags = np.unique(kmags)
        unique_sk    = np.zeros(len(unique_kmags))
        for iukmag in range(len(unique_kmags)):
            kmag    = unique_kmags[iukmag]
            idx2avg = np.where(kmags==kmag)
            unique_sk[iukmag] = np.mean(sk_arr[idx2avg])
        # end for iukmag

        # plot S(k)
        fig,ax = plt.subplots(1,1)
        ax.set_xlabel('k')
        ax.set_ylabel('S(k)')
        ax.plot(unique_kmags,unique_sk)
        plt.show()
    # end if plot

    if plot_gr:
        import matplotlib.pyplot as plt

        # analyze g(r)
        gr_df = pd.read_json(prefix+'_gofr.json')
        mydf = gr_df[gr_df['name'] == gr_name]
        rmax = mydf['cutoff'][1]
        dr   = mydf['delta'][1] # never used, save for debug

        myy = mydf['value'][1]
        myx = np.linspace(0,rmax,len(myy))

        # plot g(r)
        fig,ax = plt.subplots(1,1)
        ax.set_xlabel('r (bohr)')
        ax.set_ylabel(gr_name+'(r)')
        ax.plot(myx,myy)
        plt.show()
    # end if plot_gr
# end __main__
