#!/usr/bin/env python
import h5py
import numpy as np

def grab_stat_entries(stat_file_name,name):
    """ stat_file_name: .stat.h5 file name, name: ex. gofr """
    stat_file = h5py.File(stat_file_name)
    
    data = []
    # find entry and extract data
    for estimator in stat_file.keys():
        if estimator != name:
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
    # end for
    if len(data) == 0:
        raise RuntimeError('no entry named %s was found in %s'%(name,stat_file_name))
    # end if
    return data
# end def grab_stat_entries

def sorted_configs(fconf):
    fp = h5py.File(fconf)

    all_conf_keys = [key for key in fp.keys() if key.startswith('walkers') and key!='walkers']
    nblock = len(all_conf_keys) # number of blocks recorded
    nwalker,natom,ndim = fp[all_conf_keys[0]].value.shape

    configs = np.zeros([nwalker*nblock,natom,ndim])
    iconf = 0
    conf_idx = [int(key.strip('walkers')) for key in all_conf_keys]
    srt_idx = np.argsort(conf_idx)
    sorted_keys = np.array(all_conf_keys)[srt_idx] # place blocks in order
    for key in sorted_keys:
        iblock = int( key.strip('walkers') )
        configs[iconf:iconf+nwalker,:,:] = fp[key].value.astype(float)
        iconf += nwalker
    # end for
    mapping = '#' + '\n#'.join(sorted_keys)

    return configs,(nblock,nwalker,natom,ndim),mapping
# end def

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

def show_sk(sk_json):
    import matplotlib.pyplot as plt

    # analyze S(k)
    sk_df = pd.read_json(sk_json)
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
# end def show_sk

def show_gr(gr_json,gr_name):
    import matplotlib.pyplot as plt
    # analyze g(r)
    gr_df = pd.read_json(gr_json)
    mydf = gr_df[gr_df['name'] == gr_name]
    if len(mydf) != 1:
        raise RuntimeError('found %d entries with name %s, expected 1' %(len(mydf),gr_name))
    # end if
    rmax = mydf['cutoff'].iloc[0]
    dr   = mydf['delta'].iloc[0] # never used, save for debug

    myy = mydf['value'].iloc[0]
    myx = np.linspace(0,rmax,len(myy))

    # plot g(r)
    fig,ax = plt.subplots(1,1)
    ax.set_xlabel('r (bohr)')
    ax.set_ylabel(gr_name+'(r)')
    ax.plot(myx,myy)
    plt.show()
# end def show_gr


import pandas as pd
def raw_stats(stat_files,observable):
    """ given a list of files, save all raw data related to observable """

    data = []
    for fname in stat_files:
        mydf = pd.DataFrame( grab_stat_entries(fname,observable) )
        data.append(mydf)
    # end for
    if len(data) == 0:
        raise NotImplementedError('no observable %s found in %s etc.' % (observable,stat_files[0]))
    # end if
    df = pd.concat(data).reset_index(drop=True)
    
    return df
# end def raw_stats

def avg_twists(raw_df,meta_cols,drop_columns=[],skip_array=False):
    """ raw_df must at least contain the 'name','value_mean', meta_cols columns """
    
    err_avg = lambda x:np.sqrt(np.mean(x**2.))/np.sqrt(len(x))

    # check name of observable
    myname_list = raw_df['name'].unique()
    if len(myname_list) !=1:
        raise RuntimeError('expect one observable, found %s'%str(myname_list))
    # end if
    myname = myname_list[0] # save name

    # make twist-averaged database
    data = []
    '''
    for name in mat_cols:
        values = np.mean(raw_df[name].values,axis=0)
        errors = np.apply_along_axis(err_avg,0,raw_df[name.replace('_mean','_error'))
        data.append(values)
        data.append(errors)
    # end for name
    avg_df = pd.DataFrame(data,index=mat_cols).T
    '''
    values = np.mean(raw_df['value_mean'].values,axis=0)
    errors = np.apply_along_axis(err_avg,0,raw_df['value_error'].values)
    avg_df = pd.DataFrame([[values],[errors]],index=['value_mean','value_error']).T

    avg_df['name']  = myname
    avg_df['nblock'] = [raw_df['nblock'].values]
    for mcol in meta_cols:
        meta_data = raw_df.iloc[0][mcol]
        if type(meta_data) is np.ndarray:
            avg_df[mcol] = [meta_data]
        else:
            avg_df[mcol] = meta_data
        # end if
    # end for mcol

    return avg_df
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

def corr(trace):
    """ calculate the autocorrelation of a trace of scalar data
    pre:  trace should be a 1D iterable array of floating point numbers
    post: return the autocorrelation of this trace of scalars
    """

    mu     = np.mean(trace)
    stddev = np.std(trace,ddof=1)

    correlation_time = 0.
    for k in range(1,len(trace)):
        # calculate auto_correlation
        auto_correlation = 0.0
        num = len(trace)-k
        for i in range(num):
            auto_correlation += (trace[i]-mu)*(trace[i+k]-mu)
        # end for i
        auto_correlation *= 1.0/(num*stddev**2)
        if auto_correlation > 0:
            correlation_time += auto_correlation
        else:
            break
        # end if
    # end for k

    correlation_time = 1.0 + 2.0*correlation_time
    return correlation_time

# end def corr

def mat_error(matrix_data):
    ndat   = matrix_data.shape[0] # assume row major
    stddev = np.std(matrix_data,ddof=1,axis=0)
    corrdat= np.apply_along_axis(corr,0,matrix_data)
    return stddev/np.sqrt(ndat/corrdat)
# end def

import os
def analyze_stat_h5s(stat_files,nequil,prefix,observable,warn=True):
    """ analyze a list of stat.h5 files generated by QMCPACK, 
     the default is to analyze g(r) and s(k), assuming the names
     'gofr*', 'sk*'. Analysis steps:
      1. load raw data from each stat.h5 file
      2. perform twist average on matrix entries
      3. throw out equilibration
      4. return the analyzed database
       prefix_observable.json
    """

    # !!!! hard-coded maps
    """
    mat_entries = {
    'gofr':['value','value_squared','value_mean','value_error'],
    'sk':['value','value_squared','kpoints','value_mean','value_error']
    }
    """
    meta_entries = {
    'gofr':['name','cutoff'],
    'sk':['name','kpoints'],
    'skc':['name','kpoints'],
    'latdev':['name'],
    'skinetic':['name']
    }
    
    # load raw data
    df = raw_stats(stat_files,observable)

    # check if some twists have missing blocks
    #  also throw out equilibration (i.e. slice out nequil:minb)
    df['nblock'] = df['value'].apply(lambda x:x.shape[0])

    # cutout useful part of trace [nbegin:]
    nbegin = nequil
    minb = df['nblock'].min()
    maxb = df['nblock'].max()
    if minb != maxb:
        if warn:
            print " only using %d blocks for %s from %s" % (
                minb,observable,os.path.dirname(stat_files[0]) )
        # end if
    # end if
    if nequil >= minb:
        raise RuntimeError('not enough data (%d) for %d equilibration steps' % (minb,nequil))
    # end if
    
    for col in df: # !!!! assume only column with names * need to be averaged
        if col.startswith('value'):
            df.loc[:,col] = df[col].apply(lambda x:x[nbegin:])
        # end if
    # end for

    # average over blocks
    df['value_mean']  = df['value'].apply(np.mean,axis=0)
    df['value_error'] = df['value'].apply(mat_error)

    # average over twists
    if observable.startswith('gofr'):
        meta_data = meta_entries['gofr']
    else:
        meta_data = meta_entries[observable]
    # end if

    avg_df = avg_twists(df,meta_data)

    return avg_df
        
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
    parser.add_argument('-is','--iseries',type=int,default=2,help='series index 2 -> s002')
    parser.add_argument('-eq','--nequil',type=int,default=5,help='number of equilibration steps')
    parser.add_argument('-psk','--plot_sk',action='store_true',help='plot S(k)')
    parser.add_argument('-pgr','--plot_gr',action='store_true',help='plot g(r)')
    parser.add_argument('--gr_name',type=str,default='gofr_e_0_1',help='plot g(r)')
    parser.add_argument('-skipa','--skip_analysis',action='store_true',help='skip analysis, allow plotting if the desired json file already exists')
    parser.add_argument('-pre','--prefix',type=str,default=None,help='e.g. to plot myrun_gofr.json, use -pre myrun -skipa')
    args = parser.parse_args()
    
    qmc_path= args.path
    iseries = args.iseries
    nequil  = args.nequil
    plot_sk = args.plot_sk
    plot_gr = args.plot_gr
    gr_name = args.gr_name

    # dynamic-lattice
    observables = ['latdev','gofr_e_0_0','gofr_e_0_1','gofr_e_0_2','gofr_e_1_2','gofr_e_2_2','sk','skc']
    # static-lattice
    observables = ['gofr_e_0_0','gofr_e_0_1','sk']

    # this analysis block generates 'prefix_sk.json' and 'prefix_gofr.json'
    if not args.skip_analysis:

        stat_files = find_stat_h5_by_series(qmc_path,iseries)
        nfile = len(stat_files)
        print "found %d stat.h5 files in %s" % (nfile,qmc_path)
        if (nfile<=0):
            raise RuntimeError('no file to read')
        elif (nfile==1):
            print stat_files[0] # show file if few
        # end if
        prefix = stat_files[0].split('.')[0]+'.equil%d'%nequil

        stat_file_locs = [os.path.join(qmc_path,fname) for fname in stat_files]

        data = []
        for observable in observables:
            data.append( analyze_stat_h5s(stat_file_locs,nequil,prefix,observable) )
        # end for
        df = pd.concat(data).reset_index(drop=True)

        gr_sel = df['name'].apply(lambda x:x.startswith('gofr'))
        sk_sel = df['name'].apply(lambda x:x=='sk' or x=='skc')
        other_sel = ~(gr_sel | sk_sel)
        if gr_sel.any():
            df[gr_sel].dropna(axis=1).reset_index(drop=True).to_json(prefix+'_gofr.json')
        # end if
        if sk_sel.any():
            df[sk_sel].dropna(axis=1).reset_index(drop=True).to_json(prefix+'_sk.json')
        # end if
        if other_sel.any():
            df[other_sel].dropna(axis=1).reset_index(drop=True).to_json(prefix+'_other.json')
        # end if
    # end if

    if args.skip_analysis:
        prefix = args.prefix
        if prefix is None:
            raise RuntimeError('must provide prefix if analysis is skipped')
        # end if
    # end if

    if plot_sk:
        sk_json = prefix+'_sk.json'
        show_sk(sk_json)
    # end if plot

    if plot_gr:
        gr_json = prefix+'_gofr.json'
        show_gr(gr_json,gr_name=gr_name)
    # end if plot_gr
# end __main__
