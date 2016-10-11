import numpy as np
import pandas as pd

class NamingScheme:
    def __init__(self):
        self.subfix_mean  = "_mean"
        self.subfix_error = "_error"
        self.special_columns = set(['AcceptRatio', 'BlockCPU', 'BlockWeight'
      , 'Efficiency', 'TotalSamples', 'TotalTime', 'Variance', 'LocalEnergy'])
    # end def __init__
# end class

# global class is a little better than a bunch of global varibales?
nscheme = NamingScheme()

def add_subfix(name):
    return name + nscheme.subfix_mean, name + nscheme.subfix_error
# end def

def add_subfixes(columns):
    columnv  = [name + nscheme.subfix_mean  for name in columns]
    columne  = [name + nscheme.subfix_error for name in columns]
    return columnv,columne
# end def add_subfixes

def find_observable_names(all_columns):
    
    special_columns = nscheme.special_columns
    
    observable_names = []
    for col in all_columns:
        if col.endswith(nscheme.subfix_mean):
            col_name = col.replace(nscheme.subfix_mean,"")
            if col_name not in special_columns:
                observable_names.append( col_name )
            # end if
        # end if
    # end for
    
    return observable_names
# end def find_observable_names

import scipy.optimize as op
def ts_extrap(mydf,method='linear',ts_name='timestep'
    ,val_name='LocalEnergy_mean',err_name='LocalEnergy_error',estimate_error=True):
    ts = mydf[ts_name].values
    values = mydf[val_name].values
    errors = mydf[err_name].values
    
    if method == 'linear':
        def model(ts,a,b):
            return a*ts+b
        # end def
    else:
        raise NotImplementedError('unknown extrapolation method '+method)
    # end if method
    
    popt,pcov = op.curve_fit(model,ts,values,sigma=errors,absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    fitted_model = lambda ts:model(ts,*popt)
    
    val0 = fitted_model(0)
    err0 = 0.0
    
    if estimate_error:
        if method == 'linear':
            err0 = perr[1]
        else:
            raise NotImplementedError('how to estimate error with '+method)
        # end if method
    # end if
    
    return pd.Series({val_name:val0,err_name:err0})
# end def

def get_better_observables(one_vmc,some_dmcs,linear=True):
    """ extrapolate mixed estimators using a vmc run
      linear extrapolation: 2*DMC-VMC
      non-linear extrap.: DMC^2/VMC"""
    
    # get VMC observable
    assert len(one_vmc)==1, "must give one and only one vmc run, given " + str(len(one_vmc))
    
    extrap = some_dmcs.copy()
    extrap["method"] = "extrap"
    
    col_names = find_observable_names(some_dmcs.columns)
    
    for col_name in col_names:
        
        mean_name  = col_name + nscheme.subfix_mean
        error_name = col_name + nscheme.subfix_error
        
        # get VMC observable
        trial  = one_vmc.loc[0,mean_name]
        triale = one_vmc.loc[0,error_name]
        
        for idx in some_dmcs.index:
            
            # get DMC observable
            mixed  = some_dmcs.loc[idx,mean_name]
            mixede = some_dmcs.loc[idx,error_name]
            
            if linear:
                extrap.loc[idx,mean_name]  = 2*mixed - trial
            else:
                extrap.loc[idx,mean_name]  = mixed**2./trial
            # end if

            extrap.loc[idx,error_name] = np.sqrt(triale**2.+mixede**2.)
        # end for col_name
        
    # end for idx
    
    return extrap
# end def

def process_dmc_data_frame(df):
    """ df should be the output of pd.DataFrame( nexus_addon.scalars_from_input(qmcpack_input) ) """
    vmc = df[ df["method"] == "vmc" ]
    dmc = df[ df["method"] == "dmc" ]
    extrap = get_better_observables(vmc,dmc)
    new_df = pd.concat([df,extrap]).reset_index().drop("index",axis=1)
    try:
        new_df = pd.concat([new_df,new_df["settings"].apply(pd.Series)],axis=1).drop("settings",axis=1)
    except:
        pass
    # end try
    return new_df
# end def process_dmc_data_frame

def sum_columns(cols_to_sum,df):
    columnv, columne = add_subfixes(cols_to_sum)
    new_mean  = df[columnv].apply(np.sum,axis=1)
    new_error = df[columne].apply(lambda arr:np.sqrt(np.sum(arr**2.)),axis=1)
    return new_mean,new_error
# end def

def sub_columns(cols_to_sub,df):
    columnv, columne = add_subfixes(cols_to_sub)
    leftv,rightv = columnv
    new_mean  = df[leftv] - df[rightv]
    new_error = df[columne].apply(lambda arr:np.sqrt(np.sum(arr**2.)),axis=1)
    return new_mean,new_error
# end def

def mult_columns(cols_to_mult,df):
    columnv, columne = add_subfixes(cols_to_mult)
    new_mean  = df[columnv].apply(np.prod,axis=1)
    new_error = df[columne].apply(lambda arr:np.sqrt(np.sum(arr**2.)),axis=1)
    return new_mean,new_error
# end def

def div_columns(cols_to_div,df):
    columnv, columne = add_subfixes(cols_to_div)
    numv,denv = columnv
    nume,dene = columne
    ratio_mean  = abs(df[numv]/df[denv])
    ratio_error = np.sqrt(df[nume]**2. + df[dene]**2.)
    return ratio_mean,ratio_error
# end def

def sub_rows(irow_label,jrow_label,cols_to_sub,df):
    """ return a new row with df.iloc[irow] - df.iloc[jrow] to df"""
    columnv, columne = add_subfixes(cols_to_sub)
    new_row = df.loc[irow_label].copy()
    for icol in range(len(cols_to_sub)):
        new_row[columnv[icol]] = df.loc[irow_label,columnv[icol]] - df.loc[jrow_label,columnv[icol]]
        new_row[columne[icol]] = np.sqrt( df.loc[irow_label,columne[icol]]**2 + df.loc[jrow_label,columne[icol]]**2 )
    # end for i

    return new_row
# end def

def add_sum_columns(new_name,cols_to_sum,df):
    new_mean, new_error = sum_columns(cols_to_sum,df)
    df[new_name+nscheme.subfix_mean]  = new_mean
    df[new_name+nscheme.subfix_error] = new_error
# end def

def display_mean_error(mean,error):
    """ return string rep. of mean +- error in readable format eg. 7.21(2) """
    
    # determine the last digit above error
    if error < 1e-16:
        return "%f(0)" % mean
    # end if
    guess = int( round( np.log10(error) ) )
    for digit in range(guess+3,guess-3,-1):
        if error > 10**digit:
            break
        # end if
    # end for
    digit = -digit
    
    err_str = str( int( round(error,digit) * 10**digit ) )
    return str( round(mean,digit) ) + "(" + err_str + ")"

# end def display_mean_error

def add_display_columns(interest,df):
    """ add mean(error) string columns for every quantity of interest """
    for col in interest:
        df[col] = df[col+nscheme.subfix_mean].copy()
        for idx in df.index:
            mean  = df.loc[idx,col + nscheme.subfix_mean]
            error = df.loc[idx,col + nscheme.subfix_error]
            if np.isnan(mean) or np.isnan(error):
                df.loc[idx,col] = np.nan
            else:
                df.loc[idx,col] = display_mean_error(mean,error)
            # end if
        # end for idx
    # end for col
# end def add_display_columns

def vmc_dmc_extrap_table(static_json,dynamic_json):
    # read databases
    static   = pd.read_json(static_json)
    dynamic  = pd.read_json(dynamic_json)
    static["lattice"]  = "static"
    dynamic["lattice"] = "dynamic"

    # concatenate VMC, DMC and extrap runs
    table = pd.DataFrame()
    for mydf in [static,dynamic]:
        vmc    = mydf.loc[(mydf["method"]=="vmc")]
        dmc    = mydf.loc[(mydf["method"]=="dmc") & (mydf["iqmc"]==mydf["iqmc"].max())]
        extrap = mydf.loc[(mydf["method"]=="extrap") & (mydf["iqmc"]==mydf["iqmc"].max())]
        table  = table.append([vmc,dmc,extrap])
    # end for mydf
    table = table.reset_index().drop("index",axis=1)
    
    return table
# end def

def static_dynamic_energetics_table(vde_table
    , natom, interest = ["LocalEnergy","Te","Tp","Vee","Vep","Vpp"]):

    # select columns of interest
    columnv, columne = add_subfixes(interest)

    # per atom
    vde_table.loc[:,columnv+columne] = vde_table[columnv+columne]/natom
    vde_table = vde_table.loc[:,["method","lattice"]+columnv+columne]

    # no proton kinetic for static lattice calculations
    if "Tp" in interest:
        vde_table.loc[vde_table["lattice"]=="static","Tp_mean"] = 0
    # end if

    # subtract rows, !!!! hard-coded indices 
    dmc_dyn_stat    = sub_rows(4,1,interest,vde_table)
    extrap_dyn_stat = sub_rows(5,2,interest,vde_table)

    for entry in [dmc_dyn_stat,extrap_dyn_stat]:
        entry["lattice"] = "dyn.-stat."
        vde_table = vde_table.append(entry)
    # end for
    
    vde_table = vde_table.reset_index().drop("index",axis=1)
    
    return vde_table
# end def
