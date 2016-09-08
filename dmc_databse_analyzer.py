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

def add_subfix(columns):
    columnv  = [name + nscheme.subfix_mean  for name in columns]
    columne  = [name + nscheme.subfix_error for name in columns]
    return columnv,columne
# end def add_subfix

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
    new_df = pd.concat([new_df,new_df["settings"].apply(pd.Series)],axis=1).drop("settings",axis=1)
    return new_df
# end def process_dmc_data_frame

def sum_columns(cols_to_sum,df):
    columnv, columne = add_subfix(cols_to_sum)
    new_mean  = df[columnv].apply(np.sum,axis=1)
    new_error = df[columne].apply(lambda arr:np.sqrt(np.sum(arr**2.)),axis=1)
    return new_mean,new_error
# end def

def add_sum_columns(new_name,cols_to_sum,df):
    new_mean, new_error = sum_columns(cols_to_sum,df)
    df[new_name+nscheme.subfix_mean]  = new_mean
    df[new_name+nscheme.subfix_error] = new_error
# end def

