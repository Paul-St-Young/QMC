import numpy as np
import pandas as pd

def find_observable_names(all_columns):
    
    special_columns = set(['AcceptRatio', 'BlockCPU', 'BlockWeight'
      , 'Efficiency', 'TotalSamples', 'TotalTime', 'Variance'])
    
    observable_names = []
    for col in all_columns:
        if col.endswith("_mean"):
            col_name = col.replace("_mean","")
            if col_name not in special_columns:
                observable_names.append( col_name )
            # end if
        # end if
    # end for
    
    return observable_names
# end def find_observable_names

def get_better_observables(one_vmc,some_dmcs):
    
    # get VMC observable
    assert len(one_vmc)==1, "must give one and only one vmc run, given " + str(len(one_vmc))
    
    extrap = some_dmcs.copy()
    extrap["method"] = "extrap"
    
    col_names = find_observable_names(some_dmcs.columns)
    
    for col_name in col_names:
        
        mean_name  = col_name+"_mean"
        error_name = col_name+"_error"
        
        # get VMC observable
        trial  = one_vmc.loc[0,mean_name]
        triale = one_vmc.loc[0,error_name]
        
        for idx in some_dmcs.index:
            
            # get DMC observable
            mixed  = some_dmcs.loc[idx,mean_name]
            mixede = some_dmcs.loc[idx,error_name]
            
            extrap.loc[idx,mean_name]  = mixed**2./trial
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


