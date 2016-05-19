import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class scalar_analyzer:
    # a collection of useful functions for analyzing a trace file
    #  the trace file can contain many columns

    def __init__(self):
        self.dataframe = pd.DataFrame()
    # end def __init__

    def trace_components(self,scalar_file,header_offset=1):
        # written to read the scalar.dat file outputed by QMCPACK
        #  put all columns of data into a DataFrame and return it
        # should be applicable to any trace file of similar format
        # input: location of the scalar file
        # output: pandas DataFrame containing all columns of data
        
        # initialize empty DataFrame
        df = pd.DataFrame()
        
        # read the column names
        with open(scalar_file) as f:
            header = f.readline()
        # end with 
        
        # if scalar file is written to be numpy-readable then
        #  header shoud look like: " # Index LocalEnergy ..."
        header = header.split()[header_offset:] # remove "#"
        
        # read all columns of data (numpy uses row major)
        data = np.loadtxt(scalar_file).T
        
        for i in range(len(header)):
            # save each column of data under its header
            df[header[i]] = data[i]
        # end for i
        
        # save DataFrame object
        self.dataframe = df
        # return DataFrame object
        return df
    # end def energy_components
    
    def plot_components(self,components,df):
        # plot components from a data frame together

        fig = plt.figure()
        ax  = fig.add_subplot(111)

        for component in components:
            if component not in df.columns:
                print component + " is not in the DataFrame with columns:"
                print df.columns.values
                return
            # end if
            ax.plot(df[component],label=component)
        # end for

        ax.legend()
        fig.tight_layout()
        return fig,ax
    # end def plot_components

    def corr(self,myg):
        # calculate autocorrelation of a trace myg
        g = np.array(myg)
        mu=g.mean()
        s=g.std()

        sumR=0.0
        for k in range(1,len(g)):
            R=0.0

            # calculate autocorrelation
            for t in range(len(g)-k):
                R+=(g[t]-mu)*(g[t+k]-mu)
            #end for t
            R/=(len(g)-k)*s**2.

            # accumulate until R<=0
            if R>0:
                sumR+=R
            else:
                break
            #end if
        #end for k
        
        return 1+2.*sumR
    #end def corr

    def error(self,trace):
        return trace.std()/np.sqrt( len(trace)/self.corr(trace) )
    # end def error

# end class scalar_analyzer
