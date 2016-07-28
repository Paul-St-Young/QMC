import matplotlib.pyplot as plt

class daVinci:
    
    def __init__(self):
        
        self.fig = None
        self.ax  = None
        
    # end def
    
    def single_init(self,xlabel,ylabel,label_font):
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        ax.set_xlabel(xlabel,fontsize=label_font)
        ax.set_ylabel(ylabel,fontsize=label_font)
        return fig,ax
    # end def single_init
    
    def plot(self,x,y,xlabel="x",ylabel="y",label_font=14):
        
        fig,ax = self.single_init(xlabel,ylabel,label_font)
        ax.plot(x,y)
        fig.tight_layout()
        return fig,ax
    # end def plot
    
    def contour(self,x,y,values,xlabel="x",ylabel="y",label_font=14):
        fig,ax = self.single_init(xlabel,ylabel,label_font)
        ax.contour(x,y,values)
        fig.tight_layout()
        return fig,ax
    # end def contour
    
# end class daVinci
