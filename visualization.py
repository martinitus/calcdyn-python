from averagedata import RyRData, AverageConcentrationData
import matplotlib.pyplot as plt
import numpy as np

class Visualization(RyRData,AverageConcentrationData):
    def __init__(self,path):
        super(Visualization, self).__init__(path)
        
    def intervallplot(self,i):
        fig = plt.figure(figsize=(12,3));
        ax = fig.gca();
        l=ax.plot(self.time(),self.noch(),'r-',drawstyle='steps-post');ax.set_xlabel('time [s]');ax.set_ylabel('\# of open channels');ax.xaxis.grid(color='gray', linestyle='dashed');
        ax.set_ylim([0,10]);ax.set_xlim([30*i,30*(i+1)]);
        axs=axes([0.62,0.25,0.27,0.25],axisbg='w');axs.tick_params(axis='x', labelsize=6);axs.tick_params(axis='y', labelsize=6)
        l=axs.plot(self.time(),self.noch(),'r-',drawstyle='steps-post');axs.set_xlabel('time [s]',fontsize=8);axs.set_ylabel('\# of open channels',fontsize=8);
        axs.set_ylim([0,10]);axs.set_xlim([i*30.+0.05,i*30+0.2]);axs.xaxis.grid(color='gray', linestyle='dashed',);
        fig.legend((),[],'upper right',title = str(i*30) + "s to " + str((i+1)*30) +"s", fancybox=False);
    

    
    def ampoverview(self,normed = False, cumulative = True, stacked = True, histtype = 'stepfilled'):
        colors = ['cyan','pink','orange','gray','yellow','blue','cyan','black','purple','pink']
        labels = [str(c) for c in range(self.channels,-1,-1)]
        amps   = [self.amplitudes(lambda range: self.noch().any(lambda x: x>=c,range) and self.noch().all(lambda x: x<c+1,range)) for c in range(self.channels,-1,-1)]
        labels = [labels[i] for i in range(self.channels,-1,-1) if len(amps[i]) > 0]
        amps   = [a for a in amps if len(a) > 0]
        colors = colors[0:len(amps)]
        plt.figure();
	
        plt.hist(amps,histtype=histtype,stacked = stacked, color = colors, bins= np.arange(0,0.2,0.002),cumulative=cumulative,normed=normed,label=labels);
        if not normed:
            plt.ylim([0,self.relativefluorescense().tmax()/30]);
        else:
            plt.ylim([0,self.channels+1]);
        plt.legend(loc='upper left')
        
    
    def overview(self,timerange = None, save = None,fluolim=None):
        fig = plt.figure(figsize=[12,6],dpi=160)
        fig.text(.5, .95, self.path, horizontalalignment='center') 
        fig.subplots_adjust(hspace = 0)
        nochax = fig.add_subplot(611);nochax.tick_params(axis='y', labelsize=8)
        cycaax = fig.add_subplot(612);cycaax.tick_params(axis='y', labelsize=8)
        ercaax = fig.add_subplot(613);ercaax.tick_params(axis='y', labelsize=8)
        endoax = fig.add_subplot(614);endoax.tick_params(axis='y', labelsize=8)
        if  not self.nodye:
            dyeax  = fig.add_subplot(615);dyeax.tick_params(axis='y', labelsize=8)      
            fluoax = fig.add_subplot(616);fluoax.tick_params(axis='y', labelsize=8)      
            
        self.noch().plot(nochax,c = 'red',trange = timerange,yrange=[0,self.channels+1])
        self.cyca().plot(cycaax,c='green',trange = timerange,yscale = 'log')
        self.erca().plot(ercaax,c='blue',trange = timerange)
        self.endo().plot(endoax,c='y',trange = timerange)
    
        if not self.nodye:
            self.dye().plot(dyeax,c='orange',trange = timerange)
            self.relativefluorescense().plot(fluoax,c='orange',trange = timerange)        
    
        if(save != None ):
            fig.savefig(save,dpi = 120)