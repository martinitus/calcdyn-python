import ConfigParser
import StochasticModel
import ryanodine
import numpy
# from timeline import TimeLine

class ConfigData(object):
    def __init__(self, path):
        #print 'Initializing SimulationData....'
        self.path = path
        try:
            self.config = ConfigParser.RawConfigParser()
            self.config.read(path + "/parameters.txt")
            self.config.set('dye', 'fmin', 1.) # set fluorescence parameters, since they are not yet there...
            self.config.set('dye', 'fmax', 30.)
        except Exception as detail:
            print "Warning, failed to read parameters file\n",detail
        
        self.nodye = path.find('nodye')!=-1
        self.channels = 4

class AverageConcentrationData(ConfigData):
    def __init__(self, path):
        super(AverageConcentrationData, self).__init__(path)
        #print 'Initializing AverageConcentrationData....'
        try:
            self.concentrations = numpy.genfromtxt(path + "/avg_out.csv",\
                     dtype=[('t', '<f8'), ('cyca', '<f8'), ('endo', '<f8'), ('dye', '<f8'), ('erca', '<f8')])
            
            self.__cyca = TimeLine(self.concentrations['t'], self.concentrations['cyca'],desc = "free cytosolic calcium", ylabel = "[Ca$^{2+}$]$_{\mathrm{CY}}$", yunit = "$\mu M$")
            self.__erca = TimeLine(self.concentrations['t'], self.concentrations['erca'],desc = "ER calcium", ylabel = "[Ca$^{2+}$]$_{\mathrm{ER}}$", yunit = "$\mu M$")
            self.__endo = TimeLine(self.concentrations['t'], self.concentrations['endo'],desc = "bound endogeneous buffer", ylabel = "b$_{\mathrm{e}}$", yunit = "$\mu M$")
            self.__dye  = TimeLine(self.concentrations['t'], self.concentrations['dye'],desc = "bound dye buffer", ylabel = "b$_{\mathrm{d}}$", yunit = "$\mu M$")
        except Exception as detail:
            print "Warning, could not load average concentration files, using approximation from channel recordings\n",detail
        self.amps = None
            
    def endo(self):
        return self.__endo;
    def dye(self):
        return self.__dye;
    def cyca(self):
        return self.__cyca;
    def erca(self):
        return self.__erca;
    def totalions(self):
        return self.erca()*0.1 + self.cyca()*0.4+self.dye()*0.4+self.endo()*0.4
    def totalfluorescence(self):
        Bd     = self.config.getfloat('dye','B');
        fmin   =  self.config.getfloat('dye','fmin');
        fmax   =  self.config.getfloat('dye','fmax');
        data   = self.concentrations['dye']*fmax + (Bd-self.concentrations['dye'])*fmin
        return TimeLine(self.concentrations['t'],data,desc = 'total fluorescense', ylabel = 'F')
    def relativefluorescense(self):
        Bd     = self.config.getfloat('dye','B');
        c0     = self.config.getfloat('cytosol','c0');
        kminus = self.config.getfloat('dye','kminus');
        kplus  = self.config.getfloat('dye','kplus');
        resting =  Bd * c0 / (kminus / kplus + c0);
        fmin   =  self.config.getfloat('dye','fmin');
        fmax   =  self.config.getfloat('dye','fmax');
        data   = (self.concentrations['dye']*fmax + (Bd-self.concentrations['dye'])*fmin) / (resting*fmax + (Bd-resting)*fmin) -1;
        return TimeLine(self.concentrations['t'], data, desc = 'change of fluorescence', ylabel = '$\Delta$F/F')
    
    def amplitudes(self,condition = lambda x: True):
        rftl = self.relativefluorescense()
        intervalls = [[t+0.1,t+0.1+0.07] for t in numpy.arange(0,rftl.tmax()-10,30) if condition([t+0.1,t+0.1+0.07])]
        return [rftl.integrate(range = i) for i in intervalls]
    
class RyRData(ConfigData):
    def __init__(self, path):
        super(RyRData, self).__init__(path)
        #print 'Initializing RyRData....'
        try:
            self.data = numpy.genfromtxt(path + "/transitions.csv",\
                dtype=[('t', '<f8'), ('f1', '<i8'), ('f2', '<i8'), ('f3', '|S12'), ('noch', '<i8'), ('nocl', '<i8')]+\
                [('s'+str(i), '<i8') for i in range(self.channels)]+[('cy'+str(i), '<f8') for i in range(self.channels)]+[('er'+str(i), '<f8') for i in range(self.channels)])
        
            # Generate views
            self.__state  = self.data[['s'+str(i) for i in range(self.channels)]].view(numpy.int).reshape((self.data.size, self.channels))
            self.__cyca   = self.data[['cy'+str(i) for i in range(self.channels)]].view(numpy.float).reshape((self.data.size, self.channels))
            self.__erca   = self.data[['er'+str(i) for i in range(self.channels)]].view(numpy.float).reshape((self.data.size, self.channels))
            self.__time   = self.data['t']
            #self.__noch   = self.data['noch']
            self.__noch   = TimeLine(numpy.array([0]+[t for t in self.data['t']]),numpy.array([0]+[c for c in self.data['noch']]),interpolationorder='zero',desc='\# of open channels', ylabel = '\# of open channels')
            
            self.__meankd   = [numpy.mean(map(RyanodineModel.Kd,self.__erca[row][(self.__state[row] == 1)])) for row in range(self.data.shape[0])]
            self.__meankopen= [numpy.mean(map(RyanodineModel.kopen,self.__cyca[row][(self.__state[row] == 1)],self.__erca[row][(self.__state[row] == 1)])) for row in range(self.data.shape[0])]
            
        except Exception as detail:
            print "Warning, could not load transition file\n",detail
    
    # return a view to the calcium columns of the recarray as normal 2d numpy array
    def states(self):
        return self.__states;
    #TODO: move this to channel subclass
    #def cyca(self):
    #    return self.__cyca;
    #def erca(self):
    #    return self.__erca;
    def time(self):
        return self.__time;
    def noch(self):
        return self.__noch;
    def meantimes(self):
        return self.__meantimes;
    def meankd(self):
        return self.__meankd;
    def meankopen(self):
        return self.__meankopen;
        
        
#import matplotlib.figure

# Default plot for spatial calcium profile
class SpatialFigure(matplotlib.figure.Figure):	
	def __init__(self, data, figsize=None, dpi=None, facecolor=None, edgecolor=None, frameon=True, aoi = None, resolution = 100):
		super(SpatialFigure, self).__init__(figsize, dpi, facecolor, edgecolor, frameon)
		
		if aoi == None:
			xmin = min(data.channels[:]['x'])
			xmax = max(data.channels[:]['x'])
			width = xmax-xmin
			xmin = xmin - width * 0.2
			xmax = xmax + width * 0.2
		
			ymin = min(data.channels[:]['y'])
			ymax = max(data.channels[:]['y'])
			height = ymax-ymin
			ymin = ymin - height * 0.2
			ymax = ymax + height * 0.2
		else:
			xmin,xmax,ymin,ymax = aoi[:]
		
		xi = numpy.linspace(xmin,xmax,resolution)
		yi = numpy.linspace(ymin,ymax,resolution)
		
		# grid the data.	
		frame = data.spatial.frame(self.frame_time)
		
		data = self.calcium_data.data[frame,:]

		zi = scipy.interpolate.griddata((self.calcium_data.nodes[:,0],self.calcium_data.nodes[:,1]), data, (xi[None,:], yi[:,None]), method='cubic')

		# plot countour lines
		#~ cont2 = self.spatial.contourf(y1, z1, data, norm=plt.colors.LogNorm(lev[0],lev[len(lev)-1])) 
		from matplotlib import colors 
		from matplotlib.colors import LogNorm
		
		#~ self.contour = self.spatial.contourf(xi,yi,zi,15,norm=LogNorm())
		self.contour = self.spatial.contourf(xi,yi,zi,self.contourlevels,norm=LogNorm())
				
		#~ plot open and closed channels as X and O 
		oc = self.state_data.locations(self.frame_time, 'open')
		ic = self.state_data.locations(self.frame_time, 'inhibited')
		rc = self.state_data.locations(self.frame_time, 'resting')
		ac = self.state_data.locations(self.frame_time, 'active')
		scato = self.spatial.scatter(oc[:]['x'],oc[:]['y'],marker='o',c='black',s=100)
		scati = self.spatial.scatter(ic[:]['x'],ic[:]['y'],marker='x',c='black',s=100)
		scatr = self.spatial.scatter(rc[:]['x'],rc[:]['y'],marker='s',c='black',s=100)
		scata = self.spatial.scatter(ac[:]['x'],ac[:]['y'],marker='^',c='black',s=100)	
		n     = self.spatial.scatter([self.calcium_data.nodes[self.node,0]],[self.calcium_data.nodes[self.node,1]],marker='h',c='black',s=100)	
