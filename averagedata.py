import ConfigParser
import StochasticModel
import RyanodineModel
import numpy
from timeline import TimeLine

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
            self.concentrations = np.genfromtxt(path + "/avg_out.csv",\
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

    
class RyRData(ConfigData):
    def __init__(self, path):
        super(RyRData, self).__init__(path)
        #print 'Initializing RyRData....'
        try:
            self.data = np.genfromtxt(path + "/transitions.csv",\
                dtype=[('t', '<f8'), ('f1', '<i8'), ('f2', '<i8'), ('f3', '|S12'), ('noch', '<i8'), ('nocl', '<i8')]+\
                [('s'+str(i), '<i8') for i in range(self.channels)]+[('cy'+str(i), '<f8') for i in range(self.channels)]+[('er'+str(i), '<f8') for i in range(self.channels)])
        
            # Generate views
            self.__state  = self.data[['s'+str(i) for i in range(self.channels)]].view(np.int).reshape((self.data.size, self.channels))
            self.__cyca   = self.data[['cy'+str(i) for i in range(self.channels)]].view(np.float).reshape((self.data.size, self.channels))
            self.__erca   = self.data[['er'+str(i) for i in range(self.channels)]].view(np.float).reshape((self.data.size, self.channels))
            self.__time   = self.data['t']
            #self.__noch   = self.data['noch']
            self.__noch   = TimeLine(np.array([0]+[t for t in self.data['t']]),np.array([0]+[c for c in self.data['noch']]),interpolationorder='zero',desc='\# of open channels', ylabel = '\# of open channels')
            
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
    
class SimulationData(RyRData,AverageConcentrationData):
    def __init__(self, path):
        super(SimulationData, self).__init__(path)
        #print 'Finished'