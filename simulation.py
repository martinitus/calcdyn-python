#import CalciumData
#import ModelData

import ConfigParser
import numpy
import math

from scipy import integrate

from timeline import TimeLine
from modeldata import EventData
from binarydata import Domain
from fakedata import FakeDomain
from membrane import Membrane

class Channel(object):
	''' 
		infos:       a tuple containing channel id, channel location and cluster id
		transitions: a list of channel transitions, the channel will extract its own transitions from the global transition list
	'''
	def __init__(self,infos,eventdata,simulation):
		#print infos
		self.__location = infos['location']
		self.__index = infos['id']
		self.__radius = infos['radius']
		self.__cluster = infos['cluster']
		self.__eventdata = eventdata
		self.__simulation = simulation
		
		self.__state   = None;
	
	def state(self):
		if self.__state == None:
			i     = self.__index
			data  = self.__eventdata._data
			
			initialstate = data['states'][0][i]
			finalstate   = [data['states'][-1][i]]
			
			initialtime  = self.__simulation.tmin()
			finaltime    = self.__simulation.tmax()
			
			frames  = data[data['chid'] == i]['t']
			states  = data[data['chid'] == i]['states'][:,i]
			
			#insert initial and final state
			states  = numpy.insert(states,0, initialstate,axis = 0)
			states  = numpy.append(states,   finalstate,  axis = 0)
			
			frames  = numpy.insert(frames,0,initialtime,axis = 0)
			frames  = numpy.append(frames,  finaltime)
			
			self.__state = TimeLine(frames,states,interpolationorder = 'zero')
		return self.__state

	def radius(self):
		return self.__radius
		
	def location(self):
		return self.__location
		
	def index(self):
		return self.__index
		
	def cluster(self):
		return self.__cluster
		
	def open(self,t):
		return self.__eventdata.state('open')['condition'](self.state().at(t))
		
	def fluxcoefficient(self):
		if self.__simulation.config.has_option('membrane','pc'):
			return self.__simulation.config.getfloat('membrane','pc')
		else:
			ic   = self.__simulation.config.getfloat('membrane','ic')
			r    = self.__radius
			erc0 = self.__simulation.config.getfloat('er','c0')
			cyc0 = self.__simulation.config.getfloat('cytosol','c0')
			e    =  1.602176E-19
			Na   =  6.022141E23
			return ic * 1.E18 / (2 * r * r * math.pi * e * Na * (erc0 - cyc0))
		
	'''def current(self,t):
		ersnapshot = self.__simulation.domain('er')['calcium'].snapshot(t)
		cysnapshot = self.__simulation.domain('cytosol')['calcium'].snapshot(t)
		x,y=self.__location[0:2]
		fa = 96485 # faraday number
		foo = integrate.quad(lambda r: 2*math.pi*r*(ersnapshot([[x+r,y]])[0] - cysnapshot([[x+r,y]])[0]),0,self.__radius)
		return foo[0]*(1E-6*1E-24)*2*fa*self.fluxcoefficient()*1E12'''
		
	def current_integrate(self):
		x,y=self.__location[0:2]
		fa = 96485 # faraday number
		n  = len(self.__simulation.domain('cytosol')['calcium'].frames)
		foo = numpy.ndarray(n,dtype=numpy.float32)
		#print 'steps:', n
		for i in range(n):
		#	print i
			t = self.__simulation.domain('cytosol')['calcium'].frames[i]
			ersnapshot = self.__simulation.domain('er')['calcium'].snapshot(t)
			cysnapshot = self.__simulation.domain('cytosol')['calcium'].snapshot(t)
			if self.open(t):
				foo[i] = integrate.quad(lambda r: 2*math.pi*r*(ersnapshot([[x+r,y]])[0] - cysnapshot([[x+r,y]])[0]),0,self.__radius)[0]
			else:
				foo[i] = 0
		foo = foo*(1E-6*1E-24)*2*fa*self.fluxcoefficient()*1E12
		return self.__simulation.domain('cytosol')['calcium'].frames, foo
		
	def current(self):
		x,y= self.__location[0:2]
		r  = self.__radius
		fa = 96485 # faraday number
		t1,cy = self.__simulation.domain('cytosol')['calcium'].evolution(x,y)
		t2,er = self.__simulation.domain('er')['calcium'].evolution(x,y)
		
		state = self.state()
		
		mask = [self.__simulation._events.state('open')['condition'](state(ti)) for ti in t1]
		current = r*r*math.pi*2*fa*self.fluxcoefficient()*(1E-6*1E-24)*1E12 * (er-cy)*mask
		
		return t1,current
		
	def calcium(self):
		return self.__simulation.domain('cytosol')['calcium'].evolution(self.__location[0:2])
		
		
	'''def current_fast(self):	
		x,y    = self.__location[0:2]
		frames = self.__simulation.domain('cytosol')['calcium'].frames
		data   = 
		cur    = numpy.ndarray(len(frames),dtype=numpy.float32)
		for i in range(len(frames)):
			interpolator = scipy.interpolate.LinearNDInterpolator(nodes, data[7221,:],fill_value = numpy.nan)	
		datas         = self.data.swapaxes(0,1)
		
		# create time domain interpolator
		timeinterpol  = scipy.interpolate.interp1d(self.frames,datas,copy = False, fill_value = numpy.nan)
		
		# calculate scattered data values for time t
		datat         = timeinterpol(t)
		
		#create interpolator for spatial coordinates	
		return scipy.interpolate.LinearNDInterpolator(self.nodes, datat,fill_value = numpy.nan)		'''


class Cluster(object):
	def __init__(self, id , channels):
		self._channels = channels
		
	def channels(self):
		return self._channels
	
	def channel(self,i):
		return self._channels[i]

# hold the combination of event data, spatial data and channel data
class Simulation(object):
	def __init__(self, path):
		self.path = path
		self.config = ConfigParser.RawConfigParser()
		self.config.read(path + "/parameters.txt")
		
		# the amount of contained channels
		self._channelcount = self.config.getint('ChannelSetup','channels')
		self._clustercount = self.config.getint('ChannelSetup','clusters')
		
		# the event data (ModelData)
		modelname = self.config.get('Meta','channelmodel')
		
		self._events    = EventData(path, modelname, self._channelcount)
		
		# the spatial data
		self._spatial  = {}
		for domain in self.domains():
			self._spatial[domain] = Domain(path, domain, components = self.config.get(domain,'components').split(','))
		
		if not self._spatial.has_key('er'):
			self._spatial['er'] = FakeDomain('er',self)
			
		# the channels 
		channeldata = numpy.genfromtxt(path + '/channels.csv', dtype=[('id', int), ('cluster', int), ('location', float, (3)), ('radius', float)])
		
		# since genfromtxt returns 1d array for single line files we need to reshape
		if channeldata.ndim == 0:
			channeldata = [channeldata]
			
		self._channels = [Channel(line, self._events, self) for line in channeldata]
			
		# the cluster data 
		self._clusters = [Cluster(i, [channel for channel in self._channels if channel.cluster() == i]) for i in range(self._clustercount)]
		
		# the membrane object
		self.__membrane = Membrane(self)
			
	def domains(self):
		return self.property('Meta','domains').split(',')
		
	def domain(self,domain):
		return self._spatial[domain]
		
	def property(self, section, item):
		return self.config.get(section,item)
		
	def membrane(self):
		return self.__membrane
	
	def events(self):
		return self._events	
	
	def channels(self):
		return self._channels
		
	def channel(self,i):
		return self._channels[i]
		
	def clusters(self):
		return self._clusters
		
	def cluster(self, i):
		return self._clusters[i]
	
	def channelcount(self):
		return len(self._channels)
		
	def tmin(self):
		return self.domain('cytosol')['calcium'].tmin()
		
	def tmax(self):
		return self.domain('cytosol')['calcium'].tmax()