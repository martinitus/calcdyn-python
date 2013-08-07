#import CalciumData
#import ModelData

import ConfigParser
import numpy

from timeline import TimeLine
from modeldata import EventData
from binarydata import Domain

import dyk

channelmodels = {}
channelmodels['DYKModel'] = dyk
#channelmodels['RyRModel'] = dyk


class Channel(object):
	''' 
		infos:       a tuple containing channel id, channel location and cluster id
		transitions: a list of channel transitions, the channel will extract its own transitions from the global transition list
	'''
	def __init__(self,infos,eventdata):
		#print infos
		self.__location = infos['location']
		self.__index = infos['id']
		self.__radius = infos['radius']
		self.__cluster = infos['cluster']
		self.__eventdata = eventdata
		
		self.__state   = None;
	
	def state(self):
		if self.__state == None:
			i     = self.__index
			data  = self.__eventdata._data
			
			frames  = data[data['chid'] == i]['t']
			states  = data[data['chid'] == i]['states'][:,i]
			
			#insert initial and final state
			states  = numpy.insert(states,0,data['states'][0][i],axis = 0)
			states  = numpy.append(states,  [data['states'][-1][i]],axis = 0)
			
			frames  = numpy.insert(frames,0,data['t'][0],axis = 0)
			frames  = numpy.append(frames,  data['t'][-1])
			
			self.__state = TimeLine(frames,states,interpolationorder = 'zero')
		return self.__state
		
	def open(self,t):
		return self.state()(t)[6]>=3
		
#	def calcium(self):
#		return self.__calcium
		
	def location(self):
		return self.__location
		
	def index(self):
		return self.__index
		
	def cluster(self):
		return self.__cluster


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
		
		self._model = channelmodels[modelname]
		self._events    = EventData(path,types = self._model.types(self._channelcount), states = self._model.states, channels = self._channelcount)
		
		# the channels 
		self._channels = [Channel(line, self._events) for line in numpy.genfromtxt(path + '/channels.csv',dtype=[('id', int), ('cluster', int), ('location', float, (3)), ('radius', float)])]
		
		# the spatial data
		self._spatial  = {}
		for domain in self.domains():
			print self.config.get(domain,'components').split(',')
			self._spatial[domain] = Domain(path, domain, components = self.config.get(domain,'components').split(','))
			
		# the cluster data 
		self._clusters = [Cluster(i, [channel for channel in self._channels if channel.cluster() == i]) for i in range(self._clustercount)]
			
	def domains(self):
		return self.property('Meta','domains').split(',')
		
	def domain(self,domain):
		return self._spatial[domain]
		
	def property(self, section, item):
		return self.config.get(section,item)
	
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
		return self._events.tmin()
		
	def tmin(self):
		return self._events.tmax()
