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
from channel import Channel
from cluster import Cluster

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
		self._clusters = [Cluster(i, [channel for channel in self._channels if channel.cluster() == i],self._events,self) for i in range(self._clustercount)]
		
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