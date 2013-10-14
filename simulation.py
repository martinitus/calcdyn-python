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
		
		self._events    = EventData(path)
		
		try:
			# the spatial data
			self._spatial  = {}
			for domain in self.domains():
				self._spatial[domain] = Domain(path, domain, components = self.config.get(domain,'components').split(','))
		
			if not self._spatial.has_key('er'):
				self._spatial['er'] = FakeDomain('er',self)
		except:
			print "Warning: could not load spatial data"
			
		# the channels 
		channeldata = numpy.genfromtxt(path + '/channels.csv', dtype=[('id', int), ('cluster', int), ('location', float, (3)), ('radius', float)])
		
		# since genfromtxt returns 1d array for single line files we need to reshape
		if channeldata.ndim == 0:
			channeldata = [channeldata]

		
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
	
	def channelcount(self):
		return len(self._channels)
		
	def tmin(self):
		return self.domain('cytosol')['calcium'].tmin()
		
	def tmax(self):
		return self.domain('cytosol')['calcium'].tmax()