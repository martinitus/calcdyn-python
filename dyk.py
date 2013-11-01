import os
import numpy
import csv
#import ModelData
import scipy.integrate
import scipy.optimize
import math
import matplotlib.figure
import StringIO
import matplotlib.pyplot as plt
import matplotlib

from ConfigParser import NoOptionError
from IPython.core.pylabtools import print_figure

#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | channels x state | channels x calciumlevel
def types(channels):
	return [('t', float), ('chid', int), ('clid', int), ('tr','>S18'), ('noch',int),('nocl',int), ('states',int,(channels,8)), ('cy-calcium',float,(channels)), ('er-calcium',float,(channels))]

	ct = '<i'+str(int(numpy.floor(numpy.log10(channels)+1)))
	return [('t', float), ('chid', ct), ('clid', ct), ('tr','>S18'), ('noch',ct),('nocl',ct), ('states','<i1',(channels,8)), ('cy-calcium',float,(channels)), ('er-calcium',float,(channels))]

X000 = 0
X001 = 1
X010 = 2
X100 = 3
X011 = 4
X101 = 5
X110 = 6
X111 = 7

class Rates(object):
	pass
	
ip3   = None
rates = Rates()

def loadrates(config):
	global rates
	rates.a1 = config.getfloat('DYKModel','a1')
	rates.a2 = config.getfloat('DYKModel','a2')
	rates.a3 = config.getfloat('DYKModel','a3')
	rates.a4 = config.getfloat('DYKModel','a4')
	rates.a5 = config.getfloat('DYKModel','a5')
	
	try:
		rates.b1 = config.getfloat('DYKModel','b1')
		rates.b2 = config.getfloat('DYKModel','b2')
		rates.b3 = config.getfloat('DYKModel','b3')
		rates.b4 = config.getfloat('DYKModel','b4')
		rates.b5 = config.getfloat('DYKModel','b5')
	except NoOptionError:
		rates.b1 = rates.a1 * config.getfloat('DYKModel','dc1')
		rates.b2 = rates.a2 * config.getfloat('DYKModel','dc2')
		rates.b3 = rates.a3 * config.getfloat('DYKModel','dc3')
		rates.b4 = rates.a4 * config.getfloat('DYKModel','dc4')
		rates.b5 = rates.a5 * config.getfloat('DYKModel','dc5')
		
	global ip3
	ip3 = config.getfloat('DYKModel','ip3')


class Propensities(object):
	def __init__(self, state, cc):
		
		assert(sum(state) == 4)
		self._units = 4
		
		p = numpy.ndarray(shape = (8,8))
		p[:,:] = -1
		
		global rates
		global ip3
		
		p[X000][X010] = state[X000] * rates.a5 * cc;  # to X010 activating calcium bind
		p[X000][X001] = state[X000] * rates.a4 * cc;  # to X001 inhibiting calcium bind
		p[X000][X100] = state[X000] * rates.a1 * ip3; # to X100 ip3 bind

		p[X001][X011] = state[X001] * rates.a5 * cc;  # to X011 activating calcium bind
		p[X001][X000] = state[X001] * rates.b4;       # to X000 inhibiting calcium loss
		p[X001][X101] = state[X001] * rates.a3 * ip3; # to X101 ip3 bind

		p[X010][X000] = state[X010] * rates.b5;       # to X000 activating calcium loss
		p[X010][X011] = state[X010] * rates.a4 * cc;  # to X011 inhibiting calcium bind
		p[X010][X110] = state[X010] * rates.a1 * ip3; # to X110 ip3 bind

		p[X100][X110] = state[X100] * rates.a5 * cc;  # to X110 activating calcium bind
		p[X100][X101] = state[X100] * rates.a2 * cc;  # to X101 inhibiting calcium bind
		p[X100][X000] = state[X100] * rates.b1;       # to X000 ip3 loss

		p[X011][X001] = state[X011] * rates.b5;       # to X001 activating calcium loss
		p[X011][X010] = state[X011] * rates.b4;       # to X010 inhibiting calcium loss
		p[X011][X111] = state[X011] * rates.a3 * ip3; # to X111 ip3 bind

		p[X101][X111] = state[X101] * rates.a5 * cc;  # to X111 activating calcium bind
		p[X101][X100] = state[X101] * rates.b2;       # to X100 inhibiting calcium loss
		p[X101][X001] = state[X101] * rates.b3;       # to X001 ip3 loss

		p[X110][X100] = state[X110] * rates.b5;       # to X100 activating calcium loss
		p[X110][X111] = state[X110] * rates.a2 * cc;  # to X111 inhibiting calcium bind
		p[X110][X010] = state[X110] * rates.b1;       # to X010 ip3 loss

		p[X111][X101] = state[X111] * rates.b5;       # to X101 activating calcium loss
		p[X111][X110] = state[X111] * rates.b2;       # to X110 inhibiting calcium loss
		p[X111][X011] = state[X111] * rates.b3;       # to X011 ip3 loss
		self._p = p
		
	def __call__(self, source, target, add = False):
		assert(self._p[source][target] != -1)
		if add:
			return (self._p[source][target] - self._p[target][source]) / self._units
		else:
			return (self._p[source][target] / self._units)

	def __add__(self, other):
		p = Propensities(numpy.ndarray(shape = 8),0,0)
		p._p = self._p + other._p
		p._units = self._units + other._units
		return p
	

class Transition(object):
	def __init__(self, su1, su2):
		self.su1 = su1
		self.su2 = su2

	# return the flux between the two states
	# for net flux from su1 to su2 result is positive
	# for net flux from su2 to su1 result is negative
	def flux(self):
		return self.su1.flux(self.su2) - self.su2.flux(self.su1)


class SubUnit(object):
	
	def __init__(self, act, inhib, ip3):
		pass
	
	
	# return list of connected subunits (always 3 items)
	def targets(self):
		return []
	
	# return the flux to given target (always >=0)
	def flux(self, target):
		return 0
	
	# the negative parts of the state equation
	def outflux(self):
		return sum([self.flux(su) for su in self.targets()])
	
	# the positive parts of the state equation
	def influx(self):
		return sum([su.flux(self) for su in self.targets()])
	
	# return the total time derivative of the state counter
	def change(self):
		return self.influx() - self.outflux()


states = {'open': {'name': 'open',      'condition': lambda x: x[X110]>=3,'marker': 'o'},\
	      'closed': {'name': 'closed',    'condition': lambda x: x[X110]<3 , 'marker': '+'},\
	      'inhibited': {'name': 'inhibited', 'condition': lambda x: x[X101]+x[X111] + x[X001]+x[X011]>=2, 'marker':'x'},\
	      'noip3':{'name': 'noip3',     'condition': lambda x: x[X000]+x[X001] + x[X011]+x[X010]>=2, 'marker':'*'}}


'''
class Data(ModelData.Data):
	def __init__(self, path, channels):
		
		# define the set of data types to import from the csv
		
		
		#~  load csv file and set up channel locations in base class
		super(Data, self).__init__(path, types)

	def open(self):
		return self.observe(lambda x: x['noch'],desc = 'open')
		
	def _repr_svg_(self):
		fig,ax=plt.subplots(figsize=(10,3))
		
		ax.grid(True)
		ax.axhline(0, color='black', lw=2)
		ax.set_ylabel('channels')
		ax.set_xlabel('time [s]')
	
		#~ Plot the state evolution to axes object			
		self.open().plot(     ax,  c='red',  lw=2)
		ax.legend(loc=2)
		
		data = print_figure(fig,'svg')
		plt.close(fig)
		return data.decode('utf-8')
'''