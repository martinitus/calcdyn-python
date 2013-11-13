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

def types_ascii(channels):
	return [('t', numpy.double), ('chid', numpy.int32), ('clid', numpy.int32), ('tr','>S18'), ('noch',numpy.int32),('nocl',numpy.int32), ('states',numpy.byte,(channels,8))] #''

def types_binary(channels):
	return [('t', numpy.double), ('chid', numpy.int32), ('clid', numpy.int32), ('noch',numpy.int32),('nocl',numpy.int32), ('states',numpy.byte,(channels,8))]


X000 = 0
X001 = 1
X010 = 2
X100 = 3
X011 = 4
X101 = 5
X110 = 6
X111 = 7

# return the number of open channels for each row in data
def noch_single(data):
	return (data['states'][:,X110]>=3).sum(axis = 1)
		
def available_single(data):
	withip3    = data['states'][:,[X100,X110,X111,X101]].sum(axis=1)
	# the number of totally available channels (i.e. channels that have enough ip3 bound)
	return (withip3>=3).sum()
	
def withip3_single(data):
	return available_single(data)

def inhibited_single(data):
	return data['states'][:,[X101,X111]].sum(axis=1)


# return the number of open channels for each row in data
def noch(data):
	tmp = data['states'] if hasattr(data,'states') else data;
	return (tmp[:,:,X110]>=3).sum(axis = 1)
	
# returns number of channels that have more then two subunits with ip3 bound
def available(data):
	tmp = data['states'] if hasattr(data,'states') else data;
	return (tmp[:,:,[X100,X110,X111,X101]].sum(axis=2) > 2).sum(axis = 1)
	
# return the number of channels that have 3 or more subunits in active or open state (i.e. that have ip3 and activating calcium bound) but are NOT open yet
def active(data):
	tmp = data['states'] if hasattr(data,'states') else data;
	return (numpy.logical_and(tmp[:,:,[X100,X110]].sum(axis=2) >= 3, tmp[:,:,X110] < 3 )).sum(axis = 1)
	
# returns number of subunits that have no ip3 bound
def noip3(data):
	tmp = data['states'] if hasattr(data,'states') else data;
	return tmp[:,:,[X000,X010,X011,X001]].sum(axis=1).sum(axis = 1)

# return the number of subunits that have ip3 bound
def withip3(data):
	tmp = data['states'] if hasattr(data,'states') else data;
	return tmp[:,:,[X100,X110,X111,X101]].sum(axis=1).sum(axis = 1)

# return number of subunits that are inhibited, but still have ip3 bound
def inhibited(data):
	tmp = data['states'] if hasattr(data,'states') else data;
	return tmp[:,:,[X101,X111]].sum(axis=1).sum(axis = 1)

class Rates(object):
	# the K_i are identical to the d_i or dc_i. they are all defined as b_i/a_i
	def __init__(self,config):
		self.a1 = config.getfloat('DYKModel','a1')
		self.a2 = config.getfloat('DYKModel','a2')
		self.a3 = config.getfloat('DYKModel','a3')
		self.a4 = config.getfloat('DYKModel','a4')
		self.a5 = config.getfloat('DYKModel','a5')
		
		try:
			self.b1 = config.getfloat('DYKModel','b1'); self.d1 = self.b1/self.a1;
			self.b2 = config.getfloat('DYKModel','b2'); self.d2 = self.b2/self.a2;
			self.b3 = config.getfloat('DYKModel','b3'); self.d3 = self.b3/self.a3;
			self.b4 = config.getfloat('DYKModel','b4'); self.d4 = self.b4/self.a4;
			self.b5 = config.getfloat('DYKModel','b5'); self.d5 = self.b5/self.a5;
		except NoOptionError:
			self.d1 = config.getfloat('DYKModel','dc1'); self.b1 = self.d1*self.a1;
			self.d2 = config.getfloat('DYKModel','dc2'); self.b2 = self.d2*self.a2;
			self.d3 = config.getfloat('DYKModel','dc3'); self.b3 = self.d3*self.a3;
			self.d4 = config.getfloat('DYKModel','dc4'); self.b4 = self.d4*self.a4;
			self.d5 = config.getfloat('DYKModel','dc5'); self.b5 = self.d5*self.a5;
			
	def __repr__(self):
		fs = "a{i:0d} = {a:6.2f}, b{i:0d} = {b:8.5f} => d{i:0d} = {d:6.3f}"
		return fs.format(i=1,a=self.a1,b=self.b1,d=self.d1)+"\n"\
			  +fs.format(i=2,a=self.a2,b=self.b2,d=self.d2)+"\n"\
			  +fs.format(i=3,a=self.a3,b=self.b3,d=self.d3)+"\n"\
			  +fs.format(i=4,a=self.a4,b=self.b4,d=self.d4)+"\n"\
			  +fs.format(i=5,a=self.a5,b=self.b5,d=self.d5)
		
class SteadyState(Rates):
	def __init__(self, config):
		super(SteadyState, self).__init__(config)
	
	def Z(self,ca,ip3 = None):
		ip3 = self.ip3 if ip3 == None else ip3
		return 1 + ca/self.d4 + ca/self.d5 + ca*ca/(self.d4*self.d5) + ip3/self.d1 + ip3*ca/(self.d1*self.d2) + ip3*ca/(self.d1*self.d5) + ip3*ca*ca/(self.d1*self.d2*self.d5)
	
	def w110(self,ca,ip3=None):
		ip3 = self.ip3 if ip3 == None else ip3
		return ip3*ca/(self.d1*self.d5*self.Z(ca,ip3))
		
	def Popen(self,ca,ip3 = None):
		ip3 = self.ip3 if ip3 == None else ip3
		w110 = self.w110(ca,ip3)
		return numpy.power(w110,4) + 4* numpy.power(w110,3)*(1-w110)

class Propensities(object):
	def __init__(self, state, cc):
		
		assert(sum(state) == 4)
		self._units = 4
		
		p = numpy.ndarray(shape = (8,8))
		p[:,:] = -1
		
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