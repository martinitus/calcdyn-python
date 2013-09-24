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

from IPython.core.pylabtools import print_figure

#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | channels x state | channels x calciumlevel
def types(channels):
	return [('t', float), ('chid', int), ('clid', int), ('tr','>S18'), ('noch',int),('nocl',int), ('states',int,(channels,8)), ('cy-calcium',float,(channels)), ('er-calcium',float,(channels))]

X000 = 0
X001 = 1
X010 = 2
X100 = 3
X011 = 4
X101 = 5
X110 = 6
X111 = 7


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