
import numpy
import itertools

#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | channels x state | channels x calciumlevel

def types_binary(channels):
    return [('t', numpy.double), ('chid', numpy.int32), ('clid', numpy.int32), ('noch',numpy.int32),('nocl',numpy.int32), ('states',numpy.bool,(channels))] #''
    

def loadtxt(filename,channels):
    tmp = numpy.loadtxt(filename, dtype = [('t', numpy.double), ('chid', numpy.int32), ('clid', numpy.int32),
                                           ('tr','>S8'), ('noch',numpy.int32),('nocl',numpy.int32),
                                           ('states','>S6',(channels))])
                                   
    return numpy.fromiter(itertools.izip(tmp['t'],tmp['chid'],tmp['clid'],tmp['noch'],tmp['nocl'],[[s == 'open' for s in f] for f in tmp['states']]),
                                     dtype = types_binary(channels))

def state_type():
    return ('>S6')

def open(state):
    return state == 'open'


def types(channels):
	return [('t', float), ('chid', int), ('clid', int), ('tr','>S18'), ('noch',int),('nocl',int), ('states',int,(channels,1)), ('cy-calcium',float,(channels)), ('er-calcium',float,(channels))]

states = {'open': {'name': 'open',      'condition': lambda x: x[0]==1,  'marker': 'o'},\
	      'closed': {'name': 'closed',    'condition': lambda x: x[0]==0 , 'marker': 'x'}}


def overview(ax1,data):
    t=data.events()['t']
    noch=data.events()['noch']
    ax1.plot(t,noch,lw=1,c='black',label = '# open',drawstyle = 'steps-post')
    ax1.legend(loc='upper left')