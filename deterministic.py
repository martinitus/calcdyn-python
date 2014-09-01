
import numpy
import itertools

#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | channels x state | channels x calciumlevel

def types_binary(channels,clusters = 1):
    channels_per_cluster  = channels / clusters
    return [('t', numpy.double), ('chid', numpy.int32), ('clid', numpy.int32), ('noch',numpy.int32),('nocl',numpy.int32), ('states',numpy.bool,(clusters,channels_per_cluster))]
    

def loadtxt(filename,channels):
    tmp = numpy.loadtxt(filename, dtype = [('t', numpy.double), ('chid', numpy.int32), ('clid', numpy.int32),
                                           ('tr','>S8'), ('noch',numpy.int32),('nocl',numpy.int32),
                                           ('states',bool,(1,channels))])
                                   
    return numpy.fromiter(itertools.izip(
                                tmp['t'],tmp['chid'],tmp['clid'],tmp['noch'],tmp['nocl'],
                                tmp['states']), dtype = types_binary(channels))

def state_type():
    return (bool)

def open(data):
    '''
        If state is a onedimensional array, returns True or false,
        If state is multidimensional array of states, returns an array 
        with dimension reduced by one, masking the open channels as True, closed as False.
    '''
    print data
    print isinstance(data, numpy.recarray)
    print data.dtype
    if isinstance(data, numpy.recarray):
        tmp = data['states']
    else:
        tmp = data
    return tmp


def types(channels):
    return [('t', float), ('chid', int), ('clid', int), ('tr','>S18'), ('noch',int),('nocl',int), ('states',int,(channels,1)), ('cy-calcium',float,(channels)), ('er-calcium',float,(channels))]

states = {'open': {'name': 'open',      'condition': lambda x: x[0]==1,  'marker': 'o'},\
        'closed': {'name': 'closed',    'condition': lambda x: x[0]==0 , 'marker': 'x'}}


def overview(ax1,data):
    t=data.events()['t']
    noch=data.events()['noch']
    ax1.plot(t,noch,lw=1,c='black',label = '# open',drawstyle = 'steps-post')
    ax1.legend(loc='upper left')