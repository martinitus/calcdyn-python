
#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | channels x state | channels x calciumlevel

def types(channels):
	return [('t', float), ('chid', int), ('clid', int), ('tr','>S18'), ('noch',int),('nocl',int), ('states',int,(channels,1)), ('cy-calcium',float,(channels)), ('er-calcium',float,(channels))]

states = {'open': {'name': 'open',      'condition': lambda x: x[0]==1,  'marker': 'o'},\
	      'closed': {'name': 'closed',    'condition': lambda x: x[0]==0 , 'marker': 'x'}}

