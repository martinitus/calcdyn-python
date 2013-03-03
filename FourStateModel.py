import StochasticModel
import numpy

class Reaction(StochasticModel.Reaction):
	def __init__(self, channel, target, environment):
		self.channel = channel		
		self.environment = environment
		
		self.__source = channel.state
		self.__target = target
		self.mtt = None
		if(self.__source == 'open'):
			if(self.__target == 'inhibited'):
				self.mtt = Channel.TOI
			elif(self.__target == 'active'):
				self.mtt = Channel.TOA
		elif(self.__source == 'resting'):
			if(self.__target == 'active'):
				self.mtt = Channel.TRA
			elif(self.__target == 'inhibited'):
				self.mtt = Channel.TRI
		elif(self.__source == 'inhibited'):
			if(self.__target == 'open'):
				self.mtt = Channel.TIO
			elif(self.__target == 'resting'):
				self.mtt = Channel.TIR
		elif(self.__source == 'active'):
			if(self.__target == 'open'):
				self.mtt = Channel.TAO
			elif(self.__target == 'resting'):
				self.mtt = Channel.TAR
		
		assert(self.mtt != None)
				
	def rate(self,t):
		return  1000./self.mtt(self.environment.calcium(self.channel, t))
		
	def source(self):
		return self.__source
		
	def target(self):
		return self.__target
	
	def react(self):
		self.channel.state=self.__target
	
	def reagent(self):
		return self.channel
	
	def __repr__(self):
		return self.source() + '->' + self.target()
		
class Channel(object):	
	
	k01 = 0.0162      # 1/(microM   ms)
	k12 = 0.027       # 1/(microM^2 ms) 
	k23 = 2.1651E-4   # 1/(microM^3 ms)
	k34 = 1.0         # 1/(microM^4 ms)
	k45 = 3.5935E-8   # 1/(microM^5 ms)

	tauO = 30         # ms

	k01h = 1.4127E-3  # 1/(microM   ms)
	k12h = 1.0        # 1/(microM^2 ms)
	k23h = 1.0        # 1/(microM^3 ms)
	k34h = 1.0        # 1/(microM^4 ms)
	k45h = 5.6297E-7  # 1/(microM^5 ms)

	CA = 0.484        # microM
	CO = 0.238        # microM
	CI = 6.5          # microM
				
	def __init__(self, id, state):
		self.state = state
		self.id = id
		
	def __repr__(self):
		return "[" + str(self.id) + ", " + self.state + "]"
	
	# return a list of propensity objects describing the propensities for the different transitions
	def reactions(self, environment):
		if   self.state == 'open':
			return [Reaction(self,'inhibited',environment),Reaction(self,'active',environment)]
		elif self.state == 'active':
			return [Reaction(self,'open',environment),     Reaction(self,'resting',environment)]
		elif self.state == 'resting':
			return [Reaction(self,'active',environment),   Reaction(self,'inhibited',environment)]
		elif self.state == 'inhibited':
			return [Reaction(self,'resting',environment),  Reaction(self,'open',environment)]
		raise Exception("Invalid state: " + self.state)
		
	# return the mean transition times for the different reactions in milliseconds
	@staticmethod
	def TOI(C):
		return numpy.power(C/Channel.CO,2) * (1./(Channel.k45 * numpy.power(C,5)) + 1./(Channel.k34 * numpy.power(C,4) ) + 1./(Channel.k23 * numpy.power(C,3)))
	@staticmethod
	def TIO(C):
		return Channel.TOI(C) * numpy.power(C/Channel.CI,5) * numpy.power(Channel.CO/C,2)
	@staticmethod
	def TRI(C):
		return 1./(Channel.k01h * C) + 1./(Channel.k12h * numpy.power(C,2)) + 1./(Channel.k23h * numpy.power(C,3)) + 1./(Channel.k34h * numpy.power(C,4)) + 1./(Channel.k45h * numpy.power(C,5))	
	@staticmethod
	def TIR(C):
		return  numpy.power(C/Channel.CI,5) * Channel.TRI(C)	
	@staticmethod
	def TOA(C):
		#if C.__class__.__name__ == 'ndarray':
		#	return numpy.ones(C.shape[0],dtype=numpy.float32)*Channel.tauO
		#else:
		return Channel.tauO
	@staticmethod
	def TAO(C):
		#if C.__class__.__name__ == 'ndarray':
		#	return numpy.ones(C.shape[0],dtype=numpy.float32)*Channel.tauO * numpy.power(Channel.CO/Channel.CA,2)
		#else:
		return Channel.tauO * numpy.power(Channel.CO/Channel.CA,2)
	@staticmethod
	def TRA(C):
		return 1./(Channel.k12 * numpy.power(C,2)) + 1./(Channel.k01 * C)
	@staticmethod
	def TAR(C):
		return numpy.power(C/Channel.CA,2) * Channel.TRA(C)


