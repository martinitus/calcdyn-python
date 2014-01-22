import numpy
import scipy.integrate
import scipy.optimize
import math


class Reaction(object):
	def react(self):
		raise NotImplementedError("Please Implement this method")
	def rate(self,t):
		raise NotImplementedError("Please Implement this method")
	def source(self):
		raise NotImplementedError("Please Implement this method")
	def target(self):
		raise NotImplementedError("Please Implement this method")
	def reagent(self):
		raise NotImplementedError("Please Implement this method")

class Reagent(object):			
	def reactions(self,environment):
		raise NotImplementedError("Please Implement this method")

# a group of stochastic channels, the group object provides functionality for callback and stochastic channel evolution
class ChannelGroup(object):
	
	# environment has to provide the determinisitic evolution of all the variables the stochastic implementation depends on
	def __init__(self, channels, transition_callbacks = [], xi = None, environment = None):
		self.callbacks   = transition_callbacks
		self.environment = environment
		self.channels    = channels
		self.reactions   = [reaction for channel in self.channels for reaction in channel.reactions(environment)]
		self.tstoch   = 0.0
		self.xi       = (xi if xi != None else numpy.random.uniform())
		
		# accumulation variables
		#self.cdf0  = 0 # cdf at tstoch
		self.f0    = 0 # f at tstoch
	
	# pd for any reaction
	def rho(self,t):		
		return self.a(t) * self.p0(t)
		
	# probability that no reaction has occured up to t
	def p0(self,t):
		return numpy.exp(-self.f(t))
	
	# sum over all rates
	def a(self,t):
		return sum([reaction.rate(t) for reaction in self.reactions])
	
	# sum over all f
	def f(self,t):
		assert(t>=self.tstoch)
		dt  = t-self.tstoch
		#if(dt>1E-4):
		#	inc = scipy.integrate.quad(self.a,self.tstoch,t)
		#else:
		tm  = self.tstoch + dt/2
		inc = [dt*self.a(tm)]
		#print inc
		return self.f0 + inc[0]
		
	# cdf for any reaction
	def cdf(self,t):
		return 1-self.p0(t)
		
	# drive the set of channels by an external environment
	def applytimestep(self, t, tau):
		assert(self.environment)
		assert(self.tstoch - t <= 1E-8)
		
		while(self.tstoch < t + tau):
			
			# update set of possible reactions
			self.reactions = [reaction for channel in self.channels for reaction in channel.reactions(self.environment)]
            
			# no stochastic event in time step
			if(self.cdf(t+tau) < self.xi):
				# update accumulation variables and stochastic time
				self.f0     = self.f(t+tau)
				#self.cdf0   = self.cdf(t+tau)
				self.tstoch = t + tau;
				
				# return approximation for next event time
				#return (xi-g)/propensities.transition();
				return
			
			# stochastic event in time step
			else:
				# compute event time according to self.cdf(tevent) == self.xi
				tevent = scipy.optimize.brentq(lambda t: self.cdf(t)-self.xi,self.tstoch,t+tau)
				assert(tevent < t+tau)
				# determine next random number
				self.xi          = numpy.random.uniform()
				# reset accumulated variables
				self.f0   = 0
				#self.cdf0 = 0
				# update stochastic time
				self.tstoch      = tevent;
				
				# choose and perform the reaction that should happen
				reaction = self.choosereaction(tevent)
				reaction.react()
				
				# call all callbacks
				for callback in self.callbacks:
					callback(tevent,reaction)
			
		return #(xi-g)/ Propensities(channelsetup, simulation, tstoch).transition();
	
	# choose a reaction governed by the probability distribution defined by the individual reaction rates
	def choosereaction(self,t):
		# the total reaction rate
		a0 = self.a(t)
		
		# pick a random number in [0,a0]
		rnd = a0*numpy.random.uniform()
		
		# accumulation
		accum = 0
		
		for reaction in self.reactions:
			a = reaction.rate(t)
			if accum < rnd and rnd < accum + a:
				return reaction
			else:
				accum = accum + a
		
		raise "Never should reach here..."
