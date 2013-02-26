import os
import numpy
import csv
import ModelData
import scipy.integrate
import scipy.optimize
import math
import matplotlib.figure
import StringIO

class StochasticChannel(object):	
	
	def __init__(self, state, calcium, transition_callbacks = []):
		self.state = state
		self.transition_callbacks = transition_callbacks
		
		self.__reactions = Reactions(self)
		
		self.__lasteventtime = 0		
		
		# the history of calcium evolution since the last transition
		self.__calcium = TimeLine([0],[calcium])
		
		# propensity integration
		self.__cdf  = TimeLine([0],[0]); # cumulative density function
		
		# the state history of the channel
		self.__history = TimeLine([0],[state])
		
		# integration limit
		self.xi = numpy.random.uniform()


# a group of stochastic channels, the group object provides functionality for callback and stochastic channel evolution
class StochasticChannelGroup(object):
	
	def __init__(self, channels):
		self.channels = channels
		self.tstoch   = 0
		
	# drive the set of channels by an external environment
	def applytimestep(self, environment, t, tau):

		while(tstoch < t + tau):
      
			# always recalculate a0, since otherwise a possible changes of the propensities due to
			# changes of the calcium levels would be neglected in the stochastic time development

			# time difference between end of deterministic time step and stochastic time
			dt = t + tau - tstoch
			# mean time between stochastic time and end of deterministic time step
			tm = tstoch + 0.5*dt

			# for debug reasons...
			assert(t0 < tm && tm < t0+tau);
			if(not(t0 < tm and tm < t0+tau)) raise "Invalid time step!";

			# get propensities for mean time
			Propensities propensities(channelsetup, simulation, tm);
			i = propensities.total();
			// no stochastic event in time step
			if(g + dt*propensities.total() < xi):
				# increase accumulated "propensity" weighted with the time difference
				g = g  + dt*propensities.total();

				# update stochastic time
				tstoch = t0 + tau;
				#            std::cout << "g=" << g << " tstoch="<< tstoch<< " xi=" << xi << " prop=" << propensities.total() << std::endl;
				# return approximation for next event time
				return (xi-g)/propensities.transition();
			
			# stochastic event in time step
			else:
			
				# compute event time according to g + dt*a0 == xi
				tevent = (xi-g)/propensities.total() + tstoch;
				# determine next random number
				xi          = -std::log(getNonZeroRandomNumber());
				# reset accumulated propensities
				g           = 0.;
				# update stochastic time
				tstoch      = tevent;

				# recalculate propensities for event time
				Propensities propensities(channelsetup, simulation, tevent);

				# determine type and channel for next event
				std::pair<unsigned, Transition> transition = propensities.pick();

				# extract the channel model
				typename Traits::Model& model = channelsetup.getChannel(transition.first).getModel();

				# do the channel transition
				model.transition(tevent, transition.second);
			
		return (xi-g)/ Propensities(channelsetup, simulation, tstoch).transition();
		
	# drive this channel by an external signal (which itself might again depend on the modelstate)
	def drive(self, timeline, t0, tend, callbacks = []):
		# the integrand
		def integrand(tauprime, tt = None, propensities = None):
			if tauprime <= 0:
				return 0
			if tt+tauprime > tend:
				tauprime = tend-tt
			accum    = 0.
			exponent = 0.
			for a in propensities:
				prop  = a(tt+tauprime)
				propi = a.definitIntegral(tt,tauprime)
				assert(prop  >= 0)
				assert(propi >= 0)
				exponent = exponent + propi					
				accum    = accum    + prop
				#~ print 'prop',prop,'propi',propi
			return math.exp(-exponent) * accum
		
		# the equation which needs to be solved for the next transition time tau
		#~ def equation(tau, t = None, u = None, propensities = None):
			#~ #tau = tau[0]
			#~ if tau <= 0:
				#~ return tau - u
			#~ if t+tau > tend:
				#~ tau = tend-t
				#~ return tau
				#~ 
			#~ assert(t != None)
			#~ assert(u >= 0)
			#~ assert(propensities != None)			
					#~ 
			#~ integral = scipy.integrate.quad(lambda x: integrand(x,tt=t,propensities=propensities)[0], 0, tau)
			#~ if(integral[1] > 1E-9):
				#~ print 'Warning: Integrationerror:', integral[1]
			#~ print 'u=',u,'tau=',tau,'threshold=',(integral[0] - xi)
			#~ 
			#~ return integral[0] - u
			
		def solve(u,t,propensities,start,end):
			integral = scipy.integrate.quad(lambda x: integrand(x,tt=t,propensities=propensities), start, end)[0]
			#~ print 'u=',u,'t=',t,'start=',start,'end=',end,'i=',integral
			
			# we already found the solution
			if abs(integral - u) < 1E-9:
				return end
			
			# change of sign between start and end
			if integral > u:
				return scipy.optimize.brentq(
					lambda tau: scipy.integrate.quad(lambda x: integrand(x,tt=t,propensities=propensities), start, tau)[0] - u
					,start,end)
					
			# we reached end of timeline
			if end >= tend:
				raise Exception("tend reached, no more events!")
						
			# no change of sign, extend search range
			if integral < u:
				newend = end + min(tend-t,(end - start)*2)
				# keep searching for change of sign on the right hand side
				return solve(u - integral,t,propensities,start = end, end = newend)
						
			
			raise Exception("Should not reach here...")
			
			
		while (t0 < tend):
			# uniform random number
			u = numpy.random.rand(1)[0]
			#~ print 'u',u
			# get the next transition time
			props = self.propensities(timeline)
			
			# approximate next expected event time
			guess = -math.log(1-u)/(props[0](t0)+props[1](t0))
			# solve the equation
			#~ result = scipy.optimize.root(lambda tau: equation(tau,t=t0,xi=u,propensities=props),guess)
			#~ result,info,ierr,msg = scipy.optimize.fsolve(lambda tau: equation(tau,t=t0,xi=u,propensities=props),(guess),fprime = lambda tau: integrand(tau,tt=t0,propensities=props),full_output=True,maxfev = 1000)
			#result,info,ierr,msg = scipy.optimize.fsolve(lambda tau: equation(tau,t=t0,xi=u,propensities=props),(guess),full_output=True,maxfev = 1000)
			#result = scipy.optimize.brentq(lambda tau: equation(tau,t=t0,xi=u,propensities=props),0,tend-t0)
			#~ print 'guess:',guess
			result = solve(u,t0,props,0,guess)
			
			#~ print result
			#if ierr!=1:
			#	print 'Warning: fsolve did not converge: u:',u,'result:',result[0],'guess:',guess,'msg:',msg,'val:',info['fvec']
			
			tevent = t0 + result
			approx = t0 + guess
			
			calcium = timeline(tevent)
			
			source = self.state
			self.transition(calcium,tevent)
			target = self.state
						
			# call all callbacks
			for callback in callbacks:
				callback(t0, tevent, source, target, calcium, approx)
				
			# update loop variable
			t0 = tevent
	
	# return the mean transition times for the different reactions in milliseconds
	@staticmethod
	def TOI(C):
		return numpy.power(C/Model.CO,2) * (1./(Model.k45 * numpy.power(C,5)) + 1./(Model.k34 * numpy.power(C,4) ) + 1./(Model.k23 * numpy.power(C,3)))
	@staticmethod
	def TIO(C):
		return Model.TOI(C) * numpy.power(C/Model.CI,5) * numpy.power(Model.CO/C,2)
	@staticmethod
	def TRI(C):
		return 1./(Model.k01h * C) + 1./(Model.k12h * numpy.power(C,2)) + 1./(Model.k23h * numpy.power(C,3)) + 1./(Model.k34h * numpy.power(C,4)) + 1./(Model.k45h * numpy.power(C,5))	
	@staticmethod
	def TIR(C):
		return  numpy.power(C/Model.CI,5) * Model.TRI(C)	
	@staticmethod
	def TOA(C):
		if C.__class__.__name__ == 'ndarray':
			return numpy.ones(C.shape[0],dtype=numpy.float32)*Model.tauO
		else:
			return Model.tauO
	@staticmethod
	def TAO(C):
		if C.__class__.__name__ == 'ndarray':
			return numpy.ones(C.shape[0],dtype=numpy.float32)*Model.tauO * numpy.power(Model.CO/Model.CA,2)
		else:
			return Model.tauO * numpy.power(Model.CO/Model.CA,2)
	@staticmethod
	def TRA(C):
		return 1./(Model.k12 * numpy.power(C,2)) + 1./(Model.k01 * C)
	@staticmethod
	def TAR(C):
		return numpy.power(C/Model.CA,2) * Model.TRA(C)




#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | channels x state | channels x calciumlevel

#~ internally we drop some columns...
#~ time | {channel states} | {channel calciumlevels}



class Data(ModelData.Data):
	def __init__(self, path):
		#~  load csv file and set up channel locations in base class
		super(Data, self).__init__(path)
		
		self.states = {}
		self.states['resting']   = 0
		self.states['active']    = 1
		self.states['open']      = 2
		self.states['inhibited'] = 3
		#~ print self.channels()

	def resting(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['resting']).sum(1)
		
	def active(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['active']).sum(1)
		
	def inhibited(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['inhibited']).sum(1)
				
	def open(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['open']).sum(1)	
	
	def closed(self):
		return self.channelcount() - self.open();
		
	def _repr_svg_(self):
		foo = StateTimeLine(self)
		imgdata = StringIO.StringIO()
		foo.savefig(imgdata, format='svg')
		imgdata.seek(0)  # rewind the data
		return imgdata.getvalue()
	
	#~ return coordinates of channels in given state for given time
	def locations(self, time, state = 'open'):
		s = self.states[state]
		f = self.frame(time)
				
		#~ get the indices of channels in the given state
		indices = numpy.where(super(Data,self).states(time = time) == s)		
		return self.channels()[indices]
			
	# return an array of puffs
	# puffs[p,0] puff start
	# puffs[p,1] puff end
	# puffs[p,2] puff duration
	# puffs[p,3] puff peak open count
	def puffs(self,tolerance = 0):
		active    = False
		puffs = numpy.zeros([0,4],dtype = numpy.float32)
		puff     = 0

		for line in self.data:
			# a puff just started
			if not active and line[8] >= 1:				
				# we really have a new puff
				if puff == 0 or line[0] - puffs[-1,1] > tolerance:					
					puffs.resize([puff+1,4]) #resize output array
					puff = puff + 1          # increase puff counter
					puffs[-1,0] = line[0]    # store open time
					puffs[-1,3] = 1          # set puff peak level
					print 'puff started: t=',puffs[-1,0]
				# the last puff just continued within tolerance
				else:
					pass
					
			# a puff might terminate 
			elif active and line[8] == 0:
				puffs[-1,1] = line[0]                  # store or overwrite possible puff end time
				puffs[-1,2] = line[0] - puffs[-1,1]    # store or overwrite possible puff duration
				print 'puff terminated: t=', puffs[-1,1],'duration=',puffs[-1,1]-puffs[-1,0], 'peak=',puffs[-1,2]
					
			# save present active state for next event
			active = line[8] >= 1
			
			# update puff peak level
			if active:
				puffs[-1,3]  = max(puffs[-1,3], line[7])
		
		return puffs
	
	
	def transitionlist(self, channel):
		states = self.data[:,[0,1+channel]]
		#~ calculate the indices of the channels transitions
		transitions = numpy.where((numpy.roll(states[:,1], 1) - states[:,1]) != 0)
		print states[transitions]
		#~ array([0, 0, 1, 0, 0, 1, 0])
	
	
	def stateintervals(self, channel, statename,frames = False):
		#~ state = self.states[statename]
		
		#~ def cond_resting(data):
			#~ return data[:,1+channel] == 0
		#~ def cond_active(data):
			#~ return data[:,1+channel] == 1
		#~ def cond_open(data):
			#~ return data[:,1+channel] == 2
		#~ def cond_inhibited(data):
			#~ return data[:,1+channel] == 3
		#~ print "requested state",statename,"for channel",channel
		def cond(data, c, s):
			#~ print "requested state",s,"for channel",c
			return data[:,1+c] == s
		
		return self.intervalls(lambda data: cond(data,c=channel,s=self.states[statename]), frames)
		
		
		#~ if statename == 'resting':
			#~ return self.intervalls(cond_resting,frames)
		#~ elif statename == 'active':
			#~ return self.intervalls(cond_active,frames)
		#~ elif statename == 'open':
			#~ return self.intervalls(cond_open,frames)
		#~ elif statename == 'inhibited':
			#~ return self.intervalls(cond_inhibited,frames)			
		#~ else:
			#~ raise "Invalid State:" + statename
		
		#times  = self.data[:,0]
		#states = self.data[:,1+channel]
				
		#~ calculate the indices of the channels transitions
		#transitions = numpy.where((numpy.roll(states, 1) - states) != 0)[0]
		
		#~ filter out all transitions different than the requested state name,
		#~ and ignore state[0] and state[-1] because they will have corrupt time.
		#~ we need to add 1 since state[0] is ignored and all indices hence are shifted to the left
		#requested   = numpy.where(states[transitions][1:-1] == state)[0] + 1
		
		#~ create array with start end end times
		#result = numpy.array([times[transitions][requested], times[numpy.roll(transitions, -1)][requested]]).transpose()
		
		#return result
		
	def statedurationdistribution(self, statename):
		result = numpy.zeros(0,dtype = numpy.float32)
		
		for channel in range(self.channelcount()):			
			intervalls = self.stateintervals(channel,statename)
			#~ append this channels durations to the return array
			result = numpy.append(result,intervalls[:,1]-intervalls[:,0])
				
		#~ return the distribution in miliseconds		
		return result * 1000		
