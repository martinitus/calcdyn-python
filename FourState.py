import os
import numpy
import csv
import ModelData
import scipy.integrate
import scipy.optimize
import math
import matplotlib.figure

import StringIO

# Default plot for StateTimeLine
class StateTimeLine(matplotlib.figure.Figure):	
	def __init__(self, data, figsize=None, dpi=None, facecolor=None, edgecolor=None, frameon=True):
		super(StateTimeLine, self).__init__(figsize, dpi, facecolor, edgecolor, frameon)
		self.data  = data
			
		ax = self.add_subplot(1,1,1)
		ax.grid(True)
		ax.axhline(0, color='black', lw=2)
		ax.axis([1, self.data.tmax(), 0, 11])
		ax.set_ylabel('channels')
		ax.set_xlabel('time [s]')
	
		#~ Plot the state evolution to axes object	
		a,= ax.plot(data.times(),data.active(),   c='green',drawstyle='steps-post',lw=2)
		o,= ax.plot(data.times(),data.open(),     c='red',  drawstyle='steps-post',lw=2)
		i,= ax.plot(data.times(),data.inhibited(),c='cyan', drawstyle='steps-post',lw=2)	
		r,= ax.plot(data.times(),data.resting(),  c='blue', drawstyle='steps-post',lw=2)

		ax.legend([r, a, o, i], ["resting","active","open","inhibited"], loc=2)
		matplotlib.backends.backend_agg.FigureCanvasAgg(self)

# callable object for stochastic evolution providing the propensity depending on time
class Propensity(object):
	# transition should be a function taking calciumlevel as argument and returning the mean transition time
	def __init__(self, calciume_timeline, mean_transition_time_dependency):
		self.timeline= calciume_timeline
		self.mttd    = mean_transition_time_dependency
		
	# return the propensity for time t
	def __call__(self,t):
		# since mttd returns the mean transition time in milliseconds we need to convert to a rate in seconds
		#~ print self.timeline(t)
		#~ ca = self.timeline(t)
		#~ print self.mttd(ca)
		#~ mt = self.mttd(ca)
		#~ print 1000./mt
		#~ return 1
		return 1000./self.mttd(self.timeline(t))
	
	# return the value of the primitive of the propensity for time t, i.e the definite integral from tmin to t assuming a(tmin) = 0
	#def primitive(self,t):
	#	return scipy.integrate.quad(self,self.timeline.tmin(),t)
	
	# return the definit integral of a(t) from t to t+tau
	def definitIntegral(self,t,tau):
		#~ return tau
		return scipy.integrate.quad(self,t,t+tau)[0]
		
class Model(object):	
	
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
				
	def __init__(self,state, transition_callbacks = []):
		self.state = state
		self.transition_callbacks = transition_callbacks
	
	# return a list of propensity objects describing the propensities for the different transitions
	def propensities(self,timeline):
		if   self.state == 'open':
			return [Propensity(timeline,Model.TOI),Propensity(timeline,Model.TOA)]
		elif self.state == 'active':
			return [Propensity(timeline,Model.TAO),Propensity(timeline,Model.TAR)]
		elif self.state == 'resting':
			return [Propensity(timeline,Model.TRA),Propensity(timeline,Model.TRI)]
		elif self.state == 'inhibited':
			return [Propensity(timeline,Model.TIR),Propensity(timeline,Model.TIO)]
		raise Exception("Invalid state: " + self.state)
		
	# perform a random transition according to propensities governed by the given calcium level
	def transition(self, calcium, time):
		# pick random uniform number
		u  = numpy.random.rand(1)[0]
		
		#store old state for the callback
		oldstate = self.state
		
		if self.state == 'open':
			p1 = 1000. / Model.TOI(calcium)
			p2 = 1000. / Model.TOA(calcium)
			if u * (p1+p2) <= p1:
				self.state = 'inhibited'
			else:
				self.state = 'active'				
		elif self.state == 'active':
			p1 = 1000. / Model.TAO(calcium)
			p2 = 1000. / Model.TAR(calcium)
			if u * (p1+p2) <= p1:
				self.state = 'open'
			else:
				self.state = 'resting'				
		elif self.state == 'resting':
			p1 = 1000. / Model.TRA(calcium)
			p2 = 1000. / Model.TRI(calcium)
			if u * (p1+p2) <= p1:
				self.state = 'active'
			else:
				self.state = 'inhibited'			
		elif self.state == 'inhibited':
			p1 = 1000. / Model.TIR(calcium)
			p2 = 1000. / Model.TIO(calcium)
			#~ print 'p1',p1,'p2',p2,'u*(p1+p2)',u*(p1+p2)
			if u * (p1+p2) <= p1:
				self.state = 'resting'
			else:
				self.state = 'open'					
		else:
			raise Exception("Invalid state: " + self.state)
		
		# call all callbacks
		for callback in self.transition_callbacks:
				callback(time,oldstate,self.state)
			
		return (time,oldstate,self.state)
		
	def addTransitionCallback(self,callback):
		self.transition_callbacks.append(callback)			
		
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
