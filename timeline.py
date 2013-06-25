import scipy.interpolate
import StringIO
import matplotlib.figure
import matplotlib.pyplot as plt
import copy
import numpy

from IPython.core.pylabtools import print_figure
from IPython.display import Image, SVG, Math

# provide a continuous time evolution of a discrete variable
class TimeLine(object):
	def __init__(self,frames,data,t0 = None,tend = None, desc="", ylabel="",yunit="",interpolationorder = 'linear'):
		self.frames = frames
		self.data   = data
		self.ylabel = ylabel
		self.yunit  = yunit
		self.desc   = desc
		self.interpolationorder = interpolationorder
		
		self.t0   = frames[0] if  t0   == None else t0
		self.tend = frames[-1] if tend == None else tend
				
		self.minframe = self.frame(self.tmin())
		self.maxframe = self.frame(self.tmax())
		assert(frames[0]<=self.t0)
		assert(self.tend     <=frames[-1])
		self.interp = scipy.interpolate.interp1d(self.frames,self.data,copy = False,kind = interpolationorder)
		
		
	def __call__(self,t):
		return self.interp(t)
			
	def subrange(self, subrange):
		return TimeLine(self.frames,self.data,t0 = subrange[0],tend = subrange[1],desc =self.desc, ylabel = self.ylabel, yunit = self.yunit)
	
	def tmin(self):
		return self.t0
	
	def tmax(self):
		return self.tend	
		
	def _repr_svg_(self):
		fig,ax=plt.subplots(figsize=(10,3))
		#ax = fig.add_subplot(111)
		self.plot(ax)
		data = print_figure(fig,'svg')
		plt.close(fig)
		return data.decode('utf-8')
		#output = StringIO.StringIO()		
		#fig.savefig(output,format='svg')
		#output.seek(0)
		#return output.getvalue()
		
	def plot(self, axes, trange = None, yrange = None,c='black',lw = 0.5,yscale = 'linear'):
		trange = [self.tmin(),self.tmax()] if trange == None else trange
	
		
		axes.set_xlim(trange)
		axes.set_ylabel(self.ylabel + ("" if self.yunit == "" else " ["+self.yunit+"]"))
		axes.set_xlabel('time [s]')
		axes.plot(self.frames,self.data, c=c,lw=lw,label = self.desc, drawstyle = 'steps-post' if self.interpolationorder == 'zero' else 'default')
		axes.set_ylim(yrange)
		axes.set_yscale(yscale)
		#axes.set_ylim(yrange)
		
	def integrate(self, range = None):
		f0 = self.frame(range[0])+1
		f1 = self.frame(range[1])
		dt0 = self.frames[f0]-range[0]
		dt1 = range[1] - self.frames[f1]
		assert(dt0 >= 0)
		assert(dt1 >= 0)
		assert(range[0]+dt0/2 >=self.tmin())
		assert(range[1]+dt1/2 <=self.tmax())
		return numpy.trapz(self.data[f0:f1],self.frames[f0:f1]) + dt0*self(range[0]+dt0/2) + dt1*self(range[1]+dt1/2)
	
	#~ return true if the given condition is true during the whole period
	def all(self,condition, range = None):
		if range == None:
			return condition(self.data).all()
		else:
			f0 = self.frame(range[0]) + 1
			f1 = self.frame(range[1])
			return condition(self.data[f0:f1]).all() and condition(self(range[0])) and condition(self(range[1]))
	
	#~ return true if the given condition is true for any timepoint in the given period
	def any(self,condition, range = None):
		if range == None:
			return condition(self.data).any()
		else:
			f0 = self.frame(range[0]) + 1
			f1 = self.frame(range[1])
			return condition(self.data[f0:f1]).any() or condition(self(range[0])) or condition(self(range[1]))
		
	#~ return the smallest frame number bevore the given time	
	def frame(self, time, fmin = 0, fmax = -1):
		if fmax == -1:
			fmax = len(self.frames)-1
			
		#~ print 'time',time,'fmin',fmin,'fmax',fmax
		#abort recursion
		if fmax - fmin == 1:
			return fmin
			
		f = fmin + (fmax-fmin)/2
		#search in left hand side
		if self.frames[f] > time:
			#~ print 'lhs'
			return self.frame(time,fmin,f)
		#search in right hand side
		if self.frames[f] <= time:
			#~ print 'rhs'
			return self.frame(time,f,fmax)
		
		raise Exception("should never be reached...")
	
