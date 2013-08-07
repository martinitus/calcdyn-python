import os
import subprocess
import time
import shutil
import termcolor

from CalciumData import *

class EGTA:
	kmplus=6
	kmminus=1
	name = 'EGTA'

class BAPTA:
	kmplus=600
	kmminus=100	
	name = 'BAPTA'

class Parameters(dict):
			
	def write(self,path):
		if not os.path.exists(self['gridfile']):
			raise Exception('Cannot find gridfile: ' + self['gridfile'])
		
		if os.path.exists(os.path.join(path,'parameters.txt')):
			os.remove(os.path.join(path,'parameters.txt'))
				
		parameters = open(os.path.join(path,'parameters.txt'),'w')
			
		parameters.write(";this is an automatically created parameters file\n")
			
		parameters.write("[Time]'\n")
		parameters.write("simulation_intervall   = " + self["simulation_intervall"] + "\n")
		parameters.write("tau_channel_transition = 1E-8\n")
		parameters.write("time_ip3_max           = 0.0\n")
		parameters.write("channel_open_time      = 1E-7\n")
		parameters.write("channel_open_tau_restriction = 3.33E-8\n")		
		parameters.write("\n")
		
		parameters.write("[Solver]\n")
		parameters.write("type = BCGS_JACOBI\n")
		parameters.write("tolerance = 1E-3\n")
		parameters.write("reduction = 3E-5\n")
		parameters.write("use_last_k = false\n")
		parameters.write("amg_ssor_iterations = 20\n")
		parameters.write("jacobi_iterations   = 200\n")
		parameters.write("ilu0_iterations     = 200\n")
		parameters.write("verbosity = 0\n")
		parameters.write("factor = 0.90\n")
		parameters.write("facmax = 2.5\n")
		parameters.write("facmin = 0.1\n")
		parameters.write("taumax = 0.05\n")
		parameters.write("\n")
		parameters.write("[Parameters]\n")
		parameters.write("D0	      = 223E+6\n")
		parameters.write("D1          = 95.0E+6\n")
		parameters.write("D2          = 0.0\n")
		parameters.write("D3          = 20.0E+6\n")
		parameters.write("ER_CONC     = 700.0\n")
		parameters.write("c0          = 0.02\n")
		parameters.write("ip3         = 0.07\n")
		parameters.write("ip3max      = 0.07\n")
		parameters.write("Pp          = 40000\n")
		parameters.write("Pc          = " + self["Pc"] + "\n")
		parameters.write("Kd2         = 0.04\n")
		parameters.write("\n")
		parameters.write("kmplus      = " + str(self["buffer"].kmplus) + "\n")
		parameters.write("kmminus     = " + str(self["buffer"].kmminus) + "\n")
		parameters.write("ksplus      = 50.0\n")
		parameters.write("ksminus     = 100.0\n")
		parameters.write("kdplus      = 150\n")
		parameters.write("kdminus     = 300\n")
		parameters.write("\n")
		parameters.write("Bm          = " + self["Bm"] + "\n")
		parameters.write("Bs          = " + self["Bs"] + "\n")
		parameters.write("Bd          = 25.0\n")
		parameters.write("\n")					
		parameters.write("[DYKModel]\n")
		parameters.write("a1 = " + self["a1"] + "\n")
		parameters.write("a2  = 0.02\n")
		parameters.write("a3  = 0.2\n")
		parameters.write("a4  = 0.1\n")
		parameters.write("a5  = 100\n")
		parameters.write("dc1 = 0.001\n")
		parameters.write("dc2 = 78\n")
		parameters.write("dc3 = 2.0\n")
		parameters.write("dc4 = 0.039\n")
		parameters.write("dc5 = 0.25\n")					
		parameters.write("\n")
		parameters.write("[ChannelSetup]\n")
		parameters.write("center_x           = "+ str(self["center_x"]) + "\n")
		parameters.write("center_y           = "+ str(self["center_y"]) + "\n")
		parameters.write("channel_distance   = "+ str(self["channel_distance"]) + "\n")
		parameters.write("channel_radius     = 2.5\n")
		parameters.write("\n")
		parameters.write("[Grid]\n")
		parameters.write("filename  = " + self["gridfile"] + "\n")
		parameters.write("heap_size = 2000\n")
		parameters.write("\n")
		parameters.write("[Resume]\n")		
		parameters.write("resume = " + self["resume"] + "\n")		
		parameters.write("\n")
		parameters.write("[Output]\n")
		parameters.write("path        = ./\n")
		parameters.write("bin2d_out   = true\n")
		parameters.write("bin3d_out   = false\n")
		parameters.write("vtk_out     = true\n")
		parameters.write("vtk_interval= 1\n")
		parameters.write("vtk_4comp   = true\n")
		parameters.write("\n")
		parameters.write("[Verbosity]\n")
		parameters.write("config_set  = true\n")
		

'''
	The SimulationStatus class can provide fast informations about the current frame of the simulation directory given
'''
class SimulationStatus(object):
	def __init__(self, path):
		
		# find a datafile to open
		rank = ''		
		if os.path.exists(path + '/calcium.bin') and os.path.exists(path + '/coordinates.bin'):	
			rank = ''
		elif os.path.exists(path + '/calcium_rank_0.bin') and os.path.exists(path + '/coordinates_rank_0.bin'):
			rank = '_rank_0'
		
		try:
			# read the coordinates:				
			coordfile = open(path + '/coordinates' + rank + '.bin')		
			self.nodes = numpy.fromfile(coordfile, dtype=numpy.float32)	
			self.nodes = numpy.reshape(self.nodes,(self.nodes.size / 2,2))		
			coordfile.close()
					
			# open datafile and determine number of frames
			self.datafile = open(path + '/calcium' + rank + '.bin')
			self.frames = os.path.getsize(path + '/calcium' + rank + '.bin')/(self.nodes.shape[0]+1)/4
		#print 'Found', self.nodes.shape[0], 'nodes and', self.frames, 'frames in', 'calcium' + rank + '.bin'
		except:
			self.nodes = numpy.zeros((0,2), dtype = numpy.float32)
			self.frames = 0
		
	def data(self,f):
		if f >= self.frames or f < 0:
			raise Exception("frame " + str(f) + "not available")	
		self.datafile.seek(4*(self.nodes.shape[0]+1)*f + 4)		
		return numpy.fromfile(self.datafile, dtype=numpy.float32, count=self.nodes.shape[0])
				
	def time(self,f):
		if f >= self.frames or f < 0:
			raise Exception("frame " + str(f) + "not available")	
		self.datafile.seek(4*(self.nodes.shape[0]+1)*f)
		return numpy.fromfile(self.datafile, dtype=numpy.float32, count=1)[0]
		
'''
		The Simulation class provides basic informations about the status of a simulation directory
'''
class Simulation(object):
	def __init__(self,__path):
		self._path = __path
		
	# check if the simulation has been started at some point in the past
	def isStarted(self):
		try:			
			stage = self.currentStage()# if there is no stage this will fail
			return True
		except:
			return False
	
	# check if the simulation is currently running
	def isRunning(self):
		try:
			stage = self.currentStage()
			modtime = os.path.getmtime(os.path.join(self._path, stage + '/cout_rank_0.txt'))
			#print 'modime=',modtime
			#print 'time=', time.time()
			#print 'Time since last modification:', (time.time() - modtime)
			return time.time() - modtime < 120
		# if currentStage() fails, there is no simulationdata yet, hence it is not running
		except:
			return False
			
	
	def isFinished(self):
		try:
			stage = self.currentStage()
			status = SimulationStatus(os.path.join(self._path, stage))
			return status.time(status.frames-1) >= 200
		except:
			return False
	
	def printStatus(self, tail = False):
		path  = str(os.path.relpath(self._path))
		frame = 0
		time  = 0		
		if self.isStarted():
			stage = self.currentStage()
			
			status = SimulationStatus(os.path.join(self._path, stage))
			if status.frames > 0:
				time = status.time(status.frames-1)
			frame  = status.frames
			
			path = path + '/' + stage
			if self.isRunning():
				status = 'RUNNING'
				color  = 'blue'
			elif self.isFinished():
				status = 'FINISHED'
				color  = 'green'
			else:
				color  = 'yellow'
				status = 'WAITING'				
		else:
			color  = 'red'
			status = 'NOT STARTED'
		
		print termcolor.colored('{path:<60}: {status:<12} t = {time:>10}   frame = {frame:>5}'.format(path = path, status = status, frame = frame, time = time),color)
		
		if tail and self.isStarted():
			self.tailLast()
			
	def tailLast(self):
		stage = self.currentStage()
		subprocess.call(['tail', '-n 3', os.path.join(self._path,stage + '/cout_rank_0.txt')])				
		
	# get the folder the simulation last wrote to, or is currently writing to
	def currentStage(self):
		# check for running simulation		
		if os.path.exists(os.path.join(self._path,'1')):
			return '1'	
		
		# check for resume folders
		resume = 0
		while os.path.exists(os.path.join(self._path, 'resume' + str(resume + 1))):
			resume += 1
			
		if resume > 0:
			return 'resume' + str(resume)
		
		# check for start folder
		if os.path.exists(os.path.join(self._path, 'start')):
			return 'start'				
			
		raise Exception('Simulation ' + self._path + ' not yet started')
		
	#~ return the directory name for the next resume stage
	def nextStage(self):
		# check for start folder
		if not os.path.exists(os.path.join(self._path, 'start')):
			return 'start'			
				
		# check for resume folders
		resume = 1
		while os.path.exists(os.path.join(self._path, 'resume' + str(resume))):
			resume += 1
			
		return 'resume' + str(resume)			
		
	def cleanup(self):
		if self.isRunning():
			print str(os.path.relpath(self._path)), 'is currently running, skipping cleanup...'
			return
			
		if not self.isStarted():
			return
						
		stage = self.currentStage()		
		status = SimulationStatus(os.path.join(self._path, stage))
		
		if status.frames < 500:
			#~ if raw_input(os.path.relpath(os.path.join(self._path, stage)) + ' only contains ' +str(status.frames) + ' frames, delete? (y/n)') == 'y':
			print os.path.relpath(os.path.join(self._path, stage)),'only contains', status.frames,'frames!'
				#~ os.system('rm -rf ' + os.path.relpath(os.path.join(self._path, stage)))
				#~ subprocess.call(['rm', '-rf', os.path.relpath(os.path.join(self._path, stage))])				
				
		elif stage == '1':
			if raw_input('rename ' + os.path.join(self._path, stage) + ' to ' + self.nextStage() + '? (y/n)') == 'y':
				os.rename(os.path.join(self._path, stage),os.path.join(self._path, self.nextStage()))
			
			
		



'''
	A StartableSimulation object provides all necessary tools to run a new simulation, or to resume an old one
'''
class StartableSimulation(Simulation):
	def __init__(self, __path, __executable, __desc):
		# call base class constructor
		super(StartableSimulation, self).__init__(__path)
		self._executable = __executable
		self._description = __desc
	
	# creat status.bin file for the most recent simulation status
	def createLastStatusFile(self):		
		stage = self.currentStage()			
		status = 0
		while os.path.exists(os.path.join(self._path, stage + '/status_' + str(status+500) + '.bin')):
			status += 500				
				
		print ' Creating status symlink for', stage + '/status_' + str(status) + '.bin'
		
		try:
			os.remove(os.path.join(self._path,'status.bin'))
		except:
			pass
			
		os.symlink(os.path.join(self._path, stage, 'status_' + str(status) + '.bin'), os.path.join(self._path,'status.bin'))
		
	def createExecutableSoftLink(self):
		try:
			os.remove(os.path.join(self._path,'executable'))
		except:
			pass
			
		os.symlink(self._executable, os.path.join(self._path, 'executable'))
		
	def prepare(self, __parameters, verbose = True):
		if self.isRunning() or self.isFinished():
			return
			
		if not os.path.exists(self._path):
			os.makedirs(self._path)
		
		__parameters["resume"] = str(self.isStarted()).lower()
		__parameters.write(self._path)
		
		self.createExecutableSoftLink()
		
		if self.isStarted():
			if self.currentStage() == '1':
				print ' Renaming /1 to /' + self.nextStage()
				os.rename(os.path.join(self._path, '1'), os.path.join(self._path, self.nextStage()))
				
			self.createLastStatusFile()	
		
	def run(self, __parameters, niceness = 19, verbose = True):
		if self.isRunning():
			print 'Warning:', self._description, 'is already runnning, skipping....'
			return
			
		self.prepare(__parameters,verbose)		
		
		# cd to simulation directory
		_olddir = os.path.abspath(os.curdir)
		os.chdir(self._path)
		
		# run the job		
		subprocess.call(['nice', '-n',str(niceness),'./executable'])				
		
		# cd back to old directory
		os.chdir(_olddir)
	
	def queue(self, __parameters, verbose = True):
		if self.isRunning():
			print 'Warning:', self._description, 'is already runnning, skipping....'
			return
			
		self.prepare(__parameters, verbose)		
		# print sime info
		if self.isStarted():
			print ' Resuming',  os.path.relpath(self._path)
		else:
			print ' Starting', os.path.relpath(self._path)
			
		#submit the job	
		cwd = os.path.abspath(os.getcwd())
		os.chdir(self._path)				
		subprocess.call(['qsub', '-N', self._description[0:14], os.path.join(cwd, 'submit.sh')])				
		os.chdir(cwd)
		
		# node type selection is done in the submit.sh script
		#qsub -N icd60'_'$Bs'_'$Bm'_'$Pc -l nodes=1:ppn=1:X3300 $f/submit.sh
