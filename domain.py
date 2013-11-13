import types
		
class Domain(dict):
	def __init__(self, name, simulation):
		super(Domain, self).__init__()
		# the name of the domain
		self.__name       = name	
		self.__simulation = simulation	
			
	def components(self):
		return self.keys();
	
	def name(self):
		return self.__name
		
	def add_component(self, name, component):
		self[name] = component
		self.__dict__[name] = component
