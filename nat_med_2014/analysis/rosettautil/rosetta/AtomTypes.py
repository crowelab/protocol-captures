import os

class AtomType:
	def __init__(self,name,atom):
		self.attributes = set()
		self.name = name
		self.atom = None
		self.lk_properties = {}
		
	def add_lk_property(self,property_name,value):
		self.lk_properties[property_name] = float(value)
	
	def add_attribute(self,attribute):
		self.attributes.add(attribute)
	
	def get_property(self,property):
		return self.lk_properties[property]
		
	def has_attribute(self,attribute):
		return attribute in self.attributes
	

class AtomTypeSet:
	def __init__(self,database_path,tag):
		
		self.tag = tag
		self.atom_type_list = {}
		
		atom_path = database_path+"/chemical/atom_type_sets/"+tag+"/atom_properties.txt"
		property_lines = open(atom_path,'r').readlines()
		
		#first line is a header
		header = property_lines[0].split()
		
		for line in property_lines[1:]:
			line = line.split()
			if len(line) <= 2 or line[0][0] == "#":
				continue

			last_field = 0
			comment_exists = False
			#find the location of the comment
			for index,field in enumerate(line):
				if field[0] == "#":
					last_field = index+1
					comment_exists = True
					break
			#if there is no comment read every field
			if not comment_exists:
				last_field = len(line)

			new_atom_type = AtomType(line[0],line[1])
			#loop through non-comment fields, reading parameters and attributes
			for index,field in enumerate(line[2:last_field]):
				if index+2 < len(header):
					description = header[index+2]
					new_atom_type.add_lk_property(description,field)
				else:
					new_atom_type.add_attribute(field)
			
			self.atom_type_list[line[0]] = new_atom_type
	
	def get_atom_type(self,atom_type_name):
		return self.atom_type_list[atom_type_name]
		

		