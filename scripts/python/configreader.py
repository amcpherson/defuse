import re
import sys

class reader:
	def __init__(self, file):

		re_keyval = re.compile("^\s*([^=\s]+)\s*=\s*(.*)$")
		
		self.values = {}
		for line in file:
			if line[0:1] == "#":
				continue
			
			m = re_keyval.match(line)
			
			if not m:
				continue
			
			key = m.group(1)
			value = m.group(2)
			
			self.values[key] = value
			
		re_valref = re.compile("\$\(([^)]+)\)")
		
		for key, value in self.values.iteritems():
			while re_valref.match(value):
				m = re_valref.match(value)
				other_key = m.group(1)
				try:
					other_value = self.values[other_key]
				except KeyError:
					print "Error: no value for %s in config file" % other_key
					sys.exit(1)
					
				re_valrepl = re.compile("\$\(" + other_key + "\)")
				value = re_valrepl.sub(other_value, value)
				
			self.values[key] = value

		for key, value in self.values.iteritems():
			print key
			print value

				
	def has_value(self, key):
		return key in self.values
		
	def get_value(self, key):
		try:
			return self.values[key]
		except KeyError:
			print "Error: no value for %s in config file" % other_key
			sys.exit(1)
			
	def get_values(self, key):
		value = get_value(self, key)
		return set(value.split(","))
	
