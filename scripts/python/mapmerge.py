
class delegator(object):
	def delegate(self, jobs):
		for job in jobs:
			job()


class jobcontainer(object):
	def __init__(self, callable, input):
		self.callable = callable
		self.input = input
	def __call__(self):
		self.output = self.callable(self.input)
	
def mapmerge(input, split, mapper, output, deleg):
	splitjob = jobcontainer(input.split, split)
	deleg.delegate([splitjob])
	splitinputs = splitjob.output
	mapjobs = []
	for splitinput in splitinputs:
		mapjobs.append(jobcontainer(mapper, splitinput))
	if (len(mapjobs) == 0):
		return
	deleg.delegate(mapjobs)
	splitoutputs = []
	for mapjob in mapjobs:
		splitoutputs.append(mapjob.output)
	deleg.delegate([jobcontainer(output.merge, splitoutputs)])

class testinput(object):
	def __init__(self, data):
		self.data = data
	def split(self, split):
		print split
		splitdat = []
		for dat in self.data:
			splitdat.append(testinput([dat]))
		return splitdat
			
class testmapper(object):
	def __call__(self, input):
		a = []
		for inp in input.data:
			a.append(inp * 2)
			a.append(inp * 3)
		return a

class testoutput(object):
	def merge(self, data):
		self.merged = []
		for dat in data:
			self.merged.append(dat)

import csv

def splitfastq(self, fastqfilename, split):
	splitprefix = fastqfilename + ".split.%.fastq"
	splitlistfilename = fastq2filename + ".split.list"
	runcmd("split_fastq.pl %s %d %s > %s" % fastqfilename, split["reads_per_file"], splitprefix, splitlistfilename)
	splitlistfile = open(splitlistfilename, "rb")
	return splitlistfile.readlines()

class fastqpair(object):
	def __init__(self, fastq1filename, fastq2filename):
		self.fastq1filename = fastq1filename
		self.fastq2filename = fastq2filename
	def split(self, split):
		split1list = splitfastq(fastq1filename, split)
		split2list = splitfastq(fastq2filename, split)
		if (len(split1list) != len(split2list)):
			raise
		splitfastqpairs = []
		for splitindex in range(len(split1list)):
			splitfastqpairs.append(fastqpair(split1list[splitindex],split2list[splitindex))
		return splitfastqpairs

if __name__ == "__main__":
	input = testinput(range(10))
	mapper = testmapper()
	output = testoutput()
	deleg = delegator()
	mapmerge(input, 1, mapper, output, deleg)
	
	for m in output.merged:
		print m
	