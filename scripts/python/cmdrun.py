import os
import time
import re
import threading
import socket
import subprocess
import sys


class job(object):

	def __init__(self, command, inputs, outputs, name, script, log, mem, timeout):

		self.command = command
		self.inputs = inputs
		self.outputs = outputs
		
		self.name = name
		self.script = script
		self.log = log
		self.mem = mem
		self.timeout = timeout

		self.running = False
		self.submitted = False
		self.finished = False
		self.succeeded = False
		
		self.messages = []
		self.submit_messages = []

		# Rename all output files
		self.tmp_to_output = {}
		self.tmp_outputs = {}
		for id, output in self.outputs.items():
			tmp = "%s.tmp" % output
			self.tmp_to_output[tmp] = output
			self.tmp_outputs[id] = tmp
		
		# Interpolate command
		arguments = dict(self.inputs)
		arguments.update(self.tmp_outputs)
		self.full_command = " ".join(self.command) % arguments
		
	def __del__(self):
		self.kill()
		for filename in self.tmp_outputs.values():
			if os.path.exists(filename):
				os.remove(filename)

	def runrequired(self, remarks=[]):
		if len(self.inputs.values()) == 0:
			remarks.append("No input file specified")
		if len(self.outputs.values()) == 0:
			remarks.append("No output file specified")
		if len(self.inputs.values()) == 0 or len(self.outputs.values()) == 0:
			return True
		if not uptodate(self.inputs.values(), self.outputs.values(), remarks):
			return True
		return False
					
	def init(self):
		# Remove previously generated log
		if os.path.exists(self.log):
			os.remove(self.log)
		
		# Write script to execute job
		script_file = open(self.script, 'w')
		script_file.write("echo -e Running on $HOSTNAME\n");
		script_file.write(self.full_command + "\n");
		script_file.write("echo -e \"Finished\\nReturn codes: ${PIPESTATUS[*]}\"\n");
		script_file.close()

	def submit(self, name, script, out, mem):
		raise Exception("Error: job subclass must implement submit function")
		
	def finished(self):
		raise Exception("Error: job subclass must implement finished function")

	def kill(self):
		raise Exception("Error: job subclass must implement kill function")
		
	def run(self):
		# Return if finished
		if self.finished:
			return
		
		# Submit the job using the submit function, taking at least 1 second
		if not self.running:
			self.submit()
			time.sleep(1)
			self.running = True

		# Poll the job and update if status changed
		self.poll()
		
		# Check if the job finished
		if not self.finished:
			return
		
		# Take at least 1 second
		time.sleep(1)

		# Check if the job failed to submit
		if not self.submitted:
			self.messages.append("Command was not properly submitted/executed")
			self.messages.extend(self.submit_messages)
			return
			
		# Wait for job output to appear
		if not waitfiles([self.log],self.timeout):
			self.messages.append("Timed out waiting for %s to appear" % self.log)
			return
		
		# Get the pipestatus result from the bottom of the log file
		pipestatus = str(subprocess.Popen(["tail", "-1", self.log], stdout=subprocess.PIPE).communicate()[0])

		# Check for any failure in the pipestatus and system return code
		if pipestatus.find("Return codes") < 0 or re.search('[1-9]', pipestatus):
			self.messages.append("Job command with nonzero return code")
			return

		# Wait at least filetimout seconds for outputs to appear
		missing = []
		if not waitfiles(self.tmp_outputs.values(), self.timeout, missing):
			self.messages.append("Timed out waiting for output files:")
			self.messages.extend(missing)
			return

		# Check outputs are up to date
		remarks = []
		if len(self.inputs.values()) > 0 and len(self.tmp_outputs.values()) > 0 and not uptodate(self.inputs.values(), self.tmp_outputs.values(), remarks):
			self.messages.append("Job did not update output files:")
			self.messages.extend(remarks)
			return
			
		# Move the results to the correct filenames
		for tmp, output in self.tmp_to_output.items():
			os.rename(tmp, output)
			
		# Job succeeded if it got here
		self.tmp_outputs = {}
		self.succeeded = True
	
	def explain(self, file):
		file.write("Reason:\n")
		for message in self.messages:
			file.write("\t" + message + "\n")

	def writelog(self, file):
		if os.path.exists(self.log):
			file.write("Job Output:\n")
			job_log_file = open(self.log, 'r')
			line = job_log_file.readline()
			while line != "":
				next_line = job_log_file.readline()
				if next_line == "" and self.submitted:
					file.write(line)
				else:
					file.write("\t" + line)
				line = next_line
	
	def cleanup(self):
		os.remove(self.script)
		os.remove(self.log)
		

class direct_job(job):
	
	def submit(self):
		try:
			self.log_file = open(self.log, 'w')
			self.process = subprocess.Popen("bash " + self.script, shell=True, stdout=self.log_file , stderr=subprocess.STDOUT)
			self.submitted = True
		except OSError as e:
			self.submit_messages.append("Execution failed: " + e)
			self.submitted = False

	def poll(self):
		if self.process.poll() != None:
			self.finished = True
			self.log_file.close()
			if self.process.returncode < 0:
				self.submit_messages.append("Child was terminated by signal %d" % -self.process.returncode)
				self.submitted = False
			elif self.process.returncode > 0:
				self.submit_messages.append("Child returned %d" % self.process.returncode)
				self.submitted = False

	def kill(self):
		try:
			if self.process.poll() == None:
				self.process.kill()
		except AttributeError:
			pass

class sge_job(job):
	
	def submit(self):
		# Check qsub exists on this machine
		qsub_exists = subprocess.call(["which", "qsub"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if qsub_exists != 0:
			self.submit_messages.append("Error: no qsub on this machine")

		# Do the command with qsub and write the stdout and stderr to the log file
		try:
			qsub_args = "qsub -sync y -notify -b y -j y -S /bin/bash -q fusions.q"
			qsub_args += " -o " + self.log
			qsub_args += " -N " + self.name
			qsub_args += " -l mem_free=" + self.mem + "G";
			qsub_args += " 'bash " + self.script + "'"
			
			self.process = subprocess.Popen(qsub_args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			self.submitted = True
				
		except OSError as e:
			self.submit_messages.append("qsub execution failed: " + e)
			self.submitted = False
	
	def poll(self):
		if self.process.poll() != None:
			qsub_output = str(self.process.stdout.read(),encoding='utf8').rstrip("\n")
			qsub_messages = qsub_output.split("\n")
			self.finished = True
			if self.process.returncode < 0:
				self.submit_messages.append("qsub was terminated by signal %d" % -self.process.returncode)
				self.submit_messages.extend(qsub_messages)
				self.submitted = False
			elif self.process.returncode > 0:
				self.submit_messages.append("qsub returned %d" % self.process.returncode)
				self.submit_messages.extend(qsub_messages)
				self.submitted = False

	def kill(self):
		try:
			if self.process.poll() == None:
				self.process.kill()
		except AttributeError:
			pass

def waitfiles(filenames, filetimeout, missing=[]):
	start = time.time()
	while time.time() - start < filetimeout:
		del missing[:]
		for filename in filenames:
			if not os.path.exists(filename):
				missing.append(filename)
		if len(missing) == 0:
			return True
		else:
			time.sleep(1)
	return False

def uptodate(inputs, outputs, remarks=[]):
	in_times = []
	for input in inputs:
		if not os.path.exists(input):
			raise Exception("Error: Input %s not found" % input)
		in_times.append(os.path.getmtime(input))
	out_times = []
	for output in outputs:
		mod_time = 0
		if os.path.exists(output):
			mod_time = os.path.getmtime(output)
		out_times.append(mod_time)
		if mod_time == 0:
			remarks.append("%s missing" % output)
		elif mod_time < max(in_times):
			remarks.append("%s out of date" % output)
	if max(in_times) <= min(out_times) and min(out_times) != 0:
		return True
	return False

class cmdrun(object):

	def __init__(self, name, workdir, log_filename, job_type="direct", job_mem = "2"):
		self.name = name
		self.workdir = workdir
		self.prefix = workdir + "/" + name
		self.log_filename = log_filename
		self.max_parallel = 200
		self.file_timeout = 100
		self.jobmem = job_mem
		self.jobs = []
		if not os.path.exists(self.workdir):
			os.makedirs(self.workdir)
		if job_type == "direct":
			self.job_type = direct_job
		elif job_type == "sge":
			self.job_type = sge_job
		else:
			sys.stderr.write("Error: Invalid job type\n")
			sys.exit(1)
			
	def run(self, command, inputs, outputs):
		self.padd(command, inputs, outputs);
		self.prun()

	def padd(self, command, inputs, outputs):

		# Create job name, and script and log filenames
		job_number = len(self.jobs)
		name = "%s.%s" % (self.name, job_number)
		script = "%s.%s.sh" % (self.prefix, job_number)
		log = "%s.%s.log" % (self.prefix, job_number)

		job = self.job_type(command, inputs, outputs, name, script, log, self.jobmem, self.file_timeout)
		
		self.jobs.append(job)
		
	def prun(self):
		
		# Ensure the log file is closed on exit
		with open(self.log_filename, 'a') as log_file:
			
			# Do nothing if there were no jobs
			if len(self.jobs) == 0:
				return
	
			# If maximum parallel jobs is 0, run all concurrently
			max_parallel = self.max_parallel
			if max_parallel == 0:
				max_parallel = len(self.joblist)
	
			# Starting message
			log_file.write("Starting %d %s command(s) on %s" % (len(self.jobs), self.name, socket.gethostname()))
	
			# Start a timer
			start = time.time()
	
			# Start all jobs
			jobs = []
			for job_number, job in enumerate(self.jobs):

				# Initialize the job
				job.init()
				
				# Check if running the job is required
				runrequired_remarks = []
				if job.runrequired(runrequired_remarks):
					log_file.write("Starting " + self.name + " command:\n")
					log_file.write("\t" + job.full_command + "\n")
					log_file.write("Reasons:\n")
					for runrequired_remark in runrequired_remarks:
						log_file.write("\t" + runrequired_remark + "\n")
				else:
					log_file.write("Skipping " + self.name + " command with up to date outputs:\n")
					log_file.write("\t" + job.full_command + "\n")
					continue

				# Add to list of jobs to run
				jobs.append(job)
			
			# Reset the job list
			self.jobs = []
			
			# Check if anything was started
			if len(jobs) == 0:
				return
			
			# Catch interrupt
			try:
				# Run all jobs
				running = [] 
				queued = list(jobs)
				while len(queued) > 0 or len(running) > 0:
					while len(queued) > 0 and len(running) < max_parallel:
						running.append(queued.pop())
					new_running = []
					for job in running:
						job.run()
						if not job.finished:
							new_running.append(job)
					running = new_running
					time.sleep(1)
			except KeyboardInterrupt:
				sys.stderr.write("User interrupted jobs\n")
				sys.exit(1)
			
			# Check all jobs that were running
			num_failed = 0
			for job in jobs:
				# Track number of failed jobs
				if not job.succeeded:
					num_failed += 1
					log_file.write("Failure for " + self.name + " command:\n")
					log_file.write("\t" + job.full_command + "\n")
					job.explain(log_file)
				else:
					log_file.write("Success for " + self.name + " command:\n")
					log_file.write("\t" + job.full_command + "\n")
				
				# Write info to log
				if os.path.exists(job.log):
					job.writelog(log_file)
				
				# Write to stderr if job failed
				if not job.succeeded:
					sys.stderr.write("Failure for " + self.name + " command:\n")
					sys.stderr.write("\t" + job.full_command + "\n")
					job.explain(sys.stderr)
					job.writelog(sys.stderr)
				
				# Cleanup files if job succeeded
				if job.succeeded:
					job.cleanup()
					
			# End timer
			time_taken = time.time() - start
		
			# Check for success
			if num_failed != 0:
				log_file.write("%d commands failed after %d seconds\n" % (num_failed, time_taken));
				sys.stderr.write("%d commands failed after %d seconds\n" % (num_failed, time_taken));
				sys.exit(1)
			else:
				log_file.write("%d commands succeeded after %d seconds\n" % (num_failed, time_taken));
	
