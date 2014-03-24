
	# Assign interrupt handler to ensure correct shutdown on interrupt
	$SIG{INT} = \&job_int_handler;
	
	bless($self);

import os
import time

class cmdrun:
	def __init__(self):
		self.name = "anonymous";
		self.workdir = ".";
		self.submitter = submitter_direct;
		self.maxparallel = 200;
		self.filetimeout = 100;
		self.jobmem = "2G";
		self.tempfiles = [];
	def logfilename(self):
		return self.workdir + "/" + self.name + ".log";
	def writelog(self, *messages, prepend=""):
		logfilename = self.logfilename()
		logfile = open(logfilename, 'w')
		for message in messages:
			logfile.write(prepend)
			logfile.write(message + "\n")
		logfile.close()
	def catjoblog(self, filename):
		if not os.path.exists(filename):
			return
		retcode = subprocess.call("head --lines=-1 %s | perl -pe 's/^/\\t/' >> %s" % filename, self.logfilename())
		if retcode != 0:
			raise Exception("Error: Unable to cat to log %s" % self.logfilename())
	def waitfiles(self, filenames, missing=[]):
		start = time.time()
		while time.time() - start < self.filetimeout:
			del missing[:]
			for filename in filenames:
				if os.path.exists(filename):
					missing.append(filename)
			if len(missing) == 0:
				return True
			else:
				time.sleep(1)
		return False
	def uptodate(self, inargs, outargs, remarks=[]):
		intimes = []
		for arg in inargs:
			if not os.path.exists(arg):
				raise Exception("Error: Input %s not found" % arg)
			intimes.append(os.path.getmtime(arg))
		outtimes = []
		for arg in outargs:
			modtime = 0
			if os.path.exists(arg):
				modtime = os.path.getmtime(arg)
			outtimes.append(modtime)
			if modtime == 0:
				remarks.append("%s missing" % arg)
			elif modtime < max(intimes):
				remarks.append("%s out of date" % arg)
		if max(intimes) < min(outtimes) and min(outtimes) != 0:
			return True
		return False
	def runrequired(self, inargs, outargs, remarks=[]):
		if len(inargs) == 0:
			remarks.append("No input file specified")
		if len(outargs) == 0:
			remarks.append("No output file specified")
		if len(inargs) == 0 or len(outargs) == 0:
			return True
		if not self.uptodate(inargs, outargs, remarks):
			return True
		return False
	def run(self, cmd, inargs, outargs):
		self.padd(cmd, inargs, outargs);
		self.prun()
	def padd(self, cmd, inargs, outargs):
		self.joblist.append(cmd, inargs, outargs)


	def prun(self):

		if len(self.joblist) == 0:
			return
		numjobs = len(self.joblist)

		maxparallel = self.maxparallel
		if maxparallel == 0:
			maxparallel = numjobs

		self.writelog("Starting %d %s command(s) on %s" % numjobs, self.name, socket.gethostname())

		start = time.time()

		running = []
		jobs_running = 0
		job_number = 0
%		pid_to_cmd = {}
%		pid_to_out = {}
%		pid_to_retcode = {}
%		pid_to_pipe = {}
		running_something = False

		while len(self.joblist) > 0:

			while jobs_running < maxparallel and len(self.joblist) > 0:

				jobinfo = self.joblist.pop()
				job_number++

				job_cmd = jobinfo[0]
				in_args = jobinfo[1]
				out_args = jobinfo[2]

				job_name = "%s.%s" % self.name, job_number
				job_script = "%s.%s.sh" % self.prefix, job_number
				job_out = "%s.%s.log" % self.prefix, job_number

				# Remove previously generated log
				os.remove(job_out)

				# Check if running the job is required
				runrequired_remarks = []
				if self.runrequired(in_args, out_args, runrequired_remarks):
					self.writelog("Starting $name command:")
					self.writelog(job_cmd, prefix="\t")
					self.writelog("Reasons:")
					self.writelog(runrequired_remarks, prefix="\t")
				else:
					self.writelog("Skipping $name command with up to date outputs:")
					self.writelog(job_cmd, prefix="\t")
					continue

				# At least one thing is running
				running_something = True

				# Create a mapping between temp and output files
				argtotmp = {}
				tmptoarg = {}
				for arg in out_args:
					tmparg = "%s.tmp" % arg
					argtotmp[arg] = tmparg
					tmptoarg[tmparg] = arg

				# Add inputs to command
				inargnum = 1;
				for arg in in_args:
					placeholder = "#<%d" % inargnum
					job_cmd.replace(placeholder, arg)
					inargnum++

				# Add array inputs to command
				allin = " ".join(in_args);
				job_cmd.replace("#<A", allin)
		
				# Add outputs to command
				# Keep track of which outputs are on the command line
				# Add temporary files to static list
				# Add non command line files to a static list
				outgenerated = []
				cmdlineoutarg = {}
				outargnum = 1
				for arg in out_args
					tmpoutarg = argtotmp[arg]
					placeholder = "#>".outargnum
					if job_cmd.replace(placeholder, tmpoutarg):
						cmdlineoutarg[tmpoutarg] = 1
						outgenerated.append(tmpoutarg)
					else:
						outgenerated.append(arg)
					outargnum++

				# Check command
				if re.search("#>[0-9]", job_cmd) or re.search("#<[0-9]", job_cmd):
					raise Exception("Error: cmdrunner syntax\n%s\n" % job_cmd)

			# Fork a new process using open
			my $job_pipe;
			my $job_pid = open $job_pipe, "-|";
			die "Error: Unable to get a new process id\n" if not defined $job_pid;
			
			# Check if process is child
			if ($job_pid == 0)
			{
				# Add output files to the list of temp files generated by this job
				push @{$self->{tempfiles}}, @outgenerated;				

				# Write script to execute job
				open SCR, ">".$job_script or die "Error: Unable to open job script $job_script\n";
				print SCR "source ~/.bashrc\n";
				print SCR "echo -e Running on \$HOSTNAME\n";
				print SCR $job_cmd."\n";
				print SCR "echo -e Return codes: \${PIPESTATUS[*]}\n";
				close SCR;

				# Do the command using the run subroutine provided
				my $sysretcode = &{$submitter}($job_name,$job_script,$job_out,$self);

				# Take at least 1 second
				sleep 1;

				# Return nonzero if the job failed
				if ($sysretcode != 0)
				{
					print "Command was not properly submitted/executed\n";
					exit(1);
				}

				# Wait for job output to appear
				if (not waitfiles($filetimeout,[$job_out]))
				{
					print "Timed out waiting for $job_out to appear\n";
					exit(1);
				}

				# Get the pipestatus result from the bottom of the log file
				my $pipestatus = `tail -1 $job_out`;
				chomp($pipestatus);

				# Check for any failure in the pipestatus and system return code
				if (not $pipestatus =~ /Return codes/ or $pipestatus =~ /[1-9]/)
				{
					print "Job command with nonzero return code\n";
					exit(1);
				}

				# Wait at least filetimout seconds for outputs to appear
				my @outmissing;
				if (not waitfiles($filetimeout,\@outgenerated,\@outmissing))
				{
					print join "\n", ("Timed out waiting for output files:", @outmissing, "");
					exit(1);
				}

				# Check outputs are up to date
				my @remarks = ();
				if (scalar @{$in_args} > 0 and scalar @outgenerated > 0 and not $self->uptodate($in_args,\@outgenerated,\@remarks))
				{
					print join "\n", ("Job did not update output files:", @remarks, "");
					exit(1);
				}

				# Move the results to the correct filenames
				foreach my $tmparg (keys %tmptoarg)
				{
					next if not defined $cmdlineoutarg{$tmparg};
					die "Error: $tmparg was not found\n" if not -e $tmparg;
					rename $tmparg, $tmptoarg{$tmparg};
				}
	
				# Clear temp file list, they have all been renamed to their correct filenames
				@{$self->{tempfiles}} = ();
				
				exit;
			}
			else
			{
				# Update job numbers and store job info
				$jobs_running++;
				$pid_to_cmd{$job_pid} = $job_cmd;
				$pid_to_out{$job_pid} = $job_out;
				$pid_to_pipe{$job_pid} = $job_pipe;
			}
		}

		# Dont wait for nothing
		next if $jobs_running == 0;

		# Wait for next job	
		while ($jobs_running >= $maxparallel)
		{
			# Wait for next job	
			my $pid = wait();
		
			# Store status
			$pid_to_retcode{$pid} = $?;
			
			# One less job running
			$jobs_running--;
		}
	}

	# Wait for remaining jobs
	while ($jobs_running > 0)
	{
		# Wait for next job	
		my $pid = wait();
	
		# Store status
		$pid_to_retcode{$pid} = $?;
		
		# One less job running
		$jobs_running--;
	}
	
	# Return if everything was up to date
	return if not $running_something;
	
	# Iterate through jobs that were running to check results
	my $num_failed = 0;
	foreach my $pid (keys %pid_to_cmd)
	{
		my $job_cmd = $pid_to_cmd{$pid};
		my $job_out = $pid_to_out{$pid};
		my $job_retcode = $pid_to_retcode{$pid};
		my $job_pipe = $pid_to_pipe{$pid};
		my $job_failed = $job_retcode;
		my $job_status = $job_retcode >> 8;
		
		# Retrieve job message
		my @job_message = <$job_pipe>;
		chomp(@job_message);
	
		# Keep track of the number failing
		$num_failed++ if $job_failed;

		# Append job output to log file
		if ($job_failed)
		{
			$self->writelog("", "Failure for $name command:");
			$self->writelog("\t", $job_cmd);
			$self->writelog("", "Reason:");
			$self->writelog("\t", @job_message);
		}
		else
		{
			$self->writelog("", "Success for $name command:");
			$self->writelog("\t", $job_cmd);
		}
		
		# Write info to log
		if (-e $job_out)
		{
			my $pipestatus = `tail -1 $job_out`;
			chomp($pipestatus);
			$self->writelog("", $pipestatus);

			$self->writelog("", "Job output:");
			$self->catjoblog($job_out);
		}

		# Write to stderr if job failed
		if ($job_failed)
		{
			print STDERR "Failure for $name command:\n";
			print STDERR "\t$job_cmd\n";

			print STDERR "Reason:\n";

			foreach my $job_message_line (@job_message)
			{
				print STDERR "\t$job_message_line\n";
			}

			if (-e $job_out)
			{
				if ($job_status == 2)
				{
					print STDERR "Job output:\n";
					system "cat $job_out | perl -pe 's/^/\\t/' 1>&2";
				}
				else
				{
					print STDERR "Job output:\n";
					system "head --lines=-1 $job_out | perl -pe 's/^/\\t/' 1>&2";
					system "tail -1 $job_out 1>&2";
				}
			}
		}
	}

	# End timer
	my $end = time();
	my $timetaken = $end - $start;

	# Check for success
	if ($num_failed != 0)
	{
		$self->writelog("", "$num_failed commands failed after $timetaken seconds");

		if ($num_failed > 1)
		{
			die "$num_failed commands failed after $timetaken seconds\n";
		}
		else
		{
			die "Commands failed after $timetaken seconds\n";
		}
	}
	
	# Report time to log
	$self->writelog("", "Completed in $timetaken seconds");
}





sub DESTROY
{
	my $self = shift;
	if (scalar @{$self->{tempfiles}} > 0)
	{
		$self->writelog("", "Job interrupted");
		unlink @{$self->{tempfiles}};
	}
}

sub job_int_handler
{
	print "Job interrupted by user\n";
	exit(1);
}

sub submitter_direct
{
	my $job_name = shift;
	my $job_script = shift;
	my $job_out = shift;
	my $cmdrunner = shift;

	# Do the command and write the stdout and stderr to the log file
	my $sysretcode = system "bash $job_script > $job_out 2>&1";

	return $sysretcode;
}

sub submitter_sge_cluster
{
	my $job_name = shift;
	my $job_script = shift;
	my $job_out = shift;
	my $cmdrunner = shift;

	my $job_mem = $cmdrunner->{jobmem};
	my $qsub_commands = "-q fusions.q -l mem_free=$job_mem";
	
	# Check qsub exists on this machine
	my $qsub_exists = system "which qsub > /dev/null 2>&1";
	print STDERR "Error: no qsub on this machine\n" and return 1 if $qsub_exists != 0;
	
	# Do the command with qsub and write the stdout and stderr to the log file
	my $sysretcode = system "qsub -sync y -notify -b y -j y -o $job_out -N $job_name $qsub_commands 'bash $job_script' > /dev/null 2>&1";

	return $sysretcode;
}

sub writelog
{
	my $self = shift;
	my $prepend = shift;
	my @messages = @_;
	
	my $logfilename = $self->{logfilename};

	open LOG, ">>".$logfilename or die "Error: Unable to write to log file $logfilename: $!\n";
	foreach my $message (@messages)
	{
		print LOG "$prepend";
		print LOG $message."\n";
	}
	close LOG;
}

sub catjoblog
{
	my $self = shift;
	my $filename = shift;
	my $logfilename = $self->{logfilename};

	return if not -e $filename;

	my $cat_result = system "head --lines=-1 $filename | perl -pe 's/^/\\t/' >> $logfilename";

	$cat_result == 0 or die "Error: Unable to append to log\n";
}

sub waitfiles
{
	my $filetimeout = shift;
	my $filenames = shift;
	my $missing = shift;

	# Optional out argment
	$missing = [] if not defined $missing;
	
	# Reset missing list
	@{$missing} = ();

	# Start timer
	my $start = time();

	# Wait at least $filetimout seconds for outputs to appear
	while (time() - $start < $filetimeout)
	{
		# Reset missing list
		@{$missing} = ();

		foreach my $filename (@{$filenames})
		{
			if (not -e $filename)
			{
				push @{$missing}, $filename;
			}
		}

		if (scalar @{$missing} == 0)
		{
			return 1;
		}
		else
		{
			sleep 1;
		}
	}

	return 0;
}

sub uptodate
{
	my $self = shift;
	my $inargs = shift;
	my $outargs = shift;
	my $remarks = shift;

	# Optional argument
	$remarks = [] if not defined $remarks;

	# Reset remarks list
	@{$remarks} = ();

	# Check the arguments
	die "Error: Non-empty array ref expected for inargs argument to uptodate\n" if ref($inargs) ne 'ARRAY' or scalar @{$inargs} == 0;
	die "Error: Non-empty array ref expected for outargs argument to uptodate\n" if ref($outargs) ne 'ARRAY' or scalar @{$outargs} == 0;

	# Find modification times of inputs
	# Check all inputs exist
	my @intimes;
	foreach my $arg (@{$inargs})
	{
		die "Error: Input $arg not found\n" if not -e $arg;

		my $modtime = 0;
		$modtime = stat($arg)->mtime;
		push @intimes, $modtime;
	}
	
	# Find modification times of outputs
	my @outtimes;
	foreach my $arg (@{$outargs})
	{
		my $modtime = 0;
		$modtime = stat($arg)->mtime if -e $arg;
		push @outtimes, $modtime;

		if ($modtime == 0)
		{
			push @{$remarks}, $arg." missing";
		}
		elsif (scalar @intimes > 0 and $modtime < max(@intimes))
		{
			push @{$remarks}, $arg." out of date";
		}
	}

	# Return 1 if up to date
	return 1 if max(@intimes) < min(@outtimes) and min(@outtimes) != 0;
	
	# Return 0 if not up to date or missing
	return 0;
}

sub runrequired
{
	my $self = shift;
	my $inargs = shift;
	my $outargs = shift;
	my $remarks = shift;
		
	# Check the arguments
	die "Error: Array ref expected for inargs argument to uptodate\n" if ref($inargs) ne 'ARRAY';
	die "Error: Array ref expected for outargs argument to uptodate\n" if ref($outargs) ne 'ARRAY';

	# Add remark for no inputs
	if (scalar @{$inargs} == 0)
	{
		push @{$remarks}, "No input file specified";
	}
	
	# Add remark for no outputs
	if (scalar @{$outargs} == 0)
	{
		push @{$remarks}, "No output file specified";
	}

	# Always run commands if no inputs or outputs were specified
	return 1 if scalar @{$inargs} == 0 or scalar @{$outargs} == 0;

	# Otherwise, run commands if something is not up to date
	return 1 if not $self->uptodate($inargs,$outargs,$remarks);

	return 0;
}

sub run
{
	my $self = shift;
	my $cmd = shift;
	my $inargs = shift;
	my $outargs = shift;

	$self->padd($cmd,$inargs,$outargs);
	$self->prun();
}

sub padd
{
	my $self = shift;
	my $cmd = shift;
	my $inargs = shift;
	my $outargs = shift;
	
	push @{$self->{joblist}}, [$cmd,$inargs,$outargs];
}

sub prun
{
	my $self = shift;
	
	# Name of these jobs
	my $name = $self->{name};

	# Overridable subroutine for running scripts
	my $submitter = $self->{submitter};

	# Retreive list of jobs and clear
	my @joblist = @{$self->{joblist}};
	@{$self->{joblist}} = ();
	my $num_jobs = scalar @joblist;
	
	# Dont bother doing anything for zero jobs
	return if $num_jobs == 0;
	
	# Prefix for script filenames
	my $prefix = $self->{prefix};

	# Log filename
	my $logfilename = $self->{logfilename};
	
	# Timout for waiting for files
	my $filetimeout = $self->{filetimeout};

	# Maximum number of parallel commands, 0 for unlimited
	my $maxparallel = $self->{maxparallel};
	$maxparallel = scalar @joblist if $maxparallel == 0;

	# Write parallel starting to the log
	my $host = hostname;
	$self->writelog("", "Starting $num_jobs $name command(s) on $host");
	
	# Start timer
	my $start = time();
		
	# Keep track of running commands
	my @running;
	
	# Iterate through the list of commands
	my $jobs_running = 0;
	my $job_number = 0;
	my %pid_to_cmd;
	my %pid_to_out;
	my %pid_to_retcode;
	my %pid_to_pipe;
	my $running_something = 0;
	while (scalar @joblist > 0)
	{
		while ($jobs_running < $maxparallel and scalar @joblist > 0)
		{
			my $jobinfo = pop @joblist;
			$job_number++;

			my $job_cmd = $jobinfo->[0];
			my $in_args = $jobinfo->[1];
			my $out_args = $jobinfo->[2];

			my $job_name = $name.".".$job_number;
			my $job_script = $prefix.".".$job_number.".sh";
			my $job_out = $prefix.".".$job_number.".log";

			# Remove previously generated log
			unlink $job_out;

			# Check the arguments
			die "Error: Array ref expected as inargs for command\n$job_cmd\n" if ref($in_args) ne 'ARRAY';
			die "Error: Array ref expected as outargs for command\n$job_cmd\n" if ref($out_args) ne 'ARRAY';
		
			# Check if running the job is required
			my @runrequired_remarks;
			if ($self->runrequired($in_args,$out_args,\@runrequired_remarks))
			{
				$self->writelog("", "Starting $name command:");
				$self->writelog("\t", $job_cmd);
				$self->writelog("", "Reasons:");
				$self->writelog("\t", @runrequired_remarks);
			}
			else
			{
				$self->writelog("", "Skipping $name command with up to date outputs:");
				$self->writelog("\t", $job_cmd);
				next;			
			}

			# At least one thing is running
			$running_something = 1;
	
			# Create a mapping between temp and output files
			my %argtotmp;
			my %tmptoarg;
			foreach my $arg (@{$out_args})
			{
				my $tmparg = $arg.".tmp";
				$argtotmp{$arg} = $tmparg;
				$tmptoarg{$tmparg} = $arg;
			}
	
			# Add inputs to command
			my $inargnum = 1;
			foreach my $arg (@{$in_args})
			{
				my $placeholder = "#<".$inargnum;
				$job_cmd =~ s/$placeholder/$arg/g;
				$inargnum++;
			}

			# Add array inputs to command
			my $allin = join(" ", @{$in_args});
			$job_cmd =~ s/#<A/$allin/g;
		
			# Add outputs to command
			# Keep track of which outputs are on the command line
			# Add temporary files to static list
			# Add non command line files to a static list
			my @outgenerated;
			my %cmdlineoutarg;
			my $outargnum = 1;
			foreach my $arg (@{$out_args})
			{
				my $tmpoutarg = $argtotmp{$arg};
				my $placeholder = "#>".$outargnum;
				if ($job_cmd =~ s/$placeholder/$tmpoutarg/g)
				{
					$cmdlineoutarg{$tmpoutarg} = 1;
					push @outgenerated, $tmpoutarg;
				}
				else
				{
					push @outgenerated, $arg;
				}
	
				$outargnum++;
			}

			# Check command
			die "Error: cmdrunner syntax\n$job_cmd\n" if $job_cmd =~ /#>[0-9]/ or $job_cmd =~ /#<[0-9]/;

			# Fork a new process using open
			my $job_pipe;
			my $job_pid = open $job_pipe, "-|";
			die "Error: Unable to get a new process id\n" if not defined $job_pid;
			
			# Check if process is child
			if ($job_pid == 0)
			{
				# Add output files to the list of temp files generated by this job
				push @{$self->{tempfiles}}, @outgenerated;				

				# Write script to execute job
				open SCR, ">".$job_script or die "Error: Unable to open job script $job_script\n";
				print SCR "source ~/.bashrc\n";
				print SCR "echo -e Running on \$HOSTNAME\n";
				print SCR $job_cmd."\n";
				print SCR "echo -e Return codes: \${PIPESTATUS[*]}\n";
				close SCR;

				# Do the command using the run subroutine provided
				my $sysretcode = &{$submitter}($job_name,$job_script,$job_out,$self);

				# Take at least 1 second
				sleep 1;

				# Return nonzero if the job failed
				if ($sysretcode != 0)
				{
					print "Command was not properly submitted/executed\n";
					exit(1);
				}

				# Wait for job output to appear
				if (not waitfiles($filetimeout,[$job_out]))
				{
					print "Timed out waiting for $job_out to appear\n";
					exit(1);
				}

				# Get the pipestatus result from the bottom of the log file
				my $pipestatus = `tail -1 $job_out`;
				chomp($pipestatus);

				# Check for any failure in the pipestatus and system return code
				if (not $pipestatus =~ /Return codes/ or $pipestatus =~ /[1-9]/)
				{
					print "Job command with nonzero return code\n";
					exit(1);
				}

				# Wait at least filetimout seconds for outputs to appear
				my @outmissing;
				if (not waitfiles($filetimeout,\@outgenerated,\@outmissing))
				{
					print join "\n", ("Timed out waiting for output files:", @outmissing, "");
					exit(1);
				}

				# Check outputs are up to date
				my @remarks = ();
				if (scalar @{$in_args} > 0 and scalar @outgenerated > 0 and not $self->uptodate($in_args,\@outgenerated,\@remarks))
				{
					print join "\n", ("Job did not update output files:", @remarks, "");
					exit(1);
				}

				# Move the results to the correct filenames
				foreach my $tmparg (keys %tmptoarg)
				{
					next if not defined $cmdlineoutarg{$tmparg};
					die "Error: $tmparg was not found\n" if not -e $tmparg;
					rename $tmparg, $tmptoarg{$tmparg};
				}
	
				# Clear temp file list, they have all been renamed to their correct filenames
				@{$self->{tempfiles}} = ();
				
				exit;
			}
			else
			{
				# Update job numbers and store job info
				$jobs_running++;
				$pid_to_cmd{$job_pid} = $job_cmd;
				$pid_to_out{$job_pid} = $job_out;
				$pid_to_pipe{$job_pid} = $job_pipe;
			}
		}

		# Dont wait for nothing
		next if $jobs_running == 0;

		# Wait for next job	
		while ($jobs_running >= $maxparallel)
		{
			# Wait for next job	
			my $pid = wait();
		
			# Store status
			$pid_to_retcode{$pid} = $?;
			
			# One less job running
			$jobs_running--;
		}
	}

	# Wait for remaining jobs
	while ($jobs_running > 0)
	{
		# Wait for next job	
		my $pid = wait();
	
		# Store status
		$pid_to_retcode{$pid} = $?;
		
		# One less job running
		$jobs_running--;
	}
	
	# Return if everything was up to date
	return if not $running_something;
	
	# Iterate through jobs that were running to check results
	my $num_failed = 0;
	foreach my $pid (keys %pid_to_cmd)
	{
		my $job_cmd = $pid_to_cmd{$pid};
		my $job_out = $pid_to_out{$pid};
		my $job_retcode = $pid_to_retcode{$pid};
		my $job_pipe = $pid_to_pipe{$pid};
		my $job_failed = $job_retcode;
		my $job_status = $job_retcode >> 8;
		
		# Retrieve job message
		my @job_message = <$job_pipe>;
		chomp(@job_message);
	
		# Keep track of the number failing
		$num_failed++ if $job_failed;

		# Append job output to log file
		if ($job_failed)
		{
			$self->writelog("", "Failure for $name command:");
			$self->writelog("\t", $job_cmd);
			$self->writelog("", "Reason:");
			$self->writelog("\t", @job_message);
		}
		else
		{
			$self->writelog("", "Success for $name command:");
			$self->writelog("\t", $job_cmd);
		}
		
		# Write info to log
		if (-e $job_out)
		{
			my $pipestatus = `tail -1 $job_out`;
			chomp($pipestatus);
			$self->writelog("", $pipestatus);

			$self->writelog("", "Job output:");
			$self->catjoblog($job_out);
		}

		# Write to stderr if job failed
		if ($job_failed)
		{
			print STDERR "Failure for $name command:\n";
			print STDERR "\t$job_cmd\n";

			print STDERR "Reason:\n";

			foreach my $job_message_line (@job_message)
			{
				print STDERR "\t$job_message_line\n";
			}

			if (-e $job_out)
			{
				if ($job_status == 2)
				{
					print STDERR "Job output:\n";
					system "cat $job_out | perl -pe 's/^/\\t/' 1>&2";
				}
				else
				{
					print STDERR "Job output:\n";
					system "head --lines=-1 $job_out | perl -pe 's/^/\\t/' 1>&2";
					system "tail -1 $job_out 1>&2";
				}
			}
		}
	}

	# End timer
	my $end = time();
	my $timetaken = $end - $start;

	# Check for success
	if ($num_failed != 0)
	{
		$self->writelog("", "$num_failed commands failed after $timetaken seconds");

		if ($num_failed > 1)
		{
			die "$num_failed commands failed after $timetaken seconds\n";
		}
		else
		{
			die "Commands failed after $timetaken seconds\n";
		}
	}
	
	# Report time to log
	$self->writelog("", "Completed in $timetaken seconds");
}

sub replaceifdifferent
{
	my $self = shift;
	my $new_filename = shift;
	my $prev_filename = shift;
	
	die "Error: File $new_filename does not exist\n" if not defined $new_filename;

	if (not -e $prev_filename)
	{
		rename $new_filename, $prev_filename;
	}
	else
	{
		my $diff_result = `diff -q $new_filename $prev_filename`;
		
		if ($diff_result =~ /differ/)
		{
			rename $new_filename, $prev_filename;
		}
	}
}

1;

