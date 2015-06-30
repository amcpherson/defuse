package cmdrunner;

use strict;
use warnings FATAL => 'all';
use List::Util qw[min max];
use File::stat;
use Sys::Hostname;

sub new
{
	my $class = shift;
	
	my $self = {};
	
	$self->{name} = "cmdrunner";
	$self->{prefix} = "cmdrunner";
	$self->{logfilename} = $self->{prefix}.".log";
	$self->{submitter} = \&submitter_direct;
	$self->{qsub_params} = "";
	$self->{maxparallel} = 200;
	$self->{filetimeout} = 100;
	$self->{jobmem} = "2000000000";

	bless($self);

	return $self;
}

sub name
{
	my $self = shift;
	$self->{name} = shift;
}

sub prefix
{
	my $self = shift;
	$self->{prefix} = shift;
	$self->{logfilename} = $self->{prefix}.".log";
}

sub submitter
{
	my $self = shift;
	my $submitter = shift;
	if ($submitter eq "direct")
	{
		$self->{submitter} = \&submitter_direct
	}
	elsif ($submitter eq "sge")
	{
		$self->{submitter} = \&submitter_sge_cluster
	}
	elsif ($submitter eq "pbs")
	{
		$self->{submitter} = \&submitter_pbs_cluster
	}
	elsif($submitter eq "lsf")
	{
		$self->{submitter} = \&submitter_lsf_cluster
	}
	else
	{
		die "Error: Unrecognized submitter\nSupported submitters: direct, sge, pbs\n";
	}
}

sub qsub_params
{
	my $self = shift;
	$self->{qsub_params} = shift;
}

sub maxparallel
{
	my $self = shift;
	$self->{maxparallel} = shift;
}

sub filetimeout
{
	my $self = shift;
	$self->{filetimeout} = shift;
}

sub jobmem
{
	my $self = shift;
	$self->{jobmem} = shift;
}

sub interp_qsub_params
{
	my $self = shift;

	my $job_mem = $self->{jobmem};
	my $job_mem_kb = ($job_mem / 1000); 		# Kb
	my $job_mem_mb = ($job_mem / 1000000); 		# Mb
	my $job_mem_gb = ($job_mem / 1000000000); 	# Gb

	my $qsub_params = eval($self->{qsub_params});

	return $qsub_params;
}

sub DESTROY
{
	my $self = shift;
	
	# Check if this is a job instance or main instance
	if (defined $self->{isjob})
	{
		# Remove job temp files
		if (defined $self->{jobtempfiles})
		{
			foreach my $tempfile (@{$self->{jobtempfiles}})
			{
				unlink $tempfile;
			}
		}
	}
	else
	{
		# Remove all temp files
		if (defined $self->{alltempfiles})
		{
			foreach my $job_pid (%{$self->{alltempfiles}})
			{
				foreach my $tempfile (@{$self->{alltempfiles}{$job_pid}})
				{
					unlink $tempfile;
				}
			}
		}
	}
}

sub job_int_handler
{
	print "Jobs interrupted by user\n";
	exit(1);
}

sub create_job_bash_script
{
	my $job_cmd = shift;
	my $job_script = shift;

	# Write script to execute job
	open SCR, ">".$job_script or die "Error: Unable to open job script $job_script\n";
	print SCR "if [ -f ~/.bashrc ]\n";
	print SCR "then\n";
	print SCR "source ~/.bashrc\n";
	print SCR "fi\n";
	print SCR "echo -e Running on \$HOSTNAME\n";
	print SCR "time ".$job_cmd."\n";
	print SCR "echo -e Return codes: \${PIPESTATUS[*]}\n";
	close SCR;
}

sub submitter_direct
{
	my $job_cmd = shift;
	my $job_name = shift;
	my $job_script = shift;
	my $job_out = shift;
	my $cmdrunner = shift;
	
	create_job_bash_script($job_cmd, $job_script);

	# Do the command and write the stdout and stderr to the log file
	my $sysretcode = system "bash $job_script > $job_out 2>&1";

	return ($sysretcode,"");
}

sub check_qsub
{
	# Check qsub exists on this machine
	my $qsub_exists = system "which qsub > /dev/null 2>&1";

	return 0 if $qsub_exists != 0;
	
	return 1;
}

sub submitter_sge_cluster
{
	my $job_cmd = shift;
	my $job_name = shift;
	my $job_script = shift;
	my $job_out = shift;
	my $cmdrunner = shift;
	
	return (1,"No qsub on this machine") if not check_qsub();

	create_job_bash_script($job_cmd, $job_script);

	my $qsub_params = $cmdrunner->interp_qsub_params();
	
	# Do the command with qsub and write the stdout and stderr to the log file
	#  '-sync y' means block until finished
	#  '-notify' means warn the script with a signal before killing.  on oursystem the interrupt signal is sent, but that had to be configured properly in the sge cluster config.  basically this allows alignment jobs to clean up after themselves by removing temporary files, which is important if you are writing to a local disk and dont want to have to go and clean up each nodes disk after a cancel.
	#  '-b y' means run as a binary rather than in a bash shell, and from memory is required for the above to work
	#  '-j y' means join std out and std in
	#  '-o $job_out' sends std out and std in to the job_out file
	#  '-N' just specifies the job name for qstat
	#  '$qsub_commands' are things like queue requests or resource requests.
	
	my $qsub_output = `qsub -sync y -notify -b y -j y -o $job_out -N $job_name $qsub_params -S /bin/bash 'bash $job_script' 2>&1`;
	my $sysretcode = $? >> 8;
	
	return ($sysretcode,$qsub_output);
}


sub check_bsub
{
	# Check bsub exists on this machine
	my $bsub_exists = system "which bsub > /dev/null 2>&1";
	print STDERR "Error: no bsub on this machine\n" and return 0 if $bsub_exists != 0;
	
	return 1;
}

sub submitter_lsf_cluster
{
	my $job_cmd = shift;
	my $job_name = shift;
	my $job_script = shift;
	my $job_out = shift;
	my $cmdrunner = shift;
	
	return 1 if not check_bsub();

	create_job_bash_script($job_cmd, $job_script);

	my $qsub_params = $cmdrunner->interp_qsub_params();
	
	# Do the command with bsub and write the stdout and stderr to the log file
	#  '-K' means block until finished
	#  '-N' prevents extra output after the exit codes (which was interferring with the sysretcode parsing downstream)
	# '-oo $job_out' sends std out and std in to the job_out file
	# '-J' specifies the job name
	# '$qsub_params' are things like queue requests or resource requests.

	my $cmd = "bsub -N -K -oo $job_out -J $job_name $qsub_params 'bash $job_script' > /dev/null 2>&1";
	#print STDERR "\n\nAttempting job submission as follows:\n\n$cmd\n\n";
	my $sysretcode = system $cmd;
	#print STDERR "\n\nDEBUG sysretcode: $sysretcode\n\n";

	return $sysretcode;
}


sub get_pbs_jobstate
{
	my $job_id = shift;
	
	my @job_stat = `qstat $job_id`;
	my @job_info = split /\s+/, $job_stat[$#job_stat];
	my $job_state = $job_info[4];
	
	return $job_state;
}

sub submitter_pbs_cluster
{
	my $job_cmd = shift;
	my $job_name = shift;
	my $job_script = shift;
	my $job_out = shift;
	my $cmdrunner = shift;

	return (1,"No qsub on this machine") if not check_qsub();

	create_job_bash_script($job_cmd, $job_script);

	my $qsub_params = $cmdrunner->interp_qsub_params();

	# Do the command with qsub and write the stdout and stderr to the log file
	my $job_id = `qsub -j oe -o $job_out -N $job_name $qsub_params 'bash $job_script' > /dev/null 2>&1`;
	chomp($job_id);
	
	# Wait for the job state to be complete
	while (get_pbs_jobstate($job_id) ne "C")
	{
		sleep 10;
	}

	return (0,"");
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

sub append_joblog_filename
{
	my $self = shift;
	my $joblog = shift;
	my $filename = shift;
	my $skiplast = shift;

	return if not -e $filename;
	
	open OUT, ">>".$filename or die "Error: Unable to append job log to $filename\n";
	$self->append_joblog_file($joblog, \*OUT, $skiplast);
	close OUT;
}

sub append_joblog_file
{
	my $self = shift;
	my $joblog = shift;
	my $file = shift;
	my $skiplast = shift;

	open JLG, $joblog or die "Error: Unable to find job log $joblog\n";
	my $next_line = <JLG>;
	while (<JLG>)
	{
		print $file "\t".$next_line;
		$next_line = $_;
	}
	close JLG;
	
	print $next_line if defined $next_line and not $skiplast;
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
	return 1 if not cmdrunner::uptodate($inargs,$outargs,$remarks);

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

	# Assign interrupt handler to ensure correct shutdown on interrupt
	my $prev_int = $SIG{INT};
	$prev_int = 'DEFAULT' if not defined $prev_int;
	$SIG{INT} = \&job_int_handler;
	
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
	my %job_info;
	my $num_failed = 0;
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
		
			# Add array outputs to command
			my $allout = join(" ", @{$out_args});
			$job_cmd =~ s/#>A/$allout/g;
		
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
				# Mark this instance of cmdrunner as a job
				$self->{isjob} = 1;
				
				# Job output files to be removed if this job fails
				$self->{jobtempfiles} = [@outgenerated];

				# Do the command using the run subroutine provided
				my ($sysretcode,$message) = &{$submitter}($job_cmd,$job_name,$job_script,$job_out,$self);

				# Take at least 1 second
				sleep 1;

				# Return nonzero if the job failed
				if ($sysretcode != 0)
				{
					print "Command was not properly submitted/executed\n".$message;
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
				if (scalar @{$in_args} > 0 and scalar @outgenerated > 0 and not cmdrunner::uptodate($in_args,\@outgenerated,\@remarks))
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
	
				# Dont remove job output files on exit
				$self->{jobtempfiles} = [];

				exit;
			}
			else
			{
				# Update job numbers and store job info
				$jobs_running++;
				$job_info{$job_pid}{cmd} = $job_cmd;
				$job_info{$job_pid}{out} = $job_out;
				$job_info{$job_pid}{pipe} = $job_pipe;
				
				# Add job output files to the list of all temp files
				$self->{alltempfiles}{$job_pid} = [@outgenerated];
			}

			# Wait 1 second before submitting next job
			sleep 1;
		}

		# Dont wait for nothing
		next if $jobs_running == 0;

		# Wait for next job	
		while ($jobs_running >= $maxparallel)
		{
			# Wait for next job	
			my $job_pid = wait();

			# Store status
			$job_info{$job_pid}{retcode} = $?;
			
			# Finalize the job and add to failed count
			$num_failed += $self->finalizejob($job_pid, \%job_info);
			
			# One less job running
			$jobs_running--;
		}
	}

	# Wait for remaining jobs
	while ($jobs_running > 0)
	{
		# Wait for next job	
		my $job_pid = wait();
		
		# Store status
		$job_info{$job_pid}{retcode} = $?;
		
		# Finalize the job and add to failed count
		$num_failed += $self->finalizejob($job_pid, \%job_info);			
		
		# One less job running
		$jobs_running--;
	}
	
	# Return if everything was up to date
	if (not $running_something)
	{
		$SIG{INT} = $prev_int;
		return;
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

	# Return interrupt to previous state
	$SIG{INT} = $prev_int;
}

sub finalizejob
{
	my $self = shift;
	my $job_pid = shift;
	my $job_info = shift;

	my $name = $self->{name};
	
	my $job_cmd = $job_info->{$job_pid}{cmd};
	my $job_out = $job_info->{$job_pid}{out};
	my $job_retcode = $job_info->{$job_pid}{retcode};
	my $job_pipe = $job_info->{$job_pid}{pipe};
	my $job_failed = $job_retcode;
	my $job_status = $job_retcode >> 8;
	
	# Retrieve job message
	my @job_message = <$job_pipe>;
	chomp(@job_message);
	close $job_pipe;

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
		$self->append_joblog_filename($job_out, $self->{logfilename}, 1);
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
				$self->append_joblog_file($job_out, \*STDERR, 0);
			}
			else
			{
				print STDERR "Job output:\n";
				$self->append_joblog_file($job_out, \*STDERR, 1);
				system "tail -1 $job_out 1>&2";
			}
		}
		
		# Make sure job outputs were deleted on failure
		foreach my $tempfile (@{$self->{alltempfiles}{$job_pid}})
		{
			unlink $tempfile;
		}
	}
	
	# Clear job outputs from list of all temp files
	delete $self->{alltempfiles}{$job_pid};
	
	if ($job_failed)
	{
		return 1;
	}
	
	return 0;
}

sub uptodate
{
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
		next if not -e $arg;

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

	return 1 if scalar @intimes == 0;

	# Return 1 if up to date
	return 1 if max(@intimes) <= min(@outtimes) and min(@outtimes) != 0;
	
	# Return 0 if not up to date or missing
	return 0;
}

sub replaceifdifferent
{
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

