#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];


sub usage_exit
{
	print("Usage: $0 [options]\n");
	print("Run the a set of jobs on a multicore machine.\n");
	print("  -h, --help         Displays this information\n");
	print("  -b, --batch        Batch Filename\n");
	print("  -t, --threads      Threads to use\n");
	print("  -d, --delay        Delay in seconds\n");
	exit;
}

my $help;
my $batch_filename;
my $num_threads;
my $delay = 0;

GetOptions
(
    'help'         => \$help,
    'batch=s'      => \$batch_filename,
    'threads=i'    => \$num_threads,
    'delay=i'      => \$delay,
);

not defined $help or usage_exit();

defined $batch_filename or usage_exit();
defined $num_threads or usage_exit();
defined $delay or usage_exit();

open BAT, $batch_filename or die "Error: Unable to open $batch_filename: $!\n";
my @jobs = <BAT>;
close BAT;
chomp(@jobs);
@jobs = reverse(@jobs);

$batch_filename = abs_path($batch_filename);
my $batch_jobname_prefix = basename($batch_filename);

# Run all jobs
my $last_job_start_time;
my $jobs_to_start = scalar @jobs;
print "[runmulticore] Starting $jobs_to_start jobs\n";
my $jobs_running = 0;
my $job_number = 0;
my %pid_to_job;
my %pid_to_out;
my %num_to_pid;
my %pid_to_status;
while ($jobs_to_start > 0)
{
	while ($jobs_running < $num_threads and $jobs_to_start > 0)
	{
		my $job = pop @jobs;
		$job_number++;

		my $job_name = $batch_jobname_prefix.".".$job_number;
		my $job_script = $batch_filename.".".$job_number.".sh";
		my $job_out = $batch_filename.".".$job_number.".log";

		# Ensure a sufficient delay
		if (defined $last_job_start_time)
		{
			my $sleep_time = $delay - (time - $last_job_start_time);
			if ($sleep_time > 0)
			{
				print "[runmulticore] Pausing for $sleep_time seconds\n";
				sleep $sleep_time;
			}
		}
		$last_job_start_time = time;

		# Fork a new process for the job
		my $pid = fork();
		die "Unable to get a new process id\n" if not defined $pid;

		# Check if process is child
		if ($pid == 0)
		{
			# Write script to execute job
			open SCR, ">".$job_script or die "Error: Unable to open job script $job_script\n";
			print SCR $job."\n";
			print SCR "echo -e [runmulticore] Return codes: \${PIPESTATUS[*]}\n";
			close SCR;

			# Do the command and write the stdout and stderr to the log file
			my $sysretcode = system "bash $job_script &> $job_out";
			print "bash $job_script 2>&1 > $job_out\n";

			# Propogate failure return code
			if ($sysretcode != 0)
			{
				exit(-1);
			}

			exit;
		}
		else
		{
			# Update job numbers and store job info
			print "[runmulticore] Job $pid started\n";
			$jobs_running++;
			$jobs_to_start--;
			$num_to_pid{$job_number} = $pid;
			$pid_to_job{$pid} = $job;
			$pid_to_out{$pid} = $job_out;
		}
	}
	
	# Wait for next job	
	my $pid = wait();

	# Store status
	$pid_to_status{$pid} = $?;
	
	$jobs_running--;
}

# Wait for remaining jobs
while ($jobs_running)
{
	# Wait for next job	
	my $pid = wait();

	# Store status
	$pid_to_status{$pid} = $?;
	
	$jobs_running--;
}

# Determine which jobs failed
my @failed_jobs;
foreach my $num (sort {$a <=> $b} keys %num_to_pid)
{
	my $pid = $num_to_pid{$num};
	my $job = $pid_to_job{$pid};
	my $job_out = $pid_to_out{$pid};
	my $status = $pid_to_status{$pid};

	print "[runmulticore] $job\n";

	my $failed = 0;
	if (-e $job_out)
	{
		# Cat output file to standard out
		my $cat_result = system "cat $job_out";

		# Get the pipestatus result from the bottom of the log file
		my $pipestatus = `tail -1 $job_out`;
		chomp($pipestatus);

		# Check for any failure in the pipestatus and system return code
		$failed = 1 if not $pipestatus =~ /Return codes/;
		$failed = 1 if $pipestatus =~ /[1-9]/;
		$failed = 1 if $cat_result != 0;
	}

	# Failed if job status was failure
	$failed = 1 if $status != 0;

	if ($failed != 0)
	{
		push @failed_jobs, $job
	}
}

# Output jobs that failed
my $num_failed = scalar @failed_jobs;
if ($num_failed > 0)
{
	print "[runmulticore] The following $num_failed jobs failed:\n";
	foreach my $failed_job (@failed_jobs)
	{
		print "[runmulticore] $failed_job\n";
	}
	exit(1);
}
else
{
	print "[runmulticore] Finished succesfully\n";
	exit;
}


