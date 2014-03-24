#!/usr/bin/perl

use strict;
use warnings;
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

my $last_job_start_time;
my $jobs_to_start = scalar @jobs;
print "Starting $jobs_to_start jobs\n";
my $jobs_running = 0;
my %pid_to_job;
my @failed_jobs;
while ($jobs_to_start > 0)
{
	while ($jobs_running < $num_threads and $jobs_to_start > 0)
	{
		my $job = pop @jobs;

		if (defined $last_job_start_time)
		{
			my $sleep_time = $delay - (time - $last_job_start_time);
			if ($sleep_time > 0)
			{
				print "Pausing for $sleep_time seconds\n";
				sleep $sleep_time;
			}
		}

		$last_job_start_time = time;

		my $pid = fork();
		
		die "Unable to get a new process id\n" if not defined $pid;

		if ($pid == 0)
		{
			my $return_value = system $job;

			if ($return_value != 0)
			{
				exit(-1);
			}

			exit;
		}
		else
		{
			print "Job $pid started\n";
			$jobs_running++;
			$jobs_to_start--;
			$pid_to_job{$pid} = $job;
		}
	}
	
	my $finished_pid = wait();

	my $status = "success";
	if ($? > 0)
	{
		$status = "failure";
		push @failed_jobs, $pid_to_job{$finished_pid};
	}
	
	$jobs_running--;
	
	print "Job $finished_pid finished with status $status\n";
}

while ($jobs_running)
{
	my $finished_pid = wait();

	my $status = "success";
	if ($? > 0)
	{
		$status = "failure";
		push @failed_jobs, $pid_to_job{$finished_pid};
	}
	
	$jobs_running--;

	print "Job $finished_pid finished with status $status\n";
}

my $num_failed = scalar @failed_jobs;
if ($num_failed > 0)
{
	print "The following $num_failed jobs failed:\n";
	foreach my $failed_job (@failed_jobs)
	{
		print $failed_job."\n";
	}
	exit(1);
}
else
{
	print "Finished succesfully\n";
	exit;
}
