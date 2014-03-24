#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];


my $scriptpath = abs_path($0);

sub usage_exit
{
	print("Usage: $0 [options]\n");
	print("Run the a set of jobs on a multicore machine.\n");
	print("  -h, --help         Displays this information\n");
	print("  -s, --subhost      Submission Host\n");
	print("  -b, --batch        Batch Filename\n");
	print("  -q, --qsub         Qsub Commands\n");
	exit;
}

my $help;
my $submit_host;
my $batch_filename;
my $qsub_commands = "";
my $called_flag;

GetOptions
(
    'help'         => \$help,
    'subhost=s'    => \$submit_host,
    'batch=s'      => \$batch_filename,
    'qsub=s'       => \$qsub_commands,
    'called'       => \$called_flag,
);

not defined $help or usage_exit();

defined $batch_filename or usage_exit();

$batch_filename = abs_path($batch_filename);
my $batch_jobname_prefix = basename($batch_filename);

my $executed_host = `hostname`;
chomp($executed_host);
if ($executed_host ne $submit_host and not defined $called_flag)
{
	my $return_value = system "ssh -t $submit_host \"source ~/.bashrc; $scriptpath -s $submit_host -b $batch_filename -q '$qsub_commands' -c\"";

	if ($return_value == 0)
	{
		exit;
	}
	else
	{
		exit(1);
	}
}

open BAT, $batch_filename or die "Error: Unable to open $batch_filename: $!\n";
my @jobs = <BAT>;
close BAT;
chomp(@jobs);
@jobs = reverse(@jobs);

my $jobs_to_start = scalar @jobs;
print "Starting $jobs_to_start jobs\n";
my $jobs_running = 0;
my $job_number = 0;
my %pid_to_job;
while ($jobs_to_start > 0)
{
	my $job = pop @jobs;
	$job_number++;

	my $job_name = $batch_jobname_prefix.".".$job_number;
	my $job_out = $batch_filename.".".$job_number.".log";

	sleep 1;

	my $pid = fork();
	
	die "Unable to get a new process id\n" if not defined $pid;

	if ($pid == 0)
	{
		my $qsub_job = "qsub -sync y -notify -b y -j y -o $job_out -N $job_name $qsub_commands '$job'";
		my $return_value = system $qsub_job;

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

my @failed_jobs;
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
