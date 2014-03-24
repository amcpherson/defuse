#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];

use lib dirname($0);
use cmdrunner;

my $runname = shift;

$runname = "testouter" if not defined $runname;

my $runner = cmdrunner->new();
$runner->runparams(qsub_commands => "-q fusions.q -l mem_free=1G");
#$runner->runsub(\&cmdrunner::run_sge_cluster);
$runner->name($runname);
$runner->prefix($runname);
$runner->maxparallel(2);
$runner->filetimeout(10);

sub run_blah
{
	my $job_name = shift;
	my $job_script = shift;
	my $job_out = shift;
	my $run_params = shift;
	
	# Do the command and write the stdout and stderr to the log file
	my $sysretcode = system "bash $job_script > $job_out 2>&1";
	
	print $run_params->{qsub_commands}."\n";

	print $job_script."\n";

	system "ls /share/data4/amcpherson/fusions/transgen/test1.tmp";
	
	return $sysretcode;
}

#$runner->runsub(\&run_blah);

$runner->padd("echo asdf1 | true", [], []);
$runner->padd("echo asdf2 | false", [], []);
$runner->padd("echo asdf3 | true", [], []);
$runner->padd("true", [], []);
$runner->padd("touch #>1", [], [abs_path("test1")]);
$runner->padd("true", [], [abs_path("test2")]);
$runner->prun();

