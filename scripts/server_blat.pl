#!/usr/bin/perl

use strict;
use warnings;
use POSIX ":sys_wait_h";

my $parent_pid = $$;

my $target_filename = shift;
my $query_filename = shift;
my $output_filename = shift;

my $blat_args = "";
while (my $blat_arg = shift)
{
	$blat_args .= $blat_arg." ";
}

sub server_up
{
	my $port = shift;
	my $server_return_code = system "gfServer status localhost $port 2> /dev/null > /dev/null";
	return $server_return_code == 0;
}

my $random_port = int(rand(2000)) + 1025;
while (server_up($random_port))
{
	$random_port = int(rand(2000)) + 1025;
}

my $server_command = "gfServer start localhost $random_port -stepSize=12 -tileSize=12 -minMatch=4 $target_filename";

my $pid = fork();
if ($pid == 0)
{
	my $server_pid = open(SERVER, "$server_command |") or die;# "Couldn't launch '$server_command': $! / $?";

	local $SIG{USR1} = sub { kill 'INT', $server_pid; exit; };

	<SERVER>;

	exit 1;
}

while (not server_up($random_port))
{
	my $waited_pid = waitpid($pid, WNOHANG);
	if ($waited_pid > 0)
	{
		die "Error: Unable to run server using command\n$server_command\n";
	}

	sleep 1;
}

print "Starting alignment\n";

my $client_return_code = system "gfClient -nohead $blat_args localhost $random_port / $query_filename $output_filename";

print "Finished alignment\n";

kill 'USR1', $pid;
waitpid $pid, 0;

exit 1 if $client_return_code != 0;

