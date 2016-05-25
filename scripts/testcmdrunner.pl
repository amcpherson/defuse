#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;

use FindBin;
use lib "$FindBin::RealBin";
use cmdrunner;

my $runname = shift;

$runname = "testouter" if not defined $runname;

my $runner = cmdrunner->new();
$runner->name("delme");
$runner->prefix("delme");

$runner->padd("touch #>1 | sleep 9 | false", [], ["test1"]);
$runner->prun();

