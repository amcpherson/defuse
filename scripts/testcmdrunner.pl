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
$runner->name("delme");
$runner->prefix("delme");

$runner->padd("touch #>1 | sleep 9 | false", [], [abs_path("test1")]);
$runner->prun();

