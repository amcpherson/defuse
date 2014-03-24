#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use File::Basename;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Create circos data\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -f, --fusions   Fusion results filename\n";
push @usage, "  -c, --cnv       CNV segments filename\n";
push @usage, "  -p, --prefix    Prefix for data filenames\n";

my $help;
my $results_filename;
my $cnv_filename;
my $prefix;

GetOptions
(
	'help'        => \$help,
	'fusions=s'   => \$results_filename,
	'cnv=s'       => \$cnv_filename,
	'prefix=s'    => \$prefix,
);

not defined $help or usage() and exit;

defined $results_filename or die @usage;
defined $prefix or die @usage;

die "Error: prefix cannot be 'circos'\n" if $prefix eq "circos";

my $script_dir = dirname($0);
my $template_conf = $script_dir."/circos.conf";
my $convert_results_script = $script_dir."/results_to_circos_links.pl";
my $convert_cnv_script = $script_dir."/cnv_to_circos_points.pl";
my $filter_equal_script = $script_dir."/filter_equal.pl";

my $svg_filename = $prefix.".svg";
my $link_1_filename = $prefix.".1.link";
my $link_2_filename = $prefix.".2.link";
my $conf_filename = $prefix.".circos.conf";

my $convert_1_result = system "cat $results_filename | $filter_equal_script orf Y | $convert_results_script > $link_1_filename";
$convert_1_result == 0 or die "Error: Unable to run command\ncat $results_filename | $filter_equal_script orf Y | $convert_results_script > $link_1_filename\n";

my $convert_2_result = system "cat $results_filename | $filter_equal_script orf N | $convert_results_script > $link_2_filename";
$convert_2_result == 0 or die "Error: Unable to run command\ncat $results_filename | $filter_equal_script orf Y | $convert_results_script > $link_2_filename\n";

my %cnv_levels_filenames;
$cnv_levels_filenames{"1,2"} = $prefix.".somatic.loss.txt";
$cnv_levels_filenames{"3"} = $prefix.".neut.txt";
$cnv_levels_filenames{"4,5,6"} = $prefix.".somatic.gain.txt";
$cnv_levels_filenames{"7,8"} = $prefix.".germline.loss.txt";
$cnv_levels_filenames{"9,10,11"} = $prefix.".germline.gain.txt";

foreach my $cnv_level_list (keys %cnv_levels_filenames)
{
	my $cnv_levels_filename = $cnv_levels_filenames{$cnv_level_list};
	if (defined $cnv_filename)
	{
		my $command = "$convert_cnv_script $cnv_level_list < $cnv_filename > $cnv_levels_filename";
		my $convert_result = system $command;
		$convert_result == 0 or die "Error: Unable to run $command\n";
	}
	else
	{
		system "touch $cnv_levels_filename";
	}
}

open TP, $template_conf or die "Error: Unable to open $template_conf: $!\n";
open CFG, ">".$conf_filename or die "Error: Unable to write to $conf_filename: $!\n";
while (<TP>)
{
	my $svg_dir = dirname($svg_filename);
	my $svg_base = basename($svg_filename);

	s/defuse\.1\.link/$link_1_filename/;
	s/defuse\.2\.link/$link_2_filename/;
	s/defuse\.svg/$svg_base/;
	s/outputdirectory/$svg_dir/;

	s/somatic\.loss\.txt/$prefix.somatic.loss.txt/;
	s/neut\.txt/$prefix.neut.txt/;
	s/somatic\.gain\.txt/$prefix.somatic.gain.txt/;
	s/germline\.loss\.txt/$prefix.germline.loss.txt/;
	s/germline\.gain\.txt/$prefix.germline.gain.txt/;

	print CFG $_;
}
close TP;
close CFG;

