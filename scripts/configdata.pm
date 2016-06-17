package configdata;

use strict;
use warnings FATAL => 'all';

sub new
{
	my $self = {};
	
	bless($self);
	
	return $self;
}

sub read
{
	my $self = shift;
	my $config_filename = shift;
	my $dataset_directory = shift;
	my $source_directory = shift;

	my %config_values;

	$config_values{"dataset_directory"} = $dataset_directory;
	$config_values{"source_directory"} = $source_directory;

	open CFG, $config_filename or die "Error: Unable to open $config_filename\n";
	while (<CFG>)
	{
		chomp;
		/^\s*([^=\s]+)\s*=\s*(.*)$/;

		my $key = $1;
		my $value = $2;
		
		next if not defined $key;
		next if not defined $value;
	
		$config_values{$key} = $value;
	}
	close CFG;
	
	foreach my $key (keys %config_values)
	{
		while ($config_values{$key} =~ /\$\(([^)]+)\)/)
		{
			my $other_key = $1;
			
			if (not defined $config_values{$other_key})
			{
				die "Error: no value for $other_key in config file $config_filename\n";
			}
			
			$config_values{$key} =~ s/\$\($other_key\)/$config_values{$other_key}/;
		}
	}
	
	$self->{"config_values"} = \%config_values;
	$self->{"config_filename"} = $config_filename;
}

sub has_value
{
	my $self = shift;
	my $key = shift;

	my $config_values = $self->{"config_values"};
	my $config_filename = $self->{"config_filename"};

	defined $config_values and defined $config_filename or die "Error: config not read\n";
	
	return defined $config_values->{$key};
}

sub get_value
{
	my $self = shift;
	my $key = shift;

	my $config_values = $self->{"config_values"};
	my $config_filename = $self->{"config_filename"};

	defined $config_values and defined $config_filename or die "Error: config not read\n";
	
	if (not defined $config_values->{$key})
	{
		die "Error: no value for $key in config file $config_filename\n";
	}
	
	return $config_values->{$key};
}

sub get_list
{
	my $self = shift;
	my $key = shift;
	
	my $index = 1;
	my $key_index = $key.$index;
	my @values;
	while ($self->has_value($key_index))
	{
		push @values, $self->get_value($key_index);
		$index++;
		$key_index = $key.$index;
	}
	
	return @values;
}

sub get_hash
{
	my $self = shift;
	my $key = shift;

	my $config_values = $self->{"config_values"};
	my $config_filename = $self->{"config_filename"};
	
	defined $config_values and defined $config_filename or die "Error: config not read\n";
	
	if (not defined $config_values->{$key})
	{
		die "Error: no value for $key in config file $config_filename\n";
	}

	my %values;
	foreach my $value (split /,/, $config_values->{$key})
	{
		$value =~ s/\s//g;
		$values{$value} = 1;
	}
	
	return %values;
}

1;
