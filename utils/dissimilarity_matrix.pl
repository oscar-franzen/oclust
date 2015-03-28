#!/usr/bin/perl

use strict 'vars';

die("Usage: <*.reads_renamed.txt>") if (@ARGV == 0);

my %targets;
my %matrix;
my @all_ids;

open(fh, @ARGV[0]) or die("Cannot open @ARGV[0]");

while (my $line = <fh>) {
	chomp($line);

	$line =~ /^>(\S+) (\d+)/;
	my ($read_id, $new_id) = ($1,$2);

	$read_id =~ s/,/_/;

	if (@ARGV[0] =~ /pacbio/) {
		$read_id =~ s/\//_/g;
	}

	$targets{$read_id}=$new_id;

	push(@all_ids, $new_id);
}

close(fh);

my $first_row = " ";

foreach my $id (@all_ids) {
	$first_row .= "$id ";
}

chop($first_row);
print("$first_row\n");

my $count = 0;

for (my $i=1; $i<@ARGV; $i++) {
	my $file = @ARGV[$i];

	open(fh, $file);

	$count ++ ;

	#print(STDERR "$count\n");

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+) (\S+) (\S+) /;
		my ($read1, $read2, $identity) = ($1,$2,$3);

		if ($targets{$read1} != "" && $targets{$read2} ne "") {
			my $new_id_read1 = $targets{$read1};
			my $new_id_read2 = $targets{$read2};

			$matrix{$new_id_read1}{$new_id_read2} = $identity;
		}
	}

	close(fh);
}

# Print the matrix
foreach my $id1 (@all_ids) {
	my $row_line = "$id1 ";
	
	foreach my $id2 (@all_ids) {
		$row_line .= $matrix{$id1}{$id2} . " ";

		if ($matrix{$id1}{$id2} eq "") {
			die("Cannot find $id1 and $id2");
		}
	}

	chop($row_line);
	print("$row_line\n");
}
