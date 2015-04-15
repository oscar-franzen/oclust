#!/usr/bin/perl

use strict 'vars';

die("Usage: <fastq file>") if (@ARGV != 1);

open(fh, @ARGV[0]) or die("Cannot open @ARGV[0]\n");

while (my $line1 = <fh>) {
	my $line2 = <fh>;
	my $line3 = <fh>;
	my $line4 = <fh>;

	chomp($line1);
	chomp($line2);
	chomp($line3);
	chomp($line4);

	$line1 =~ /^\@(\S+)/;
	my $id = $1;

	print(">$id\n$line2\n");
}

close(fh);