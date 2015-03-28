#!/usr/bin/perl

##################################################################################
#    Written by Oscar Franzen, <p.oscar.franzen@gmail.com>, Mt. Sinai, NY, USA.  #
##################################################################################

use strict 'vars';
use Cwd 'abs_path';

my $cwd = abs_path($0);

$cwd =~ /^(\S+\/)utils/;
$cwd = $1;

sub workaround {
	my $p = abs_path($0);
	$p =~ /^(\S+\/)utils/;
	$p = $1;

	my $path = $p . "modules";
	return $path;
}

use lib workaround();
use Bio::AlignIO;

open(fh, @ARGV[0]) or die("Cannot open @ARGV[0]");

my @ids;
my %id_targets;

while (my $line = <fh>) {
	chomp($line);

	if ($line =~ /^\# [12]: /) {
		$line =~ /\s(\S+)$/;
		
		push(@ids, $1);
		$id_targets{$1} = 1;
	}
}

close(fh);

# read in an alignment from the EMBOSS program water
my $in = new Bio::AlignIO(-format => 'emboss', -file => @ARGV[0]);

my $alignment_count = 0;

while( my $aln = $in->next_aln ) {
	my @foo;

	# % identity calculated as in http://drive5.com/usearch/manual/id_threshold.html
	# i.e., identity = (number of identities)/(number of columns)
	# i.s., the definition used in BLAST and USEARCH

	# An "identity" is an alignment column containing two identical letters. With global alignments, terminal gaps are discarded and do not count towards the number of columns; internal gaps are included.

	foreach my $seq ($aln->each_seq()) {
		push(@foo, $seq->seq);
	}

	# Remove terminal gaps
	my $str1 = @foo[0];
	my $str2 = @foo[1];

	$str1 = uc($str1);
	$str2 = uc($str2);

	if ($str1 =~ /^\-/) {
		$str1 =~ /^(-+)(\S+)/;
		my $gaps_str = $1;
		$str1 = $2;
		$str2 = substr($str2, length($gaps_str), length($str2));
	}

	if ($str2 =~ /^\-/) {
		$str2 =~ /^(-+)(\S+)/;
		my $gaps_str = $1;
		$str2 = $2;
		$str1 = substr($str1, length($gaps_str), length($str1));
	}

	if ($str1 =~ /\-+$/) {
		$str1 =~ /(\-+)$/;
		my $gaps_str = $1;

		$str1 = substr($str1, 0, length($str1) - length($gaps_str));
		$str2 = substr($str2, 0, length($str2) - length($gaps_str));
	}

	if ($str2 =~ /\-+$/) {
		$str2 =~ /(\-+)$/;
		my $gaps_str = $1;

		$str2 = substr($str2, 0, length($str2) - length($gaps_str));
		$str1 = substr($str1, 0, length($str1) - length($gaps_str));
	}

	# Delete single base indels
	my @positions_to_delete;

	while ($str2 =~ m/[ATGC](-)[ATGC]/g) {
		my $p = pos($str2);
		push(@positions_to_delete, $p - 2);
	}

	for (my $i=@positions_to_delete-1; $i>=0; $i--) {
		my $item = @positions_to_delete[$i];
		
		substr($str1, $item, 1) = "";
		substr($str2, $item, 1) = "";
	}

	my @positions_to_delete;

	while ($str1 =~ m/[ATGC](-)[ATGC]/g) {
		my $p = pos($str1);
		push(@positions_to_delete, $p - 2);
	}

	for (my $i=@positions_to_delete-1; $i>=0; $i--) {
		my $item = @positions_to_delete[$i];
		
		substr($str1, $item, 1) = "";
		substr($str2, $item, 1) = "";
	}

	# Delete indels larger than 4
	my @positions_to_delete;

	while ($str1 =~ m/[ATGC](-{5,})[ATGC]/g) {
		my $del_start = pos($str1) - length($1);
		push(@positions_to_delete, { del_start => $del_start, del_len => length($1) });
	}

	 for (my $i=@positions_to_delete-1; $i>=0; $i--) {
		my $item = @positions_to_delete[$i];
		
		substr($str1, $item->{del_start} - 1, $item->{del_len}) = "";
		substr($str2, $item->{del_start} - 1, $item->{del_len}) = "";
	}

	my @positions_to_delete;

	while ($str2 =~ m/[ATGC](-{5,})[ATGC]/g) {
		my $del_start = pos($str2) - length($1);
		push(@positions_to_delete, { del_start => $del_start, del_len => length($1) });
	}

	 for (my $i=@positions_to_delete-1; $i>=0; $i--) {
		my $item = @positions_to_delete[$i];
		
		substr($str1, $item->{del_start} - 1, $item->{del_len}) = "";
		substr($str2, $item->{del_start} - 1, $item->{del_len}) = "";
	}

	# Calc identity
	my $id1 = @ids[$alignment_count];
	my $id2 = @ids[$alignment_count + 1];

	my @t1 = split(//, $str1);
	my @t2 = split(//, $str2);

	my $differences = 0;
	my $total = 0;

	for (my $i=0; $i<@t1; $i++) {
		if (@t1[$i] ne "N" && @t2[$i] ne "N") {

			# There are some IUPAC-encoded nucleotides in the miseq data - IGNORE!
			if (@t1[$i] =~ /[ATGC\-]/ && @t2[$i] =~ /[ATGC\-]/) {
				if (@t1[$i] ne @t2[$i]) {
					$differences ++ ;
				}

				$total ++ ;
			}
		}
	}

	# # Count number of insertions >= 2 nt
	# my @count_seq_1 = ($str1 =~ m/-{2,}/g);
	# my @count_seq_2 = ($str2 =~ m/-{2,}/g);

	# Treat each gap as one difference no matter how long it is
	# $differences += int(@count_seq_1) + int(@count_seq_2);

	my $fraction_differences = $differences / $total;
	$fraction_differences = sprintf("%.5f", $fraction_differences);

	print("$id1 $id2 $fraction_differences $differences $total\n");

	$alignment_count += 2 ;
}
