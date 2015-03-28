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
use Bio::SeqIO;

my $is = Bio::SeqIO->new(-file => @ARGV[0], -format => 'fasta');
my $counter = 0;
my $db = @ARGV[1];
my $out_needle = @ARGV[2];
my $tmp_dir = @ARGV[3];
my $job_id = @ARGV[4];

while (my $obj = $is->next_seq()) {
	print(STDERR "$counter\n");

	@ARGV[0] =~ /.+\/(\S+)$/;
	my $fn = $1;
	my $temp_file = $tmp_dir . "/". $fn . "." . $counter;

	print("$temp_file\n");

	open(fh_out, ">" . $temp_file);
	print(fh_out ">" . $obj->display_id . "\n");
	print(fh_out $obj->seq . "\n");
	close(fh_out);

	$ENV{EMBOSS_ACDROOT} = $cwd . "/bin/acd";

	my $cmd = $cwd."bin/emboss.needle.static.linux.x86-64 -datafile $cwd/bin/EDNAFULL -asequence $temp_file -bsequence $db -auto -stdout > $out_needle";

	`$cmd`;

	my $out_identity = $out_needle . ".$counter" . ".pi";
	my $cmd2 = $cwd . "utils/perc_identity.pl $out_needle > $out_identity";

	`$cmd2`;

	my $cmd3 = "rm $out_needle $temp_file";
	`$cmd3`;

	$counter ++ ;
}
