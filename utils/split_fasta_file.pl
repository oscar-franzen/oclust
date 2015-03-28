#!/usr/bin/perl

use strict 'vars';
use Cwd 'abs_path';

sub workaround {
    my $p = abs_path($0);
    $p =~ /^(\S+\/)utils$/;
    $p = $1;

    my $path = $p . "modules";
    return $path;
}

use lib workaround();
use Bio::SeqIO;

use Bio::SeqIO;

### LOCALS
my $infile    = $ARGV[0];
my $size    = $ARGV[1];
my $i        = 0;
my $j        = 1;

### OBJECT INSTANTIATION
my $in    = Bio::SeqIO->new(
            -file    => $infile,
            -format    => 'Fasta',
            );
my $out    = Bio::SeqIO->new(
            -file    => ">$infile"."_$j.fasta",
            -format    => 'Fasta',
            );

### PROCESS FILES
while ( my $seq = $in->next_seq() ) {

    $i++;

    if ($i == $size) {
        $j++;
        $out = Bio::SeqIO->new(
                    -file    => ">$infile"."_$j.fasta",
                    -format    => 'Fasta');
        $i = 0;
        }

    $out->write_seq($seq);

    }
