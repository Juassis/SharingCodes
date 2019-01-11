#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

my $usage = "\nperl extract_sequence_by_length.pl file.fasta cutoff_length output_name\n";

my $file = shift or die(print "$usage\n");
my $length = shift or die(print "$usage\n");
my $output = shift or die(print "$usage\n");

my $seqin  = Bio::SeqIO->new(-file => "$file", -format => "fasta");

my $seqout = Bio::SeqIO->new(-file => ">$output", -format => "fasta");
while(my $seq = $seqin->next_seq) {
  if($seq->length >= $length) {
    $seqout->write_seq($seq);
  }
}
