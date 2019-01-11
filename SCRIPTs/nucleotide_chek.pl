#!/usr/bin/env perl
#count Start and Stop Codons in nucleotides sequence
#
use strict;
use warnings;
use Bio::SeqIO;
my ($f1) = @ARGV;
my $in = Bio::SeqIO->new(-file => "<$f1");
my $out = Bio::SeqIO->new(-file => ">$f1.INTERNAL_STOPS.fa", -format=>"fasta");
printf("%s\n", join("\t", qw/ID FirstAA Length StartCheck StopCheck InternalStops/));
while(my $s = $in->next_seq){
	my $id = $s->primary_id;
	my $alpha = $s->alphabet;
	my $aa = $s->seq;
	my $start = ($aa =~ /^ATG/) || '0';
	my $end = ($aa =~ /([TAA][TAG][TGA])$/) || '0';
	$aa =~ s/^.//;
	$aa =~ s/.$//;
	my @stops =($aa =~ /(\.|\*)/g);
	my $stop = scalar @stops;
	printf("%s\t%s\t%d\tSTART=%s\tEND=%s\tINTERNAL=%s\n", $id, $alpha || '.', length($aa), $start, $end, $stop || "0");
	$out->write_seq($s) if $stop;
}
