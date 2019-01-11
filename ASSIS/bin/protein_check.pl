#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

my ($f1) = @ARGV;
my $in = Bio::SeqIO->new(-file => "<$f1");

printf("%s\n", join("\t", qw/ID	Length StartCheck StopCheck InternalStops/));
while(my $s = $in->next_seq){
	my $id = $s->primary_id;
	my $alpha = $s->alphabet;
	my $aa = $s->seq;

	my $start = ($aa =~ /^M/) || 0;
	my $end = ($aa =~ /([^A-Z])$/) || 0;

	my @stops = ($aa =~ /(\.|\*)/g);
	my $inner_stop; #in-frame stop codons

	if($end){
		$inner_stop = scalar @stops - 1;
	}else{
		$inner_stop = scalar @stops;
	}

	printf("%s\t%d\t%s\t%s\t%s\n", $id, $s->length, $start, $end, $inner_stop);
}
