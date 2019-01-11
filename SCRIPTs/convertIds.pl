#!/usr/bin/perl -w
# ex: perl convertIds.pl Thr.fasta listaoficialIDS out
use strict;

my $Fasta = shift; 
my $Ids = shift;
my $OutFile = shift;

my @array_ids=();

open(IDS, $Ids); 

while(<IDS>){
	chomp;
	push(@array_ids, $_);
}

close(IDS);

#for(@array_ids){
	#print "$_\n";
#}


open(FASTA, $Fasta);
open(OUT, ">".$OutFile);


while(<FASTA>){
	if(/^>/){
		my $id = shift(@array_ids);
		print OUT "$id\n";
	}else{
		print OUT $_;
	}
}

close(OUT);
close(FASTA);
