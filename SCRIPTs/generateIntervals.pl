use strict;
use Bio::SeqIO;

my $file = shift;

validateParameters();

my $db = Bio::SeqIO->new(-file=>$file, -format=>"fasta");

my $countGlobal = 0;
my $start = 0;
my $end = 0;

open(OUT, ">$file"."_positions");
while(my $seq = $db->next_seq){
	my $count = $seq->length;
	$start = $countGlobal+1;
	$end = $start+$count-1;
	$countGlobal += $count;
	print OUT $seq->id,"\t$start\t$end\n";
}
close(OUT);



sub validateParameters{
        my $allExists = 1;
        my $haveAllParameters=1;

        unless(length $file > 0){
                $haveAllParameters = 0;
        }

        unless($haveAllParameters){
                print <<FOO;
Usage:
        perl generateIntervals.pl fasta_file_name
FOO
                exit 0;
        }

        unless (-e $file) {
                print "$file não existe\n";
                $allExists = 0;
        }

        exit 0 unless($allExists);


}

