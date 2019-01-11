use strict;
use Bio::SeqIO;

my $file = shift;

my $db = Bio::SeqIO->new(-file=>$file, -format=>"fasta");

my $countGlobal = 0;

open(OUT, ">$file"."_length");
while(my $seq = $db->next_seq){
        my $count = $seq->length;
        print OUT $seq ->id, "\t$count\n"
}
close(OUT);

