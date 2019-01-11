#!/usr/bin/perl
######################################################################
#######################################################################
# Copyright 2011 Fundação Oswaldo Cruz
# Author: Eric Aguiar and Juliana Assis
# Alteração: Inserção maiores que 500pb - 01/10/12
# template.pl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
#
# template.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with template.pl (file: COPYING).
# If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
#######################################################################


use Getopt::Long;
use Bio::SeqIO;

my $usage = "

$0 -i <input file>
$0 -h

-i <input file>		: Fasta Input file
-h			: Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.

my $inputFile;
my $outputFile;
my %hash;
GetOptions ("i=s" => \$inputFile,
			"h!" => \$help);

if ($help) {
	die $usage;
}

if (not(defined($inputFile))) {
	die "\nGive an input file name",$usage;
}


our $sum=0;
my $div;
our $nContigs;
my $nMaior100b=0;
my $totBaseMaior100b=0;
my $nMaior200b=0;
my $totBaseMaior200b=0;
my $nMaior500b=0;
my $totBaseMaior500b=0;
my $nMaior1k=0;
my $totBaseMaior1k=0;
my @v;
my $j=0;
readFile();
our $average = $sum / $nContigs;
calc();

sub readFile{

	#print "\nLoading Contigs... \n";
	my $seqio_object = Bio::SeqIO->new(-file => $inputFile);
	while($seq_object = $seqio_object->next_seq){

		$nContigs++;
		$hash{$seq_object->id}=length($seq_object->seq);
		$totBase+=length($seq_object->seq);
		$hash{$seq_object->id}=length($seq_object->seq);
		$v[$j]=length($seq_object->seq);
		$j++;

		if (length($seq_object->seq)> 100){
			$nMaior100b++;
			$totBaseMaior100b+=length($seq_object->seq);
				if (length($seq_object->seq)> 200){
                        		$nMaior200b++;
                        		$totBaseMaior200b+=length($seq_object->seq);
                                        if (length($seq_object->seq)> 500){
                        		$nMaior500b++;
                        		$totBaseMaior500b+=length($seq_object->seq);
						if (length($seq_object->seq)> 1000){
                                        	$nMaior1k++;
                                        	$totBaseMaior1k+=length($seq_object->seq);
                                               }
					}
				}
		}
		$sum+=length($seq_object->seq);
	}
	#print "Close file\n";
	$div = $sum /2;
	#print "Sum: $sum\n Sum /2 :$div\n";
}

sub calcMediana{
	@v2= sort{$a <=> $b} @v;
	$pm=$nContigs / 2;
	$parImpar=$nContigs %2;
#	print "@v2\n $pm -$nContigs -  $parImpar\n";
	if ($nContigs%2 == 0 ){
#	print "mediana:".$v2[$pm]."\n";
		return $v2[$pm];
}	else{
		$pm = int($pm) +1;
		return $v2[$pm];
	}


}
sub hashValueAscendingNum {
   $hash{$a} <=> $hash{$b};
}
sub hashValueDescendingNum {
   $hash{$b} <=> $hash{$a};
}

sub calc{

#print "\nGRADES IN DESCENDING NUMERIC ORDER:\n";
my $value=0;
$average = $sum / $nContigs;
my $largest;
foreach $key (sort hashValueDescendingNum (keys(%hash))) {
  #print "\t$hash{$key} \t\t $key\n";
   if ($value==0){
      $largest = $hash{$key};
   }
  # print "$hash{$key} \n";
   $value+=$hash{$key};
   if ($value > $div){
	   my $perc100b=($nMaior100b*100)/$nContigs;
	   my $perc200b=($nMaior200b*100)/$nContigs;
           my $perc500b=($nMaior500b*100)/$nContigs;
	   my $perc1k=($nMaior1k*100)/$nContigs;
	   print "Execution time (hours)\t?\nMMU(GB)\t?\nNumber of contigs\t".$nContigs."\nNumber of bases in all contigs\t".$totBase."\nN50\t".$hash{$key}."\nLongest contig\t".$largest."\nMedian contig size\t".calcMediana()."\nNumber of contigs > 100 b\t".$nMaior100b."\nTotal of bases into contigs > 100 b\t".$totBaseMaior100b."\n% contigs > 100 b\t".$perc100b."\nNumber of contigs > 200 b\t".$nMaior200b."\nTotal of bases into contigs > 200 b\t".$totBaseMaior200b."\n% contigs > 200 b\t".$perc200b."\nNumber of contigs > 500 b\t".$nMaior500b."\nTotal of bases into contigs > 500 b\t".$totBaseMaior500b."\n% contigs > 500 b\t".$perc500b."\nNumber of contigs > 1 kb\t".$nMaior1k."\nTotal of bases into contigs > 1 kb\t".$totBaseMaior1k."\n% contigs > 1 kb\t".$perc1k."\n";
   		die;
   }


}

}
