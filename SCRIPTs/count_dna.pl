#!/usr/bin/perl
# While executing this script it asks for the file name of the DNA sequence. 
#24/04/2013 - Juliana Assis - Alterado cntagem de Ns.
##########################################################################

print "\n\n\t\#################### FREQUENCY OF NUCLIOTIDE #################### \n\n";

use strict;

my $countA = 0;
my $countC = 0;
my $countG = 0;
my $countT = 0;
my $countN = 0;
my $line;
my $sequence;
my $nt;
my $dna_filename;
my $dna;
my @dna;

print "PLEASE ENTER THE FILENAME OF THE DNA SEQUENCE:=";
chomp($dna_filename=<STDIN>);

open(DNAFILE,$dna_filename) or die "unable to open the file";
@dna=<DNAFILE>;
close DNAFILE;

foreach $line (@dna) {

# discard blank line
if ($line =~ /^\s*$/) {
next;

# discard comment line
} elsif($line =~ /^\s*#/) {
next;

# discard fasta header line
} elsif($line =~ /^>/) {
next;

# keep line, add to sequence string
} else {
$sequence .= $line;
}
}

# remove non-sequence data (in this case, whitespace) from $sequence string
$sequence =~ s/\s//g;
@dna=split("",$sequence); #splits string into an array
#print " \nThe original DNA file is:\n$sequence \n";

while(@dna){
$nt = shift (@dna);
if($nt =~/[A]/ig){
$countA++;
}
if($nt=~/[T]/ig){
$countT++;
}
if($nt=~/[C]/ig){
$countC++;
}
if($nt=~/[G]/ig){
$countG++;
}
if($nt=~/[N]/ig){
$countN++;
}
}

print "\nAdeninas:".$countA."\n";
print "Citosinas:".$countC."\n";
print "Guaninas:".$countG."\n";
print "Timinas:".$countT."\n";
print "Ns:".$countN."\n";

                          
