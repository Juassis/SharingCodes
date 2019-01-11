#!/usr/bin/perl
#
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my $PATTERN_TRANSCRIPT = "\^.*ID=(transcript:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id


my $gff_file;
my $accepted_features;

GetOptions(
  'i=s'    => \$gff_file,
	'a=s'    => \$accepted_features
);
$| = 1;    # Flush output

# my @suffixlist = ("fasta", "fna", "fa", "fas");
# my $new_prefix = basename($file_fasta,@suffixlist);
my %hash_features = map{$_ => 1} split(/,/, $accepted_features);
read_gff_file();
#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
#########################################################################
sub read_gff_file{
# primeiro a gente abre o arquivo

  open(IN, $gff_file);
  # my %hash_mrna_by_scaffolds = ();
  while(<IN>){
    chomp;

    unless(/^#/ || /^\s*$/){  # To avoid lines which start with '#'
        my @cols = split /\t/;

        my $last_column = $cols[-1];
        my @array_last_col = split(";", $last_column);
        my $bedvalues = $array_last_col[-1];
        $bedvalues =~ s/bedtools=//g;

        my $type = $cols[2];

        if (exists $hash_features{$type}) {
          my $feature_id = get_feature_id($last_column, $PATTERN_TRANSCRIPT);
          print($feature_id, "\t", $bedvalues, "\n");
        }
    } # end unless(/^#/){
  }# end while(<IN>){

    # print(Dumper(%hash_gff));

  close(IN);
}

#########################################################################
# Extract the gene id based on the global variable pattern
#########################################################################
sub get_feature_id{
	my $attribute_text = shift;
  my $pattern = shift;
	my $regex = qr/$pattern/;
	if($attribute_text =~ $regex){
		$attribute_text =~ s/$regex/$2/g;
	}else{
    $attribute_text = "";
  }
	return $attribute_text;
}
