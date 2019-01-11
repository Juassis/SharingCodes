#!/usr/bin/env perl

# GFF
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

# Some constants to easily change the columns index if necessary
use constant{
	CHROMOSOME => 0,
	SOURCE => 1,
	FEATURE => 2, # the column from the GFF that represents the feature (cdna, mrna, gene, etc)
	POS_START => 3, # the column that represents the start position of the feature
	POS_END => 4, # the column that represents the end position of the feature
	SCORE => 5,
	STRAND => 6, # the column from the GFF that represents the dna strand (plus or minus),
	PHASE => 7,
	ATTRIBUTES => 8,
	GENE_START => 's', # this variable is the start key of the hash with the genes
	GENE_END => 'e',# this variable is the end key of the hash with the genes
	GENE_STRAND => 't', # this variable is the strand key of the hash with the genes
	GENE_PARENT => 'p', # this variable is the strand key of the hash with the genes
};

#########################################################################
## Global variables to use in this script
#########################################################################

my $file;
my $accepted_features;
my $file_output;

my %hash_mrna_by_scaffolds;
my %hash_start_end_cds = ();


my $PATTERN_CDS = "\^.*ID=(CDS:)?([a-zA-Z0-9\-_]+(\.\\d+)?)(\.cds)?;?.*\$"; # pattern to extract the gene id
my $PATTERN_TRANSCRIPT = "\^.*ID=(transcript:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id
my $PATTERN_PARENT_TRANSCRIPT = "\^.*Parent=(transcript:)?([a-zA-Z0-9\-_\.]+);?.*\$"; # pattern to extract the gene id

my %hash_features = ();

###########################################################################
#receiving parameters
###########################################################################
GetOptions(
  'i=s'    => \$file,
	'a=s'    => \$accepted_features
);

$| = 1;

main();


###########################################################################
sub main {
	%hash_features = map{$_ => 1} split(/,/, $accepted_features);

	read_gff_file(); # Calling the subroutine to read the GFF file
	# print(Dumper(%hash_start_end_cds));
	sort_cds(); # Calling the subroutine to sort the CDS in each gene
	geting_overlap_in_gene();

}

#########################################################################
# This subroutine  sorts the features from all genes to keep the features in
# crescent order.
#########################################################################
sub sort_cds{
	my @scaffolds = keys %hash_mrna_by_scaffolds;
	for(@scaffolds){
		my $scaffold = $_;
		my $array_mrnas = $hash_mrna_by_scaffolds{$scaffold};
		my @array_ordered = sort {$a->{GENE_START} <=> $b->{GENE_START} || $a->{GENE_END} <=> $b->{GENE_END} }  @$array_mrnas;
		for(@array_ordered){
			my $hash_mrna = $_;
			if($hash_mrna->{"GENE_START"} < $hash_start_end_cds{$hash_mrna->{"ID"}}{"START"}){
				$hash_mrna->{"GENE_START"} = $hash_start_end_cds{$hash_mrna->{"ID"}}{"START"};
			}
			if($hash_mrna->{"GENE_END"} > $hash_start_end_cds{$hash_mrna->{"ID"}}{"END"}){
				$hash_mrna->{"GENE_END"} = $hash_start_end_cds{$hash_mrna->{"ID"}}{"END"};
			}
		}
		$hash_mrna_by_scaffolds{$scaffold} = \@array_ordered;
	} # end for(@scaffolds)

}# end sort_cds
#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
#########################################################################
sub read_gff_file{
# primeiro a gente abre o arquivo

  open(IN, $file);
  # my %hash_mrna_by_scaffolds = ();
  while(<IN>){
    chomp;

    unless(/^#/ || /^\s*$/){  # To avoid lines which start with '#'
        my $line = $_;
        my @cols = split /\t/;

        my $gene_name = $cols[-1]; # The gene name is at last column
        my $scaffold = $cols[0]; # The scaffold name is the first column
        my $strand = $cols [STRAND];

        my $gene_id = get_feature_id($gene_name, $PATTERN_TRANSCRIPT);
				my $parent_id = get_feature_id($gene_name, $PATTERN_PARENT_TRANSCRIPT);

        my $feature = $cols[FEATURE];
				my $attributes = $cols[ATTRIBUTES];

        # just in case: when the start of the gene is greater than the end position we invert
        if($cols [POS_START] > $cols [POS_END]){
          my $aux = $cols [POS_START];
          $cols [POS_START] = $cols [POS_END];
          $cols [POS_END] = $aux;
        }

        if(exists $hash_features{$feature}){
          my %hash_mrna = ();
          $hash_mrna{ID} = $gene_id;
          $hash_mrna{STRAND} = $strand;
          $hash_mrna{GENE_START} = $cols [POS_START];
          $hash_mrna{GENE_END} = $cols [POS_END];
					$hash_mrna{GENE_PARENT} = $parent_id;

          if (not exists $hash_mrna_by_scaffolds {$scaffold}){
            my @array_mrnas = ();
            push @array_mrnas, \%hash_mrna;
            $hash_mrna_by_scaffolds{$scaffold} = \@array_mrnas;
          }else{
            my $array_mrnas = $hash_mrna_by_scaffolds {$scaffold};
            push @$array_mrnas, \%hash_mrna;
          }
        }elsif(lc($feature) eq "cds"){
					# my $cds_id = get_feature_id($gene_name, $PATTERN_CDS);
					if(exists $hash_start_end_cds{$parent_id}){
						if($cols [POS_START] < $hash_start_end_cds{$parent_id}{"START"}){
							$hash_start_end_cds{$parent_id}{"START"} = $cols [POS_START] ;
						}
						if($cols [POS_END] > $hash_start_end_cds{$parent_id}{"END"}){
							$hash_start_end_cds{$parent_id}{"END"} = $cols [POS_END];
						}
					}else{
						$hash_start_end_cds{$parent_id}{"START"} = $cols [POS_START];
						$hash_start_end_cds{$parent_id}{"END"} = $cols [POS_END];
					}
				}
    } # end unless(/^#/){
  }# end while(<IN>){
  close(IN);
}

#########################################################################
# Sub overlap in genes
#########################################################################
sub geting_overlap_in_gene{
	print("ID\tOVERLAP_WITH\n");
	while(my ($scaffold, $array_mrna) = each(%hash_mrna_by_scaffolds)){
		for(my $i = 0; $i<scalar @{$array_mrna}-1; $i++){
			my $gene_i = $array_mrna->[$i];
			for(my $j = $i+1; $j<scalar @{$array_mrna}; $j++){
				my $gene_j = $array_mrna->[$j];
				if($gene_i->{"GENE_END"} >= $gene_j->{"GENE_START"}){
					if($gene_i->{"GENE_PARENT"} ne $gene_j->{"GENE_PARENT"}){
						print($gene_i->{"ID"}, "\t", $gene_j->{"ID"}, "\n");
					}
				}
			}
		}#end for(my $i = 0; $i< scalar @{$array_mrna}; $i++){
	}#end while(my ($scaffold, $array_mrna) = each(%hash_mrna_by_scaffolds)){
}

########################################################################
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
