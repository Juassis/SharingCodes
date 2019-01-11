#!/usr/bin/perl
#
use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use File::Basename;

my $PATTERN_CDS = "\^.*ID=(CDS:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id
my $PATTERN_TRANSCRIPT = "\^.*ID=(transcript:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id
my $PATTERN_PARENT_TRANSCRIPT = "\^.*Parent=(transcript:)?([a-zA-Z0-9\-_\.]+);?.*\$"; # pattern to extract the gene id

my $file_fasta = shift or die "Usage: perl $0 Genome.fasta Annotation.gff OutputDir\n";
my $gff_file = shift or die "Usage: perl $0 Genome.fasta Annotation.gff OutputDir\n";
my $output_dir = shift or die "Usage: perl $0 Genome.fasta Annotation.gff OutputDir\n";

my %hash_transcript = ();

$| = 1;    # Flush output

# my @suffixlist = ("fasta", "fna", "fa", "fas");
# my $new_prefix = basename($file_fasta,@suffixlist);

my $new_prefix = $file_fasta;
$new_prefix =~ s/^(.+\/)+//g;
# $new_prefix=~ s/\.(fasta|fna|fa|fas)$//g;

my $outfile_pep = Bio::SeqIO->new( -format => 'fasta', -file => ">$output_dir/$new_prefix" );

###### Output type description ######
# cds - translated sequence (starting with ATG and ending with a stop codon included)
# protein - cds translated (includes a * as the stop codon)

### First, index the genome
my $db = Bio::DB::Fasta->new($file_fasta);

print("File ", $file_fasta," parsed\n");
my %hash_gff = ();


read_gff_file();
sort_gff();

sub sort_gff{
  my @scaffolds = keys %hash_gff;
  for(@scaffolds){
    my $scaffold = $_;
    my @array_genes = sort {$hash_gff{$scaffold}{$a}{"MIN"} <=> $hash_gff{$scaffold}{$b}{"MIN"}} keys $hash_gff{$scaffold};
  	for(@array_genes){
      my $parent_id = $_;

        my @array_cds = sort {$a->{"START"} <=> $b->{"START"}} @{$hash_gff{$scaffold}{$parent_id}{"CDS"}};
        my $mergedCDS_seq = "";
        my $strand = "+";
        for(@array_cds){
          $strand = $_->{"STRAND"};
          my $cds_seq = $db->seq( $scaffold, $_->{"START"}, $_->{"END"} );
          $mergedCDS_seq .= $cds_seq;
          # print($_->{"TYPE"}, "\n");
        }
        my $output_cds = Bio::Seq->new(
            -seq        => $mergedCDS_seq,
            -id         => $parent_id,
            -display_id => $parent_id,
            -alphabet   => 'dna',
        );
        if ($strand eq "-") {
            $output_cds = $output_cds->revcom();
        }
        #translate CDS to peptide for protein sequence
        my $output_pep = $output_cds->translate();
        #write to file
        if (length($mergedCDS_seq) != 0) {
              $outfile_pep->write_seq($output_pep);
        }
      # if(exists $hash_random{$parent_id}){
      # } # end if if(exists $hash_random{$parent_id}){
      # else{
      #   print("ERROR: CDS without a valid mRNA: $parent_id\n");
      # }

    } # end 	for(@array_genes){
  }# end for(@scaffolds){
}

#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
#########################################################################
sub read_gff_file{
# primeiro a gente abre o arquivo

  open(IN, $gff_file);
  # my %hash_mrna_by_scaffolds = ();
  my $transcript_id = "";
  while(<IN>){
    chomp;

    unless(/^#/ || /^\s*$/){  # To avoid lines which start with '#'
        my $line = $_;
        my @cols = split /\t/;

        my $gene_name = $cols[-1]; # The gene name is at last column
        my $scaffold = $cols[0]; # The scaffold name is the first column
        my $type = $cols[2];
        my $gene_start = $cols[3];
        my $gene_end = $cols[4];
        my $strand = $cols [6];

				my $parent_id = get_feature_id($gene_name, $PATTERN_PARENT_TRANSCRIPT);

        if (lc($type) eq 'cds') {
          my $feature_id = get_feature_id($gene_name, $PATTERN_CDS);

          my %hash_cds = ();
          $hash_cds{"START"} = int($gene_start);
          $hash_cds{"END"} = int($gene_end);
          $hash_cds{"STRAND"} = $strand;
          $hash_cds{"ID"} = $feature_id;
          $hash_cds{"TYPE"} = $transcript_id;

          if(exists $hash_gff{$scaffold}{$parent_id}){
            my $array_cds = $hash_gff{$scaffold}{$parent_id}{"CDS"};
            if($gene_start < $hash_gff{$scaffold}{$parent_id}{"MIN"}){
              $hash_gff{$scaffold}{$parent_id}{"MIN"} = $gene_start;
            }
            push(@$array_cds, \%hash_cds);

          }else{
            $hash_gff{$scaffold}{$parent_id}{"MIN"} = $gene_start;
            my @array_cds = (\%hash_cds);
            $hash_gff{$scaffold}{$parent_id}{"CDS"} = \@array_cds;
          }
        }else{
          my $feature_id = get_feature_id($gene_name, $PATTERN_TRANSCRIPT);
          if($feature_id ne ""){
            $transcript_id = $type;
          }
          # if(exists $hash_transcript{$feature_id}){
          #   # print("ERROR: ")
          # }else{
          #   $hash_transcript{$feature_id} = 1;
          # }
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
