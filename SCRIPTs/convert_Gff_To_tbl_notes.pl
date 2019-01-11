use strict;
no warnings;
use Getopt::Long;

# Some constants to easily change our columns index if necessary
use constant{
	FEATURE => 2,
	POS_START => 3,
	POS_END => 4,
	STRAND => 6,
	NOTE => -1,
	FEATURES => 'n',
	GENE_START => 's',
	GENE_END => 'e',
	GENE_STRAND => 't'
};

#########################################################################
# First we will declare our global variables to use in this script
#########################################################################
my $file_input_gff; # gff file to convert into the tbl
my $file_input_notes; # file with two columns with: 1st column: Gene name; 2nd column: Note to put into the tbl file
my $file_output; # Name of the output file
my $help; # to control if the user want help of the script

my %hash_scaffolds=(); # hash to store all scaffolds in the gff input
my @array_scaffolds_order=(); # array to store the order of the scaffold read from the file

my %hash_genes=(); # hash with all informations about features, genes and scaffolds
my %hash_notes = (); # hashtable with the relation GENE_NAME -> NOTES from the notes file 

#########################################################################
# Second we will parse the parameters provided by the user using the
# function GetOptions from the package Getopt::Long
#########################################################################
GetOptions(
	'i=s'    => \$file_input_gff,
	'n=s'    => \$file_input_notes,
	'o=s'    => \$file_output,
	'help|?' => \$help
);

#########################################################################
# Our script is divided into subroutines. The main() subroutine call the
# another subroutines. Tha main subroutine is just to organize our script
#########################################################################
main();

#########################################################################
# The main subroutine
#########################################################################
sub main {
	# We validate the parameters provided by the user since all parameters are mandatories
	validate_parameters();
	
	print "Reading $file_input_notes...\n";
	read_note_file(); # Calling the subroutine to read the notes file
	print "Reading and organising $file_input_gff...\n";
	read_gff_file(); # Calling the subroutine to read the GFF file
	sort_cds(); # Calling the subroutine to sort the CDS in each gene
	print "Generating $file_output...\n";
	generate_tbl(); # Calling the subroutine to export our output file
	print "Output generated [$file_output]...\n";
	print "Process finished!\n";

}

#########################################################################
# Reads a two columns file where the first column is the key of the hashtable
# and the second column is the note to be used
#########################################################################
sub read_note_file{
	open(IN, $file_input_notes);
	
	while(<IN>){
		chomp;
		unless(/^#/){ # To avoid lines which start with '#'
			unless(/^\s*$/){ # This is to avoid lines with no relevant content
				my @cols = split /\t/;
				$hash_notes{$cols[0]} = $cols[1]; # The key of our hash is the first column (Gene name)
												  # Tha value of the hash is the note
			}
		}
	}
	
	close(IN);
}

#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
#########################################################################
sub read_gff_file{
	
	open(IN, $file_input_gff);
	
	# HtScaffold0001  GeneMark.hmm    start_codon     654     656     .       +       .       codon_start=1;GFF Attribute="gene_id "HtGene0001"";GFF Attribute="transcript_id """;GFF Attribute="protein_id "HtProtein0001""
	
	while(<IN>){
		chomp;
		unless(/^#/){  # To avoid lines which start with '#'
			unless(/^\s*$/){  # This is to avoid lines with no relevant content
				my $line = $_;
				my @cols = split /\t/;
				my $gene_name = $cols[-1]; # The gene name is at last column
				my $scaffold = $cols[0]; # The scaffold name is the first column
				$gene_name =~ s/^.+gene_id "(\w+)"";GFF.+/$1/g; # We remove the garbage content at the last column
				
				# This is to control the names of the scaffolds and the genes inside each scaffold
				if(exists $hash_scaffolds{$scaffold}){
					my $array = $hash_scaffolds{$scaffold};
					unless($gene_name ~~ @$array){
						push(@$array, $gene_name)
					}
					
				}else{
					my @array = ($gene_name);
					$hash_scaffolds{$scaffold} = \@array;
					
					push(@array_scaffolds_order, $scaffold);
				}
				
				# For each line of the GFF we organize the info the levels: scaffold -> gene name -> [FEATURES]
				# For each line we have another hash with the information about start, end, strand and feature
				# Our final hash we will have this structure:
				# %{hash key=scaffold name} = @{array with informations about genes inside a scaffold} 
				# For each gene we will have another hash with the structure:
				# %{hash keys: gene name and FEATURES} = @{array with all features from this gene}
				# %{hash keys: gene name and GENE_START} = the start position of this gene
				# %{hash keys: gene name and GENE_END} = the end position of this gene 
				# %{hash keys: gene name and GENE_STRAND} = the strand of this gene
				# %{hash keys: gene name and NOTE} = the note of this gene
				if(exists $hash_genes{$scaffold}{$gene_name}{FEATURES}){
					my $array_genes = $hash_genes{$scaffold}{$gene_name}{FEATURES};
					my $hash_line = get_hash_from_line($line);
					
					push(@$array_genes, $hash_line);
					
					# we consider the gene start position as the minor position of all features
					if($hash_genes{$scaffold}{$gene_name}{GENE_START} > $hash_line->{POS_START}){
						$hash_genes{$scaffold}{$gene_name}{GENE_START} = $hash_line->{POS_START};
					}
					# we consider the gene end position as the major position of all features
					if($hash_genes{$scaffold}{$gene_name}{GENE_END} < $hash_line->{POS_END}){
						$hash_genes{$scaffold}{$gene_name}{GENE_END} = $hash_line->{POS_END};
					}
					
				}else{
					my $hash_line = get_hash_from_line($line);
					my @array_genes = ($hash_line);
					
					$hash_genes{$scaffold}{$gene_name}{FEATURES} = \@array_genes;
					$hash_genes{$scaffold}{$gene_name}{GENE_START} = $hash_line->{POS_START};
					$hash_genes{$scaffold}{$gene_name}{GENE_END} = $hash_line->{POS_END};
					$hash_genes{$scaffold}{$gene_name}{GENE_STRAND} = $hash_line->{STRAND};
					$hash_genes{$scaffold}{$gene_name}{NOTE} = $hash_notes{$gene_name};
				}
			}
		}
	}
	
	close(IN);
	
}

#########################################################################
# This subroutine  sorts the features from all genes to keep the features in
# crescent order
#########################################################################
sub sort_cds{
	for(@array_scaffolds_order){
		my $scaffold = $_;
		my $array_genes_order = $hash_scaffolds{$scaffold};
		for(@$array_genes_order){
			if($hash_genes{$scaffold}{$_}{GENE_STRAND} eq "+"){
				my @array_ordered = sort {$a->{POS_START} <=> $b->{POS_START} || $a->{POS_END} <=> $b->{POS_END} }  @{$hash_genes{$scaffold}{$_}{FEATURES}};
				$hash_genes{$scaffold}{$_}{FEATURES} = \@array_ordered;	
			}else{
				my @array_ordered = sort {$a->{POS_END} <=> $b->{POS_END} || $a->{POS_START} <=> $b->{POS_START} }  @{$hash_genes{$scaffold}{$_}{FEATURES}};
				$hash_genes{$scaffold}{$_}{FEATURES} = \@array_ordered;
			}	
		}
	}
}

#########################################################################
# This subroutine export the output in tbl format 
#########################################################################
sub generate_tbl{
	for(@array_scaffolds_order){
		open(OUT, ">".$file_output."/".$_.".tbl");
		my $scaffold = $_;
		print OUT ">Features $scaffold\n";
		my $array_genes_order = $hash_scaffolds{$scaffold};
		for(@$array_genes_order){
			my $gene_name = $_;
			if($hash_genes{$scaffold}{$gene_name}{GENE_STRAND} eq "+"){
				print OUT sprintf("%d\t%d\tgene\n", $hash_genes{$scaffold}{$gene_name}{GENE_START}, $hash_genes{$scaffold}{$gene_name}{GENE_END});	
			}else{
				print OUT sprintf("%d\t%d\tgene\n", $hash_genes{$scaffold}{$gene_name}{GENE_END}, $hash_genes{$scaffold}{$gene_name}{GENE_START});
			}
			print OUT sprintf("\t\t\tlocus_tag\t%s\n", $gene_name);
			for(@{$hash_genes{$scaffold}{$gene_name}{FEATURES}}){
				my $feature_info = $_;
				if($feature_info->{FEATURE} eq "CDS"){
					if($feature_info->{STRAND} eq "+"){
						print OUT sprintf("%d\t%d\tCDS\n", $feature_info->{POS_START}, $feature_info->{POS_END});
					}else{
						print OUT sprintf("%d\t%d\tCDS\n", $feature_info->{POS_END}, $feature_info->{POS_START});
					}
				}
			}
			print OUT sprintf("\t\t\tnote\t%s\n", $hash_genes{$scaffold}{$gene_name}{NOTE});
		}
		close(OUT);
	}
}

#########################################################################
# Splits a line from gff file and transform into a hash
#########################################################################
sub get_hash_from_line{
	my $line = shift;
	my %hash = ();
	my @cols = split(/\t/, $line);
	$hash{FEATURE} = $cols[FEATURE];
	$hash{POS_START} = $cols[POS_START];
	$hash{POS_END} = $cols[POS_END];
	$hash{STRAND} = $cols[STRAND];
	
	return \%hash;
}


#########################################################################
# Just a report to see strange things at the gff
#########################################################################
sub print_report{
	
	my $count_total = 0;
	my $count_good = 0;
	my $count_two_start = 0;
	my $count_two_stop = 0;
	my $count_only_cds = 0;
	my $count_zero_cds = 0;
	my $count_one_start_zero_stop = 0;
	my $count_one_stop_zero_start = 0;
	for(@array_scaffolds_order){
		my $scaffold = $_;
		my $array_genes_order = $hash_scaffolds{$scaffold};
		for(@$array_genes_order){
			my $gene_name = $_;
			my $has_start = 0;
			my $has_stop = 0;
			my $has_cds = 0;
			for(@{$hash_genes{$scaffold}{$gene_name}{FEATURES}}){
				my $feature_info = $_;
				my @cols = split /\t/;
				if($feature_info->{FEATURE} eq "start_codon"){
					$has_start++;
				}elsif($feature_info->{FEATURE} eq "stop_codon"){
					$has_stop++;
				}elsif($feature_info->{FEATURE} eq "CDS"){
					$has_cds++;
				}
			}
			if($has_start == 1 and $has_stop == 1 and $has_cds >=1){
				$count_good++;
			}
			if($has_start >=2){
				$count_two_start++;
			}
			if($has_stop >=2){
				$count_two_stop++;
			}
			if($has_start == 0 and $has_stop == 0 and $has_cds >=1){
				$count_only_cds++;
			}
			
			if($has_start == 1 and $has_stop == 1 and $has_cds == 0){
				$count_zero_cds++;
			}
			
			if($has_start == 1 and $has_stop == 0  and $has_cds >=1){
				$count_one_start_zero_stop++;
			}
			if($has_start == 0 and $has_stop == 1 and $has_cds >=1){
				$count_one_stop_zero_start++;
			}
			
			$count_total++;
		}
	}
	print "Good: $count_good - $count_total\n";
	print "Two start: $count_two_start\n";
	print "Two stop: $count_two_stop\n";
	print "Only CDS: $count_only_cds\n";
	print "Zero CDS: $count_zero_cds\n";
	print "One Start - Zero Stop: $count_one_start_zero_stop\n";
	print "Zero Start - One Stop: $count_one_stop_zero_start\n";
	
	
}

#########################################################################
# Subroutine to validate the parameters
#########################################################################
sub validate_parameters {
	my $allExists  = 1;
	my $fileExists = 1;

	if ( defined $help ) {
		print usage();
		exit 0;
	}

	unless ( defined $file_input_gff ) {
		$allExists = 0;
	}
	
	unless ( defined $file_input_notes ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output ) {
		$allExists = 0;
	}

	unless ( defined $file_output ) {
		$allExists = 0;
	}

	if ($allExists) {
		unless ( -f $file_input_gff ) {
			print STDERR "$file_input_gff doesn't exists.\n";
			$fileExists = 0;
		}
		unless ( -f $file_input_notes ) {
			print STDERR "$file_input_notes doesn't exists.\n";
			$fileExists = 0;
		}
	}
	else {
		print usage();
		exit 0;
	}

	unless ($fileExists) {
		print STDERR "Program execution aborted.\n";
		exit 0;
	}

}

sub usage {
	my $usage = <<FOO;
	
Usage:
	perl $0 -i input_gff -n input_notes -o output_dir
FOO
	return $usage;

}
