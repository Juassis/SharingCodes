# You may distribute this module under the same terms as perl itself
# Cobertura das reads IsoSeq nos genes completos - Caso 1

#########################################################################
## First we will declare our global variables to use in this script
#########################################################################
use strict;
use Getopt::Long;
use Data::Dumper;

my $file_rnaseq;    # input file
my $file_gff;       # input file

my $file_output;
my $help;

my @posicao_rnaseq = ();
my @posicao_gene   = ();

my @genes_with_transcripts     = ();
my %hash_gene_with_transcripts = ();

###########################################################################
#receiving parameters
###########################################################################
GetOptions(
	's=s' => \$file_rnaseq,
	'g=s' => \$file_gff,

	'o=s'    => \$file_output,
	'help|?' => \$help
);

$| = 1;

main();
##################################################################################
# Our script is divided into subroutines. The main() subroutine call the
# another subroutines. Tha main subroutine is just to organize our script
##################################################################################
sub main {

	#validating parameters
	validate_parameters();

	print "Reading rnaseq file...\n";
	read_rnaseq();    #reading rnaseq input file and generating output

	print "Reading GFF file...\n";
	read_gff();

	print "Crossing information...\n";
	search_transcripts_per_gene();

	print "Exporting genes...\n";
	export_genes();

}

############################################################################
# Fazer a leitura do arquivo BED-BAM
############################################################################
sub read_rnaseq {

	# opening the input file
	open( IN, $file_rnaseq );

	my $count     = 0;
	my $wrap_line = 0;
	my $line;

	while (<IN>) {
		chomp;
		unless (/^(\s|\t)*$/) {
			unless (/^#/) {

				#   0   1       2       3       4       5       6
				#CHROM  Began     END      READ_Name     ALT     STRAND

				my @columns = split /\t/;

				my $id_scaffold = $columns[0];

				#               my $id_read = $columns [3];
				#print $id_read; pegando informacao da read
				my $posicaobegin = $columns[1];
				my $posicaoend   = $columns[2];

				if ( $posicaobegin > $posicaoend ) {
					my $aux = $posicaoend;
					$posicaoend   = $posicaobegin;
					$posicaobegin = $aux;
				}

				my %hash_rnaseq = ( 'BEGIN' => $posicaobegin, 'END' => $posicaoend );

				push( @posicao_rnaseq, \%hash_rnaseq );

			}

		}
	}

	@posicao_rnaseq = sort { $a->{'BEGIN'} <=> $b->{'BEGIN'} || $a->{'END'} <=> $b->{'END'} } @posicao_rnaseq;
}

###############################################
#GFF file
##############################################
sub read_gff {

	open( IN, $file_gff );

	my $count     = 0;
	my $wrap_line = 0;
	my $line;

	while (<IN>) {
		chomp;
		unless (/^(\s|\t)*$/) {
			unless (/^#/) {

				#0  1   2   3   4   5   6   7   8
				#ID Score   Feature Start   End Score   Strand  Phase   Attributes
				my @columns = split /\t/;

				my $id_scaffold = $columns[0];

				my $feature = $columns[2];

				#Alterado para pegar CDS

				#Arquivo add.utrs.gff
				#Parece o mais certo
				#eq "gene" eq "CDS"

				if ( lc($feature) eq "gene" ) {
					my $gene = $columns[8];

					$gene =~ s/ID=(.+);Name=.*$/$1/g;
					my $begin = $columns[3];
					my $end   = $columns[4];

					if ( $begin > $end ) {
						my $aux = $end;
						$end   = $begin;
						$begin = $aux;
					}
					my %hash_gene = ( 'ID' => $gene, 'BEGIN' => $begin, 'END' => $end );

					push( @posicao_gene, \%hash_gene );

				}

			}

		}
	}

	@posicao_gene = sort { $a->{'BEGIN'} <=> $b->{'BEGIN'} } @posicao_gene;
}

#########################################################################
#Transcripts per gene
#########################################################################
sub search_transcripts_per_gene {

	for ( my $i = 0; $i < scalar @posicao_gene; $i++ ) {

		my $gene    = $posicao_gene[$i];
		my $gene_id = $posicao_gene[$i]->{'ID'};
		for ( my $j = 0; $j < scalar @posicao_rnaseq; $j++ ) {
			my $transcript = $posicao_rnaseq[$j];
			if ( $transcript->{'BEGIN'} >= $gene->{'BEGIN'} and $transcript->{'END'} <= $gene->{'END'} ) {
				if ( exists $hash_gene_with_transcripts{$gene_id} ) {
					my $array = $hash_gene_with_transcripts{$gene_id};
					push( @$array, $transcript );
				}
				else {
					my @array = ($transcript);
					$hash_gene_with_transcripts{$gene_id} = \@array;
					$hash_gene_with_transcripts{'POS'}{$gene_id} = $i;
					push( @genes_with_transcripts, $gene_id );
				}
			}
		}
	}
}

#########################################################################
#Gerando output
#########################################################################
sub export_genes {
	open(OUT, ">".$file_output);
	for (@genes_with_transcripts) {
		my $gene_id              = $_;
		my $gene                 = $posicao_gene[ $hash_gene_with_transcripts{'POS'}{$gene_id} ];
		my $gene_begin           = $gene->{'BEGIN'};
		my $gene_end             = $gene->{'END'};
		my $gene_length          = $gene_end - $gene_begin + 1;
		my @array_gene_positions = (1) x $gene_length;

		my $array_transcripts = $hash_gene_with_transcripts{$gene_id};
		for (@$array_transcripts) {
			my $transcript = $_;
			change_gene_array( \@array_gene_positions, $gene, $transcript );
		}
		my $coverage = calculate_coverage(\@array_gene_positions);
		if($coverage == 100){
			print OUT $gene_id."\n";
		}
	}
	close(OUT);
}
#########################################################################

sub calculate_coverage{
	my $array_gene_positions = shift;
	my $count_uncovered = 0;
	my $array_length = scalar @$array_gene_positions;
	for(my $i = 0; $i < $array_length; $i++){
		if($array_gene_positions->[$i] == 1){
			$count_uncovered++;
		}
	}
	my $percentage_coverage = (($array_length-$count_uncovered)*100) / $array_length;
	return $percentage_coverage;
}

#########################################################################
sub change_gene_array {
	my $array_gene_positions = shift;
	my $gene                 = shift;
	my $transcript           = shift;

	my $gene_begin  = $gene->{'BEGIN'};
	my $transcript_begin = $transcript->{'BEGIN'};
	my $transcript_end = $transcript->{'END'};

	my $begin =  $transcript_begin - $gene_begin;
	my $end = $transcript_end - $gene_begin;
	for(my $i = $begin; $i<=$end; $i++){
		$array_gene_positions->[$i] = 0;
	}
}

###################################################################################
# validando parametros
##################################################################################
sub validate_parameters {
	my $allExists  = 1;
	my $fileExists = 1;

	if ( defined $help ) {
		print usage();
		exit 0;
	}

	unless ( defined $file_rnaseq ) {
		$allExists = 0;
	}

	unless ( defined $file_output ) {
		$allExists = 0;
	}

	unless ( defined $file_gff ) {
		$allExists = 0;
	}

	if ($allExists) {
		unless ( -e $file_rnaseq ) {
			print STDERR "$file_rnaseq doesn't exists.\n";
			$fileExists = 0;
		}
		unless ( -e $file_gff ) {
			print STDERR "$file_gff doesn't exists.\n";
			$fileExists = 0;
		}
	}

	unless ($allExists) {
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
  perl $0 -s rnaseq_file -g gff_file -o output
FOO
	return $usage;
}
