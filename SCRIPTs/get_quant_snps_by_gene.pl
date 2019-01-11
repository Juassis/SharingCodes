# Copyright Center for Excelence in Bioinformatics < www dot cebio dot org>
#
# You may distribute this module under the same terms as perl itself

use strict;
use Getopt::Long;

my $file_samtools; # input file
my $file_gff; # input file

my $file_output; #ID	QUANT

my $count_lines_samtools = 0;
my $count_lines_gff = 0;
my %hash_count = ();

my %hash_scaffolds_samtools = ();
my @ids_samtools=();
my %hash_scaffolds_gff = ();

my %hash_results=();
my @ids_results=();


my $help;

#receiving parameters
GetOptions('s=s' => \$file_samtools,
			'g=s'	=> \$file_gff,
			'o=s' => \$file_output,
			'help|?' => \$help
);
$| = 1;

main();

sub main{
	#validating parameters
	validate_parameters();

	print "Preprocessing files...\n";
	pre_processing(); # counting quantity of lines in input file

	print "Reading Samtools file...\n";
	read_samtools(); #reading samtools input file and generating output

	print "Reading GFF file...\n";
	read_gff();

	print "Crossing information...\n";
	count_snps();

	print "Exporting results...\n";
	export_results();

	print "\nFile generated: [$file_output]...\n";

	print "\nProcess finished!\n";

}

# Essa função é responsável por contas a quantidade de linhas do arquivo para mostrar
# A porcentagem de progresso de leitura dos arquivos
sub pre_processing{
	my $exclude = 0;

	$count_lines_samtools = `wc -l $file_samtools`;
	$exclude = `grep -cP '^#' $file_samtools`;
	$count_lines_samtools -= $exclude;

	$count_lines_gff = `wc -l $file_gff`;
	$exclude = `grep -cP '^#' $file_gff`;
	$count_lines_gff -= $exclude;
}

# Função responsável por exportar os resultados
sub export_results{
	open(OUT, ">".$file_output);

	# o array @ids_results possui os ids na ordem de leitura do arquivo
	for(@ids_results){
		my $id_scaffold = $_;

		# pra cada scaffold vamos imprimir os genes e a quantidade de snps associados
		my $hash_genes = $hash_results{$id_scaffold};

		my @ids_genes = sort {$hash_genes->{$b} <=> $hash_genes->{$a}} keys %{$hash_genes}; # vamos ordenar pela quantidade de snps por scaffolds

		for(@ids_genes){
			print OUT "$id_scaffold\t$_\t".$hash_genes->{$_}."\n";
		}
	}

	close(OUT);
}


# Função responsável por contar a quantidade de snps por gene
sub count_snps{

	# a variavel @ids_samtools armazena os identificadores dos scaffolds provenientes do arquivo do samtools
	for(@ids_samtools){
		my $id_scaffold = $_;
		my @positions = sort {$a <=> $b} @{$hash_scaffolds_samtools{$id_scaffold}};

		if(exists $hash_scaffolds_gff{$id_scaffold}){

			my @genes = sort {$a->{'BEGIN'}<=>$b->{'BEGIN'}} @{$hash_scaffolds_gff{$id_scaffold}};

			my $gene = shift(@genes);
			my $gene_id = $gene->{'ID'};
			my $begin = $gene->{'BEGIN'};
			my $end = $gene->{'END'};

			for(@positions){ # pra cada posição identificada como snp proveniente do samtools
				my $position = $_;

				if($position < $begin){ # se a posição do snp é menor do que a posição inicial do primeiro gene então ignoramos esse snp e vamos para o próximo
					next;
				}else{
					while($position > $end){ # se a posição do snp é maior do que a posição final do gene então ignoramos o gene e vamos para o próximo
						if(scalar @genes > 0){ # se a quantidade de genes para o scaffold atual for maior do que 0
							$gene = shift(@genes);
							$gene_id = $gene->{'ID'};
							$begin = $gene->{'BEGIN'};
							$end = $gene->{'END'};
						}else{ # se não existir mais genes no scaffold atual saimos do laço
							last;
						}
					}

					if($position >= $begin and $position <= $end){ # se o snp encontra-se entre a posição inicial e final do gene então ele nos interessa

						if(exists $hash_results{$id_scaffold}){ # se o scaffold já possui outro gene com snp

							my $hash_genes = $hash_results{$id_scaffold}; # recuperamos o hash de genes do scaffold atual

							if(exists $hash_genes->{$gene_id}){ # se o gene já possui algum snp apenas incrementamos
								$hash_genes->{$gene_id}++;
							}else{ # se o gene não possui nenhum snp atribuimos o valor de 1
								$hash_genes->{$gene_id} = 1;
							}

						}else{# se o scaffold não possui gene com snp
							my %hash_genes = ($gene_id => 1); # o gene passa a ter 1 snp

							$hash_results{$id_scaffold} = \%hash_genes; # armezenamos esse gene no hash de resultados

							push(@ids_results, $id_scaffold); # esse array só serve pra armazenar a ordem dos ids visto que hash não guarda a ordem das chaves

						} # end if exists $hash_results{$id_scaffold}
					}#end if $position >= $begin and $position <= $end
				} # end if $position < $begin
			} # end for(@positions)
		} # end if exists $hash_scaffolds_gff{$id_scaffold}
	} # end for(@ids_samtools)
}


sub read_samtools{
	# opening the input file
	open(IN, $file_samtools);

	my $count = 0;
	my $wrap_line = 0;
	my $line;

	while(<IN>){
		chomp;
		unless(/^(\s|\t)*$/){
			unless(/^#/){
				#	0	1		2		3		4		5		6		7		8		9								10								11								12								13								14								15								16								17								18								19								20
				#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sort_accepted_hits_AGcat.bam    sort_accepted_hits_BUCcat.bam   sort_accepted_hits_CNScat.bam   sort_accepted_hits_DGHPcat.bam  sort_accepted_hits_FOOTcat.bam  sort_accepted_hits_HAPOcat.bam  sort_accepted_hits_KIDcat.bam   sort_accepted_hits_MANcat.bam   sort_accepted_hits_OVOcat.bam   sort_accepted_hits_SALcat.bam   sort_accepted_hits_STOcat.bam   sort_accepted_hits_TRGcat.bam

				my @columns = split /\t/ ;

				my $id_scaffold = $columns[0];
				my $posicao = $columns[1];

				if(exists $hash_scaffolds_samtools{$id_scaffold}){
					my $array = $hash_scaffolds_samtools{$id_scaffold};

					push(@$array, $posicao);

					$hash_scaffolds_samtools{$id_scaffold} = $array;

				}else{
					my @array = ($posicao);

					$hash_scaffolds_samtools{$id_scaffold} = \@array;
					push(@ids_samtools, $id_scaffold);
				}


				# the following lines is just to calculate and printing the percentage of progress
				my $percent_lines = int(($count++ * 100)/$count_lines_samtools);
				unless(exists $hash_count{$percent_lines}){
					print sprintf(" | %2d%%", $percent_lines);
					$wrap_line++;
					$hash_count{$percent_lines} = $percent_lines;
					if($percent_lines > 0 and $wrap_line % 10 == 0){
						print "\n";
					}
				}
			}
		}else{
			$count_lines_samtools--;
		}
	}
	print "\n";

	close(IN);

	%hash_count = ();

}

sub read_gff{

	open(IN, $file_gff);

	my $count = 0;
	my $wrap_line = 0;
	my $line;

	while(<IN>){
		chomp;
		unless(/^(\s|\t)*$/){
			unless(/^#/){

				my @columns = split /\t/ ;

				my $id_scaffold = $columns[0];

				my $feature = $columns[2];

				if(lc($feature) eq "gene"){
					my $gene = $columns[8];

					$gene =~ s/ID=(.+);.*$/$1/g;

					my $begin = $columns[3];
					my $end = $columns[4];

					if($begin > $end){
						my $aux = $end;
						$end = $begin;
						$begin = $aux;
					}

					my %hash_gene = ('ID'=>$gene, 'BEGIN'=>$begin, 'END'=>$end);

					if(exists $hash_scaffolds_gff{$id_scaffold}){

						my $array = $hash_scaffolds_gff{$id_scaffold};

						push(@$array, \%hash_gene);

						$hash_scaffolds_gff{$id_scaffold} = $array;

					}else{
						my @array=(\%hash_gene);

						$hash_scaffolds_gff{$id_scaffold} = \@array;
					}


				}

				# the following lines is just to calculate and printing the percentage of progress
				my $percent_lines = int(($count++ * 100)/$count_lines_gff);
				unless(exists $hash_count{$percent_lines}){
					print sprintf(" | %2d%%", $percent_lines);
					$wrap_line++;
					$hash_count{$percent_lines} = $percent_lines;
					if($percent_lines > 0 and $wrap_line % 10 == 0){
						print "\n";
					}
				}
			}
		}else{
			$count_lines_gff--;
		}
	}
	print "\n";

	close(IN);
	%hash_count = ();

}


sub validate_parameters{
	my $allExists = 1;
	my $fileExists = 1;

	if ( defined $help ) {
		print usage();
		exit 0;
	}

	unless ( defined $file_samtools ) {
		$allExists = 0;
	}

	unless ( defined $file_output ) {
		$allExists = 0;
	}

	unless ( defined $file_gff ) {
		$allExists = 0;
	}



	if ($allExists) {
		unless ( -e $file_samtools ) {
			print STDERR "$file_samtools doesn't exists.\n";
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


sub usage{
	my $usage = <<FOO;
Usage:
	perl $0 -s samtools_file -g gff_file -o output
FOO
	return $usage;

}
