# GFF e Orthofinder
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

# Some constants to easily change the columns index if necessary
use constant{
	PRINT_DEV => 0,
	PRINT_LOG => 1, # to show our development log: 1=show; 2=no show
	PRINT_SCORE => 0, # to show our development log: 1=show; 2=no show
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
	MAX_ERRORS => 2, # the number maximum of errors that we accept without breaking the block
	MAX_WINDOW => 5, # the number of genes in the sliding window
	MIN_SCAFFOLD_LENGTH => 10 # a ju vai estudar esse tamanho mínimo aceitável
};

#########################################################################
## Global variables to use in this script
#########################################################################

my @array_input_gff; # each item in this array is one GFF file
my $file_output; # the output file with the results
my $help;

my $file_orthofinder;    # input file from the orthofinder
my %hash_orth_group_by_gene = (); # chave: gene_id; valor: grupo de ortologo
my %hash_genes_by_orth_group = (); # chave: orthologous group; valor: array de genes

my %hash_orth_group_by_gene_uniq = (); # key: gene_id; value: ortholog group for single genes
my %hash_orth_group_by_gene_mult = (); # key: gene_id; value: ortholog group for multigenes family

my @array_scaffolds_by_species = (); # array with hashes with key: scaffold; value: array of genes
my @array_all_gff = (); # array with all gff file to report

my %hash_score_genes = ();

my $pattern = "\^.*ID=(transcript:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id
my $pattern_parent_name = "\^.*Parent=([^;]+);?.*\$"; # pattern to extract the gene id

###########################################################################
#receiving parameters
###########################################################################
GetOptions(
	's=s' => \$file_orthofinder,
  'i=s'    => \@array_input_gff,
	'o=s'    => \$file_output,
);

$| = 1;

main();
##################################################################################
# The script is divided into subroutines. The main() subroutine call the
# another subroutines. Tha main subroutine is just to organize our script
##################################################################################
sub main {

	# print "Reading orthofinder file...\n";
	read_orthofinder();    #reading orthofinder input file and mRNA rating output
	# o arquio do orthofinder é lido e o hash %hash_orth_group_by_gene é populado
	# sendo a chave do hash: id do gene; e o valor: o grupo de ortologo

	# # print "Reading GFF files...\n";
  read_gff_file(); # Calling the subroutine to read the GFF file @array_scaffolds_by_species

	# for(@array_scaffolds_by_species){
	# 	print keys %{$_};
	# 	print("\n");
	# }

	# # # o arquivo gff é lido e o @array_scaffolds_by_species é populado
	# # cada elemento do array é um hash onde a chave é o cromossomo e o valor é um arrays
	# # de genes com posição de início e fim
  sort_cds(); # Calling the subroutine to sort the CDS in each gene
	#
	compare_orthologous_by_synteny();
	#
	# if(PRINT_SCORE){
	# 	print(Dumper(%hash_score_genes));
	# }

}

sub compare_orthologous_by_synteny {
	# na linha abaixo vamos pegar o hash referente ao gff da espécie referência (coluna 1 do orthofinder)
	my $hash_gff1 = $array_scaffolds_by_species[0];

	# aqui pegamos todos os scaffolds da referência
	my @array_scaffolds_gff1 = keys %$hash_gff1;
	my %hash_synteny_error = ();
	for(@array_scaffolds_gff1){
		my $scaffold_reference = $_;
		# print($scaffold_reference, "\n");
		# next;
		my $array_genes_by_sccafold = $hash_gff1->{$_};
		if(scalar @$array_genes_by_sccafold > MIN_SCAFFOLD_LENGTH){
			my $idx_uniq_reference = 0;
			my $flag_no_og = 0;
			my $nth_gene = $array_genes_by_sccafold->[$idx_uniq_reference];

			# print STDERR $scaffold_reference, " - ";
			# print STDERR scalar @$array_genes_by_sccafold, " - ";

			while(not exists $hash_orth_group_by_gene_uniq{$nth_gene->{"ID"}}){
				# we need to think what to do with the multigenic families
				if($idx_uniq_reference <= (scalar @$array_genes_by_sccafold-2) ){
					$nth_gene = $array_genes_by_sccafold->[++$idx_uniq_reference];
				}else{
						print("No orthologous group for all genes in $scaffold_reference\n");
						$flag_no_og=1;
						last;
				}
			}
			if($flag_no_og){
				next;
			}
			# print($scaffold_reference, "\n");
			# print("ele saiu do while\n");
			# print STDERR $i, "\n";
			my $orth_group = $hash_orth_group_by_gene_uniq{$nth_gene->{"ID"}};
			my $array_genes = $hash_genes_by_orth_group{$orth_group};

			my $number_species = scalar @$array_genes;
			# print($scaffold_reference, "\n");
			# next;
			for(my $n = 1; $n<$number_species; $n++){
				my $second_gene = $array_genes->[$n];
				# gene: PCOAH_000435300.1
				# chromosome: CP016250.1
				my $hash_gff_n = $array_scaffolds_by_species[$n];
				my @array_scaffolds_gff_n = keys %$hash_gff_n;
				my $idx_homologous = 0;
				my $scaffold_idx = 0;

				# if($scaffold_reference eq "Schisto_mansoni.Chr_6"){
				# 	print($orth_group, " - ");
				# 	print(@$array_genes, "\n");
				#
				# 	my $scaffold = "";
				# 	for(my $sc_idx = 0; $sc_idx < scalar @array_scaffolds_gff_n; $sc_idx++){
				# 		$scaffold = $array_scaffolds_gff_n[$sc_idx];
				# 		print($scaffold, "\n");
				# 		my $array_genes_by_sccafold_n = $hash_gff_n->{$scaffold};
				#
				# 		for(my $idx = 0; $idx<scalar @$array_genes_by_sccafold_n; $idx++){
				# 			my $gene = $array_genes_by_sccafold_n->[$idx];
				# 			# print Data::Dumper::qquote($second_gene);
				# 			# print Data::Dumper::qquote($gene->{"ID"});
				#
				# 			if($gene->{"ID"} eq $second_gene){
				# 				print("entrou satanás\n");
				# 				$idx_homologous = $idx;
				# 				$scaffold_idx = $sc_idx;
				# 				$sc_idx = scalar @array_scaffolds_gff_n;
				# 				last;
				# 			}
				# 		}
				# 		# g has the index of the homologous gen
				# 	} # end @array_scaffolds_gff_n
				# 	# if($scaffold_reference eq "Schisto_mansoni.Chr_6"){
				# 	# 	print($scaffold, "\n");
				# 	# }
				# 	# if($scaffold_reference eq "Schisto_mansoni.Chr_6"){
				# 	# 	print($scaffold, "\n");
				# 	# 	print($scaffold_idx, "\n");
				# 	# 	print scalar @$array_genes_by_sccafold_n;
				# 	# 	print("\n");
				# 	#
				# 	# }
				# 	print($scaffold, "\n");
				# 	print("sai demonio\n");
				# 	return(0);
				# }

				my $scaffold = "";
				for(my $sc_idx = 0; $sc_idx < scalar @array_scaffolds_gff_n; $sc_idx++){
					$scaffold = $array_scaffolds_gff_n[$sc_idx];
					my $array_genes_by_sccafold_n = $hash_gff_n->{$scaffold};

					for(my $idx = 0; $idx<scalar @$array_genes_by_sccafold_n; $idx++){
						my $gene = $array_genes_by_sccafold_n->[$idx];
						# print Data::Dumper::qquote($second_gene);
						# print Data::Dumper::qquote($gene->{"ID"});
						if($gene->{"ID"} eq $second_gene){
							$idx_homologous = $idx;
							$scaffold_idx = $sc_idx;
							$sc_idx = scalar @array_scaffolds_gff_n;
							last;
						}
					}
					# g has the index of the homologous gen
				} # end @array_scaffolds_gff_n

				# g has the index of the ;homologous gene in the species n
				# i has the index of the gene in the species reference
				# $scaffold_idx has the index of the scaffold
				my $array_genes_by_sccafold_n = $hash_gff_n->{$scaffold};

				#PCOAH_000435300.1

				# numero de genes do scaffold na referencia
				my $number_genes_reference = scalar @$array_genes_by_sccafold;

				# montamos uma pilha com todos os genes do scaffold a ser comparado
				my @array_queue_reference = @{$array_genes_by_sccafold}[($idx_uniq_reference)..($number_genes_reference-1)];

				# numero de genes do scaffold no homologo
				my $number_genes_orth = scalar @$array_genes_by_sccafold_n;

				# montamos uma pilha com os genes no scaffold do homologo
				my @array_queue_homologous = @{$array_genes_by_sccafold_n}[($idx_homologous)..($number_genes_orth-1)];

				my @array_queue_homologous_copy = @array_queue_homologous;
				my $minor_queue = ($number_genes_reference < $number_genes_orth) ? $number_genes_reference : $number_genes_orth;

				my @array_queue_errors_ref_nogh = (); # pilha de genes que não possuem grupo de ortologos

				my @array_queue_errors_ref_genes = (); # pilha de erros de genes que não quebram o bloco
				my @array_queue_errors_orth_genes = (); # pilha de erros de genes que não quebram o bloco

				# my @array_queue_errors_ref_category = ();
				my $flag_pop_reference = 1;
				my $flag_pop_homologous = 1;
				my $gene_of_reference;
				my $gene_homologous;

				# ideia inicial colocar num hash todos os genes com erros e IGUAIS
				# aqueles que não estiverem no hash, provavelmente sao genes que nao
				# quebram o bloco
				my %hash_report = ();

				while($minor_queue > 0){
					# print($minor_queue, "\n");
					if($flag_pop_reference){
						$gene_of_reference = shift(@array_queue_reference);
						# print("ref - ", $gene_of_reference->{"ID"}, "\n");
					}
					if($flag_pop_homologous){
						$gene_homologous = shift(@array_queue_homologous);
						# print("orth - ", $gene_homologous->{"ID"}, "\n");
					}

					my $flag_orth_reference_exists = exists $hash_orth_group_by_gene{$gene_of_reference->{"ID"}};
					my $flag_orth_homologous_exists = exists $hash_orth_group_by_gene{$gene_homologous->{"ID"}};

					if($flag_orth_reference_exists && $flag_orth_homologous_exists){
						my $orth_group_reference = $hash_orth_group_by_gene{$gene_of_reference->{"ID"}};
						my $orth_group_homologous = $hash_orth_group_by_gene{$gene_homologous->{"ID"}};
						if(PRINT_DEV){
							print(sprintf("%s\t%s\t%s\t%s\n", $gene_of_reference->{"ID"}, $orth_group_reference,$orth_group_homologous, $gene_homologous->{"ID"}));
						}
						if($orth_group_reference eq $orth_group_homologous){
							$hash_score_genes{$gene_of_reference->{"ID"}}{"SCORE"} = 5;
							$hash_score_genes{$gene_of_reference->{"ID"}}{"CATEGORY"} = "EQUAL";
							$hash_score_genes{$gene_of_reference->{"ID"}}{"GENE_HOMO"} = $gene_homologous->{"ID"}; #ELENAO

							if(scalar @array_queue_errors_ref_genes > 0){
								for(@array_queue_errors_ref_genes){
									if(exists $hash_synteny_error{$_->{"ID"}}){
										$hash_score_genes{$_->{"ID"}}{"SCORE"} = 2;
										$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "ERROR_NO_BREAK_FROM_SYNTENY";
										$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "ERROR_NO_BREAK_FROM_SYNTENY";
									}else{
										$hash_score_genes{$_->{"ID"}}{"SCORE"} = 4;
										$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "ERROR_NO_BREAK";
										$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "ERROR_NO_BREAK";
									}
								}
							}

							if(scalar @array_queue_errors_ref_nogh > 0){
								for(@array_queue_errors_ref_nogh){
									$hash_score_genes{$_->{"ID"}}{"SCORE"} = 3;
									$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "NO_GH";
									$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "NO_GH";
								}
							}



							if(PRINT_LOG){
								print($gene_of_reference->{"ID"}, " - ", $gene_homologous->{"ID"});
								print(" - SAO IGUAISSSSSS AEWWWWWW\n");
								if(scalar @array_queue_errors_ref_genes > 0){
									for(@array_queue_errors_ref_genes){
										print($_->{"ID"}, " - ERRO SEM QUEBRA DE BLOCO\n");
									}
									# for(@array_queue_errors_orth_genes){
									# 	print($_->{"ID"}, " - ERRO SEM QUEBRA DE BLOCO\n");
									# }
								}
								if(scalar @array_queue_errors_ref_nogh > 0){
									for(@array_queue_errors_ref_nogh){
										print($_->{"ID"}, " - NO GH\n");
									}
								}
							}#PKNOH_9999999997-t35_1

							@array_queue_errors_ref_genes = ();
							@array_queue_errors_orth_genes = ();
							@array_queue_errors_ref_nogh = ();

							$flag_pop_reference = 1;
							$flag_pop_homologous = 1;
						}else{ # end if($orth_group_reference eq $orth_group_homologous){
							my @array_queue_one;

							if(scalar @array_queue_reference > MAX_WINDOW){
								# @array_queue_one = @array_queue_reference[0..(MAX_WINDOW-1)];
								for(my $idx = 0; ($idx < scalar @array_queue_reference) && (scalar @array_queue_one < MAX_WINDOW);$idx++){
									my $gene = $array_queue_reference[$idx];
									if(exists $hash_orth_group_by_gene{$gene->{"ID"}}){
										push(@array_queue_one, $gene);
									}
								}

							}else{
								# @array_queue_one = @array_queue_reference;
								for(my $idx = 0; $idx < scalar @array_queue_reference; $idx++){
									my $gene = $array_queue_reference[$idx];
									if(exists $hash_orth_group_by_gene{$gene->{"ID"}}){
										push(@array_queue_one, $gene);
									}
								}
							}

							my @array_queue_second;
							if(scalar @array_queue_homologous > MAX_WINDOW){
								# @array_queue_second = @array_queue_homologous[0..(MAX_WINDOW-1)];
								for(my $idx = 0; ($idx < scalar @array_queue_homologous) && (scalar @array_queue_second < MAX_WINDOW);$idx++){
									my $new_gene = $array_queue_homologous[$idx];

									if(exists $hash_orth_group_by_gene{$new_gene->{"ID"}}){
										push(@array_queue_second, $new_gene);
									}
								}
							}else{
								# @array_queue_second = @array_queue_homologous;
								for(my $idx = 0; $idx < scalar @array_queue_homologous; $idx++){
									my $gene = $array_queue_homologous[$idx];
									if(exists $hash_orth_group_by_gene{$gene->{"ID"}}){
										push(@array_queue_second, $gene);
									}
								}
							}
							if(scalar @array_queue_one > 0 && scalar @array_queue_second > 0){
								my @array_comparison = ($gene_of_reference);

								my $ref_to_orth = compare_from_one_to_sec(\@array_comparison, \@array_queue_second);
								@array_comparison = ($gene_homologous);

								my $orth_to_ref = compare_from_one_to_sec(\@array_comparison, \@array_queue_one);

								if($ref_to_orth & !$orth_to_ref){
									$flag_pop_reference = 0;
									$flag_pop_homologous = 1;
								}elsif(!$ref_to_orth & $orth_to_ref){
									push(@array_queue_errors_ref_genes, $gene_of_reference);
									push(@array_queue_errors_orth_genes, $gene_homologous);
									$flag_pop_reference = 1;
									$flag_pop_homologous = 0;
								}elsif(!$ref_to_orth & !$orth_to_ref){
									push(@array_queue_errors_ref_genes, $gene_of_reference);
									push(@array_queue_errors_orth_genes, $gene_homologous);

									$flag_pop_reference = 1;
									$flag_pop_homologous = 1;
								} # end if($ref_to_orth & !$orth_to_ref){
								else{
									@array_queue_one = ($gene_of_reference, @array_queue_one);
									@array_queue_second = ($gene_homologous, @array_queue_second);

									my $array_compare_two_lists = compare_two_lists(\@array_queue_one, \@array_queue_second);

									if(scalar @$array_compare_two_lists == 0){
										# se os arrays comparados acima retorna 0 significa que existe uma inversão
										# em que os 5 genes da referência possuem grupo de ortologos iguais aos 5 grupos de ortologos dos genes dos homologos

										$hash_score_genes{$gene_of_reference->{"ID"}}{"SCORE"} = 4;
										$hash_score_genes{$gene_of_reference->{"ID"}}{"CATEGORY"} = "SINTENY_CONSERVED_BLOCK";
										$hash_score_genes{$gene_of_reference->{"ID"}}{"GENE_HOMO"} = $gene_homologous->{"ID"};

										if(PRINT_LOG){
											print("POSSUI BLOCO DE SINTENIA IGUAL\n");
											print($gene_of_reference->{"ID"}, " - ", $gene_homologous->{"ID"}, "\n");
										}

										for(my $idx = 0; $idx < (scalar @array_queue_one - 1 ); $idx++){
											$gene_of_reference = shift(@array_queue_reference);
											$gene_homologous = shift(@array_queue_homologous);

											$hash_score_genes{$gene_of_reference->{"ID"}}{"SCORE"} = 4;
											$hash_score_genes{$gene_of_reference->{"ID"}}{"CATEGORY"} = "SYNTENY_CONSERVED_BLOCK";
											$hash_score_genes{$gene_of_reference->{"ID"}}{"GENE_HOMO"} = $gene_homologous->{"ID"};

											if(PRINT_LOG){
												print($gene_of_reference->{"ID"}, " - ", $gene_homologous->{"ID"}, "\n");
											}
										}

										if(PRINT_LOG){
											if(scalar @array_queue_errors_ref_genes > 0){
												for(@array_queue_errors_ref_genes){
													print($_->{"ID"}, " - ERROS SEM QUEBRA DE BLOCO\n");
												}
												# for(@array_queue_errors_orth_genes){
												# 	print($_->{"ID"}, " - ERROS SEM QUEBRA DE BLOCO\n");
												# }
											}
											if(scalar @array_queue_errors_ref_nogh > 0){
												for(@array_queue_errors_ref_nogh){
													print($_->{"ID"}, " - NO GH\n");
												}
											}
										}

										if(scalar @array_queue_errors_ref_genes > 0){
											for(@array_queue_errors_ref_genes){
												if(exists $hash_synteny_error{$_->{"ID"}}){
													$hash_score_genes{$_->{"ID"}}{"SCORE"} = 2;
													$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "ERROR_NO_BREAK_FROM_SYNTENY";
													$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "ERROR_NO_BREAK_FROM_SYNTENY";
												}else{
													$hash_score_genes{$_->{"ID"}}{"SCORE"} = 4;
													$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "ERROR_NO_BREAK";
													$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "ERROR_NO_BREAK";
												}
											}
										}

										if(scalar @array_queue_errors_ref_nogh > 0){
											for(@array_queue_errors_ref_nogh){
												$hash_score_genes{$_->{"ID"}}{"SCORE"} = 3;
												$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "NO_GH";
												$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "NO_GH";
											}
										}

										@array_queue_errors_ref_genes = ();
										@array_queue_errors_orth_genes = ();
										@array_queue_errors_ref_nogh = ();
										# para verificar quais genes que sao ignorados pq só acontecem uma vez
										# a gente pode verificar nesse momento antes de apagar a pilha de ERROS
										# se ela possui algum elemento
										$flag_pop_reference = 1;
										$flag_pop_homologous = 1;
									}else{ # end if(scalar @$array_compare_two_lists == 0){
										if(PRINT_LOG){
											print($gene_of_reference->{"ID"}, " - ", $gene_homologous->{"ID"});
											print(" - BLOCO SEM SINTENIA\n");
											# talvez imprimir todo o bloco com problema de sintenia em um arquivo separado
										}
										# we create a hash with all the problems from reference
										# this will be useful to score the problematic genes inside a block without synteny
										for(@$array_compare_two_lists){
											$hash_synteny_error{$_} = 1;
										}
										$flag_pop_reference = 0;
										$flag_pop_homologous = 1;
										# aqui não existe sintenia perfeita entre os cinco genes comparados de cada organismo
										# precisa pensar em uma estrategia para tratar estes casos
									}


								}
								my $flag_full_queue_different = scalar @array_queue_errors_ref_genes > MAX_ERRORS;
								if($flag_full_queue_different){
									for(@array_queue_errors_ref_genes){
										if(exists $hash_synteny_error{$_->{"ID"}}){
											$hash_score_genes{$_->{"ID"}}{"SCORE"} = 1;
											$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "ERROR_BREAK_FROM_SYNTENY";
											$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "ERROR_BREAK_FROM_SYNTENY";
										}else{
											$hash_score_genes{$_->{"ID"}}{"SCORE"} = 1;
											$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "ERROR_BREAK";
											$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "ERROR_BREAK";
										}
									}
									for(@array_queue_errors_ref_nogh){
										$hash_score_genes{$_->{"ID"}}{"SCORE"} = 3;
										$hash_score_genes{$_->{"ID"}}{"CATEGORY"} = "NO_GH";
										$hash_score_genes{$_->{"ID"}}{"GENE_HOMO"} = "NO_GH";
									}

									if(PRINT_LOG){
										for(@array_queue_errors_ref_genes){
											print($_->{"ID"}, " - QUEBRA DE BLOCO\n");
										}
										# for(@array_queue_errors_orth_genes){
										# 	print($_->{"ID"}, " - QUEBRA DE BLOCO\n");
										# }
										for(@array_queue_errors_ref_nogh){
											print($_->{"ID"}, " - NO GH EM QUEBRA DE BLOCO\n");
										}
										print("--------------------------------------------\n");
									}
									@array_queue_errors_ref_nogh = ();
									@array_queue_errors_ref_genes = ();
									@array_queue_errors_orth_genes = ();
								}
							}
						}
					} # end if($flag_orth_reference_exists && $flag_orth_homologous_exists)
					else{
						# print($gene_of_reference->{"ID"}, " - ", $gene_homologous->{"ID"});
						# print(" - nao tem grupo de ortologos\n");
						if(! $flag_orth_reference_exists & $flag_orth_homologous_exists){
							push(@array_queue_errors_ref_nogh, $gene_of_reference);
							$flag_pop_reference = 1;
							$flag_pop_homologous = 0;
						} elsif($flag_orth_reference_exists & !$flag_orth_homologous_exists){

							$flag_pop_reference = 0;
							$flag_pop_homologous = 1;
						}elsif(!$flag_orth_reference_exists & !$flag_orth_homologous_exists){
							push(@array_queue_errors_ref_nogh, $gene_of_reference);
							$flag_pop_reference = 1;
							$flag_pop_homologous = 1;
						}
					}

					$number_genes_reference = scalar @array_queue_reference;
					$number_genes_orth = scalar @array_queue_homologous;
					$minor_queue = ($number_genes_reference < $number_genes_orth) ? $number_genes_reference : $number_genes_orth;
				} # end while($minor_queue > 0){

			}# end for(my $n = 0; $n<$number_species; $n++){
		} # end if(scalar @$array_genes_by_sccafold > MIN_SCAFFOLD_LENGTH){
		else{
			print STDERR "Scaffold/chromosome [",$scaffold_reference,"] not used because number of genes [",(scalar @$array_genes_by_sccafold),"] is less than ",MIN_SCAFFOLD_LENGTH;
			print STDERR "\n";
		}

	} # end for @array_scaffolds_gff1

}# end function compare_orthologous_by_synteny

sub compare_two_lists{
	my $array_queue_one = shift;
	my $array_queue_second = shift;
	my %hash_gh1 = ();
	my %hash_gh2 = ();

	for(my $i = 0; $i<scalar @$array_queue_second; $i++){
		my $first_gene = $array_queue_second->[$i];
		if(exists $hash_orth_group_by_gene{$first_gene->{"ID"}}){
			$hash_gh2{$hash_orth_group_by_gene{$first_gene->{"ID"}}}= 0;
		}
	}

	my @uniques_from_ref = ();
	for(my $i = 0; $i<scalar @$array_queue_one; $i++){
		my $first_gene = $array_queue_one->[$i];
		if(exists $hash_orth_group_by_gene{$first_gene->{"ID"}}){
			my $orth_group = $hash_orth_group_by_gene{$first_gene->{"ID"}};
			if(not exists $hash_gh2{$orth_group}){
				push(@uniques_from_ref, $first_gene->{"ID"});
			}
		}

	}

	return \@uniques_from_ref;
}


# the idea of this function is compare two arrays of genes
# and verify if the first gene of the first queue has the ortolog group
# in the second queue
sub compare_from_one_to_sec{
	my $array_queue_one = shift;
	my $array_queue_second = shift;

	for(my $i = 0; $i<scalar @$array_queue_one; $i++){
			my $first_gene = $array_queue_one->[$i];
			if(exists $hash_orth_group_by_gene{$first_gene->{"ID"}}){
				my $orth_group_first = $hash_orth_group_by_gene{$first_gene->{"ID"}};

				for(my $j = 0; $j < scalar @$array_queue_second; $j++){
					my $second_gene = $array_queue_second->[$j];
					if(exists $hash_orth_group_by_gene{$second_gene->{"ID"}}){
						my $orth_group_homologous = $hash_orth_group_by_gene{$second_gene->{"ID"}};
						# print($orth_group_first, " - ", $orth_group_homologous, "\n");
						if($orth_group_first eq $orth_group_homologous){
							return 1;
						}
					}else{
						next;
					}
				}
			}else{
				next;
			}
	}
	return 0;
}

sub verify_if_error_queue_is_full{
	my $array_queue = shift;
	if(scalar @$array_queue > MAX_ERRORS){
		return 1;
	}
	return 0;

}

############################################################################
# Fazer a leitura do arquivo Orthofinder
############################################################################
sub read_orthofinder {

	# opening the input file
	open( IN, $file_orthofinder );

	my $count     = 0;
	my $wrap_line = 0;
	my $line = <IN>;
	chomp $line;
	my @columns = split /\t/, $line;
	my $nspecies = scalar  @columns -1;

	while (<IN>) {
		chomp;
				#   0   1	2	3	4	5	6	7	...
				#Orthogroup  mRNAs
				# OG0004424       PCOAH_0006314001        PKNOH_S140255000-t35_1  PVP01_1436100.1

		my @columns = split /\t/;
		my $grupo_ortologos = $columns[0];
		my $flag = 0;
		if( (scalar @columns -1 ) != $nspecies ){
			$flag =1;
		}else{
			for (my $i = 1; $i <= $nspecies; $i++) {
				my $geneN = $columns [$i];

				my @geneCount = split /,/, $geneN;
				if (scalar @geneCount  >1) {
					$flag =1;
				}
			}
		}

		my @genes_per_species = @columns[1..$nspecies];

		if ($flag ==0) {

			$hash_genes_by_orth_group{$grupo_ortologos} = \@genes_per_species;
			for (@genes_per_species){
					$hash_orth_group_by_gene{$_} = $grupo_ortologos;

					$hash_orth_group_by_gene_uniq{$_} = $grupo_ortologos;
			}
		}
		else {
			for(@genes_per_species){
				if(defined($_) && $_ ne ""){
					my @array_genes = split(", ", $_);
					for(@array_genes){
						$hash_orth_group_by_gene{$_} = $grupo_ortologos;
						# grupo de ortologos com genes multiplos estao sendo inseridos
						# no hash de grupo de ortologos de genes unicos
						$hash_orth_group_by_gene_mult{$_} = $grupo_ortologos;
					}
				}
			}
		}

	}
	close(IN);

#print Dumper (%hash_genes_multiplos);

			}

#########################################################################
# Extract the gene id based on the global variable pattern
#########################################################################
sub get_gene_id{
	my $gene_text = shift;
	my $rex = qr/$pattern/;
	if($gene_text =~ $rex){
		$gene_text =~ s/$rex/$2/g;
	}
	return $gene_text;
}

sub get_parent_gene{
	my $attribute_text = shift;
	my $regex = qr/$pattern_parent_name/;
	if($attribute_text =~ $regex){
		$attribute_text =~ s/$regex/$1/g;
	}
	return $attribute_text;
}

#########################################################################
# This subroutine  sorts the features from all genes to keep the features in
# crescent order.
#########################################################################
sub sort_cds{
	for(@array_scaffolds_by_species){
		my $hash_mrna_by_scaffolds = $_;
		my @scaffolds = keys %$hash_mrna_by_scaffolds;
		for(@scaffolds){
			my $scaffold = $_;
			my $array_mrnas = $hash_mrna_by_scaffolds->{$scaffold};
			my @array_ordered = sort {$a->{GENE_START} <=> $b->{GENE_START} || $a->{GENE_END} <=> $b->{GENE_END} }  @$array_mrnas;
			$hash_mrna_by_scaffolds->{$scaffold} = \@array_ordered;
		} # end for(@scaffolds)
	}
}# end sort_cds

#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
#########################################################################
sub read_gff_file{

	for(my $idx = 0; $idx < scalar @array_input_gff; $idx++){
		my $file_input_gff = $array_input_gff[$idx];
		# print "Reading $file_input_gff...\n";
		open(IN, $file_input_gff);
		my %hash_mrna_by_scaffolds = ();
		my %hash_parents_by_scaffolds = ();
		while(<IN>){
			chomp;
			if(/^>/){
				last;
			}
			unless(/^#/){  # To avoid lines which start with '#'
				unless(/^\s*$/){  # This is to avoid lines with no relevant content
					my $line = $_;
					my @cols = split /\t/;

					#saving lines from the gff to iterate later
					if($idx == 0){
						push(@array_all_gff, \@cols);
					}

					my $gene_name = $cols[-1]; # The gene name is at last column
					my $scaffold = $cols[0]; # The scaffold name is the first column
					my $strand = $cols [STRAND];

					my $gene_id = get_gene_id($gene_name);
					my $parent_id = get_parent_gene($gene_name);

					my $feature = $cols[FEATURE];

					# just in case: when the start of the gene is greater than the end position we invert
					if($cols [POS_START] > $cols [POS_END]){
						my $aux = $cols [POS_START];
						$cols [POS_START] = $cols [POS_END];
						$cols [POS_END] = $aux;
					}

					if(lc($feature) eq "mrna"){
						if(!exists $hash_parents_by_scaffolds{$scaffold}{$parent_id}){
							$hash_parents_by_scaffolds{$scaffold}{$parent_id} = 1;
							my %hash_mrna = ();
							$hash_mrna{ID} = $gene_id;
							$hash_mrna{STRAND} = $strand;
							$hash_mrna{GENE_START} = $cols [POS_START];
							$hash_mrna{GENE_END} = $cols [POS_END];

							if (not exists $hash_mrna_by_scaffolds {$scaffold}){
								my @array_mrnas = ();
								push @array_mrnas, \%hash_mrna;
								$hash_mrna_by_scaffolds {$scaffold} = \@array_mrnas;
							}else{
								my $array_mrnas = $hash_mrna_by_scaffolds {$scaffold};
								push @$array_mrnas, \%hash_mrna;
							}
						}# end if(!exists $hash_parents_by_scaffolds{$parent_id}){
					} # end if(lc($feature) eq "mrna"){
				}
			}
		}
		close(IN);
		push(@array_scaffolds_by_species, \%hash_mrna_by_scaffolds);
	}

}
