#!/usr/bin/perl

# Copyright Center for Excelence in Bioinformatics < www dot cebio dot org>
#
# You may distribute this module under the same terms as perl itself

use strict;
use Getopt::Long;

my $help;
my $file_input1;
my $file_input2;
my $file_output;

my $extension="";

GetOptions(
	'i1=s'    => \$file_input1,
	'i2=s'    => \$file_input2,
	'help|?' => \$help
);

validate_parameters();

main();

my %hash_file_1=();
my %hash_file_2=();

my %hash_result_intersection=();
my %hash_result_exclusives1=();
my %hash_result_exclusives2=();

my $count_intersection = 0;
my $count_only1 = 0;
my $count_only2 = 0;

# -i1 posicao_basesc_gatk -i2 posicao_basesc_mPileup

sub main{
	print "Reading file $file_input1 and $file_input2 ...\n";
	read_file();

	find_intersection();


	print "Exporting files...\n";
	export_intersection();

	export_exclusives_file1();

	export_exclusives_file2();

	print "Process finished!\n";
}

# salva o arquivo com informações de interseção entre os arquivos
sub export_intersection{
	my $output = get_base_name().".intersection";
	if($extension ne ""){
		$output .= ".$extension";
	}
	open(OUT, ">$output");

	my @keys = sort {$a cmp $b} keys %hash_result_intersection;

	for(@keys){
		my $key = $_;
		my $value = $hash_result_intersection{$_};
		my @value_sorted = sort {$a <=> $b} @$value;

		for(@value_sorted){
			print OUT "$key\t$_\n";
		}
	}

	close(OUT);

	print "# Intersection: $count_intersection\n";
}

# salva o arquivo com informações exclusivas do arquivo 1
sub export_exclusives_file1{
	my $output = get_base_name().".only1";
	my $quant_col1 = 0;
	if($extension ne ""){
		$output .= ".$extension";
	}
	open(OUT, ">$output");

	my @keys = sort {$a cmp $b} keys %hash_result_exclusives1;

	for(@keys){
		my $key = $_;
		my $value = $hash_result_exclusives1{$_};
		my @value_sorted = sort {$a <=> $b} @$value;

		for(@value_sorted){
			print OUT "$key\t$_\n";
		}
	}


	close(OUT);

	print "# Exclusives file 1: $count_only1\n";

}

# salva o arquivo com informações exclusivas do arquivo 2
sub export_exclusives_file2{
	my $output = get_base_name().".only2";
	my $quant_col2 = 0;
	if($extension ne ""){
		$output .= ".$extension";
	}
	open(OUT, ">$output");

	my @keys = sort {$a cmp $b} keys %hash_result_exclusives2;

	for(@keys){
		my $key = $_;
		my $value = $hash_result_exclusives2{$_};
		my @value_sorted = sort {$a <=> $b} @$value;

		for(@value_sorted){
			print OUT "$key\t$_\n";
		}
	}

	close(OUT);

	print "# Exclusives file 2: $count_only2\n";
}

# função responsável por ler os arquivos passados como parâmetros
sub read_file{
	open(IN, $file_input1);
	while(<IN>){
		chomp;
		my @line = split(/\t/, $_);

		if(exists $hash_file_1{$line[0]}){
			my $ref = $hash_file_1{$line[0]};
			push(@$ref, $line[1]);


		}else{
			my @array = ($line[1]);
			$hash_file_1{$line[0]} = \@array;
		}

	}
	close(IN);

	open(IN, $file_input2);


	while(<IN>){
		chomp;
		my @line = split(/\t/, $_);
		if(exists $hash_file_2{$line[0]}){
			my $ref = $hash_file_2{$line[0]};
			push(@$ref, $line[1]);


		}else{
			my @array = ($line[1]);
			$hash_file_2{$line[0]} = \@array;
		}
	}

	close(IN);
}

# função responsável por filtrar as interseções e exclusividades entre os arquivos
sub find_intersection{

	while(my ($key, $value) = each(%hash_file_1)){

		if(exists $hash_file_2{$key}){
			my $array1 = $value;
			my $array2 = $hash_file_2{$key};

			my @intersection = get_intersection($array1, $array2);

			$hash_result_intersection{$key} = \@intersection;

			my $size_a1 = scalar @$array1;
			my $size_a2 = scalar @$array2;
			my $size_intersection = scalar @intersection;

			$count_intersection += $size_intersection;

			if($size_intersection == $size_a1 and $size_a1 < $size_a2){ # array1 é subconjunto de array2
				# vou ter somente exclusivos do arquivo 2
				my @result = get_exclusive($array1, $array2);
				$hash_result_exclusives2{$key} = \@result;

				$count_only2 += scalar @result;

			}elsif($size_intersection == $size_a2 and $size_a1 > $size_a2){# array2 é subconjunto de array1
				# vou ter somente exclusivos do arquivo 1
				my @result = get_exclusive($array2, $array1);
				$hash_result_exclusives1{$key} = \@result;

				$count_only1 += scalar @result;
			}elsif($size_intersection != $size_a1 and $size_intersection != $size_a2){ # tem exclusivos em array1 e em array2

				my @results_exclusive_1 = get_exclusive($array2, $array1);
				my @results_exclusive_2 = get_exclusive($array1, $array2);

				$hash_result_exclusives1{$key} = \@results_exclusive_1;
				$hash_result_exclusives2{$key} = \@results_exclusive_2;

				$count_only1 += scalar @results_exclusive_1;
				$count_only2 += scalar @results_exclusive_2;
			}else{ # são iguais

			}

			delete $hash_file_1{$key};
			delete $hash_file_2{$key};

		}else{ # não existe chave no arquivo 2, logo ele é exclusivo do arquivo 1
			$hash_result_exclusives1{$key} = $value;
			$count_only1 += scalar @$value;

			delete $hash_file_1{$key};
		}
	}

	# por fim o que restou no arquivo 2 são exclusivos do arquivo 2, então adicionaremos
	# ao hash de exclusivos
	if(scalar keys %hash_file_2 > 0){
		%hash_result_exclusives2 = (%hash_result_exclusives2, %hash_file_2);

		while(my ($key, $value) = each(%hash_file_2)){
			$count_only2 += scalar @$value;
		}
	}
}

# retorna os dados exclusivos do segundo parâmetro em relação com o primeiro
sub get_exclusive{
	my $array1 = shift or die "Parameter not found";
	my $array2 = shift or die "Parameter not found";

	my %lookup;
	my @result=();
	@lookup{@$array1} = ();
	for(@$array2) {
	  push(@result, $_) unless exists $lookup{$_};
	}
	return @result;
}

# retorna array com a interseção de dois arrays passados como parâmetros
sub get_intersection{
	my $array1=shift or die "Parameter not found\n";
	my $array2=shift or die "Parameter not found\n";

	my @intersection =
    grep { defined }
        @{ { map { lc ,=> $_ } @$array1 } }
           { map { lc } @$array2 };

	return @intersection;
}

# retorna nome base para o nome do arquivo de saida
sub get_base_name {
	my @name = split( /\./, $file_input1 );

	my $base_name;
	if ( scalar @name >= 1 ) {
		$base_name = $name[0];
	}

	@name = split( /\./, $file_input2 );
	if ( scalar @name >= 1 ) {
		$base_name .= "_".$name[0];
	}


	return $base_name;
}

# função responsáel por validar os parâmetros de linha de comando
sub validate_parameters {
	my $allExists  = 1;
	my $fileExists = 1;

	if ( defined $help ) {
		print usage();
		exit 0;
	}

	unless ( defined $file_input1 ) {
		$allExists = 0;
	}
	unless ( defined $file_input2 ) {
		$allExists = 0;
	}


	if ($allExists) {
		unless ( -e $file_input1 ) {
			print STDERR "$file_input1 doesn't exists.\n";
			$fileExists = 0;
		}
		unless ( -e $file_input2 ) {
			print STDERR "$file_input2 doesn't exists.\n";
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
	perl $0 -i1 file1 -i2 file2
FOO
	return $usage;

}
