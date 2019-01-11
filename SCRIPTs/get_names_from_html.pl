#!/usr/bin/perl

# Copyright Center for Excelence in Bioinformatics < www dot cebio dot org>
#
# You may distribute this module under the same terms as perl itself
#
# Pay attention: This script will only work if you have 3 levels


use strict;
use Getopt::Long;
use Data::Dumper;
my $help;
my $file_input;
my $file_output;



my %hash_contar=();
my @array_niveisA=(); # para armazenar os niveis de A
my %hash_niveisB = ();


GetOptions(
	'i=s'    => \$file_input,
	'o=s'    => \$file_output,
	'help|?' => \$help
);

main();


sub main{

	validate_parameters();

	print "Reading file $file_input ...\n\n";
	read_file();

	export_file();

	print "\nProcess finished!\n";
}

sub export_file{

	open(OUT, ">$file_output");

	for(@array_niveisA){
		my $nivelA = $_;
		my $array_niveisB = $hash_niveisB{$nivelA};
		for(@$array_niveisB){
			my $nivelB = $_;

			while(my ($nivelD, $contagem) = each %{$hash_contar{$nivelA}{$nivelB}}){
				print OUT "$nivelA\t$nivelB\t$nivelD\t$contagem\n";
			}

		}
	}
# as linhas comentadas abaixo seria caso nós levassemos em consideração o C na hierarquia
#	for(@array_niveisA){
#		my $nivelA = $_;
#		my $array_niveisB = $hash_niveisB{$nivelA};
#		for(@$array_niveisB){
#			my $nivelB = $_;
#
#			while(my ($nivelC, $value) = each %{$hash_contar{$nivelA}{$nivelB}}){
#				while(my ($nivelD, $contagem) = each(%$value)){
#					print OUT "$nivelA\t$nivelB\t$nivelC\t$nivelD\t$contagem\n";
#				}
#			}
#
#		}
#	}

	close(OUT);

}

sub read_file{

	open(IN, $file_input);

	my $nivelA_id;
	my $nivelB_id;
	my $nivelC_id;
	my $nivelD_id;

	while(<IN>){
		chomp;
		unless(/^[#+!]/){
			if(/^A<b>(.+)<\/b>$/){
				$nivelA_id = $1;
				push(@array_niveisA, $nivelA_id);
			}elsif(/B\s{2}<b>(.+)<\/b>$/){
				$nivelB_id=$1;
				if(exists $hash_niveisB{$nivelA_id}){
					my $array = $hash_niveisB{$nivelA_id};
					push(@$array, $nivelB_id);
				}else{
					my @array = ($nivelB_id);
					$hash_niveisB{$nivelA_id} = \@array;
				}
			}elsif(/C\s{4}(.+)$/){
				$nivelC_id=$1;
			}elsif(/D\s{6}(.+)$/){
				my $value = $1;
				my @array = split(/;/, $value);
				$nivelD_id=$array[0];
#				if(exists $hash_contar{$nivelA_id}{$nivelB_id}{$nivelC_id}{$nivelD_id}){
#					$hash_contar{$nivelA_id}{$nivelB_id}{$nivelC_id}{$nivelD_id}++;
#				}else{
#					$hash_contar{$nivelA_id}{$nivelB_id}{$nivelC_id}{$nivelD_id}=1;
#				}
				if(exists $hash_contar{$nivelA_id}{$nivelB_id}{$nivelD_id}){
					$hash_contar{$nivelA_id}{$nivelB_id}{$nivelD_id}++;
				}else{
					$hash_contar{$nivelA_id}{$nivelB_id}{$nivelD_id}=1;
				}
			}
		}

	}

	close(IN);
}

sub validate_parameters {
	my $allExists  = 1;
	my $fileExists = 1;

	if ( defined $help ) {
		print usage();
		exit 0;
	}

	unless ( defined $file_input ) {
		$allExists = 0;
	}

	unless ( defined $file_output ) {
		$allExists = 0;
	}

	if ($allExists) {
		unless ( -f $file_input ) {
			print STDERR "$file_input doesn't exists.\n";
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
	perl $0 -i inputfile -o outputfile
FOO
	return $usage;

}
