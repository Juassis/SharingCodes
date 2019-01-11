#!/usr/bin/perl

# Copyright Center for Excelence in Bioinformatics < www dot cebio dot org>
#
# You may distribute this module under the same terms as perl itself

use strict;
use Getopt::Long;

my $help;
my $file_input;
my $file_output;

my $extension="";

GetOptions(
	'i=s'    => \$file_input,
	'help|?' => \$help
);

validate_parameters();

main();

my @array1=();
my @array2=();
my @intersection=();
my %hash_intersection=();

sub main{
	print "Reading file $file_input ...\n";
	read_file();

	find_intersection();

	my $quant_intersec = scalar @intersection;
	print "# Intersection: $quant_intersec\n";
	export_intersection();


	export_exclusives_column1();

	export_exclusives_column2();

	print "Process finished!\n";
}

sub export_intersection{
	my $output = get_base_name().".intersection";
	if($extension ne ""){
		$output .= ".$extension";
	}
	open(OUT, ">$output");
	for(@intersection){
		print OUT "$_\n";
	}
	close(OUT);
}

sub export_exclusives_column1{
	my $output = get_base_name().".only1";
	my $quant_col1 = 0;
	if($extension ne ""){
		$output .= ".$extension";
	}
	open(OUT, ">$output");
	for(@array1){
		unless(exists $hash_intersection{$_}){
			print OUT "$_\n";
			$quant_col1++;
		}
	}
	close(OUT);

	print "# Exclusives column 1: $quant_col1\n";

}

sub export_exclusives_column2{
	my $output = get_base_name().".only2";
	my $quant_col2 = 0;
	if($extension ne ""){
		$output .= ".$extension";
	}
	open(OUT, ">$output");
	for(@array2){
		unless(exists $hash_intersection{$_}){
			print OUT "$_\n";
			$quant_col2++;
		}
	}
	close(OUT);

	print "# Exclusives column 2: $quant_col2\n";
}

sub read_file{
	open(IN, $file_input);

	while(<IN>){
		chomp;
		my @line = split(/\t/, $_);

		if($line[0] ne ""){
			push(@array1, $line[0]);

		}
		if($line[1] ne ""){
			push(@array2, $line[1]);

		}
	}

	close(IN);
}

sub find_intersection{
	@intersection =
    grep { defined }
        @{ { map { lc ,=> $_ } @array1 } }
           { map { lc } @array2 };

   # colocaremos num hash para otimizar a busca
   for(@intersection){
   		$hash_intersection{$_} = 1;
   }
}

sub get_base_name {
	my @name = split( /\./, $file_input );
	my $base_name;
	if ( scalar @name > 1 ) {
		$extension = pop(@name);
		$base_name = join( ".", @name );
	}
	else {
		$base_name = $name[0];
	}

	return $base_name;
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


	if ($allExists) {
		unless ( -e $file_input ) {
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
	perl $0 -i file
FOO
	return $usage;

}
