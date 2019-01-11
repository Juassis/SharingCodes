
#!/usr/bin/perl
#
###########################################################################################################
### You may distribute this module under the same terms as perl itself
### Laura Rabelo and Juliana Assis
### GGBC Fiocruz-MG - 09_07_2015
###########################################################################################################

#use strict;
use Getopt::Long;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use List::MoreUtils qw(uniq);
use Data::Dumper;

my $help;
my $input1;
my $input2;
my $output;

GetOptions(
	'1=s'    => \$input1,
	'2=s'    => \$input2,
        'o=s'    => \$output,
	'help|?' => \$help
);


main();


########################
## MAIN 
#########################


sub main{
	validate_parameters();
	print "Reading files $input1 and $input2 \n";
	my $hash_file_a = read_files($input1);
	my $hash_file_b = read_files($input2);
	my $out = compare($hash_file_a, $hash_file_b);
	my @output = @$out;
	open(OUT, ">".$output);
	while(<OUT>){
		print "@output\n";
	}
	close(OUT);
	print "\nFile generated: $output\n";
	print "Process finished!\n";
}


########################
### VALIDATION
##########################

#verificando se todos os parâmetros foram passados
sub validate_parameters {
	my $allExists  = 1;
	my $fileExists = 1;

	if ( defined $help ) {
		print usage();
		exit 0;
	}

	unless ( defined $input1 ) {
		$allExists = 0;
	}
	unless ( defined $input2 ) {
		$allExists = 0;
	}
	unless ( defined $output ) {
		$allExists = 0;
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

#mostrando modo de uso do programa
sub usage {
	my $usage = <<FOO;
Usage:
	perl $0 -1 input1 -2 input2 -o output
	input1		txt file
	input2		txt file
	output		Name of output txt
FOO
	return $usage;

}



########################
#### READ INPUT FILES
###########################


#lê os arquivos.txt
sub read_files{
	my $files=shift;
#	print "$files\n";
	my %hash_file;
	open(IN, "<".$files);
	while(<IN>){
#		print $_;
		chomp;
		my @line = split ("\t", $_);
		$hash_file{$line[0]} = $line[0];
	}	
	close(IN);
	return \%hash_file;
#	while(my ($key, $value) = each(%hash_file)){
#		print "$key\n";
#    	}
	
}



#

sub compare{
#comparando dois arquivos por vez 
	my $hash_file_a = shift;
	my $hash_file_b = shift;
#	print Dumper (\%$hash_file_a);
	my @out;
	while(my ($key_a, $value) = each(%$hash_file_a)){
#		print $value."\n";
		if(exists $hash_file_b->{$key_a}){
			my $tmp = $key_a."\t".$hash_file_b->{$key_a}."\t1\n";
			push(@out, $tmp);
			print $tmp;
		}else{
			my $tmp = $key_a."\t0\n";
			push(@out, $tmp);
			print $tmp;
		}
	}
	return \@out;
}
