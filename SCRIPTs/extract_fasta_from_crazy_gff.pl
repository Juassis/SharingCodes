use strict;
no warnings;
use Getopt::Long;

my $file_input;
my $file_gff_output;
my $file_fasta_output;
my $help;

GetOptions(
	'i=s'    => \$file_input,
	'g=s'    => \$file_gff_output,
	'f=s'    => \$file_fasta_output,
	'help|?' => \$help
);

main();

sub main {
	validate_parameters();
		
	print "Reading $file_input...\n";
	read_file();
	
	print "Output generated [$file_fasta_output]...\n";

	print "Process finished!\n";

}

sub read_file{
	open(IN, $file_input);
	open(OUT_GFF, ">$file_gff_output");
	open(OUT_FASTA, ">$file_fasta_output");
	my $header_gff = <IN>; # ##gff-version 3
	$header_gff .= <IN>; # ##source-version geneious 6.1.8
	
	my $line = <IN>;
	my %hash_fasta_id = ();
	my $fasta_id="";
	my $sequence = "";
	print OUT_GFF $header_gff;
	do{
		if($line =~ /^##Type DNA\s+(.+)$/){
			# não faz nada
		}elsif(not $line =~ /^##/){ 
			print OUT_GFF $line;
		}else{
			if($line =~ /##DNA (.+)/){
				if($fasta_id eq ""){
					$fasta_id = $1;
				}else{
					print OUT_FASTA ">$fasta_id\n";
					print OUT_FASTA "$sequence";
					$sequence = "";
					$fasta_id = $1;
				}
			}else{
				unless($line =~ /^##end-DNA/){
					$line =~ s/^##//g;
					$sequence .= $line;	
				}
			}
		}
		
	}while($line = <IN>);
	
	print OUT_FASTA ">$fasta_id\n";
	print OUT_FASTA "$sequence";
	
	close(OUT_GFF);
	close(OUT_FASTA);
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
	
	unless ( defined $file_gff_output ) {
		$allExists = 0;
	}

	unless ( defined $file_fasta_output ) {
		$allExists = 0;
	}

	if ($allExists) {
		unless ( -f $file_input ) {
			print STDERR "$file_input doesn't exists.\n";
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
	perl $0 -i input_crazy_gff -g output_gff -f output_fasta
FOO
	return $usage;

}
