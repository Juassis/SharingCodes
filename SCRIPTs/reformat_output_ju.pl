# Copyright Center for Excelence in Bioinformatics < www dot cebio dot org>
#
# You may distribute this module under the same terms as perl itself

use strict;
use Getopt::Long;

my $file_input; # input file
my $file_output; # output file

my $help;

#receiving parameters
GetOptions('i=s' => \$file_input,
			'o=s'	=> \$file_output,
			'help|?' => \$help
);
$| = 1;

main();

sub main{
	#validating parameters
	validate_parameters();

	read_file();

	print "\nFile generated: [$file_output]...\n";

	print "\nProcess finished!\n";

}

sub read_file{

	open(IN, $file_input);
	open(OUT, ">".$file_output);

	while(<IN>){
		chomp;
		my $first_line = $_;
		$first_line =~ s/^>(.+)$/$1/g;
		my $second_line = <IN>;
		chomp($second_line);

		print OUT "$first_line\t$second_line\n";

	}
	close(OUT);
	close(IN);


}



sub validate_parameters{
	my $allExists = 1;
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


sub usage{
	my $usage = <<FOO;
Usage:
	perl $0 -i input -o output
FOO
	return $usage;

}
