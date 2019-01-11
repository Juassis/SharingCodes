
# Copyright Center for Excelence in Bioinformatics < www dot cebio dot org>
#
# You may distribute this module under the same terms as perl itself

use strict;
use Getopt::Long;

my $file_clustal; # input file

my $file_output;
my $count_lines = 0; 
my %hash_count = ();
my $total_diff=0;

my $help;

#receiving parameters
GetOptions('i=s' => \$file_clustal,
			'o=s' => \$file_output,
			'help|?' => \$help
);

main();

sub main{
	#validating parameters
	validate_parameters();
	
	
	print "Reading input [$file_clustal]...\n";
	pre_processing(); # counting quantity of lines in input file
	
	read_file();
	
	print "\nFile generated: [$file_output]...\n";
	print "Total differences: $total_diff\n";
	
	print "\nProcess finished!\n";
	
}

#######################################
sub pre_processing{
	
	$count_lines = `wc -l $file_clustal`;
	$count_lines -= 1; # removing head line
}

sub read_file{
	# opening the input file and the output file
	open(IN, $file_clustal);
	
	open(OUT, ">$file_output");
	
	my $count = 0;
	my $line;
	
	my $head = <IN>;
	my $offset = 0;
	while(<IN>){
		chomp;
		unless(/^(\s|\t)*$/){
			
			my $first_line = $_; # getting the first line in alignment, represents the sequence query
	
			$first_line =~ /\s([ACTG-])/;
			
			my $pos = $+[0]; # position of first nucleotide in sequence
			
			my $size = length($first_line); # size of all line
			
			my $second_line = <IN>; # getting the second line in alignment, represents the sequence original 
			my $third_line = <IN>; # getting the third line in alignment (* if nucleotide aligned or \s if not)
			chomp($third_line);
			
			my $third_sequence = substr($third_line, $pos-1, $size); # getting only the alignment information without space in beggining
			
			my @positions=();
		    while ($third_sequence =~ /\s/g) {
		    	my $pos =  $offset + $+[0];
		    	push(@positions, $pos);
		    }
			
			for(@positions){
				print OUT "$_\n";
			}
			
			$total_diff += (scalar @positions);
			
			my $off_calc = $size-$pos+1;
			$offset += $off_calc;
			
			
			# the following lines is just to calculate and printing the percentage of progress
			$count+=3;
			my $percent_lines = int(($count * 100)/$count_lines);		
			unless(exists $hash_count{$percent_lines}){
				print sprintf(" | %-3s%%", $percent_lines);
				$hash_count{$percent_lines} = $percent_lines;
				if($percent_lines > 0 and $count % 10 == 0){
					print "\n";
				}
			}
		}else{
			$count_lines--;
		}
	}
	print "\n";
	
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
	
	unless ( defined $file_clustal ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output ) {
		$allExists = 0;
	}
	
	if ($allExists) {
		unless ( -e $file_clustal ) {
			print STDERR "$file_clustal doesn't exists.\n";
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
	perl $0 -i clustal_file -o output_file
FOO
	return $usage; 
	
}