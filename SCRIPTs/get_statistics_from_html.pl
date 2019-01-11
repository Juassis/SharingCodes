#!/usr/bin/perl

# Copyright Center for Excelence in Bioinformatics < www dot cebio dot org>
#
# You may distribute this module under the same terms as perl itself


use strict;
use Getopt::Long;

my $help;
my $file_html;
my $file_output;


my %hash_count = ();
my @array_ids=();

GetOptions(
	'i=s'    => \$file_html,
	'o=s'    => \$file_output,
	'help|?' => \$help
);

validate_parameters();

#$file_html = "KEGG_BRITE_Enzymes_Bg.html";
#$file_html = "KEGG_BRITE_Peptidases_Bg.html";
#$file_html = "KEGG_BRITE_Transcription_factors_Bg.html";
#$file_html = "KEGG_BRITE_Protein_kinases_Bg.html";


# http://www.genome.jp/kegg-bin/get_htext?htext=q01000.keg&filedir=/tools/kaas/files/log/result/1375981368&hier=4

main();


sub main{
	print "Reading file $file_html ...\n";
	read_file_html();

	print "Exporting results...\n";
	export_results();

	print "File generated: $file_output...\n";
	print "\nProcess finished!\n";
}

sub export_results{

	open(OUT, ">$file_output");

	for(@array_ids){
		my $array = $hash_count{$_};
		print OUT "$_\t".scalar @$array."\n";
	}
	close(OUT);
}

sub read_file_html{

	open(IN, $file_html);

	my $classification;
	while(<IN>){
		chomp;
		if(/<table id="grid"/){
			my $line = <IN>;
			chomp($line);
			while(not $line =~ /<\/table/){
				my $content="";
				if($line =~ /<tr/){
					while(not $line =~ /<\/tr/){
						$content .= $line;
						$line = <IN>;
						chomp($line);
					}
					if($content =~ /<\/a><b>(\d+\.\s*)*(.+)<\/b><\/pre>/){
						$classification=$2;
						push(@array_ids, $classification);

						my @array = ();
						$hash_count{$classification} = \@array;

					}else{

						my @array_split = split(/<img/, $content);

						for(@array_split){
							if(/whiteSP\.png"\s*>\s*(.+)\s*;\s*<a/g){
								if(exists $hash_count{$classification}){
									my $array = $hash_count{$classification};

									unless($1 ~~ @$array){
										push(@$array, $1);
#										$hash_count{$classification} = $array;
									}

								} else{
									my @array = ($1);
									$hash_count{$classification} = \@array;
								}
							}
						}

					}
				}

				$line = <IN>;
				chomp($line);
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

	unless ( defined $file_html ) {
		$allExists = 0;
	}

	unless ( defined $file_output ) {
		$allExists = 0;
	}

	if ($allExists) {
		unless ( -e $file_html ) {
			print STDERR "$file_html doesn't exists.\n";
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
	perl $0 -i html_file -o output_name
FOO
	return $usage;

}
