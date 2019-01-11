#!/usr/bin/env perl

# use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

# Some constants to easily change the columns index if necessary
my %hash_params;
my $gff_reference;
my $fasta_reference;
my @fasta_files;
my @gff_files;

my $project_dir;
my $pipeline_dir;
my $params_file;
my $GFF_DIR;
my $FASTA_DIR;
my $OUTPUT_DIR;

GetOptions(
	'd=s' => \$project_dir,
  'p=s'    => \$params_file
);

$| = 1;


main();

sub main{
  $pipeline_dir = $0;
  $pipeline_dir =~ s/Pipeline\.pl//g;
  $pipeline_dir .= "bin";

  if(not defined $project_dir){
      $project_dir = ".";
  }else{
      $project_dir =~ s/\/$//g;
  }
  if(not defined $params_file){
      $params_file = "params";
  }

  print("Reading parameters...\n");
  read_params ();

  $GFF_DIR = fix_directory($hash_params{"GFF_DIR"});
  $FASTA_DIR = fix_directory($hash_params{"FASTA_DIR"});
  $OUTPUT_DIR = fix_directory($hash_params{"OUTPUT_DIR"});
  #
  print("Reading input files...\n");
  read_directory ();

  print("Extracting protein sequences...\n");
  extract_protein ();
  print("Checking stop codons...\n");
  check_stop_codon();
  print("Checking overlap...\n");
  check_overlap();
  if(lc($hash_params{"RNASEQ"}) eq "yes"){
      print("Exporting RNASeq coverage...\n");
      export_rnaseq();
  }


}

sub fix_directory {
  my $dir = shift;
  $dir =~ s/\/$//g;
  return($project_dir."/".$dir);
}

sub export_rnaseq{
  unless(-e "$OUTPUT_DIR/REPORT"){
    `mkdir -p $OUTPUT_DIR/REPORT`;
  }
  my $command = sprintf ("$pipeline_dir/export_rnaseq.pl -i $GFF_DIR/%s -a mRNA > %s/REPORT/check_rnaseq", $gff_reference, $OUTPUT_DIR);
  my $output = `$command`;
  print($output, "\n");

}

sub check_overlap{
  unless(-e "$OUTPUT_DIR/REPORT"){
    `mkdir -p $OUTPUT_DIR/REPORT`;
  }
  my $command = sprintf ("$pipeline_dir/Overlap_gff.pl -i $GFF_DIR/%s -a %s > %s/REPORT/check_overlap", $gff_reference, $hash_params{"ACCEPTED_TRANSCRIPTS"}, $OUTPUT_DIR);
  my $output = `$command`;
  print($output, "\n");


}

sub check_stop_codon{
  unless(-e "$OUTPUT_DIR/REPORT"){
    `mkdir -p $OUTPUT_DIR/REPORT`;
  }
  my $command = sprintf ("$pipeline_dir/protein_check.pl $OUTPUT_DIR/PROTEINS/%s > %s/REPORT/check_protein", $fasta_reference, $OUTPUT_DIR);
  my $output = `$command`;
  print($output, "\n");

}

sub extract_protein{
  unless(-e "$OUTPUT_DIR/PROTEINS"){
    `mkdir -p $OUTPUT_DIR/PROTEINS`;
  }

  my $command = sprintf ("$pipeline_dir/gff2fasta.pl $FASTA_DIR/%s $GFF_DIR/%s %s/PROTEINS/", $fasta_reference, $gff_reference, $OUTPUT_DIR);
  print($command);

  my $output = `$command`;
  print($output, "\n");

  for(my $i = 0; $i<scalar @fasta_files; $i++){
    my $command = sprintf ("$pipeline_dir/gff2fasta.pl $FASTA_DIR/%s $GFF_DIR/%s %s/PROTEINS/", $fasta_files[$i], $gff_files[$i], $OUTPUT_DIR);
    my $output = `$command`;
    print($output, "\n");
  }

  #gff2perl Genome.fasta Annotation.gff OutputPrefix
}


####################
sub read_directory{
  opendir(my $dh, $GFF_DIR) || die "Can't opendir $GFF_DIR: $!";
      my @gff_directory = grep { /^[^.]/ && -f "$GFF_DIR/$_" } readdir($dh);
      closedir $dh;



  for (@gff_directory) {
    my @split_name =  split (/\./, $_);
    pop (@split_name);
    my $joined_name = join (".", @split_name);
    my $fasta_name = $joined_name.".fasta";
    # print ($fasta_name, "\n");
    if($joined_name eq $hash_params {"reference_genome"} ){
      $fasta_reference = $fasta_name;
      $gff_reference = $_;
    }
    else{
      push @fasta_files, $fasta_name;
      push @gff_files, $_;
    }

  }
}

###########################################
sub read_params{
  open (IN, $project_dir."/".$params_file);
  while (<IN>) {
    chomp;
    unless(/^#/ || /^\s*$/){
      my @split_name = split (/\s*=\s*/, $_);
      $hash_params {$split_name [0]} = $split_name [1];
    }
  }

  close (IN);
}
