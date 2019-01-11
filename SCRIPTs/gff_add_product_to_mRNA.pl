#!/usr/bin/env perl
#args:
# 1: GFF annotation file
# 2: Product list from BLAST results
#output to STDOUT
# GFF annotation file with product attribute added to mRNA and CDS features
use strict;
use warnings;
use Data::Dumper;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;

my ($gfffile,$productfile) = @ARGV;

my %products;

open(FH, "<$productfile") or die "$!\n";
my @lines = <FH>;
close(FH);
foreach my $line (@lines){
  chomp $line;
  my ($id,$prod)=split(/\t/, $line);
  $prod =~ s/^[^"]*"([^"]*)"\s.*$/$1/;
  $prod =~ s/Similar to \S+\s*//i;
  $prod =~ s/\s*\(strain H\)//i;
  $prod =~ s/\s*\(Plasmo[^)]+\)"$/"/i;
  my @str = split(//,$prod);
  my $parens = 0;
  for(my $i = $#str; $i >= 0; $i--){
    if($str[$i] eq '('){
      $parens--;
    }
    elsif($str[$i] eq ')'){
      $parens++;
    }
    elsif($parens){
    }
    else{
      $prod = join("", @str[0 .. $i]);
      $prod =~ s/\s*$//;
      $prod =~ s/,\s*/ - /g; ### WATCH OUT! if this becomes " - putative" possibly tbl2asn will move 'putative' to beginning and leave a training ' - ' which is not allowed
      last;
    }
  }
  $products{$id} = $prod;
#printf("%s\t%s\n", $id, $prod);
}
  
@lines = ();
  


open(FH, "<$gfffile") or die "$!\n";

@lines = <FH>;
close(FH);
chomp @lines;

my $gff = Bio::Tools::GFF->new(-gff_version => '3');

for(my $i=0; $i<=$#lines; $i++){
  last if ($lines[$i] =~ /^##FASTA/i);
  next if $lines[$i] =~ /^#/;
  my $f = Bio::SeqFeature::Generic->new();
  $gff->from_gff_string($f, $lines[$i]);
  my $pid;
  if($f->primary_tag eq 'mRNA'){
    ($pid) = $f->get_tag_values('ID');
  }
  elsif($f->primary_tag eq 'CDS'){
    ($pid) = $f->get_tag_values('Parent');
  }
  else {
   print "$lines[$i]\n";
   next;
  }
  if($pid && $products{$pid}){
		if($pid && $products{$pid}){
			$lines[$i] =~ s/;product=[^;]+//;
    	$lines[$i] .= ";product=$products{$pid}";
		}
  }
  elsif($pid){
    print STDERR ("NO PRODUCT at $i:$pid\n");
  }
  else{
    print STDERR ("NO ID at $i\n");
  }
  print "$lines[$i]\n";
}
    
  
