#!/usr/bin/perl
 
use strict;
use warnings;
use Parallel::ForkManager;
 
## Define numero de processadores
my $MAX_PROCESSES = 78;



my $pm = new Parallel::ForkManager($MAX_PROCESSES);

while (<>) {
        
        my $pid = $pm->start and next;
        
        chdir($wd);
        
        my $input = $_;
        chomp $input;
        

        system("fastacmd -s $input -d  /sto2data-5/ju/anotaca_onca/Interpro/teste/arquivo_1_.all.maker.proteins.fasta| sed 's/lcl|//g' > $input.fasta");


        system("/usr/local/bioinformatics/interproscan/interproscan-5.7-48.0/interproscan.sh $input.fasta -goterms -iprlookup -pa -f xml -f tsv -b $input");
	system("rm $input.fasta");
	

        $pm->finish;
}

$pm->wait_all_children;

