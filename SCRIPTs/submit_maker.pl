#!/usr/bin/perl
 
use strict;
use warnings;
use Parallel::ForkManager;
 
## Define numero de processadores
my $MAX_PROCESSES = 78;

## Define o diret—rio de trabalho  *****NÌO USE O "/" NO FINAL*****
## Diret—rio de trabalho tem que conter os arquivos de controle do maker j‡ preenchidos mas com o linha "genome=" em branco
## Exemplo do maker_opts.ctl
#genome= #genome sequence (fasta file or fasta embeded in GFF3 file)
# perl submit_maker.pl lista
my $wd = "/sto2data-5/ju/anotaca_onca/data";

#Path para o arquivo do genoma a ser anotado, formatado para o BLAST
my $gf = "/sto2data-5/ju/anotaca_onca/data/Jaguar_assembly_scaffolds.fasta";

my $pm = new Parallel::ForkManager($MAX_PROCESSES);

while (<>) {
	
	my $pid = $pm->start and next;
	
	chdir($wd);
	
	my $input = $_;
	chomp $input;
	
	system("mkdir $wd/$input");
	system("cp maker* $wd/$input");
	chdir($wd."/".$input);
	system("/data/thiago/instaladores/blast-2.2.17/bin/fastacmd -s $input -d $gf | sed 's/lcl|//g' > $input.fasta");
	my $if = $wd."/".$input."/".$input.".fasta";
	$if =~ s/\//\\\//g;
	system("sed 's/^genome=/genome=$if/' maker_opts.ctl > maker_opts.ctl.new; mv maker_opts.ctl.new maker_opts.ctl");
	system("maker -fix_nucleotides");
	system("fasta_merge -d $input.maker.output/$input\_master_datastore_index.log -o $input");
	system("gff3_merge -d $input.maker.output/$input\_master_datastore_index.log -o $input.gff");	
	$pm->finish;
}

$pm->wait_all_children;
