open COLLAPSED,"Arquivo.txt";

while(<COLLAPSED>)
{
    if($_=~/ko\:(K.+)\s\w+\;\s*(.+\s*\[*.*\]*)\s*/)
    {
    #print "$1,$2";
    $KHash{$1}=$2;
    }
}

open FH,"Lista.Ko";

while(<FH>)
{
    #print $_;    
    if($_=~/(.+)\s+(.+)/ && exists($KHash{$2}))
    {
    print "$1\t$2\t$KHash{$2}";
    }
}

close(COLLAPSED);
close(FH);
