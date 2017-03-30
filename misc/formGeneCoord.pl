#!/usr/bin/perl -w

#perl ../Blueberries/formatAlign.pl --coord1 tigr4_genes.txt --coord2 19F_genes.txt --m ../3-strainAlignment/TIGR4vs.Taiwan19FpcKE.txt 

use strict;
use Getopt::Long;

our ($coord1,$coord2,$matched,$align);

GetOptions(
'coord1:s' => \$coord1,
'coord2:s' => \$coord2,
'm:s' =>\$matched,
'a:s' =>\$align
);

#genes of both aligned in small file
open (DIC,'<', $matched);
my %genePairs;
while(my $entry=<DIC>){
    chomp $entry;
    my @line=split("\t",$entry);
    if (scalar(@line)==2){
        $genePairs{$line[0]}=$line[1];
    }
}
close DIC;

#for (keys %genePairs){
#  print $_, "\t",$genePairs{$_},"\t";
#}

#First strain file
open (C1,'<', $coord1);
my %genes1;
my $dummy=<C1>;
while(my $entry=<C1>){
    chomp $entry;
    my @line=split("\t",$entry);
    my @coords=($line[1],$line[2]);
    $genes1{$line[5]}=[@coords];
    
}
close C1;


#Second strain file
open (C2,'<', $coord2);
my %genes2;
my $dummy2=<C2>;
while(my $entry=<C2>){
    chomp $entry;
    my @line=split("\t",$entry);
    my @coords=($line[1],$line[2]);
    $genes2{$line[5]}=[@coords];
}
close C2;


#Now add the coordinates to the alignment file
open (AL,'<',$align);
my %out;

while(my $entry=<AL>){
    chomp $entry;
    my @newline;
    my @line=split(",",$entry);
    
    if (((scalar @line)==4) and (defined $genes1{$line[0]}) and (defined $genes2{$line[2]})){
        push (@newline,$line[0]);
        my @sten1=@{$genes1{$line[0]}};
        foreach (@sten1){
            push(@newline,$_); #start coordinate
        }
        push (@newline,sprintf("%0.3f",$line[1])); #fitness
        
        #for corresponding strain2 gene
        push (@newline,$line[2]); #id of gene for strain2
        my @sten2=@{$genes2{$line[2]}};
        foreach (@sten2){
            push(@newline,$_); #start coordinate
        }
        push (@newline,sprintf("%0.3f",$line[3])); #fitness
        $out{$line[0]}=[@newline];
        
    }
}

close AL;

open (CSV,'>', "align_withCoords.txt");
foreach (sort keys %out){
    my @s=@{$out{$_}};
    print CSV join( ',', @s );
    print CSV "\n";
}

close CSV;







