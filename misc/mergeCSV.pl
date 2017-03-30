#!/usr/bin/perl -w

#Margaret Antonio
#16.01.12

use strict;
use warnings;
use autodie;
use File::Path;
use File::Basename;
use Getopt::Long;

#USAGE: perl ../../Blueberries/mergeCSV.pl -indir ../T4Dapto/ -o T4Dmerged.csv

#ASSIGN INPUTS TO VARIABLES

our ($infile,$indir,$outfile);
GetOptions(
'i:s' => \$infile,
'o:s' => \$outfile,
'indir:s'=>\$indir,
);

if (!$outfile){$outfile="merged.csv";}

my @files;
if ($indir){
    my $directory="$indir";
    opendir(DIR, $directory) or die "couldn't open $directory: $!\n";
    my @direct= readdir DIR;
    my $tail=".csv";
    foreach (@direct){
        if (index($_, $tail) != -1){
            $_=$indir.$_;
            push (@files,$_);
        }
    }
    closedir DIR;
}
else{
    @files=@ARGV;
}
my $num=(scalar @files);
my @unsorted;

for (my $i=0; $i<$num; $i++){   #Read files from ARGV
    print "File #",$i+1,"\t";
    
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    my $dummy=<DATA>;
    while (my $entry = <DATA>) {
        chomp $entry;
        my @line=split(",",$entry);
        push(@unsorted,\@line);   #keep track of actual insertion site position
        }
        close DATA;
        }



my @sorted = sort { $a->[0] <=> $b->[0] } @unsorted;


#PRINT TO FILE HERE

open (OUT, '>', $outfile);

print OUT (join(",",@$_),"\n") for @sorted;

close OUT;




#15.03.27
#my $files = $#ARGV; # get number of files - 1
#while (my $file = shift @ARGV) {
#   open my $fh, "<", $file;
#   <$fh> unless $files == @ARGV; # discard header unless first file
#   while (<$fh>){
#       push ; # output the rest
#}



#In terminal, use command $ perl ../script/mergeCSV.pl results/L1_2394eVI_PennG.csv results/L2_2394eVI_PennG.csv results/L3_2394eVI_PennG.csv results/L4_2394eVI_PennG.csv results/L5_2394eVI_PennG.csv results/L6_2394eVI_PennG.csv >compiled.csv

