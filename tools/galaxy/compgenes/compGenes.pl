#!/usr/bin/perl -w

#Margaret Antonio 17.01.06 without essentiality. For original, unmodified aggregate.pl output (where unique insertions are not included). 

use Getopt::Long;
use Statistics::Distributions;
use strict;
use autodie;
no warnings;

#no warnings;

#ASSIGN INPUTS TO VARIABLES USING FLAGS
our ($input1, $input2, $h, $out, $sortkey, $round, $l1, $l2, $outfile);

GetOptions(
'input1:s' => \$input1,
'input2:s' => \$input2,
'h' => \$h,
'o:s' =>\$outfile,
'r:i'=> \$round,
'l1:s'=>\$l1,
'l2:s'=>\$l2
);

sub print_usage() {
    print "\n####################################################################\n";
    
    print "compGenes: compare genes for an organism under different conditions\n\n";
    print "DESCRIPTION: Takes two aggregate.pl outputs and compares them by\n";
    print "calculating the difference in mean fitness for each gene.\n";
    print "Example: compare organism in presence of control vs antibiotic.\n";
    print "Note: For different strains/genomes, use compStrains.pl\n";
    
    print "\nUSAGE:";
    print "perl compGenes.pl -d inputs/ \n\n";
    
    print "REQUIRED:\n";
    print " -d\tDirectory containing all input files (files from\n";
    print "   \taggregate script)\n";
    print "   \tOR\n";
    print "   \tIn the command line (without a flag), input the name(s) of\n";
    print "   \ttwo files containing aggregate gene fitness values. \n\n";
    
    print "OPTIONAL:\n";
    print " -h\tPrints usage and exits program\n";
    print " -o\tOutput file for comparison data. Default: label1label2.csv\n";
    print " -r\tRound final output numbers to this number of decimals\n";
    print " -l1\tLabels for input files. Default: filenames\n";
    print "   \tTwo strings, comma separated (i.e. -l expt1,expt2).\n";
    print "   \tOrder should match file order.\n";
    
    print " \n~~~~Always check that file paths are correctly specified~~~~\n";
    print "\n##################################################################\n";

}

# Check if help needed or if improper inputs

if ($h){
    print_usage();
    exit;
}

if (! $outfile){
    $outfile="geneComp.csv"
}
$round = '%.4f';

my @files=($input1,$input2);

#GET LABELS: USE (-l) OR USE FILNEAMES AS LABELS FOR COLUMNS IN OUTPUT FILE

my @labels=($l1,$l2);

#CHECK IF REQ. VARIABLES WERE DEFINED USING FLAGS. IF NOT THEN USE DEFAULT VALUES

#$out="comp.".$labels[0].$labels[1].".csv";

#READ INPUT FILES AND PUT GENE INFO FROM EACH LINE INTO HASH %all WHERE KEY=GENE_ID AND
#VALUE=GENE INFO (MEAN FITNESS, # INSERTIONS, ETC.)

my %all;
my @header;

#OPEN TWO COMPARISON FILES, ONE PER LOOP
for (my $i=0; $i<2; $i++){
    print "File #",$i+1,"\t";
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or die "Could not open '$file'\n";
    
    #EXTRACT THE HEADER (COLUMN NAMES) OF THE FILE AND KEEP FOR OUTPUT HEADER
    #APPEND FILE NAME OR GIVEN LABEL (-l) TO COLUMN NAME SO ITS DIFFERENT FROM OTHER INPUT FILE
    
    my $head=<DATA>;
    my @temp = split("\n",$head);
    $head = $temp[0];
    
    my @cols=split(',',$head);
    @cols = @cols[0,1,2,3,4,5];
    for (my $j=0;$j<scalar @cols;$j++){
        $cols[$j]=$cols[$j].'-'.$labels[$i];
    }
    push (@header,@cols);
    while (my $entry = <DATA>) {
        chomp $entry;
        my @line=split(",",$entry);
        if (!$line[5]){
            $line[5]="NA";
        }

        @line=@line[0,1,2,3,4,5];
        my $gene=$line[0];
        chomp($gene);
        
        #PUT GENE AND INFO INTO THE HASH FOR EXISTING KEY (1ST FILE) OR CREATE NEW KEY (2ND FILE)
        
        if (!exists $all{$gene}){
            my @info;
            push (@info,@line);
            $all{$gene}=\@info;
        }
        else{
            my @info=@{$all{$gene}};
            push (@info,@line);
            my $diff=sprintf("$round",($info[1]-$info[7]));
            my $total1=$info[2];
            my $total2=$info[8];
            my $sd1=$info[3];
            my $se1=$info[4];
            my $sd2=$info[9];
            my $se2=$info[10];
            my $df=$total1+$total2-2;
            my $tdist;
            my $pval;
            
            # CHECK TO MAKE SURE ALL VARIABLES IN TDIST,PVAL CALCULATIONS ARE NUMBERS AND NO
            # ZEROS (0) IN THE DENOMINATOR
            
            if ($se1 eq "X" or $se2 eq "X" or $sd1 eq "X" or $sd2 eq "X" or $total1==0 or $total2==0 or $sd1==0 or $sd2==0 or $df<=0){
                ($tdist,$pval)=("NA","NA");
            }
            else{
                $tdist=sqrt((($diff)/(sqrt((($sd1**2)/$total1)+(($sd2**2)/$total2))))**2);
                
                $pval=Statistics::Distributions::tprob($df,$tdist);
            }
            push (@info,$diff,$df,$tdist,$pval);
            $all{$gene}=\@info;
        }
    }
    close DATA;
}

#READ ALL HASH CONTENTS INTO 2D ARRAY FOR EASY SORTING AND PRINTING TO OUT FILE
my @unsorted;

foreach my $entry (keys %all) {
    my @info=@{$all{$entry}};
    my @temp;
    push (@temp,@info);
    push (@unsorted,\@temp);
}

#SORT GENES BY PVALUE OR FITNESS DEPENDING ON WHICH FLAG WAS USED IN COMMANDLINE

$sortkey=15; #default: sort by difference of means

my @sorted = sort { $a->[$sortkey] <=> $b->[$sortkey] } @unsorted;

#ADD NEW FIELDS TO HEADER (COLUMN NAMES)
my $field="Mean".$labels[0].'.'.$labels[1];
push (@header,$field,"DOF","TDIST","PVALUE");

#PRINT TO OUT FILE
open OUT, '>',"$outfile";
print OUT (join(',',@header),"\n");
foreach (@sorted){
    my @woo=@{$_};
    print OUT join(',',@woo),"\n";
    }
close OUT;


