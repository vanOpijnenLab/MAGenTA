#!/usr/bin/perl -w

#Margaret Antonio 16.08.30

use Data::Dumper;
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Statistics::Distributions;
no warnings;

#ASSIGN INPUTS TO VARIABLES USING FLAGS
our ($indir,$help,$out,$round,$l,$sortkey);
GetOptions(
'd:s' => \$indir,
'h' => \$help,
'o:s' =>\$out,
's:i' => \$sortkey,
'r:i'=> \$round,
'l:s'=> \$l,
);

sub print_usage() {
    print "\n####################################################################\n";
    print "\n";
    print "compWindows.pl: compare regions outputted by the slidingWindow tool\n\n";
    print "DESCRIPTION: Takes two slidingWindows.csv files and compares them by\n";
    print "calculating the difference in mean fitness, the pval for each gene.\n";
    print "Example: compare control vs antibiotic, where both from same strain (genome).\n";
    print "Note: For different strains/genomes, use compStrains for a gene comparison.pl\n";
    
    print "\nUSAGE:\n";
    print "perl compGenes.pl -d inputs/ -l cond1,cond2 -r 2 -o comp-cond1cond2_date.csv\n";
    
    print "\nREQUIRED:\n";
    print " -d\tDirectory containing two slidingWindow.csv files (output of\n";
    print "   \tslidingWindow tool)\n";
    print "   \tOR\n";
    print "   \tIn the command line (without a flag), input the name(s) of\n";
    print "   \ttwo files.\n";
    
    print "\nOPTIONAL:\n";
    print " -h\tPrints usage and exits program\n";
    print " -o\tOutput file for comparison data. Default: label1label2.csv\n";
    print " -r\tRound final output numbers to this number of decimals\n";
    print " -s\tColumn number to sort by. Default: difference of means\n";
    print " -l\tLabels for input files. Default: filenames\n";
    print "   \tTwo strings, comma separated (i.e. -l expt1,expt2).\n";
    print "   \tOrder should match file order.\n";
    
    print " \n~~~~Always check that file paths are correctly specified~~~~\n";
    print "\n##################################################################\n";
}
if ($help){
    print_usage();
    exit;
}

if (!$indir and (scalar @ARGV==0)){
	print "\nERROR: Please correctly specify input files or directory\n";
    print_usage();
	print "\n";
	exit;
}


#THE @files ARRAY WILL CONTAIN INPUT FILE NAMES, EXTRACTED FROM A DIRECTORY (-indir) OR ARGV
my @files;
if ($indir){
    my $directory="$indir";
    opendir(DIR, $directory) or (print "Couldn't open $directory: $!\n" and print_usage() and exit);
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

#GET LABELS: USE (-l) OR USE FILNEAMES AS LABELS FOR COLUMNS IN OUTPUT FILE

my @labels;
if ($l){
    @labels=split(',',$l);
}
else{
    foreach (@files){
        my $filename=basename($_);
        my @temp=split(/\./,$filename);
        my $colName=$temp[0];
        print $colName,"\t";
        push (@labels,$colName);
    }
}

#INDICES 0:start	1:end	2:mutants	3:insertions	4:TA_sites	5:ratio	6:p-value
#7:average	8:variance	9:stdev	10:stderr	11:genes	12:fit-mean

#CHECK IF REQ. VARIABLES WERE DEFINED USING FLAGS. IF NOT THEN USE DEFAULT VALUES

if (!$out) {$out="comp.".$labels[0].$labels[1].".csv"}
if (!$round){$round='%.4f'}

#READ INPUT FILES AND PUT GENE INFO FROM EACH LINE INTO HASH %all WHERE KEY=GENE_ID AND
#VALUE=GENE INFO (MEAN FITNESS, # INSERTIONS, ETC.)

my %all;
my @header;
for (my $i=0; $i<2; $i++){
    print "File #",$i+1,"\t";
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or (print "Could not open '$file'\n" and print_usage() and exit);
    
    #EXTRACT THE HEADER (COLUMN NAMES) OF THE FILE AND KEEP FOR OUTPUT HEADER
    #APPEND FILE NAME OR GIVEN LABEL (-l) TO COLUMN NAME SO ITS DIFFERENT FROM OTHER INPUT FILE
    
    my $head=<DATA>;
    my @cols=split(',',$head);
    for (my $j=0;$j<scalar @cols;$j++){
        my $colname=$cols[$j];
        chomp($colname);
        $cols[$j]=$colname.'-'.$labels[$i];
        print $cols[$j],"\n";
    }
    push (@header,@cols);

    while (my $entry = <DATA>) {
        chomp $entry;
        my @line=split(",",$entry);
        if (!$line[11]){
            $line[11]="NA";
        }
        my ($start,$end,$mutants,$insertions,$TAsites,$ratio,$pval,$avg,$var,$sd,$se,$genes)=@line;
        $avg=sprintf("%.3f",$avg);

        #if (!$line[6]){
        #    $line[6]=0;
        #}
        #@line=@line[0,1,2,3,4,5,6];

        
        #PUT GENE AND INFO INTO THE HASH FOR EXISTING KEY (1ST FILE) OR CREATE NEW KEY (2ND FILE)
        my $key=$start;
        if(!exists $all{$key}){
            my @info;
            push (@info,@line);
            $all{$key}=\@info;
        }
        #Otherwise the window existed for the prior file, and now need to calcualte extra stats
        else{
            my @info=@{$all{$key}};
            my ($start2,$end2,$mutants2,$insertions2,$TAsites2,$ratio2,$pval2,$avg2,$var2,$sd2,$se2,$genes2)=@info;
            my $diff=sprintf("$round",($avg-$avg2));
            my $df=$insertions+$insertions2-2;
            my $tdist;
            my $pval;

            # CHECK TO MAKE SURE ALL VARIABLES IN TDIST,PVAL CALCULATIONS ARE NUMBERS AND NO
            # ZEROS (0) IN THE DENOMINATOR
            
            if ($se eq "NA" or $se2 eq "NA" or $sd eq "NA" or $sd2 eq "NA" or $insertions==0 or $insertions2==0 or $sd==0 or $sd2==0 or $df<=0){
                ($tdist,$pval)=(" "," ");
            }
            else{
                $tdist=sqrt((($diff)/(sqrt((($sd**2)/$insertions)+(($sd2**2)/$insertions2))))**2);
                $pval=Statistics::Distributions::tprob($df,$tdist);
            }
            push (@info,@line,$diff,$df,$tdist,$pval);
            $all{$key}=\@info;
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
if (!$sortkey){
       $sortkey=26; #default: sort by difference of means
    }
my @sorted = sort { $b->[$sortkey] <=> $a->[$sortkey] } @unsorted;

#ADD NEW FIELDS TO HEADER (COLUMN NAMES)
my $field="Mean".$labels[0].'.'.$labels[1];
push (@header,$field,"DOF","TDIST","PVALUE");

#PRINT TO OUT FILE
open OUT, '>',"$out";
print OUT (join(',',@header),"\n");
foreach (@sorted){
    my @woo=@{$_};
    print OUT join(',',@woo),"\n";
}
close OUT;


