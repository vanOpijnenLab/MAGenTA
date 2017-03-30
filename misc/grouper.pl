#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: After sliding windows have been filtered and ordered, can be grouped
#must decide to group by fitness or insertion representation (one or the other)

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use List::Util qw(sum);
use List::BinarySearch qw( :all );
use List::BinarySearch::XS;
use List::MoreUtils qw(uniq);
use feature qw/say/;
use autodie;
use Scalar::Util qw(looks_like_number);
use List::MoreUtils qw(uniq);
#no warnings 'numeric';

sub print_usage() {
    print "\n";
    print "grouper.pl: Takes slidingWindow.csv file and groups consecutive like windows\n\n";
    print "Example: Group consecutive windows by fitness. If fitness is the same then expand.\n";
    print "For different strains/genomes, use compStrains.pl\n";
    print "\nUSAGE: perl grouper.pl <options> <slidingWindow.csv>\n\n";
    print "OPTIONS:\n\n";
    print " -h\tPrints usage and quits\n";
    print " -o\tOutput file for grouped windows. Default: groupedWindows.csv\n";
    print " -s\tSort output by this index of the file (indices begin at 0).\n";
    print " -sig\t Group consecutive windows of same pvalue\n";
    print " -fit\t Group consecutive windows of same fitness value\n";
    print " -r\tRound final output numbers to this number of decimals\n";
    print "\n\n";
}


#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
# 0:start,1:end,2:fitness,3:mutant_count,4:insertions,5:TA_sites,6:ratio,7:p-value

#ASSIGN INPUTS TO VARIABLES
our ($h,$out,$sortby,$round,$sig,$key);
GetOptions(
'o:s' =>\$out,
'h' => \$h,
's:i' => \$sortby,
'r:i'=> \$round,
'sig'=>\$sig,
'key'=>\$key,

);

if ($h){
    print_usage();
    exit;
}


#Assign defaults
if (!$out) {$out="groupedWindows.csv"}
if (!$round){$round='%.2f'}
if (!$sortby){$sortby=0}
sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
sub mean {
    return sum(@_)/@_;
}
sub cleaner{
    my $line=$_[0];
    chomp($line);
    $line =~ s/\x0d{0,1}\x0a{0,1}\Z//s;
    return $line;
}

#READ INPUT FILE INTO A 2D ARRAY @windows

my $in=$ARGV[0]; #get input file from first entry in ARGV
my @unsorted;
open(DATA, '<', $in) or die "Could not open '$in' \n";
#Store first line for use as column names of output file
my $line=<DATA>;
$line=cleaner($line);
my @header=split(',',$line);
while (my $entry = <DATA>) {
	$entry=cleaner($entry);
	my @fields=split(',',$entry);
	my $start=$fields[0];
    #if ($order){
        #my $bool=looks_like_number($fields[0]);
        #if ($bool and ($fields[7] !=0)){
            #$fields[7]=sprintf($round,$fields[7]);
        #push (@unsorted,\@fields);
        #}
        
    #}
   # else{
        push (@unsorted,\@fields);
    #}
}
close DATA;

my $count=0; #number of windows added to cumm that will be used for avg calc
my $sumFit=0;
my $sumRatio=0;
my $string;


#If sort by fitness, most interesting regions have high mean differece so sort largest to smallest
my @outarray;

if ($order){
  
    #PRINT COLUMN NAMES (HEADER) TO THE OUTPUT FILE
    @header=("start","end","avg","genes");
    $string = join(',', @header);
    
    my @windows= sort {$b->[$sortby]<=>$a->[$sortby] || $a->[0]<=>$b->[0] } @unsorted;
    
    
    #Start cummulative array to begin for loop
    my $lastLine=scalar @windows;
    my @cumu=@{$windows[0]};
    my $cstart=$cumu[0];
    my $cend=$cumu[1];
    my $cavg=$cumu[7];
    my @cgenes=split(/ /,$cumu[11]);
    
    for (my $i=1;$i<$lastLine;$i++){
        my @field=@{$windows[$i]};
        my $fstart=$field[0];
        my $fend=$field[1];
        my $favg=$field[7];
        my $gene=$field[12];
        my @fgenes=split(/ /,$field[11]);
        #Either this window needs to be added or need to start new cummulative
        
        #Add field window (@field) if overlaps with cumulative window (@cumu)
    	if ($gene eq "intergenic"){
			if (($cend>=$fstart and $cstart<=$fstart) or ($cend>=$fend and $cstart<=$fend)){
				#Keep cstart as it is but change the end coordinate
				$cend=$fend;
				#Append new genes
				foreach my $gene(@fgenes){
					#Don't know there are double quotes are coming from
					$gene =~ s/"//g;
					$gene =~ s/'//g;
					push (@cgenes,$gene);
				}
			}
			#need to output this cumm region with average calcs
			else{
				@cgenes=uniq(sort @cgenes);
				my $allgenes=join(" ",@cgenes);
				my @final=($cstart,$cend,$cavg,$allgenes);
				push (@outarray,\@final);
				#Set up current entry as new cumulative
				@cumu=@field;
				$cstart=$cumu[0];
				$cend=$cumu[1];
				$cavg=$cumu[7];
				@cgenes=split (/ /,$cumu[11]);
			}
    	}
    }
}

#If sort by significance, most interesting regions are low pvals, so sort smallest to largest
else{
    
    #PRINT COLUMN NAMES (HEADER) TO THE OUTPUT FILE
    @header=("start","end","pvalue","genes");
    $string = join(',', @header);
    
    my @windows= sort {$a->[$sortby]<=>$b->[$sortby] || $a->[0]<=>$b->[0] } @unsorted;
    
    #Start cummulative array to begin for loop
    my $lastLine=scalar @windows;
    my @cumu=@{$windows[0]};
    my $cstart=$cumu[0];
    my $cend=$cumu[1];
    my $cpval=$cumu[5];
    my @cgenes=split('',$cumu[11]);
    
    for (my $i=1;$i<$lastLine;$i++){
        my @field=@{$windows[$i]};
        my $fstart=$field[0];
        my $fend=$field[1];
        my $fpval=$field[5];
        my @fgenes=split(' ',$field[11]);
        #Either this window needs to be added or need to start new cummulative
        
        #Add field window (@field) if overlaps with cumulative window (@cumu)
        if (($cend>=$fstart and $cstart<=$fstart) or ($cend>=$fend and $cstart<=$fend)){
            #Keep cstart as it is but change the end coordinate
            $cend=$fend;
            #Append new genes
            foreach my $gene(@fgenes){
                if (! grep( /^$gene$/, @cgenes )){
                    push (@cgenes,$gene);
                }
            }
        }
        #need to output this cumm region with average calcs
        else{
            @cgenes=uniq(sort @cgenes);
            my $allgenes=join('',@cgenes);
            my @final=($cstart,$cend,$cpval,$allgenes);
            push (@outarray,\@final);
            #Set up current entry as new cumulative
            @cumu=@field;
            $cstart=$cumu[0];
            $cend=$cumu[1];
            $cpval=$cumu[5];
            @cgenes=split('',$cumu[11]);
        }
    }
}

#SORT first by whether or not there are annotated genes, then by fitness/pval

@outarray= sort {$a->[3] cmp $b->[3] || $a->[2]<=>$b->[2] } @outarray;

open my $gOut, '>', "$out" or die $!;

#OUTPUT to csv file
print $gOut "$string\n";

foreach (@outarray){
    print $gOut ($_,",") for @{$_};
    print $gOut ("\n");
}

close $gOut;

