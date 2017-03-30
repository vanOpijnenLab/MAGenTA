#!/usr/bin/perl -w

#Margaret Antonio updated 15.09.15

<<<<<<< HEAD
#VERSION 16: combines version 14 and the multiple file input to array and sorting feature from version 13. Still need to: make window cutoffs...

#/Volumes/Macintosh HD/Users/margaretantonio/Documents/TvOp_Lab/Tn_SeqAnalysisScripts/windowFit.pl --ref=NC_003028b2.gbk --cutoff 15 --csv windowFit/16A.csv --step 10 --size 500 --txtg ../viewer/19A.txt --txt ../viewer/19A.txt --wig 19A.wig results/L1_2394eVI_PennG.csv results/L3_2394eVI_PennG.csv results/L4_2394eVI_PennG.csv results/L5_2394eVI_PennG.csv results/L6_2394eVI_PennG.csv
=======
#../Tn_SeqAnalysisScripts/windowFit.pl --csv windowFit/22.csv results/L6_2394eVI_PennG.csv results/L4_2394eVI_PennG.csv >runner.txt
>>>>>>> bc08bd1ae7252375e96a48e5c44e0d346a63a591

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use Text::CSV_XS;
use Bio::SeqIO;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    print "\nDescription:\n";
    print "Integrates multiple files of transposon insertion data and outputs aggregate fitness within a sliding window (specified by size and step). Can ouput files as text, csv, wig.\n";
    print "\nCommand line: windowFit.pl <OPTIONS> <REQ OUTPUT TYPE(S)> <INPUT FILE(S)>\n\n";
    print "\nRequired:\n";
    print "In the command line (without a flag), input the name(s) of the file(s) containing fitness values for individual insertion mutants.\n";
    
    print "\nOptional:\n";
    print "--size \t The size of the sliding window(default=500) \n";
    print "--step \t The window spacing (default=10) \n";
    print "--cutoff \tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    
    print "\nMust choose at least one type of output:\n";
    print "--csv \t Name of a file to enter the .csv output for sliding windows.\n";
    print "--wig\tCreate a wiggle file for viewing in a genome browser. Provide a filename. Also provide genome under --ref\n";
    print "--txt\t Output all data [start,end,W,count] into a text of bed file.\n";
    print "--txtg\t If consecutive windows have the same value, then group them into one window. Ouput into txt file or bed file.\n";
    print "--ref\tThe name of the reference genome file, in GenBank format. Needed for wig and txt file creation\n";
    
    print "--h for help\n\n";
    
}


#ASSIGN INPUTS TO VARIABLES
our ($txt,$txtg,$cutoff,$h, $wig,$ref_genome,$infile, $csv, $step, $size);
GetOptions(
'wig:s' => \$wig,
'ref:s' => \$ref_genome,
'cutoff:i'=>\$cutoff,
'in:s' => \$infile,
'csv:s'  => \$csv,
'step:i' => \$step,
'size:i' => \$size,
'txtg:s' => \$txtg,
'txt:s' => \$txt,
'h'=>\$h,
);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
# Just to test out the script opening

print "\n";
if ($h){
    print_usage();
    exit;
}

if ($csv){print "CSV output file: ", $csv,"\n";}
if ($txt){print "Text file for window data: $txt\n";}
if ($txtg){print "Text file for grouped windows: $txtg\n";}


#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step = 10 };   #set the default step size to 10
if (!$cutoff) {$cutoff=0};
print "Window size: $size\n";
print "Step value: $step\n";
print "Cutoff: $cutoff\n";

if ((!$csv) and (!$txt) and (!$txtg) and (!$wig)) {
    print "\nThere must be some kind of output file assigned (--csv, --txt, or --txtg)\n";
    print_usage();
    die;
}
#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)
print "\nStart input array ",get_time(),"\n";

my $rowCount=-1;
my $last;
my @unsorted;
my $num=$#ARGV+1;
print "\nNumber of files in csv: ", $num,"\n";

my $csvtemp=Text::CSV->new;

for (my $i=0; $i<$num; $i++){   #Read files from ARGV
    my $file=$ARGV[$i];
    open(my $data, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    $csvtemp->getline($data);
    print "\t",$file,"\n";
    while (my $line = $csvtemp->getline($data)) {
        chomp $line;
        my $w = $line->[12];
        #my $temp=$line->[0];
        #if ($temp<500){
        #    print $temp,"\t", $w, "\n";
        #}
        if (!$w){next;} # For blanks
        
        else{
            my $c1 = $line->[2];
            my $c2 = $line->[3];
            my $avg = ($c1+$c2)/2;
            if ($avg < $cutoff) {
                next;
            } # Skip cutoff genes.
            else {
                my @select=($line->[0],$line->[12]);
                my $select=\@select;
                push(@unsorted,$select);
                $rowCount++;
                $last=$select->[0];
                
            }
        }
    }
    close $data;
    
}
my @sorted = sort { $a->[0] <=> $b->[0] } @unsorted;

#Print test of array
#for (my $i=0;$i<30;$i++){
#    foreach my $element ( @{ $sorted[$i] }){print $element,"\t";}
#   print "\n";}

print "Finished input array ",get_time(),"\n";

###############################################################################################

print "\n---------This is the sliding window fitness calculation part--------\n\n";

my $index=-1;
my $marker=0;

#SUBROUTINE FOR EACH WINDOW CALCULATION
sub OneWindow{
    my $Wstart=shift @_;
    my $Wend=shift@_;
    my $Wavg=0;
    my $Wcount=0;
    my $Wsum=0;
    my @indiv;
    my $i;
    for ($i=$marker;$i<$rowCount;$i++){   #looping through whole file, begin at marker until end of file but really stops at window end
        my @fields=@{$sorted[$i]};
        my $site=$fields[0];
        #if ($fields[0]<$Wstart){  #if deleted, error shows up
        #next;
        #}
        if ($site<=$Wend){
            if ($site<($Wstart+$step)){
                $marker++;
            }
            $Wsum+=$fields[1];
            $Wcount++; #print $Wcount, "\t",$fields[0],"\t",$fields[1],"\n";
        }
        else{   #if finished with that window ($site>$Wend) then calculate window fitness
            if ($Wcount!=0){
                $Wavg=sprintf("%.2f",$Wsum/$Wcount);
                my @window=($Wstart,$Wend,$Wavg,$Wcount);
                #print @Wwindow;
                return (\@window);
            }
            else{
                return -1;  #Because count=0 (i.e. there were no insertion mutants in that window)
            }
            
        }
    }
}

print "Start calculation: ",get_time(),"\n";
my $start=1;
my $end=$size+1;
my $windowNum=0;
my @allWindows=(); #will be a 2D array containing all window info to be written into output file

#WHILE LOOP TO CALL THE ONE WINDOW SUBROUTINE FOR CALCULATIONS===INCREMENTS START AND END VALUES OF THE WINDOW

while ($end<=$last-$size){  #100,000bp-->9,950 windows--> only 8500 windows in csv because 0
    my($window)=OneWindow($start,$end);
    if ($window!=-1){
        push (@allWindows,$window);
    }
    $start=$start+$step;
    $end=$end+$step;
    #print "$end  ";
}
print "End calculation: ",get_time(),"\n\n";

#MAKE OUTPUT CSV FILE WITH WINDOW CALCULATIONS
if ($csv){
    print "Start csv ouput file creation: ",get_time(),"\n";
    my $csvBIG = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  # open in append mode
    open my $file, ">", "$csv" or die "Failed to open file";
    $csvBIG->print($file, [ "start", "end","fitness","mutant_count" ]); #header
    foreach my $winLine(@allWindows){
        $csvBIG->print($file,$winLine);
    }
    close $file;
    print "End csv ouput file creation: ",get_time(),"\n\n";
}

#MAKE WIG FILE---->later make BW--->IGV
if ($wig){
    print "Start wig file creation: ",get_time(),"\n";
    my $in = Bio::SeqIO->new(-file=>$ref_genome);
    my $refseq = $in->next_seq;
    my $refname = $refseq->id;
    open WIG, ">$wig";
    print WIG "track type=wiggle_0 name=$wig\n";
    print WIG "variableStep chrom=$refname\n";
    foreach my $wigLine(@allWindows){
        my @wigFields=@$wigLine;
        my $position=$wigFields[0];
        #while ($position<=$wigFields[1]){
        print WIG $position," ",$wigFields[2],"\n";
        #$position=$position+1;
        #}
        #print  WIG $wigFields[0]," ",$wigFields[2],"\n";
    }
    close WIG;
    print "End wig file creation: ",get_time(),"\n\n";
    print "If this wig file needs to be converted to a Big Wig, then use USCS program wigToBigWig in terminal: \n \t./wigToBigWig gview/12G.wig organism.txt BigWig/output.bw \n\n";
}

#IF GOING TO MAKE A TEXT FILE FOR BED CONVERSION TO BIGBED, NEED CHROM # IN COLUMN 0
my @cummulative;
if ($txtg or $txt){
    for my $line(@allWindows){
        unshift($line, "NC_003028");
    }
}
#IF MAKING A REGULAR TEXT FILE fields: [chrom,start,end,fitness,count]
if ($txt){
    print "Start text file creation time: ",get_time(),"\n";
    open my $TXT, '>', "$txt" or die $!;
    print $TXT (join("\t",@$_),"\n") for @allWindows;
    close $TXT;
    print "End text file creation: ",get_time(),"\n\n";
}

#IF MAKING A TEXT FILE OF GROUPED CONSECUTIVE WINDOWS WITH SAME FITNESS
if ($txtg){
    print "Start grouped txt file creation time: ",get_time(),"\n";
    open my $TXTg, '>', "$txtg" or die $!;
    for my $line(@allWindows){
        my @field=@$line;
        if (!@cummulative){
            @cummulative=@field;
            if ($cummulative[4]>1000){$cummulative[4]=1000}
        }
        else{
            if ($cummulative[3]==$field[3]){
                $cummulative[2]=$field[2];
                if ($cummulative[4]<=1000-$field[4]){
                    $cummulative[4]+=$field[4];
                }
            }
            else{
                print $TXTg ($_,"\t") for @cummulative;
                print $TXTg ("\n");
                @cummulative=@field;
            }
        }
    }
    close $TXTg;
    print "End grouped text file creation: ",get_time(),"\n\n";
    if ($txt or $txtg){
        print "\nTo make a BigBed file from this text file, rename file to .bed and use USCS program bedToBigBed in terminal \n\t\n";
    }
}





