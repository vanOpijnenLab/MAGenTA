#!/usr/bin/perl -w

#Margaret Antonio 15.10.18

#Description: Given a tab delimited text file of genome region name, start position, and end position, customWindow.pl
    #outputs Tn-Seq stats for those regions: # of mutants, # of unique insertions, # of TA sites, ratio of insertions
    #to TA sites, aggregate fitness value, p-value for essentiality


#perl ../Bluberries/customWindow.pl  --ref=NC_003028b2.gbk --fasta tigr4_genome.fasta
        #results/L1_2394eVI_PennG.csv results/L3_2394eVI_PennG.csv results/L4_2394eVI_PennG.csv results/L5_2394eVI_PennG.csv results/L6_2394eVI_PennG.csv
        #--log --custom <file>

#For custom file, accepts any csv file with start and end as first two columns.


use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use Bio::SeqIO;
use Data::Random qw(:all);
use List::Util qw(sum);
use List::BinarySearch qw( :all );
use List::BinarySearch::XS;
use List::MoreUtils qw(uniq);
use File::Path;
use File::Basename;
use feature qw/say/;
use autodie;

use Text::CSV;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    print "\nRequired:\n";
    print "In the command line (without a flag), input the name(s) of the file(s) containing fitness values for individual insertion mutants.\n";
    print "\n slidingWindow.pl <OPTIONS> <OUTPUT> <--essentials AND --genome> <INPUT CSV FILES>\n\n";
    
    print "\nOPTIONS:\n";
    print "--size \t The size of the sliding window(default=500) \n";
    print "--step \t The window spacing (default=10) \n";
    print "--csv \t Name of a file to enter the .csv output for sliding windows.\n";
    print "--cutoff \tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    print "--essentials \t Calculate genome region essentiality based on transposon insertion representation\n";
    print "--outdir\tSpecify name of new directory for all output files\n";
    print "--log\t Send all output to a log file instead of the terminal";
    print "--usage\t Print usage\n";
    print "--tan\t Max number of TA sites in each window---used for creating null distribution library (default: 100)";
    
    
    print "\nREQUIRED: Must choose at least one type of output:\n";
    print "--wig\tCreate a wiggle file for viewing in a genome browser. Provide a filename. Also provide genome under --ref\n";
    print "--txt\t Output all data [start,end,W,count] into a text of bed file.\n";
    print "--txtg\t If consecutive windows have the same value, then group them into one window. Ouput into txt file or bed file.\n";
    print "--ref\tThe name of the reference genome file, in GenBank format. Needed for wig and txt file creation\n";
    print "--custom\t Tab delimited text file for custom window regions : name | start | end \n";
    print "--indir\t Input file directory\n";

}


#ASSIGN INPUTS TO VARIABLES
our ($round,$random,$txt,$txtg,$cutoff,$wig,$infile, $csv, $step, $h, $outdir,$size,$fasta, $log, $ref_genome,$tan,$custom,$indir,$label);
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
'random:s' =>\$random,
'round:i' =>\$round,
'fasta:s' => \$fasta,
'outdir' =>\$outdir,
'log' => \$log,
'usage' => \$h,
'tan'=>\$tan,
'c:s'=>\$custom,
'indir:s'=>\$indir,
'label'=>\$label,
);

sub cleaner{
	my $line=$_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}
	
sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}


#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step=10 };   #set the default step size to 10
if (!$cutoff) {$cutoff=15};
if (!$tan){$tan=100};
if (!$round){$round='%.3f';}
if (!$outdir){
    $outdir="customWindow";
}
mkpath($outdir);
if (!$custom){ print "Must give custom file"; die;}

print "Window size: $size\n";
print "Step value: $step\n";
print "Cutoff: $cutoff\n";
# Just to test out the script opening
if ($h){
    print print_usage(),"\n";
    print "\n";
}
if ($csv){print "CSV output file: ", $csv,"\n";}
if ($txt){print "Text file for window data: $txt\n";}
if ($txtg){print "Text file for grouped windows: $txtg\n";}


if ($log){
    print "\nSending all output to log file\n";
    # redirect STDOUT to log.txt
    open (STDOUT, ">>$outdir/log.txt");
}

#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)

my $rowCount=-1;
my $last=0;
my @unsorted;
my @insertPos; #array to hold all positions of insertions. Going to use this later to match up with TA sites

print "\n---------Importing files--------\n";
print "\tStart input array ",get_time(),"\n";



#Go through each file from the commandline (ARGV array) and read each line as an array into select array if values satisfy the cutoff
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
print "\n---------Importing files--------\n";
print "\tStart input array ",get_time(),"\n";
print "\tNumber of csv files: ", $num,"\n";
my @header;

#Go through each file from the commandline (ARGV array) and read each line as an array into select array if values satisfy the cutoff
for (my $i=0; $i<$num; $i++){   #Read files from ARGV
    print "File #",$i+1,"\t";
    
    my $file=$files[$i];
    print $file,"\n";
  
    open(DATA, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    my $line=<DATA>;
    @header=split(',',$line);
    push (@header,"pval");
    
    while (my $entry = <DATA>) {
    	$entry=cleaner($entry);
		my @line=split(",",$entry);
        my $w = $line[12];
        if (!$w){next;} # For blanks
        else{
            my $c1 = $line[2];
            my $c2 = $line[3];
            my $avg = ($c1+$c2)/2;
            if ($avg > $cutoff) {
                my @select=($line[0],$line[12],$line[9]);
                my $select=\@select;
                push(@unsorted,$select);
                push(@insertPos,$line[0]);   #keep track of actual insertion site position
                $last=$select[0];
                $rowCount++;  
            }
        }
    }
    close DATA;
}


my @sorted = sort { $a->[0] <=> $b->[0] } @unsorted;

@insertPos = sort { $a <=> $b } @insertPos;
@insertPos= uniq @insertPos;

print "\n\tFinished input array ",get_time(),"\n";


print "\n---------Creating sliding windows across the genome--------\n\n";

my $index=-1;
my $marker=0;
my $totalInsert=0;
my $totalWindows=0;


#SUBROUTINE FOR EACH WINDOW CALCULATION
sub OneWindow{
    my $Wstart=shift @_;
    my $Wend=shift@_;
    my $Wavg=0;
    my $Wcount=0;
    my $insertion=0;
    my $Wsum=0;
    my $lastPos=0;
    my $i;
    my @allgenes=();

    
    for ($i=$marker;$i<$rowCount;$i++){
        my @fields=@{$sorted[$i]};
        my $gene=$fields[2];
        if ($fields[0]<$Wstart){  #if deleted, error shows up
            next;
        }
        if ($fields[0]<=$Wend){
            if ($fields[0]<($Wstart+$step)){
                $marker++;
            }
            $Wsum+=$fields[1];
            $Wcount++;
            
            $gene=~ s/"//g;
            $gene=~ s/ //g;
            $gene=~ s/'//g;
            if (!($gene =~ /^ *$/)){
                push @allgenes,$gene;
            }
            
            if ($fields[0]!=$lastPos){
                $insertion+=1;
                $lastPos=$fields[0];
            }
        }
        
        else{   #if ($fields[0]>$Wend){         #if finished with that window, then:
            if ($Wcount!=0){
                $Wavg=sprintf("$round",$Wsum/$Wcount);
            }
            else{
                $Wavg=0;
            }
            
            @allgenes= uniq (@allgenes);
			my @sortedGenes = sort { lc($a) cmp lc($b) } @allgenes;
			#foreach (@sortedGenes){
			#	print TEST "Array: ",$_,"\t";
			#}
			#print TEST "\n";
			my $wind_genes=join(" ",@sortedGenes);
			#print TEST $wind_genes,"\n";		
       
            my @window=($Wstart,$Wend,$Wavg,$Wcount,$insertion,$wind_genes);
    
            $totalWindows++;
            $totalInsert+=$insertion;
            #print @Wwindow, "\n";
            return (\@window);
 #Because count=0 (i.e. there were no insertion mutants in that window)
        }
    }
}




###############################################################################################


#WHILE LOOP TO CALL THE ONE WINDOW SUBROUTINE FOR CALCULATIONS
#Read in user specified windows. Here, use rois.txt for test file.

print "Start calculation: ",get_time(),"\n";

open (CUST, '<', $custom);
my $dummy=<CUST>;
my @allWindows;
while(my $line=<CUST>){
    $line=cleaner($line);
    my @cwind=split(",",$line);
    my $start=$cwind[0];
    my $end=$cwind[1];
    my $window=OneWindow($start,$end);
    push (@allWindows,$window);
    

}

close CUST;

print "End calculation: ",get_time(),"\n";


##############################################################################################

#ESSENTIALS: Counting the number of TA sites in the genome and whether an insertion occurred there or not

print "\n---------Assessing essentiality of genome region in each window--------\n\n";

    my @sites;

    #First read fasta file into a string
    my $seqio = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
    my $prev;
    my $total=0;
    while(my $seq = $seqio->next_seq) {
        $fasta = $seq->seq;
    }
    #Get number of "TA" sites, regardless of insertion---this is just the fasta file
    my $x="TA";
    my @c = $fasta =~ /$x/g;
    my $countTA = @c;
    
    #At what positions in the genome do TA sites occur?
    my $pos=0;
    my $countInsert=0;
    my @selector;   #use this array to hold all 1 and 0 of whether insertion happened or not.

    #Go through all TA sites identified above and see if an insertion occurs there.
    #Push results onto two arrays a 2d array with the position and a 1D array with just the 0 or 1
   
    my @unmatched;  #hold all unmatched ta sites
    my @allTAsites; #2d array to hold all occurences of TA sites in genome
    my $unmatchedCount=0;
    my $offset=0;
    my $genPos = index($fasta,'TA',$offset);

while (($genPos != -1) and ($pos!=scalar @insertPos)) { #as long as the TA site is foun
    my $res=0;
    if ($genPos>$insertPos[$pos]){
        push @unmatched,$insertPos[$pos];
        $unmatchedCount++;
        $pos++;
    }
    if ($genPos==$insertPos[$pos]){
        $res=1;
        $countInsert++;
        $pos++;
    }
    my @sites=($genPos,$res);
    push @selector,$res;
    #push the 0 or 1 onto the array @selector---going to use this to draw random sets for the null distribution
    push (@allTAsites,\@sites);
    $offset = $genPos + 1;
    $genPos = index($fasta, 'TA', $offset);
    $countTA++;
}
my $FILE1 = "$outdir/allTAsites.txt";
open (ALL_TA, ">", $FILE1);
foreach my $sit(@allTAsites){
    foreach (@$sit){
        print ALL_TA $_, "\t";
    }
    printf ALL_TA "\n";
}
close ALL_TA;
my $FILE2 = "$outdir/unmatched.txt";
unless(open UNM, ">", $FILE2){
    die "\nUnable to create $FILE2:\n$!";
}
foreach (@unmatched){
    print UNM $_, "\n";
}
close UNM;

print "\n\nTotal of unmatched insertions: $unmatchedCount\n";
print "See unmatched.txt for genome indices of unmatched sites\n";
print "See allTAsites.txt for details on all TA sites and insertions\n\n";

#print "\nTotal: $countInsert insertions in $countTA TA sites.\n";

#-----------------------------------------------------------------------------------------------------------------------------

    #Now, have an array for each TA site and if an insertion occurred there. So per site @sites=(position, 0 or 1 for insertion).
    #Next step, create null distribution of 10,000 random sets with same number of TA sites as the window and then calculate p-value
  
#SUBROUTINE FOR MAKING THE NULL DISTRIBUTION SPECIFIC TO THE WINDOW

#DISTRIBUTION'S STANDARD DEVIATION
my $N=10000;

sub mean {
	return sum(@_)/@_;
}
sub stdev{
	my @data = @{shift @_};
	my $average=shift @_;
	my $sqtotal = 0;
	foreach(@data) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / ($N-1)) ** 0.5;
	return $std;  
}



#MAKE LIBRARY OF NULL DISTRIBUTIONS:
print "Making a library of null distributions.\n For information on distribution library, see nullDist.txt\n";

my @distLib;

my $FILE3 = "$outdir/nullDist.txt";
unless(open DIST, ">", $FILE3){
	die "\nUnable to create $FILE3:\n$!";
}
printf DIST "Sites\tSize(n)\tMin\tMax\tMean\tstdev\tvar\n";

#Loop once for each distribution in the library 
#(i.e. distribution of 10,000 sites each with 35 TA sites, then a distribution of 10,000 sites each with 36 TA sites, etc)
    
for (my $sitez=1; $sitez<=$tan;$sitez++){ 
    #print "In the first for loop to make a null distribution\n";
    my @unsorted;
    my $count=0;
    my $sum=0;
    
    for (my $i=1; $i<=$N; $i++){
        my @random_set = rand_set( set => \@selector, size => $sitez);
        my $setMean=mean(@random_set);
        push (@unsorted, $setMean);
        #print "$i:\t", "$setMean\n";
        $count++; 
        $sum+=$setMean;
    }
    my @nullDist= sort { $a <=> $b } @unsorted;
    my $nullMean=sprintf("$round",($sum/$count));
    my $standev =sprintf("$round",stdev(\@nullDist, $nullMean));
    my $variance=sprintf("$round",($standev*$standev));
    my $min=sprintf("$round",$nullDist[0]);
    my $max=sprintf("$round",$nullDist[scalar @nullDist-1]);
    my $setScalar=scalar @nullDist;
    printf DIST "$sitez\t$N\t$min\t$max\t$nullMean\t$standev\t$variance\n";
    push (@distLib,\@nullDist);
    
}
close DIST;

#SUBROUTINE TO CALCULATE THE P-VALUE OF THE WINDOW INSERTIONS AGAINST THE NULL DISTRIBUTION

sub pvalue{
    
    #takes in window count average (countAvg) and number of TAsites and makes a null distribution to calculate the pvalue, which it returns
    my $countAvg=shift@_;
    my $TAsites=shift@_;
    my $N=10000;
    my $rank= binsearch_pos { $a cmp $b } $countAvg, @{$distLib[$TAsites+1]};
    my $pval=$rank/$N; #calculate pval as rank/N
    return $pval;
    
}

print "\nThe size of genome is: ", length($fasta), " bp\n";

#Now we have an array called @allTAsites which contains every TAsite position with a 0 next to it for "no insertion".
#Now just need to replace 0 with 1 if there IS and insertion at that site


    my @newWindows=();
    my $printNum=0;

print "Start p-value calculation for individual windows: ",get_time(),"\n\n";

#my $allWindows=\@allWindows;
for (my $i=0;$i<scalar @allWindows;$i++){
    my @win=@{$allWindows[$i]};
    my $starter=$win[0];
    my $ender=$win[1];
    my $length=$ender-$starter;
    my $seq = substr($fasta,$starter-1,$length);  #start-1 becase $start and $end are positions in genome starting at 1,2,3.... substr(string,start, length) needs indexes
    my $ta="TA";
    my @c = $seq =~ /$ta/g;
    my $TAsites = scalar @c;
    push(@win,$TAsites);
    my $ratio=sprintf("$round",($win[4]/$TAsites));
    push (@win,$ratio);
    my $pval=pvalue($ratio,$TAsites);
    push (@win,$pval);
    push (@newWindows,\@win);
    $printNum++;
    }

   
print "End p-value calculation: ",get_time(),"\n\n";

#print "This is the TA count: $count\nTotal genome size is: $total\n\n";

#-------------------------------------------Essentials OUTPUTS-------------------------------------

#MAKE OUTPUT CSV FILE WITH ESSENTIAL WINDOW CALCULATIONS

    print "Start csv ouput file creation: ",get_time(),"\n";
    my $csvBIG = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  
    # open in append mode
    
	open (my $FH8, ">$outdir/customWindows.csv");
    my @head=("start", "end","fitness","mutant_count","insertions","TA_sites","ratio","p-value");
    $csvBIG->print($FH8, \@head); #header
    foreach my $winLine(@newWindows){
        $csvBIG->print($FH8,$winLine);
    }
    close $FH8;
    print "End csv ouput file creation: ",get_time(),"\n\n";

my $in = Bio::SeqIO->new(-file=>$ref_genome);
    my $refseq = $in->next_seq;
    my $refname = $refseq->id;

#MAKE essentials WIG FILE---->later make BW--->IGV
 sub printwig{
    print "Start wig file creation: ",get_time(),"\n";
    my $in = Bio::SeqIO->new(-file=>$ref_genome);
    my $refseq = $in->next_seq;
    my $refname = $refseq->id;
    my $ewig="essentialWindows.wig";
    my $FILE5 ="$outdir/$ewig";
    open (eWIG, ">$FILE5");
    printf eWIG "track type=wiggle_0 name=$ewig\n";
    printf eWIG "variableStep chrom=$refname\n";
    foreach my $wigLine(@newWindows){
        my @wigFields=@$wigLine;
        my $position=$wigFields[0];
        while ($position<=$wigFields[1]){
        printf eWIG $position, " ",$wigFields[7],"\n";    #7 for pvalue, but 2 for fitness!!!!!!
        $position=$position+1;
        }
        #print  WIG $wigFields[0]," ",$wigFields[2],"\n";
    }
    close eWIG;
    print "End wig file creation: ",get_time(),"\n\n";
    print "If this wig file needs to be converted to a Big Wig, then use USCS program wigToBigWig in terminal: \n \t./wigToBigWig gview/12G.wig organism.txt BigWig/output.bw \n\n";
}
#printwig();

#IF GOING TO MAKE A TEXT FILE FOR BED CONVERSION TO BIGBED, NEED CHROM # IN COLUMN 0
my @ecummulative;
my @temp=split("$ref_genome",".");
my $refName=$temp[0];
if ($txtg or $txt){
    for my $line(@newWindows){
        unshift($line, "$refName");
    }
}
#IF MAKING A REGULAR TEXT FILE fields: [chrom,start,end,fitness,count]
#Don't think this is necessary....what would we do with a regular text file?:
if ($txt){
    print "Start text file creation time: ",get_time(),"\n";
    open my $TXT, '>', "$outdir/$txt" or die "\nUnable to create $txt:\n$!";
    print $TXT (join("\t",@$_),"\n") for @newWindows;
    close $TXT;
    print "End text file creation: ",get_time(),"\n\n";
}




