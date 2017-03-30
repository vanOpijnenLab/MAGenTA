#!/usr/bin/perl -w

#Margaret Antonio 15.08.129

use strict;
use Getopt::Long;
use Set::Scalar;
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
no warnings;
use Text::CSV;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {

    print "\n####################################################################\n";
    print "USAGE:\n";
    print "slidingWindow.pl -d inputs/ -o slidWind_output/ -f genome.fasta -r genome.gbk\n";
    
    print "\nREQUIRED:\n";
    print " -d\tDirectory containing all input files (output files from\n";
    print "   \tcalcFitness tool)\n";
    print "   \tOR\n";
    print "   \tIn the command line (without a flag), input filename(s)\n";
    print " -f\tFilename for genome sequence, in fasta format\n";
    print " -r\tFilename for genome annotation, in GenBank format\n";
    
    print "\nOPTIONAL:\n";
    print " -h\tPrint usage\n";
    print " -size\tThe size of the sliding window. Default=50\n";
    print " -step\tThe window spacing. Default=10\n";
    print " -x\tExclude values with avg. counts less than x where (c1+c2)/2<x\n";
    print "   \tDefault=15\n";
    print " -log\tSend all output to a log file instead of the terminal\n";
    print " -max\tExpected max number of TA sites in a window.\n";
    print "     \tUsed for creating null distribution library. Default=100\n";
    print " -o\tSpecify name of new directory for all output files\n";
    print "   \tDefault=sw_out/\n";
    print " -w\tDo weighted average for fitness per insertion\n";
    print " -wc\tInteger value for weight ceiling. Default=50\n";
    
    print " \n~~~~Always check that file paths are correctly specified~~~~\n";
    print "\n##################################################################\n";
    
    #print "--wig\tCreate a wiggle file for viewing in a genome browser. Provide a filename.\n";
    #print "Also provide genome under --ref\n";
    #print "--txt\t Output all data [start,end,W,count] into a text of bed file.\n";
    #print "--txtg\t If consecutive windows have the same value,\n";
    #print "then group them into one window. Ouput into txt file or bed file.\n";
    #print "--ref\tThe name of the reference genome file, in GenBank format.\n";
    #print "Needed for wig and txt file creation\n";
    #print "--c\tComma separated list of start and end coordinates for custom windows\n";
}

#ASSIGN INPUTS TO VARIABLES
our ($round,$cutoff,$step, $size, $help, $outdir,$fasta, $log, $ref_genome,$tan,$indir,$weight_ceiling,$weight,$custom);
GetOptions(

'x:i'=>\$cutoff,
'in:s' => \$infile,
'csv:s'  => \$csv,
'step:i' => \$step,
'size:i' => \$size,
'round:i' =>\$round,
'f:s' => \$fasta,
'r:s' => \$ref_genome,
'o:s' =>\$outdir,
'log' => \$log,
'h' => \$help,
'max:i'=>\$tan,
'd:s'=>\$indir,
'wc:i'=>\$weight_ceiling,
'w'=>\$weight,
'c:s'=>\$custom,
);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
# Just to test out the script opening
if ($help){
	print_usage();
	print "\n";
    exit;
}
if (!$round){$round='%.3f';}
if (!$outdir){$outdir="sw_out";}
if (!$weight_ceiling){$weight_ceiling=50;}
if (!$cutoff){$cutoff=15;}
mkpath($outdir);

if ($log){
	print "\nSending all output to log file\n";
	# redirect STDOUT to log.txt
    open (STDOUT, ">>$outdir/log.txt");
}


#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step=10 };   #set the default step size to 10
if (!$cutoff) {$cutoff=15};
if (!$tan){$tan=100};

print "Window size: $size\n";
print "Step value: $step\n";
print "Cutoff: $cutoff\n";

#if ((!$csv) and (!$txt) and (!$txtg) and (!$wig)) {
#print "\nThere must be some kind of output file assigned (--csv, --txt, or --txtg)\n";
#print_usage();
#die;
#}

# Returns mean, variance, sd, se
sub average {
    my $scoreref = shift @_;
    my @scores = @{$scoreref};
    
    my $sum=0;
    my $num=0;
    
    # Get the average.
    foreach my $w (@scores) {
        $sum += $w;
        $num++;
    }
    my $average= $sum/$num;
    my $xminusxbars = 0;
    
    # Get the variance.
    foreach my $w (@scores) {
        $xminusxbars += ($w-$average)**2;
    }
    my $variance;
    if ($num<=1){
        $variance=0.10;
    }
    else{
        $variance = sprintf($round,(1/($num-1)) * $xminusxbars);
    }
    my $sd = sprintf($round,sqrt($variance));
    my $se = sprintf($round,$sd / sqrt($num));
    
    return ($average, $variance, $sd, $se);
    
}

sub cleaner{
    my $line=$_[0];
    chomp($line);
    $line =~ s/\x0d{0,1}\x0a{0,1}\Z//s;
    return $line;
}

# Takes two parameters, both hashrefs to lists.
# 1) hashref to list of scores
# 2) hashref to list of weights, equal in length to the scores.
sub weighted_average {
    
    my $scoreref = shift @_;
    my $weightref = shift @_;
    my @scores = @{$scoreref};
    my @weights = @{$weightref};
    
    my $sum=0;
    my ($weighted_average, $weighted_variance)=(0,0);
    my $v2;
    
    # Go through once to get total, calculate V2.
    for (my $i=0; $i<@weights; $i++) {
        $sum += $weights[$i];
        $v2 += $weights[$i]**2;
    }
    if ($sum == 0) { return 0; } # No scores given?
    
    my $scor = join (' ', @scores);
    my $wght = join (' ', @weights);
    #print "Scores: $scor\n";
    #print "Weights: $wght\n";
    
    # Now calculated weighted average.
    my ($top, $bottom) = (0,0);
    for (my $i=0; $i<@weights; $i++) {
        $top += $weights[$i] * $scores[$i];
        $bottom += $weights[$i];
    }
    $weighted_average = sprintf($round,$top/$bottom);
    #print "WA: $weighted_average\n";
    
    ($top, $bottom) = (0,0);
    # Now calculate weighted sample variance.
    for (my $i=0; $i<@weights; $i++) {
        $top += ( $weights[$i] * ($scores[$i] - $weighted_average)**2);
        $bottom += $weights[$i];
    }
    $weighted_variance = sprintf($round,$top/$bottom);
    #print "WV: $weighted_variance\n";
    
    my $weighted_stdev = sprintf($round,sqrt($weighted_variance));
    my $weighted_stder = sprintf($round,$weighted_stdev / sqrt(@scores));  # / length scores.
    
    #print "$weighted_average, $weighted_variance, $weighted_stdev\n";
    return ($weighted_average, $weighted_variance, $weighted_stdev, $weighted_stder);
}

#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)

my $rowCount=-1;
my $last=0;
my @unsorted;
my @insertPos; #array to hold all positions of insertions. Going to use this later to match up with TA sites
if (!$weight_ceiling){$weight_ceiling=999999;}
my %wind_summary;

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

my $datestring = localtime();
print "Local date and time $datestring\n";

print "\n---------Importing files--------\n";
print "\tStart input array ",get_time(),"\n";
print "\tNumber of csv files: ", $num,"\n";


#Go through each file from the commandline (ARGV array) and read each line as an array
#into select array if values satisfy the cutoff
for (my $i=0; $i<$num; $i++){
	print "File #",$i+1,"\t";
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    my $dummy=<DATA>; #read and store column names in dummy variable
    while (my $entry = <DATA>) {
    	chomp $entry;
		my @line=split(",",$entry);
        my $locus = $line[9]; #gene id (SP_0000)
        my $w = $line[12]; #nW
        if (!$w){ $w=0 }   # For blanks
        my $c1 = $line[2];
        my $c2 = $line[3];
        my $avg = ($c1+$c2)/2;
        #Average counts must be greater than cutoff (minimum allowed)
        if ($avg > $cutoff) {
            my @select=($line[0],$w,$avg,$line[9]);
            my $select=\@select;
            push(@unsorted,$select);
            push(@insertPos,$line[0]);   #keep track of actual insertion site position
            $last=$select[0];
            $rowCount++;
        }
        if ($avg >= $weight_ceiling) { $avg = $weight_ceiling; } # Maximum weight
    }
    close DATA;
}

my @sorted = sort { $a->[0] <=> $b->[0] } @unsorted;
@insertPos = sort { $a <=> $b } @insertPos;
@insertPos= uniq @insertPos;

print "\n\tFinished input array ",get_time(),"\n";

###############################################################################################

print "\n---------Creating sliding windows across the genome--------\n\n";

my $index=-1;
my $marker=0;
my $totalInsert=0;
my $totalWindows=0;

open TEST,'>',"testQuotes.txt";

#SUBROUTINE FOR EACH WINDOW CALCULATION
sub OneWindow{
    my $Wstart=shift @_;
    my $Wend=shift@_;
    my $Wcount=0;
    my $insertion=0;
    my $Wsum=0;
    my $lastPos=0;
    my $i;
    my @allgenes=();
    
    for ($i=$marker;$i<$rowCount;$i++){
        my @fields=@{$sorted[$i]};
        my $pos=$fields[0];
        my $w=$fields[1];         #fitness value for single insertion
        my $avgCount=$fields[2];
        my $gene=$fields[3];
        
        if ($pos<$Wstart){  #if deleted, error shows up
            next;
        }
        if ($pos<=$Wend){
            if ($pos<($Wstart+$step)){
                $marker++;
            }
            my @empty;
            if (!$wind_summary{$Wstart}) {
                $wind_summary{$Wstart}{w} = [@empty];
                $wind_summary{$Wstart}{s} = [@empty];
                $wind_summary{$Wstart}{g} = [@empty];
            }
            # Hash of Fitness scores.
            $wind_summary{$Wstart}{w} = [@{$wind_summary{$Wstart}{w}}, $w];
            # Hash of counts used to generate those fitness scores.
            $wind_summary{$Wstart}{s} = [@{$wind_summary{$Wstart}{s}}, $avgCount];
            
            
            if (!($gene =~ /^ *$/)){
                $gene=~ s/"//g;
                $gene=~ s/ //g;
                $gene=~ s/'//g;
                push @allgenes,$gene;
            }
            $Wsum+=$w;
            $Wcount++;
            if ($pos!=$lastPos){
                $insertion+=1;
                $lastPos=$pos;
            }
        }
        #if ($fields[0]>$Wend) #finished with that window, then:
        else{
            if ($Wcount>0){
                $totalWindows++;
                $totalInsert+=$insertion;
                
                my ($average, $variance, $stdev, $stderr);
                
                if ($num <=1 ) {
                    ($average, $variance, $stdev, $stderr)=(0.10,0.10,"NA","NA");
                }
                else{
                    if (!$weight) {
                        ($average, $variance, $stdev, $stderr) = &average($wind_summary{$Wstart}{w});
                    }
                    else {
                        ($average, $variance, $stdev, $stderr)= &weighted_average($wind_summary{$Wstart}{w},$wind_summary{$Wstart}{s});
                    }
                }
                
                @allgenes= uniq (@allgenes);
                my @sortedGenes = sort { lc($a) cmp lc($b) } @allgenes;
                my $wind_genes=join(" ",@sortedGenes);
                my @window=($Wstart,$Wend,$Wcount,$insertion,$wind_genes,$average,$variance,$stdev,$stderr);
                return (\@window);
            }
            
            else{ 
            #Even if there were no insertions, still want window in file for consistent start/end
                my @window=($Wstart,$Wend,$Wcount,$insertion," ","NA","NA","NA","NA");
                return (\@window);
            }  	#Because count=0 (i.e. there were no insertion mutants in that window)
        }
    }
}

sub customWindow{

    my $Wstart=shift @_;
    my $Wend=shift@_;
    my $Wavg=0;
    my $Wcount=0;
    my $insertion=0;
    my $Wsum=0;
    my $lastPos=0;
    my @allgenes=();
    my $i;
    
    
    for ($i=$marker;$i<$rowCount;$i++){
        my @fields=@{$sorted[$i]};
        my $pos=$fields[0];
        my $w=$fields[1];         #fitness value for single insertion
        my $avgCount=$fields[2];
        my $gene=$fields[3];


        
        if ($pos<$Wstart){  #if deleted, error shows up
            next;
        }
        if ($pos<=$Wend){
            my @empty;
            if (!$wind_summary{$Wstart}) {
                $wind_summary{$Wstart}{w} = [@empty];
                $wind_summary{$Wstart}{s} = [@empty];
                $wind_summary{$Wstart}{g} = [@empty];
            }
            # Hash of fitness scores
            $wind_summary{$Wstart}{w} = [@{$wind_summary{$Wstart}{w}}, $w];
            # Hash of counts used to generate those fitness scores.
            $wind_summary{$Wstart}{s} = [@{$wind_summary{$Wstart}{s}}, $avgCount];

            if (!($gene =~ /^ *$/)){
                $gene=~ s/"//g;
                $gene=~ s/ //g;
                $gene=~ s/'//g;
                push @allgenes,$gene;
            }
            $Wsum+=$w;
            $Wcount++;
            if ($pos!=$lastPos){
                $insertion+=1;
                $lastPos=$pos;
            }
        }
        
        #if ($fields[0]>$Wend){         #if finished with that window, then:
        else{
            if ($Wcount>0){
                $totalWindows++;
                $totalInsert+=$insertion;
                
                my ($average, $variance, $stdev, $stderr);
                
                if ($num <=1 ) {
                    ($average, $variance, $stdev, $stderr)=(0.10,0.10,"NA","NA");
                }
                else{
                    if (!$weight) {
                        ($average, $variance, $stdev, $stderr) = &average($wind_summary{$Wstart}{w});
                    }
                    else {
                        ($average, $variance, $stdev, $stderr)= &weighted_average($wind_summary{$Wstart}{w},$wind_summary{$Wstart}{s});
                    }
                }
                
                @allgenes= uniq (@allgenes);
                my @sortedGenes = sort { lc($a) cmp lc($b) } @allgenes;
                my $wind_genes=join(" ",@sortedGenes);
                print TEST $wind_genes,"\n";
                my @window=($Wstart,$Wend,$Wcount,$insertion,$wind_genes,$average,$variance,$stdev,$stderr);
                return (\@window);
            }
            
            else{ #Even if there were no insertions, still want window in file for consistent start/end
                my @window=($Wstart,$Wend,$Wcount,$insertion," ","NA","NA","NA","NA");
                return (\@window);
            }  	#Because count=0 (i.e. there were no insertion mutants in that window)
        }
    }
}

print "Start calculation: ",get_time(),"\n";
my $start=1;
my $end=$size+1;
my $windowNum=0;
my @allWindows; #will be a 2D array containing all window info to be written into output file



#IF A CUSTOM LIST OF START/END COORDINATES IS GIVEN THEN CREATE WINDOWS USING IT
if ($custom){
    open (CUST, '<', $custom);
    my $dummy=<CUST>;
    while(my $line=<CUST>){ 
        $line=cleaner($line);
        my @cwind=split(",",$line);
        my $start=$cwind[0];
        my $end=$cwind[1];
        my $window=customWindow($start,$end);
        push @allWindows,$window;
    }
    close CUST;
}

else{ #if not custom list
    #WHILE LOOP TO CALL THE ONE WINDOW SUBROUTINE FOR CALCULATIONS===INCREMENTS START AND END VALUES OF THE WINDOW
    while ($end<=$last-$size){  #100,000bp-->9,950 windows--> only 8500 windows in csv because 0
        my($window)=OneWindow($start,$end);
        push (@allWindows,$window);
        $start=$start+$step;
        $end=$end+$step;
    }
    print "End calculation: ",get_time(),"\n";
}


#my $avgInsert=$totalInsert/$totalWindows;
#print "Average number of insertions for $size base pair windows: $avgInsert\n";

close TEST;

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

while (($genPos != -1) and ($pos<scalar @insertPos)) { #as long as the TA site is foun
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

#--------------------------------------------------------------------------------------------------
    
#Now, have an array for each TA site and if an insertion occurred there.
    #So per site @sites=(position, 0 or 1 for insertion).
#Next step, create null distribution of 10,000 random sets with same number of
    #TA sites as the window and then calculate p-value
  
#SUBROUTINE FOR MAKING THE NULL DISTRIBUTION SPECIFIC TO THE WINDOW

#DISTRIBUTION'S STANDARD DEVIATION
my $N=10000;

sub mean {
    if (scalar @_ ==0){
        return 0;
    }
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
printf DIST "Sites\tSize(n)\tMin\tMax\tMean\tstdev\tvar\tMinPval\n";

#Loop once for each distribution in the library 
#(i.e. distribution of 10,000 sites each with 35 TA sites, then a distribution of 10,000 sites each with 36 TA sites, etc)

sub pvalue{
    
    #takes in window count average (countAvg) and number of TAsites and makes a null distribution to calculate the pvalue, which it returns
    my $mean=shift@_;
    my $TAsites=shift@_;
    my $N=10000;
    my @specDist=@{$distLib[$TAsites-1]};
    my $rank= binsearch_pos { $a cmp $b } $mean,@specDist;
    my $i=$rank;
    while ($i<scalar(@specDist)-1 and $specDist[$i+1]==$specDist[$rank]){
        $i++;
    }
    $rank=$i;
    my $pval=$rank/$N; #calculate pval as rank/N
    return $pval;
    
}

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
    push (@distLib,\@nullDist);
    my $minp=pvalue(0,$sitez);
    printf DIST "$sitez\t$N\t$min\t$max\t$nullMean\t$standev\t$variance\t$minp\n";
}
close DIST;

#SUBROUTINE TO CALCULATE THE P-VALUE OF THE WINDOW INSERTIONS AGAINST THE NULL DISTRIBUTION    
#---------------------------------------------------------------------------------------------
    
    #Now we have an array called @allTAsites which contains every TAsite position 
    #with a 0 next to it for "no insertion".
    #Now just need to replace 0 with 1 if there IS and insertion at that site
    
my @newWindows=();
my $printNum=0;

#my $allWindows=\@allWindows;
for (my $i=0;$i<scalar @allWindows;$i++){
    my @win=@{$allWindows[$i]};
    my $starter=$win[0];
    my $ender=$win[1];
    my $insertions=$win[3];
    my $mutcount=$win[2];
    #print "num $printNum -->\tStart pos: $starter\tEnd pos: $ender\n";
    #How many TA sites are there from $genome[$start] to $genome[$end]?

    my $seq = substr($fasta,$starter-1,$size);  #start-1 becase $start and $end are positions in genome starting at 1,2,3.... substr(string,start, length) needs indexes
    my $ta="TA";
    my @c = $seq =~ /$ta/g;
    my $TAsites = scalar @c;
    my $avgInsert=sprintf("$round",$insertions/$TAsites);
    my $pval=pvalue($avgInsert,$TAsites);
    
    #reorder so things make more sense: all fitness related calcs together and all sig together
    my @new=($starter,$ender,$mutcount,$insertions,$TAsites,$avgInsert,$pval,$win[5],$win[6],$win[7],$win[8],$win[4]);
    
    push (@newWindows,\@new);
    $printNum++;
    }

   
 print "End p-value calculation: ",get_time(),"\n\n";

#print "This is the TA count: $count\nTotal genome size is: $total\n\n";

#CALCULATE THE ABSOLUTE DIFFERENCE BETWEEN REGION'S MEAN FITNESS AND AVERAGE MEAN FITNESS
### For all windows, add a field that has the difference between that window's mean fitness and the
#average mean fitness for all of the windows

my @allFits = map $_->[7], @newWindows;
my $meanFit=mean(@allFits);
print "Average fitness for all windows: ",$meanFit,"\n";
my @expWindows=();
foreach (@newWindows){
    my @entry=@{$_};
    my $mean=$entry[7];
    my $absdev=sprintf("$round",$mean-$meanFit);
    push (@entry,$absdev);
    push @expWindows,\@entry;
}

@newWindows=@expWindows;

#MAKE OUTPUT CSV FILE WITH ESSENTIAL WINDOW CALCULATIONS

print "Start csv ouput file creation: ",get_time(),"\n";

open CSV, '>', "$outdir/slidingWindows.csv" or die "Cannot open $outdir/slidingWindows.csv";
	
my @csvHeader= ("start", "end","mutants","insertions","TA_sites","ratio","p-value","average", "variance","stdev","stderr","genes","fit-mean");
my $header=join(",",@csvHeader);
print CSV $header,"\n";
foreach my $winLine(@newWindows){
    my $string=join(",",@{$winLine});
    print CSV $string, "\n";
}
close CSV;
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
    }
    close eWIG;
    print "End wig file creation: ",get_time(),"\n\n";
    print "If this wig file needs to be converted to a Big Wig, then use USCS program wigToBigWig in terminal: \n \t./wigToBigWig gview/12G.wig organism.txt BigWig/output.bw \n\n";
}


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

#IF MAKING A TEXT FILE OF GROUPED CONSECUTIVE WINDOWS WITH SAME FITNESS


#Edit: what is the 1000?....probably should be a variable

    my $etxtg="groupedWindows.txt";
    print "Start grouped txt file creation time: ",get_time(),"\n";
    open my $eTXTg, '>', "$outdir/$etxtg" or die $!;
    for my $line(@newWindows){
        my @field=@$line;
        if (!@ecummulative){
            @ecummulative=@field;
            if ($ecummulative[4]>1000){$ecummulative[4]=1000}
        }
        else{
            if ($ecummulative[3]==$field[3]){
                $ecummulative[2]=$field[2];
                if ($ecummulative[4]<=1000-$field[4]){
                    $ecummulative[4]+=$field[4];
                }
            }
            else{
                print $eTXTg ($_,"\t") for @ecummulative;
                print $eTXTg ("\n");
                @ecummulative=@field;
            }
        }
    }
    close $eTXTg;
    print "End grouped text file creation: ",get_time(),"\n\n";
    
    print "\nTo make a BigBed file from this text file, rename file to .bed and use USCS program bedToBigBed in terminal \n\t\n";



#OUTPUT REGULAR SLIDING INFORMATION WITHOUT ESSENTIALS CALCULATIONS (P-VALUES)


#MAKE OUTPUT CSV FILE WITH FITNESS WINDOW CALCULATIONS
	my $fcsv="fitWindows.csv";
    print "Start csv ouput file creation: ",get_time(),"\n";
    my $fcsvBIG = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  # open in append mode
    open my $file, ">", "$outdir/$fcsv" or die "Failed to open file";
    $fcsvBIG->print($file, [ "start", "end","fitness","mutant_count" ]); #header
    foreach my $winLine(@allWindows){
        $fcsvBIG->print($file,$winLine);
    }
    close $file;
    print "End csv ouput file creation: ",get_time(),"\n\n";

#MAKE WIG FILE---->later make BW--->IGV
    print "Start wig file creation: ",get_time(),"\n";
    my $fin = Bio::SeqIO->new(-file=>$ref_genome);
    $refseq = $fin->next_seq;
    $refname = $refseq->id;
    my $fwig="fitWindows.wig";
    open fWIG, ">$outdir/$fwig";
    print fWIG "track type=wiggle_0 name=$fwig\n";
    print fWIG "variableStep chrom=$refname\n";
    foreach my $wigLine(@allWindows){
        my @wigFields=@$wigLine;
        my $pos=$wigFields[0];
        my $fit=$wigFields[7];
        print fWIG $pos," ",$fit,"\n";
    }
    close fWIG;
    print "End wig file creation: ",get_time(),"\n\n";
    print "If this wig file needs to be converted to a Big Wig,\n";
    print " then use USCS program wigToBigWig in terminal: \n";
    print "\t./wigToBigWig gview/12G.wig organism.txt BigWig/output.bw \n\n";


#GOING TO MAKE A TEXT FILE FOR BED CONVERSION TO BIGBED, NEED CHROM # IN COLUMN 0
my @cummulative;
if ($txtg or $txt){
    for my $line(@allWindows){
        unshift($line, $refName);
    }
}
#MAKING A REGULAR TEXT FILE fields: [chrom,start,end,fitness,count]

print "Start text file creation time: ",get_time(),"\n";
open my $TXT, '>', "$outdir/fitWindows.txt" or die $!;
print $TXT ("start","end","W","mutant_count","insertion_count\n");
print $TXT (join("\t",@$_),"\n") for @allWindows;
close $TXT;
print "End text file creation: ",get_time(),"\n\n";


