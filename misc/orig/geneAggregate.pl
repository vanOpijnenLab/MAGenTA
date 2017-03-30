#!/usr/bin/perl -w

#This script (aggregate_mla.pl) is a modified version of aggregate.pl. It includes the following changes:
#1 Differentiated between boolean and argument flags. So now random argument isn't needed for -w
#2 Since we always run the following parameters for flags, made these the defaults
    # -w 1 -x 10 -l 50 -b 0
#3 Added directory input so individual files don't have to be specified in the command line

#To do: 1) fix the maxLength option

use strict;
use Getopt::Std;
use Bio::SeqIO;
use Set::Scalar;

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
use Data::Dumper;

#Example on the command line: ../Blueberries/geneAggregate.pl -i ../9-daptomycin/data/T4Dapto/results/ -o testAgg_160607.csv -w -l 50 -x 10 -g ../0-genome/NC_003028b2.gbk -f ../0-genome/NC_003028b2.fasta -n ../9-daptomycin/nullDist/nullDist_3000t4dapto.csv

sub print_usage() {

	print "\n\n***********************************************\n";

    print "\nUsage: geneAggregate.pl -d inputDirectory/ -o outfileName.csv 
    		-g genbankFile.gbk -f fastaFile.fasta -n nullLibrary.csv -w\n\n";
    print "Description: Finds average fitness and insertion representation over\n";
    print "\t\t an annotated gene.\n\n";
    print "***********************************************\n";

    print "Required:\n\n";
    print " -i\tDirectory containing input files. Make sure / is included after name\n";
    print " -g\tGenbank reference file\n";
    print "   \twith a list of genes seperated by newlines.\n";
    print " -f\tFasta file for the reference genome\n";
    print " -n\tNull distribution library created by the makeNullDist.pl tool\n";
    print "   \t(default: make library--takes a long time)\n";
    	#Currently the default for null distribution doesn't work efficiently
    
    print "Optional:\n";
    print " -o\tOutput csv file for gene aggregate data 
    		(default: geneAggOutput.csv)\n";
    print " -m\tLine separated file that has genes which should be marked\n";
    print " -x\tCutoff: Don't include fitness scores with average counts 
    		(T1+T2)/2 < x\n";
    print "   \t(default: 10)\n";
	print " -b\tBlanks: Exclude -b % of blank fitness scores (write as decimal)
			(where c2 = 0)\n";
    print "   \t(default: 0 = 0%)\n";
    print " -w\tUse weighted algorithm to calculate averages, variance, sd, se\n";
    print " -l\tWeight ceiling: maximum value to use as a weight (default: 50)\n";
    print " -r\tRound to this number of decimals (default: 3 decimal places)\n";
    print " -p\tPrint out usage and all option flags without running program\n";
    print " -N\tSize of null distributions (default=10000)\n";

    print "***********************************************\n\n";
}

#Global options
my %options;
getopts('i:g:f:n:o:m:x:b:wl:r:N:s:',\%options);
my $indir=$options{i};
my $ref = $options{g};
my $fastaFile=$options{f};
my $library=$options{n};
my $out = $options{o} || "geneAggOutput.csv";
my $marked = $options{m};
my $cutoff = $options{x} || 10;
my $blank_pc = $options{b} || 1;
my $weight_ceiling = $options{l} || 50;
my $round=$options{r} || '%.3f';
my $print=$options {p};
my $maxLength=$options{a} ||2000;
my $N=$options{N} || 10000; #Size of the null distribution
my $s=$options{s} || "TA"; #dinucleotide for insertion location


#If any of the required flags are omitted, then print usage and terminate
if ($print || !$ref || !$indir || !$fastaFile || !$library ) { 
	print_usage(); 
	exit; 
}

#Print out all input and output files. Helps with easy troubleshooting.
print "\n*** Input and output files***\n";
print "Input directory: \t\t$indir\n";
print "Outfile: \t\t$out\n";
print "Reference genome: \t\t$ref\n";
print "Fasta file: \t\t$fastaFile\n";
print "Distribution Library: \t\t$library\n";
print "Cutoff: \t\t$cutoff\n";
my $percent=$blank_pc*100;
print "Percentage of blanks to remove: \t\t",$percent,"%\n";



#Read in the fasta file so it's stored as a string
my $seqio = Bio::SeqIO->new(-file => $fastaFile, '-format' => 'Fasta');
my $fasta;
while(my $seq = $seqio->next_seq) {
	$fasta = $seq->seq;
}

#Read in input files that are in the input directory
my @files;
if ($indir){
    my $directory=$indir;
    print "Input directory: ", $directory,"\n";

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

# Returns mean, variance, sd, se
sub average {
    my $scoreref = shift @_;
    my @scores = @{$scoreref};
    my ($variance,$sd,$se);

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
    if ($num>1){
   		$variance = (1/($num-1)) * $xminusxbars;
   		$sd = sqrt($variance);
   		$se = $sd / sqrt($num);
   	 }
   	 else{
   	 	$variance=0;
   	 	$sd=0;
   	 	$se=0;
   	 	}

    return ($average, $variance, $sd, $se);

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
    if ($sum == 0) { return 0; } 
    my $scor = join (' ', @scores);
    my $wght = join (' ', @weights);

    # Calculate weighted average for fitness
    my ($top, $bottom) = (0,0);
    for (my $i=0; $i<@weights; $i++) {
        $top += $weights[$i] * $scores[$i];
        $bottom += $weights[$i];
    }
    $weighted_average = $top/$bottom;

    ($top, $bottom) = (0,0);
    
    # Calculate weighted sample variance and standard deviation and error
    for (my $i=0; $i<@weights; $i++) {
       $top += ( $weights[$i] * ($scores[$i] - $weighted_average)**2);
       $bottom += $weights[$i];
    }
    $weighted_variance = $top/$bottom;
    my $weighted_stdev = sqrt($weighted_variance);
    my $weighted_stder = $weighted_stdev / sqrt(@scores); 

    return ($weighted_average, $weighted_variance, $weighted_stdev, $weighted_stder);
}

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

#Keep track of ALL insertion (positions) for manual null distributions
my @insertPos; 

#Read in all input data files and put insertions into gene keys for data hash
my %data;
foreach my $filename (@files) {
   print $filename,"\n";
   open (IN,'<', $filename);
   my $headerLine=<IN>; #read column names and store in dummy variable
   while (my $line = <IN>) {
      chomp $line;
      my @lines = split(/,/,$line);
      my $locus = $lines[9]; #gene id (SP_0000)
      my $w = $lines[12];    #nW
      my $pos = $lines[0];
      push (@insertPos,$pos);
	  #if blank then w=0
      my $c1 = $lines[2];   #Count at T1
      my $c2 = $lines[3];   #Count at T2
      my $avg = ($c1+$c2)/2; 
       
      # Skip insertions that don't meet the cutoff for reads
      if ($avg < $cutoff) { next; } 
      if ($avg >= $weight_ceiling) { $avg = $weight_ceiling; } 
	  
	  #Start empty arrays if first insertion for this gene locus
	  my @empty;
      if (!$data{$locus}) {
        $data{$locus}{w} = [@empty];
        $data{$locus}{s} = [@empty];
        $data{$locus}{p} = [@empty];
      }
      #Otherwise, add to data already present
      
      # List of fitness scores
      $data{$locus}{w} = [@{$data{$locus}{w}}, $w];

      # List of counts used to generate those fitness scores.
      $data{$locus}{s} = [@{$data{$locus}{s}}, $avg]; 
      
      #List of insertion positions in the gene
      $data{$locus}{p} = [@{$data{$locus}{p}}, $pos];
      #later can get UNIQUE insertion positions and total of that mutant
   }
   
   ######TESTING

		#####TESTING
   close IN;
}
my $same=0;
my $diff=0;
foreach my $locus (keys %data){
	my @pos =  @{$data{$locus}{p}};
	my $first=scalar @pos;
	#print "\n Unfiltered: $locus\t", scalar @pos,"\n";
	@pos = uniq sort {$a<=>$b }@{$data{$locus}{p}};
	my $sec=scalar @pos;
	#print "\nUnique: $locus\t", scalar @pos,"\n";

	if ($first==$sec){
	$same++;
	}
	else{$diff++;}
	
	}
print "\n\nsame: $same\tdiff: $diff\n\n";
#TESTINGGGGGG
#Get unique insertion positions for each gene in the %data hash


	

#Set up the null distribution library to determine insertion representation
my @distLib;
my $FILE3 = "nullDist.txt";
unless(open DIST, ">", $FILE3){
    die "\nUnable to create $FILE3:\n$!";
}
printf DIST "Sites\tSize(n)\tMin\tMax\tMean\tstdev\tvar\tMinPval\n";

#Loop once for each distribution in the library
	#(i.e. distribution of 10,000 sites each with 35 TA sites, 
	#then a distribution of 10,000 sites each with 36 TA sites, etc)

my %distLibrary;
my @selector;
my @nullDist;
my @unsorted;
my $sum;
my $insertions;

if ($library){
	#If a library was already made using the MakeNullDist tool then read it in
	open LIB,'<',$library;
	while (<LIB>){
		chomp $_;
		my @dist=split (",",$_);
		my $key=$dist[0];
		shift @dist;
		@dist= sort { $a <=> $b } @dist;
		$distLibrary{$key}=\@dist;

	}	
	close LIB;
	
	#Print stats on each array in the distribution library
	foreach (sort {$a <=> $b} keys %distLibrary){
	    my $sites=$_;
		@nullDist=@{$distLibrary{$sites}};
		my $size=scalar @nullDist;
		#my ($nullMean, $var, $stdev, $sterr)=average(@nullDist);
		my $counter=scalar @nullDist;
		my $nullMean=sprintf("$round",mean(@nullDist));
		my $stdev =sprintf("$round",stdev(\@nullDist, $nullMean));
		my $var=sprintf("$round",($stdev*$stdev));
		my $min=sprintf("$round",$nullDist[0]);
		my $max=sprintf("$round",$nullDist[scalar @nullDist-1]);
		my $minp=pvalue(0,$sites);
		printf DIST "$sites\t$N\t$min\t$max\t$nullMean\t$stdev\t$var\t$minp\n";
	}
}

else{
	#If no input library then make distribution library
	for (my $sitez=1; $sitez<=$maxLength;$sitez++){
		#print "In the first for loop to make a null distribution\n";
	    @unsorted=();
		$insertions=0;
		$sum=0;
	
		for (my $i=1; $i<=$N; $i++){
				my @random_set = rand_set( set => \@selector, size => $sitez);
				my $setMean=mean(@random_set);
				push (@unsorted, $setMean);
			}
		@nullDist= sort { $a <=> $b } @unsorted;	
		$distLibrary{$sitez}=\@nullDist;
		my ($nullMean, $nullVar, $nullsd, $nullse)=average(@nullDist);
		my $counter=scalar @nullDist;
		$nullMean=sprintf("$round",($sum/$insertions));
		my $stdev =sprintf("$round",stdev(\@nullDist, $nullMean));
		my $variance=sprintf("$round",($stdev*$stdev));
		my $min=sprintf("$round",$nullDist[0]);
		my $max=sprintf("$round",$nullDist[scalar @nullDist-1]);
		my $minp=pvalue(0,$sitez);
		printf DIST "$sitez\t$N\t$min\t$max\t$nullMean\t$stdev\t$variance\t$minp\n";
	}
	close DIST;

}

sub pvalue{
    
    #Takes in window count average (countAvg) and number of TAsites 
    #and makes a null distribution to calculate and return the pvalue
    
    my $mean=shift@_; 	#mean= (number_of_insertions)/(number_of_TAsites)
    my $TAsites=shift@_;
    my @specDist= @{$distLibrary{$TAsites}};
    my @unsortedtemp=@specDist;
    push (@unsortedtemp, $mean);
    my @temp=sort {$a <=> $b} @unsortedtemp;
    my $rank= binsearch_pos {$a <=> $b} $mean,@temp;
    my $i=$rank;
    while ($i<scalar(@temp) and $temp[$i]==$temp[$rank]){
        $i++;
    }
    $rank=$i;
    my $pval=$rank/$N;  
    return $pval;     
}

#If a list of genes to mark was specified then read it in
my @marked;
my $marked_set = Set::Scalar->new;
if ($marked) {
   open MARKED, $marked;
   while (<MARKED>) {
      chomp $_;
      $marked_set->insert($_);
   }
   close MARKED;
}

#Store all final gene summary information in a hash
my %outHash;

#Get fasta file as string so we can look up the number of TA sites given genomic coordinates
my $in = Bio::SeqIO->new(-file=>$ref); 
my $refseq = $in->next_seq;
my @features = $refseq->get_SeqFeatures;

my @header=("locus","start","stop","length","strand","gene", "countTA", 
			"insertions","essPval","avgFit","variance","stdev","stderr",
			"total","blank_ws","removed","to_remove");
if ($marked){
	push (@header,"marked");
	}

#Go through all features (gene, rna, etc) in the genbank file
foreach my $feature (@features) {
    
    # Check if the feature is a gene
    if ($feature->primary_tag eq 'gene') {
		my $start=$feature->start;
		my $stop=$feature->end;
		my $strand=$feature->strand;
		my $length=$stop-$start;
		my @locus = $feature->get_tagset_values('locus_tag');
		my $locus = $locus[0] || '';
		my @gene = $feature->get_tagset_values('gene');
		my $gene = $gene[0] || '';
		
		#Get the number of TA sites in the gene
		my $seq=substr ($fasta, $start, $stop - $start); 
		my @ctemp = $seq =~ /$s/g; 
		my $countTA = @ctemp;
         
		#Check if we have tn-seq fitness information about that gene and 
		#Calculate average fitness, var, stdev, and sterr.
	
		#Array to hold all info for a gene that will be printed to out file
        my @summary;
		
		#For every gene summary, we need to know the locus, start, stop,
			#length, strand, geneID, countTA, insertions, essPval, avgFit, 
			#variance, stdev, stderr, blank_ws, removed
			
		my ($average, $variance, $stdev, $sterr,$removed,$total,$insertions,$to_remove) =
			("NA","NA","NA","NA",0,0,0,0);
		my ($sum,$nonblank,$avgsum,$blank_ws)=(0,0,0,0,0);
		
		if ($data{$locus}) {
			#Go through all of the fitness values that fall in the gene locus
		    #and count the blank fitness scores (w=0 because no reads at t1/t2)
			foreach my $w (@{$data{$locus}{w}}) {
				if ($w==0) {
					$blank_ws++; #blank_ws = blank fitness score
				}
			  	else {
			  		$sum += $w; 
			  		$nonblank++;
			  	}
			  	$total++;
			}
			
			#Remove blanks from scores if the option flag was specified
			#blank_pc = Percentage of blanks to remove as -b flag 
		   	if ($blank_pc !=1){
				$to_remove = int($blank_pc * $total);
           		if ($blank_ws > 0) {
           			my @wscores=@{$data{$locus}{w}};
              		for (my $i=0; $i < @wscores; $i++) {
                		if ($removed > $to_remove) {last}
                		my $w = $wscores[$i];
                		if ($w==0) {
                    		$removed++;
                    		splice( @{$data{$locus}{w}}, $i, 1);
                    		splice( @{$data{$locus}{s}}, $i, 1);
                    		$i-=1;
                		}
              		}
				}
			}
		
			#Get stats on gene: average fitness, var, stdev, sterr, etc.
			#Can't do stats if no observations (insertions)
			if ($nonblank == 0 ) {
				($average, $variance, $stdev, $sterr) = ("NA","NA","NA","NA");
			}
			else {	
			   #Do weighted or unweighted average depending on -w flag
			   if ($options{w}) {
				   ($average, $variance, $stdev, $sterr) = 
				   	weighted_average ($data{$locus}{w},$data{$locus}{s});
			   }
			   else {
				   ($average, $variance, $stdev, $sterr) = 
				   average($data{$locus}{w});
			   }
			}
		} 
		#End of routine for when we have Tn-Seq data (insertions) for gene
		
		#Insertion representation assessment
		my ($avgInsert,$pval);
		#How many UNIQUE insertion mutations are in the gene?
		my @unique=uniq @{$data{$locus}{p}};
		$insertions=scalar @unique;
		#If there are TA sites in the gene
		if ($countTA >0){
			$avgInsert=sprintf("$round",$insertions/$countTA);
			$pval=pvalue($avgInsert,$countTA);
		}	
		else{
			$avgInsert=0; 
			$pval="NA";
		}
		
		
			
		#Store all of these variables into a summary for the gene
		@summary = ($locus,$start,$stop,$length,$strand,$gene,$countTA,
			   		$insertions,$pval,$average,$variance,$stdev,$sterr,
			   		$total,$blank_ws,$removed,$to_remove); 
			   		
		#If this gene is in the list of genes to be marked, then mark it
		if ($marked && $marked_set->contains($locus)) {
			   push(@summary,"M");
		}
		
		$outHash{$locus}=\@summary;

	} 
	#end of routine for features that are genes
}

open OUT, ">",$out;
print OUT join(",",@header),"\n";
foreach my $key (sort keys %outHash){
	my @vals=@{$outHash{$key}};
	print OUT join(",",@vals);
	print OUT "\n";
	}
close OUT;
