#!/usr/bin/perl -w

#Margaret Antonio updated 16.10.03

use strict;
use Getopt::Long;
use warnings;
use Bio::SeqIO;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
	print "\n####################################################################\n";
    print "singleVal: creates a wig file with position and insertion count OR fitness\n";
    print "\nDESCRIPTION: ";
    print "Integrates multiple files of transposon insertion data and outputs\n";
    print "aggregate fitness within a sliding window (specified by size and step).\n";
    
    print "\nUSAGE:\n";
    print "perl singleVal.pl <OPTIONS> <REQ OUTPUT TYPE(S)> <INPUT FILE(S)>\n\n";
    
    print "\nREQUIRED:\n";
    print " -d\tDirectory containing all input files (files from\n";
    print "   \taggregate script)\n";
    print "   \tOR\n";
    print "   \tIn the command line (without a flag), input the name(s) of\n";
    print "   \tfiles containing gene fitness values (output of calcFit). \n\n";
    print " -x\tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 10)\n";
    print " -b\tBlanks: Exclude -b % of blank fitness scores (scores where c2 = 0) (default: 0 = 0%)\n";
    print " -w\tUse weighted algorithm to calculate averages, variance, sd, se\n";
    print " -l\tWeight ceiling: maximum value to use as a weight (default: 50)\n";
    
    print "OPTIONAL:\n";
    print " -h\tPrints usage and exits program\n";
    print " -o\tOutput file for comparison data. Default: singleVal.wig\n";
    print " -v\tString value for output: 'fit' for fitness OR 'count' for count\n";
    print "   \tDefault: fit for fitness\n";
    print " -n\tName of the reference genome, to be included in the wig header\n";
    print "   \tDefault: genome\n";

    print " \n~~~~Always check that file paths are correctly specified~~~~\n";
    print "\n##################################################################\n";
    
}


#ASSIGN INPUTS TO VARIABLES
our ($help,$ref_genome,$indir,$out,$cutoff,$blank_pc,$weight_ceiling);
GetOptions(
'h'=>\$help,
'o:s' =>\$out,
'c:i' =>\$cutoff,
'b:i' =>\$blank_pc,
'l:i' =>\$weight_ceiling,
);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
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

#ADDED BY MLA allows input to be directory---good for inputting L1-L6
my @files=@ARGV;


#SET DEFAULTS
if (!$cutoff){$cutoff=10}
if (!$blank_pc){$blank_pc=0}
if (!$weight_ceiling){$weight_ceiling=50}
if (!$out){$out="singleVal.csv"}

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
    my $variance = (1/($num-1)) * $xminusxbars;
    my $sd = sqrt($variance);
    my $se = $sd / sqrt($num);
    
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
    if ($sum == 0) { return 0; } # No scores given?
    
    my $scor = join (' ', @scores);
    my $wght = join (' ', @weights);
    
    # Now calculated weighted average.
    my ($top, $bottom) = (0,0);
    for (my $i=0; $i<@weights; $i++) {
        $top += $weights[$i] * $scores[$i];
        $bottom += $weights[$i];
    }
    $weighted_average = $top/$bottom;
    #print "WA: $weighted_average\n";
    
    ($top, $bottom) = (0,0);
    # Now calculate weighted sample variance.
    for (my $i=0; $i<@weights; $i++) {
        $top += ( $weights[$i] * ($scores[$i] - $weighted_average)**2);
        $bottom += $weights[$i];
    }
    $weighted_variance = $top/$bottom;
    
    my $weighted_stdev = sqrt($weighted_variance);
    my $weighted_stder = $weighted_stdev / sqrt(@scores);  # / length scores.
    
    return ($weighted_average, $weighted_variance, $weighted_stdev, $weighted_stder);
}

my %pos_summary;


foreach my $filename (@files) {
    print "\n",$filename,"\n";
    open IN, $filename;
    my %hash;
    while (my $line = <IN>) {
        chomp $line;
        my @lines = split(/,/,$line);
        my $pos = $lines[0];
        my $w = $lines[12];
        if ($w and $w eq 'nW') {next;}
        if (!$w) { $w = 0 } # For blanks
        my $c1 = $lines[2];
        my $c2 = $lines[3];
        my $avg = ($c1+$c2)/2; # Later: change which function to use? C1? AVG(C1+C2)?
        if ($avg < $cutoff) { next; } # Skip cutoff genes.
        if ($avg >= $weight_ceiling) { $avg = $weight_ceiling; } # Maximum weight.
        
        my @empty;
        if (!$pos_summary{$pos}) {
            $pos_summary{$pos}{w} = [@empty];
            $pos_summary{$pos}{s} = [@empty];
        }
        $pos_summary{$pos}{w} = [@{$pos_summary{$pos}{w}}, $w];  # List of Fitness scores.
        $pos_summary{$pos}{s} = [@{$pos_summary{$pos}{s}}, $avg]; # List of counts used to generate those fitness scores.
    }
    close IN;
    
}

open SUMMARY, ">",$out;


    print SUMMARY "pos,fitness,ins_count,fitness_sd,fitness_se\n";
    # Now print out summary stats.
    foreach my $key (sort {$a<=>$b} keys %pos_summary) {
        if (!$key) {next}
        my $sum=0;
        my $num=0;
        my $avgsum = 0;
        
        # Get the average.
        foreach my $w (@{$pos_summary{$key}{w}}) {
            $sum += $w;
            $num++;
        }
        my $average= $sum/$num;
        my $xminusxbars = 0;
        
        # Get the standard deviation.
        foreach my $w (@{$pos_summary{$key}{w}}) {
            $xminusxbars += ($w-$average)**2;
        }
        my ($sd, $se) = ('','');
        if ($num > 1) {
            $sd = sqrt($xminusxbars/($num-1));
            $se = $sd / sqrt($num);
        }
        
        print SUMMARY "$key,$average,$num,$sd,$se\n";

    }


close SUMMARY;
