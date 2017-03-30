#!/usr/bin/perl -w
# Count and compare transcripts at each position in a reference genome.
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::GenBank;
use Getopt::Long;
use List::Util 'shuffle';

sub print_usage() {
    
    print "Usage: ./ctrans.pl <options>\n\n";
    
    print "Required\n";
    print "--ref\t\tThe name of the reference genome file, in GenBank format.\n";
    print "--t1\t\tThe name of the bowtie mapfile from time 1.\n";
    print "--t2\t\tThe name of the bowtie mapfile from time 2.\n";
    print "--out\t\tName of a file to enter the .csv output.\n";
    print "\n";
    
    print "Optional\n";
    print "--expansion\t\tExpansion factor (default: 250)\n";
    print "--reads1\t\tThe number of reads to be used to calculate the correction factor for time 0.\n\t\t(default counted from bowtie output)\n";
    print "--reads2\t\tThe number of reads to be used to calculate the correction factor for time 6.\n\t\t(default counted from bowtie output)\n";
    print "--cutoff\t\tDiscard any positions where the average of counted transcripts at time 0 and time 6 is below this number (default 0)\n";
    print "--strand\t\tUse only the specified strand (+ or -) when counting transcripts (default: both)\n";
    print "--normalize\t\tA file that contains a list of genes that should have a fitness of 1\n";
    print "--multiply\t\tMultiply all fitness scores by a certain value (e.g., the fitness of a knockout). You should normalize the data.\n";
    print "--ef\n";
    print "--el\n";
    print "--wig\t\tCreate a wiggle file for viewing in a genome browser. Provide a filename.\n";
    print "\n";
    
}

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}

# Takes two parameters, both hashrefs to lists.
# 1) hashref to list of scores
# 2) hashref to list of weights, equal in length to the scores.
sub weighted_average {
   
   my $weightref = shift @_;
   my $scoreref = shift @_;
   my @scores = @{$scoreref};
   my @weights = @{$weightref};
   
   my $sum=0;
   my $weighted_average=0;

   #print "Scores: @scores\n";
   #print "Weights: @weights\n\n";
   
   # Go through once to get total
   foreach my $w (@weights) { $sum += $w; }
   
   if ($sum == 0) { return 0; } # No scores given?
   # Now calculated weighted average.
    for (my $i=0; $i<@weights; $i++) {
        $weighted_average += ( ($weights[$i]/$sum) * $scores[$i]);
    }
    return $weighted_average;
}

# GET ALL OPTIONS.
our ($ref_genome, $mapfile1, $mapfile2, $outfile, $reads1, $reads2, $usestrand, $wig,
    $use_averaging, $normalize, $expansion_factor, $cutoff, $multiply, $exclude_first, $exclude_last);
GetOptions('ref:s' => \$ref_genome,
        't1:s'  => \$mapfile1,
        't2:s'  => \$mapfile2,
        'out:s'  => \$outfile,
        'expansion:i' => \$expansion_factor,
        'reads1:i'   => \$reads1,
        'reads2:i'   => \$reads2,
        'cutoff:i'   => \$cutoff,
        'strand:s'   => \$usestrand,
        'normalize:s' => \$normalize,
        'multiply:f' => \$multiply,
        'ef:f' => \$exclude_first,
        'el:f' => \$exclude_last,
        'wig:s' => \$wig
);

# Check required vars and set defaults.
if (!$ref_genome or !$mapfile1 or !$mapfile2 or !$outfile) {
    print_usage();
    die; }
if (!$expansion_factor) { $expansion_factor = 250 }
if (!$cutoff) { $cutoff = 0 };

# Start!
print STDOUT "Starting: ", get_time(), "\n";

# Pull in reference genome and feature list.
# create array of features for every feature in the genome.
my $in = Bio::SeqIO->new(-file=>$ref_genome);
my $refseq = $in->next_seq;
my $refname = $refseq->id;
my @features = $refseq->get_SeqFeatures;
# Hashify features for easy lookup
my @array_of_features;
foreach my $feature (@features) {
    # Check if it's a gene.
    if ($feature->primary_tag eq 'gene') {
        
        # Get gene name
        my @values = $feature->get_tagset_values('locus_tag');
        my $gene = $values[0];
        my $start = $feature->start;
        my $end = $feature->end;
        
        if ($exclude_first) {$start = int($start + (($end-$start) * $exclude_first))}
        if ($exclude_last)  {$end   = int($end   - (($end-$start) * $exclude_last))}
        
        #print $feature->start, " ", $feature->end, "\n";
        
        my %hash = ( gene => $gene,
                     start => $start,
                     end => $end,
                     strand => $feature->strand );
        push @array_of_features, \%hash;
        
    } 
}
print STDOUT "Done generating feature lookup: ", get_time(), "\n";

# Read a MAQ mapview file.
my (@reads1, @reads2);
open IN, $mapfile1; {
    while (my $line = <IN>) {
        push @reads1, $line;
    }
} close IN;
open IN, $mapfile2; {
    while (my $line = <IN>) {
        push @reads2, $line;
    }
} close IN;

sub read_mapfile {
    my $reads = shift @_;
    my $usestrand = shift @_; # Optional
    
    my ($plustotal, $minustotal) = (0,0);
    my (%plus_counts, %minus_counts);
    
    my @reads = @{$reads};
    foreach my $line (@reads) {
        
        my @lines = split(/\t/,$line);
    
        # All we're interested in is the position that the line occurs.
        my @counts = split(/-/, $lines[0]);
        my $count = $counts[1];
        my $strand = $lines[1];
        my $pos = $lines[3];
        
        # Skip this if we only want one strand (OPTIONAL)
        if ($usestrand  && ($strand ne $usestrand)) {
            next;
        }
        
        # Store in different hash depending on strand.
        my $hashref;
        if ($strand eq '+') {
            my $length = length ($lines[4]);
            #print "$length\n";
            $pos+=($length-2); # TA site
            $plustotal+= $count;
            $hashref = \%plus_counts;
        }
        else {
            $minustotal+= $count;
            $hashref = \%minus_counts;
        }
        
        # Increase the count at that position.
        if ($$hashref{sites}) { $$hashref{sites}++ }
        else {$$hashref{sites}=1}
        
        if (!$$hashref{$pos}) {
            $$hashref{$pos} = $count;
        }else {
            $$hashref{$pos} += $count;
        }
        
        
    }
    close IN;
    
    $plus_counts{total} = $plustotal;
    $minus_counts{total} = $minustotal;
    return (\%plus_counts, \%minus_counts);
}

# Next - count the data from the mapfiles.
# Read mapfiles into counts.
my $r1ref = \@reads1;
my $r2ref = \@reads2;
my ($plus_ref_1, $minus_ref_1) = read_mapfile($r1ref, $usestrand);
print STDOUT "Read first file: ", get_time(), "\n";
my ($plus_ref_2, $minus_ref_2) = read_mapfile($r2ref, $usestrand);
print STDOUT "Read second file: ", get_time(), "\n";

my %plus_counts_1 = %{$plus_ref_1};
my %plus_counts_2 = %{$plus_ref_2};
my %minus_counts_1 = %{$minus_ref_1};
my %minus_counts_2 = %{$minus_ref_2};

my $total;

# APPLY CORRECTION FACTOR
print "Reads:\n";
print "1: + ", $plus_counts_1{'total'}, " - ", $minus_counts_1{'total'}, "\n";
print "2: + ", $plus_counts_2{'total'}, " - ", $minus_counts_2{'total'}, "\n";

print "Sites:\n";
print "1: + ", $plus_counts_1{sites}, " - ", $minus_counts_1{sites}, "\n";
print "2: + ", $plus_counts_2{sites}, " - ", $minus_counts_2{sites}, "\n";

my ($cfactor1, $cfactor2);
if (!($reads1 && $reads2)) { # Use total number of reads from original output if provided.
    $reads1 = $plus_counts_1{'total'} + $minus_counts_1{'total'};
    $reads2 = $plus_counts_2{'total'} + $minus_counts_2{'total'};
}
$total = ($reads1 + $reads2)/2;
$cfactor1 = $reads1 / $total;
$cfactor2 = $reads2 / $total;
print "Cfactor1: $cfactor1\n";
print "Cfactor2: $cfactor2\n";
printf "Using given reads to get correction factors: %.03f : %.03f\n", $cfactor1, $cfactor2; 

# Header line
open OUT, ">$outfile" or die "Could not open outfile!\n";
print OUT "position,strand,count_1,count_2,ratio,mt_freq_t1,mt_freq_t2,pop_freq_t1,pop_freq_t2,gene,D,W,nW\n";

my ($genic, $total_inserts) = (0,0);
# Now go through the counts and compare.
for(my $i=0; $i<3000000; $i++) {
    my ($c1, $c2);
    my $strand;
    
    # FIRST SET, TIME 0
    if ( ($plus_counts_1{$i}) ) {
        
        # We have reads. Combine.
        $c1 = $plus_counts_1{$i};       # + singles along
        $strand = "+/";
        if ($minus_counts_1{$i}) {
            $c1 = ($plus_counts_1{$i} + $minus_counts_1{$i});  # +,- combined, average
            $strand = "b/";
        }

    }elsif ( $minus_counts_1{$i} and !$plus_counts_1{$i}) { # - singles
        $c1 = $minus_counts_1{$i};
        $strand = "-/";
    }else {
        next;   
    }
    
    # SECOND SET, TIME 1
    if ( ($plus_counts_2{$i}) ) {
        
        # We have reads. Combine.
        $c2 = $plus_counts_2{$i};       # + singles along
        
        if ($minus_counts_2{$i}) {
            $c2 = ($plus_counts_2{$i} + $minus_counts_2{$i});  # +,- combined
            $strand .= "b";
        }else {
            $strand .= "+";
        }

    }elsif ( $minus_counts_2{$i} and !$plus_counts_2{$i}) { # - singles
        $c2 = $minus_counts_2{$i};
        $strand .= "-";
    }
    
    # CUTOFF - before normalization
    # if (($c1 + $c2)/2 < $cutoff) { next }
    
    $c1 /= $cfactor1;
        
    my $ratio;
    if ($c2) {
        $c2 /= $cfactor2;
        $ratio = $c2 / $c1;
    }
    else {
        $c2 = 0;
        $ratio = 0;
    }
    
    # CUTOFF - after normalization
    if (($c1 + $c2)/2 < $cutoff) { next }
    
    # Get mutant frequency at each time point.
    my $mt_freq_t1 = $c1 / $total; #($reads1 * $cfactor1); #$plus_counts_1{'total'};
    my $mt_freq_t2 = $c2 / $total; #($reads1 * $cfactor1); #$plus_counts_1{'total'};
        
    my $pop_freq_t1 = 1 - $mt_freq_t1;
    my $pop_freq_t2 = 1 - $mt_freq_t2;
    
    # Calculate fitness.
    # W =
    my $w = '';
    if ($mt_freq_t2 != 0) {
        #print "Equation: ln ($mt_freq_t1 * ( $expansion_factor / $mt_freq_t2) )\n";
        my $top_w = log($mt_freq_t2 * ($expansion_factor/$mt_freq_t1));
        my $bot_w = log($pop_freq_t2 * ($expansion_factor/$pop_freq_t1));
        $w = $top_w / $bot_w;
    }
    
    my $gene = '';
        
    # Pull features
    foreach my $feature (@array_of_features) {
        
        my %hash = %{$feature};
        #print $hash{'start'};
        #print "\n";
        if ($hash{'start'} <= $i && $hash{'end'} >= $i) {
            $gene = $hash{'gene'};
            $genic++;
            #print "Found: $gene\n";
            last;
        }
        
    }
    $total_inserts++;
    
    # Output everything to a comma separated value file.
    print OUT "$i,$strand,$c1,$c2,$ratio,$mt_freq_t1,$mt_freq_t2,$pop_freq_t1,$pop_freq_t2,$gene,$expansion_factor,$w,$w\n";   

}
close OUT;

print STDOUT "Done comparing mapfiles ", get_time(), "\n";
print "Genic: $genic\n";
print "Total: $total_inserts\n";

######################
# NORMALIZATION STEP
######################

if ($wig) {
    open WIG, ">$wig";
    print WIG "track type=wiggle_0 name=$wig\n";
    print WIG "variableStep chrom=$refname\n";
    
}

# Normalize if filename given.
if ($normalize) {
    
    # Get a list of the genes that should have a fitness of 1.
    open IN, $normalize;
    my $transposon_genes = '';
    while (<IN>) { chomp $_; $transposon_genes .= "$_ "; }
    close IN;
    print "Normalize genes loaded\n";#: $transposon_genes\n";
    
    # Get all of the fitness values for those genes.
    open IN, $outfile;
    my $sum=0;
    my $count=0;
    
    my @weights;
    my @scores;
    my $blank_ws=0;
    
    while (<IN>) {
        # Find all entries for those genes.
        if ($_ =~ /position/) { next; } # Header
        chomp $_;
        my @list = split(/,/, $_);
        
        if ($list[9] ne '' and $transposon_genes =~ $list[9] and $list[11]) {
            # For regular average
            my $c1 = $list[2];
            my $c2 = $list[3];
            my $avg = ($c1+$c2)/2;
            if ($c1 >= 10) { # and $c2 >= 25) {
                if ($avg >= 75) { $avg = 75; } # Maximum weight.
                # Average
                if ($list[11] eq 0) { $blank_ws++;}
                $sum += $list[11];
                $count++;
                # For weighted average
                push @weights, $avg;
                push @scores, $list[11];
                
                print $list[9], " ", $list[11], "\n";
                    
            } # Skip cutoff genes.
        }
    }
    close IN;
    
    # Remove 100% of blanks, but report what percentage of the sets were blank
    my $blank_count = 0;
    my $original_count = @scores;
    for (my $i=0; $i < @scores; $i++) {
        my $w = $scores[$i];
        #print "$w\n";
        if ($w == 0) {  # Remove blank.
            $blank_count++;
            splice( @weights, $i, 1);
            splice( @scores, $i, 1);
            $i-=1;
        }
    }
    my $pc_blank_normals = $blank_count / $original_count;
    print "$blank_count blank out of $original_count: $pc_blank_normals\n";
    
    my $average = $sum/$count;
    my $weighted_average = &weighted_average(\@weights, \@scores);
    
    print "Normalization step:\n";
    print "Regular average: $average\n";
    print "Weighted Average: $weighted_average\n";
    print "Total Insertions: $count\n";   
    
    # Next, divide every fitness score by the average.
    open IN, $outfile;
    open OUT, ">$outfile.normalized";
    my $header = <IN>;
    print OUT $header; # header
    my $old_ws=0;
    my $new_ws=0;
    my $wcount=0;
    while (<IN>) {
        chomp $_;
        my @lines = split(/,/,$_);
        if (!$lines[11]) { $lines[11] = 0; }
        
        #print "Old W: ", $lines[11], "\n";
        my $new_w = $lines[11] / $weighted_average; #$weighted_average; #$average; # # Set new W.
        if ($multiply) {$new_w *= $multiply} # Multiply by given value .. if given.
        
        if ($lines[11]> 0) {
            $old_ws += $lines[11];
            $new_ws += $new_w;
            $wcount++;
        }
        
        $lines[12] = $new_w;
        #print "New W: ", $lines[11], "\n";
        
        my $row = join(',', @lines);
        print OUT "$row\n";
        
        # Create wiggle file.
        if ($wig) {
            my $clean_new_w = sprintf("%.04f", $new_w);
            print WIG $lines[0], " ", $clean_new_w, "\n";
        }
        
    }
    close IN;
    close OUT;
    if ($wig) {close WIG}
    
    my $old_w_mean = $old_ws / $wcount;
    my $new_w_mean = $new_ws / $wcount;
    
    print "Old W Average: $old_w_mean\n";
    print "New W Average: $new_w_mean\n";
    
    # Overwrite old file with normalized file.
    my $cmd = "mv $outfile.normalized $outfile\n";
    exec $cmd;
}

