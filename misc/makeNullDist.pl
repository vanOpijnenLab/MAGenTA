#!/usr/bin/perl -w



#Margaret Antonio Start date: 1 March 2016

#Creates the null distribution used in the essential test

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


#Global options
my ($indir,@files,$max,$fasta,$round);


GetOptions(
'indir|d=s' => \$indir, #Input directory
'input|i=s' => \@files, #Input files
'max|a=i'   => \$max,   #Null distribution max
'fasta|f=s' => \$fasta, #Fasta file to get TA site positions
'round|r=s' => \$round,

);

if (!$round){ $round='%.3f';}

#Get input libraries
if ($indir){
    print "Input directory: ", $indir,"\n";
    
    opendir(DIR, $indir) or die "couldn't open $indir: $!\n";
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

#Now that files are identified, read them in and store the insertions
my @insertPos;
foreach my $filename (@files) {
    print $filename,"\n";
    open (IN,'<', $filename);
    my $headerLine=<IN>; #read the header (column names) line and store in dummy variable
    my %hash;
    while (my $line = <IN>) {
        chomp $line;
        my @lines = split(/,/,$line);
        my $pos = $lines[0];
        push (@insertPos,$pos);
    }
}

#Have array of all positions where an insertion occurred. Each insertion should be represented only
#once and should be sorted

@insertPos = sort { $a <=> $b } @insertPos;
@insertPos = uniq @insertPos;


#Functions for mean and standard deviation (used in null distribution creation and stats)
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

#

my @sites;

#First read fasta file into a string
my $seqio = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
my $prev;
my $total=0;
while(my $seq = $seqio->next_seq) {
    $fasta = $seq->seq;
}


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


#Get number of "TA" sites, regardless of insertion---this is just the fasta file
my $x="TA";
my @c = $fasta =~ /$x/g;
my $countTA = @c;

print "\nNow looking for where insertions occurred\n";

while (($genPos != -1) and ($pos!=scalar @insertPos)) { #as long as the TA site is found
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

# the @allTAsites array is a 2d array with arrays containing (position,0 or 1 depending on insertion)
print "Total number of TA sites: $countTA\n";
my $FILE1 = "allTAsites.txt";
open (ALL_TA, ">", $FILE1);

foreach (@allTAsites){
    foreach (@{$_}){
        print ALL_TA$_," ";
    }
    print ALL_TA"\n";
}

close ALL_TA;


#Save all the unmatched insertion sites
my $FILE2 = "unmatched.txt";
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


#MAKE LIBRARY OF NULL DISTRIBUTIONS:
print "Making a library of null distributions.\n For information on distribution library, see nullDist.txt\n";

my @distLib;

my $file = "nullDist.txt";
open (DIST, '>', $file) or die "\nUnable to create $file:\n$!";
my @header=("sites","size(N)","min","max","mean","stdev","var","minPval"); #######
print join(",",@header);
print "\n";

#Loop once for each distribution in the library
#(i.e. distribution of 10,000 sites each with 35 TA sites, then a distribution of 10,000 sites each with 36 TA sites, etc)


my %distLibrary;

for (my $sitez=1; $sitez<=$max;$sitez++){
    #print "In the first for loop to make a null distribution\n";
    my @unsorted=();
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
    #print "sum: ",$sum, "\t", "count: ",$count,"\n";
    
    my @nullDist= sort { $a <=> $b } @unsorted;
    my $nullMean=sprintf("$round",($sum/$count));
    my $stdev =sprintf("$round",stdev(\@nullDist, $nullMean));
    my $var=sprintf("$round",($stdev*$stdev));
    my $min=sprintf("$round",$nullDist[0]);
    my $max=sprintf("$round",$nullDist[scalar @nullDist-1]);
    my $setScalar=scalar @nullDist;
    
    $distLibrary{$sitez}=\@nullDist;
    my $minp=pvalue(0,$sitez);
    #print $sum,$count;
    #print $nullMean, "\t",$stdev;
    my @printArray=($sitez,$N,$min,$max,$nullMean,$stdev,$var,$minp);
    print $sitez;
    print "\n";
    print DIST join(",",@printArray);
    print DIST "\n";
}
close DIST;
#my $maxDist=scalar @distLib -1;

open LIB,'>',"nullDistLibrary.csv";
foreach my $key(keys %distLibrary){
    print LIB $key,",";
    my @values=@{$distLibrary{$key}};
    print LIB join( ',', @values ),"\n";
}
close LIB;



sub pvalue{
    
    #takes in window count average (countAvg) and number of TAsites and makes a null distribution to calculate the pvalue, which it returns
    my $mean=shift@_;
    my $TAsites=shift@_;;
    my $N=10000;
    my @specDist;
    @specDist=@{$distLibrary{$TAsites}};
    
    my $rank= binsearch_pos { $a cmp $b } $mean,@specDist;
    my $i=$rank;
    while ($i<scalar(@specDist)-1 and $specDist[$i+1]==$specDist[$rank]){
        $i++;
    }
    $rank=$i;
    my $pval=$rank/$N; #calculate pval as rank/N
    return $pval;
    
}



