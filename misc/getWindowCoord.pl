#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: After sliding windows have been filtered and ordered, can be grouped

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use List::Util qw(sum);

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


#USAGE from /8-essFilters
#perl ../Blueberries/grouper.pl --in apples.csv

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR


# 0:start,1:end,2:fitness,3:mutant_count,4:insertions,5:TA_sites,6:ratio,7:p-value
print "Sort by options:\n0:Start coord.\n1:End coord\n2:Fitness\n3: mutant count\n4: insertions\n5:TA_sites\n6:ratio\n7: p-value\n8: deviation from mean fitness\n";

#ASSIGN INPUTS TO VARIABLES
our ($in,$h,$size, $bracket,$step,$defEss,$out,$sortby,$round,$super);
GetOptions(
'i:s' => \$in,
'o:s' =>\$out,
'h' => \$h,
'b'=>\$bracket,
'size'=>\$size,
'step'=>\$step,
'ess'=>\$defEss,
's:i' => \$sortby,
'r:i'=> \$round,
'super'=> \$super,
);

if (!$size) {$size=500}   #set the default sliding window size to 500
if (!$step) {$step=10}
if (!$out) {$out="orderedGroup.csv"}
if (!$round){$round='%.3f'}



sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}

#GET DATA OUT OF FILE AND INTO 2D ARRAY

sub mean {
	return sum(@_)/@_;
}
my @windows;

sub cleaner{
	my $line=$_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}
	
open(DATA, '<', $in) or die "Could not open '$in' \n";
	
my $line=<DATA>;
#print "This is the line: ",$line,"stoooooop";
$line=cleaner($line); #gets rid of carriage returns (^M)
my @header=split(',',$line);

my $tick=0;
while (my $entry = <DATA>) {
	$entry=cleaner($entry);
	my @fields=split(',',$entry);
	push (@windows,\@fields);     
}
close DATA;

#my @finalOutput;
#my $lastLine=scalar @outArray;

my $count=0; #number of windows added to cumm that will be used for avg calc
my $sumFit=0;
my $sumRatio=0;

print "Start grouped txt file creation time: ",get_time(),"\n";



#What's the mean fitness value for all of the windows?
my @allFits = map $_->[2], @windows;
my $meanFit=mean(@allFits); #not sure what module this needs

#Add the absolute deviation from mean fitness to each window array (use this to sort)
my @expWindows=();
foreach (@windows){
	my @entry=@{$_};	
	my @foo=($entry[0],$entry[1]);
	push @expWindows,\@foo;	
	} 
	
#Now sort @expWindows by the abs. dev. of fitness (index 8). Reuse old @windows variable
@windows= sort {$a->[0]<=>$b->[0] } @expWindows;


open (OUT, '>',"filteredCoords.csv") or die $!;

#print the header
push(@header, "abs(diff_mean)");
my $string = join(",", @header);
#print OUT "$string\n"; 
my $i=1;
foreach (@windows){
	my $outer = join("\t", @{$_});
	#print OUT "A$i\t";
	print OUT $outer;
	print OUT "\n";
	$i++;
	}
	
close OUT;
	
