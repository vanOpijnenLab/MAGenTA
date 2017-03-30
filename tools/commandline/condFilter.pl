#!/usr/bin/perl -w

#Margaret Antonio
#Began: 15.12.26 Updated: 16.01.17
#DESCRIPTION: Filter to keep (lines) that satisfy the specified condition
#Improvements to make: allow multiple indices and conditions for filtering

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use diagnostics;

#ASSIGN INPUTS TO VARIABLES
our ($help,$out,$kIndex,$kOper,$kVal);
GetOptions(
'o:s' =>\$out,
'h' => \$help,
'ki:i' => \$kIndex,   #keep index #any number that is an index (i for integer)
'ko:s' => \$kOper,    #keep operator condition. Can be e,l,g,le,ge,ne. (s for string)
'kv:f' => \$kVal,     #keep value #any number (f for float)
);


sub print_usage() {
    print "\n";
    print "\n####################################################################\n";
    print "\n";
    print "condFilter.pl: Conditionally filter lines in a csv file\n\n";
    print "DESCRIPTION: Takes a single file and keeps only lines that satisfy condition\n";
    print "the difference in mean fitness, the pval for each gene.\n";
    print "EXAMPLE: Keep all lines in a comparison file where the difference between\n";
    print "\tmean fitnesses is greater than .15\n";
    
    print "\nUSAGE: perl condFilter.pl -ki 5 -ko .05 -kv le inputfile.csv\n\n";
    
    print "OPTIONS:\n\n";
    print " -h \t Prints usage and quits\n";
    print "    \t Input comma-separated-value file (CSV) without a flag\n";
    print " -ki\t Index of column to test keep value and keep operator.\n";
    print "    \t Indices begin at 0.\n";
    print " -ko\t String version of operator to use (e, l, g, le, ge, ne).\n";
    print " -kv\t Value to use in conditional test.\n";
    print " -o \t Output file for filtered data. Default: filt.csv\n";

    print "\n####################################################################\n";
    print "\n";
}
#IF HELP FLAG (-h) WAS SPECIFIED THEN PRINT USAGE AND QUIT
if ($help){
    print_usage();
    exit;
}

#INTERPRET KEEP CONDITION FROM STRING TO OPERATOR. AND CHECK FOR OTHER PARTS OF CONDTIONAL TEST.
#Can't directly use operator because > is interpreted as output into (> file)


if ($kOper and $kOper eq "e"){
    $kOper="==";
}
elsif ($kOper and $kOper eq "l"){
    $kOper="<";
}
elsif ($kOper and $kOper eq "g"){
    $kOper=">";
}
elsif ($kOper and $kOper eq "le"){
    $kOper="<=";
}
elsif ($kOper and $kOper eq "ge"){
    $kOper=">=";
}
elsif ($kOper and $kOper eq "ne"){
    $kOper="!=";
}
else {
    print "\nPlease enter one of the following for the keep operator (-ko): e, l, g, le, ge, ne";
    print_usage();
    exit;
}
if (!$kVal){
    print "\nPlease enter the value to use in the conditional test for filtering. -kv <value> \n\n";
    print_usage();
    exit;
}
if (!$kIndex){
    print "\nPlease enter the column index to use for filtering. -ki <index> \n\n";
    print_usage();
    exit;
}
if (!$out) {$out="filt.csv"};

sub cleaner{
    my $line=$_[0];
    chomp($line);
    $line =~ s/\x0d{0,1}\x0a{0,1}\Z//s;
    return $line;
}

#READ INPUT FILE, FILTER, AND IMMEDIATELY WRITE TO THE OUTPUT FILE
my @outArray;

my $in=$ARGV[0];
open (OUT,'>',"$out") or (print "Could not open '$out'\n" and print_usage() and exit);
open(DATA, '<', "$in") or (print "Could not open '$in'\n" and print_usage() and exit);

#REUSE the header (column names) from the input file for the output file
my $head=<DATA>;
$head=cleaner($head);
my @header=split(',',$head);
print OUT (join(",",@header),"\n");

while (my $entry = <DATA>) {
	$entry=cleaner($entry);
	my @fields=split(',',$entry);
    my $testVal=$fields[$kIndex];
    if (eval($testVal.$kOper.$kVal)){
        print OUT (join(",",@fields),"\n");
    }
}
close DATA;
close OUT;

