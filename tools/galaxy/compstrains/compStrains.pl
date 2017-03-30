#!/usr/bin/perl -w

#Margaret Antonio 16.01.13

#DESCRIPTION: Takes two aggregate.pl outputs and compares them using mean difference, pval for each
#gene. Can compare, for example, 19F in glucose and TIGR4 in glucose.
#DIFFERENT GENOMES (ie. diff. strains).
#Requires CONVERSION FILE


use Data::Dumper;
use strict;
use Getopt::Long;
#use warnings;
use File::Path;
use File::Basename;
use Statistics::Distributions;

#ASSIGN INPUTS TO VARIABLES USING FLAGS
our ($input1,$input2,$out,$sortkey,$round,$l1,$l2,$cfile);
GetOptions(
'o:s' =>\$out,
's:i' => \$sortkey,
'r:i'=> \$round,
'l1:s'=> \$l1,
'l2:s'=> \$l2,
'c:s'=> \$cfile,
'input1:s'=>\$input1,
'input2:s'=>\$input2
);


#THE @files ARRAY WILL CONTAIN INPUT FILE NAMES
my @files=($input1,$input2);

#GET LABELS:
my @labels = ($l1,$l2);

sub cleaner{
    my $line=$_[0];
    chomp($line);
    $line =~ s/\x0d{0,1}\x0a{0,1}\Z//s;
    return $line;
}


#CHECK IF REQ. VARIABLES WERE DEFINED USING FLAGS. IF NOT THEN USE DEFAULT VALUES
if (!$out) {$out="comp.".$labels[0].$labels[1].".csv"}
if (!$round){$round='%.4f'}

#OPEN INPUTTED AGGREGATE GENE FILES AND STORE THEIR CONTENTS INTO TWO HASHES
#FILE1 GOES INTO HASH %ONE AND FILE2 GOES INTO HASH %TWO.

#FILE1 OPENING ---> %one WHERE KEY:VALUE IS GENE_ID:(GENE_ID,INSERTIONS,MEAN,ETC.)
my @header;
my %one;

open (F1,'<',$files[0]);

#STORE COLUMN NAMES (FIRST LINE OF FILE1) FOR HEADER AND APPEND LABELS
my $head=<F1>; #the header in the file
my @cols=split(',',$head);
@cols=@cols[0,1,2,3,4,5]; #get rid of blank columns
for (my $j=0;$j<scalar @cols;$j++){
    $cols[$j]=cleaner($cols[$j]).'-'.$labels[0];   #mark each column name with file it comes from
}
push (@header,@cols);

while (my $line=<F1>){
    chomp $line;
    my @info=split(",",$line);
    #Only keep the first 6 columns (Ones about blanks aren't needed for comparisons)
    @info=@info[0,1,2,3,4,5];
    #Sometimes genes that don't have a gene name can't be blank, so fill with NA
    if (!$info[5]){
        $info[5]="";
    }

    $one{$info[0]}=\@info;
}
close F1;

#FILE2 OPENING ---> %two WHERE KEY:VALUE IS GENE_ID:(GENE_ID,INSERTIONS,MEAN,ETC.)

my %two;
open (F2,'<',$files[1]);

#STORE COLUMN NAMES (FIRST LINE OF FILE2) FOR HEADER AND APPEND LABELS
$head=<F2>; #the header in the file
@cols=split(',',$head);
@cols=@cols[0,1,2,3,4,5]; #get rid of blank columns
for (my $j=0;$j<scalar @cols;$j++){
    $cols[$j]=cleaner($cols[$j]).'-'.$labels[1];   #mark each column name with file it comes from
}
push (@header,@cols);

while (my $line=<F2>){
    $line = cleaner($line);
    my @info=split(",",$line);
    @info=@info[0,1,2,3,4,5];
    if (!$info[5]){
        $info[5]="";
    }

    $two{$info[0]}=\@info;
}
close F2;


#READ CONVERSION FILE INTO ARRAY.
#Conversion file must have strain 1 for file 1 in column 1 (index 0) and
    #strain 2 for file 2 in column 2 (index 1)
    #conversion file must be tab delimited with no NA fields
#If homologs (exist then take info from hashes (%one and %two) by referring to gene_id in KEY

my @all; #store all homologs in this hash
open (CONV,'<',$cfile);
while (my $line=<CONV>){
    $line = cleaner($line);
    my @genes=split("\t",$line);   #Array @genes will contain two genes (SP_0000,SPT_0000)
    if (scalar @genes==2 and (exists $one{$genes[0]}) and (exists $two{$genes[1]})){
        my @info;
        my @oneArray=@{$one{$genes[0]}};
        my @twoArray=@{$two{$genes[1]}};
        push (@info,@oneArray,@twoArray);
        # Fitness values at index 1 and 7
        my $diff=sprintf("$round",($info[1]-$info[7]));
        my $total1=$info[2];
        my $total2=$info[8];
        if (!$info[3] or $info[3] eq ""){$info[3]="NA"};
        if (!$info[4] or $info[4] eq ""){$info[4]="NA"};
        if (!$info[9] or $info[9] eq ""){$info[9]="NA"};
        if (!$info[10] or $info[10] eq ""){$info[10]="NA"};
        my $sd1=$info[3];
        my $se1=$info[4];
        my $sd2=$info[9];
        my $se2=$info[10];
        my $df=$total1+$total2-2;
        my ($tdist,$pval);
        #TDIST, PVAL calculations with fail if standard dev, error, or counts are not real numbers
        #or if 0 ends up in denominator
        if ($sd1 eq "NA" or $sd2 eq "NA" or $total1==0 or $total2==0 or $sd1==0 or $sd2==0){
            ($tdist,$pval)=("NA","NA");
        }
        else{
            $tdist=sqrt((($diff)/(sqrt((($sd1**2)/$total1)+(($sd2**2)/$total2))))**2);
            $pval=Statistics::Distributions::tprob($df,$tdist);
        }
        push (@info,$diff,$df,$tdist,$pval);
        push (@all,\@info);
    }
}
close CONV;

#SORT THE HOMOLOGS BY THE SORTKEY OR BY DEFAULT DIFFERENCE IN MEAN FITNESSES
if (!$sortkey){
    $sortkey=12; #for mean difference
}
my @sorted = sort { $b->[$sortkey] <=> $a->[$sortkey] } @all;

#FINISH THE HEADER BY ADDING COLUMN NAMES FOR MEAN-DIFF, DOF, TDIST, AND PVALUE
my $field="MeanDiff(".$labels[0].'.'.$labels[1].")";
push (@header,$field,"DOF","TDIST","PVALUE");

#PRINT MATCHED HOMOLOG INFORMATION INTO A SINGLE OUTPUT FILE
open OUT, '>',"$out";
print OUT (join(',',@header),"\n");
foreach (@sorted){
    my @comparison=@{$_};
    print OUT join(',',@comparison),"\n";
    }

close OUT;


