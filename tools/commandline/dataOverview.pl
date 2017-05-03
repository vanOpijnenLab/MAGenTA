#!/usr/bin/perl -w

#Margaret Antonio 16.08.29

#use strict;
use Getopt::Long;
use Bio::SeqIO;
use autodie;
no warnings;

#AVAILABLE OPTIONS. WILL print OUT UPON ERROR
sub print_usage() {

   print "\n###############################################################\n";
    print "dataOverview: outputs basic statistics for tn-seq library files \n\n";

    print "USAGE:\n";
    print  "perl dataOverview.pl -i inputs/ -f genome.fasta -r genome.gbk\n";
        
    print  "\nREQUIRED:\n";
    print  " -d\tDirectory containing all input files (results files from\n";
    print  "   \tcalc fitness script)\n";
    print  "   \t  OR\n";
    print  " \tIn the command line (without a flag), input the name(s) of \n";
    print  " \tthe files containing fitness values for individual \n\tinsertion mutants\n";
    print  " -f\tFilename for genome sequence, in fasta format\n";
    print  " -r\tFilename for genome annotation, in GenBank format\n";

    print  "\nOPTIONAL:\n";
    print  " -h\tprint OUT usage\n";
    print  " -c\tCutoff average(c1+c2)>c. Default: 15\n";
    print  " -o\tFilename for output. Default: standard output\n";
    print  " \n~~~~Always check that file paths are correctly specified~~~~\n";
    print  " \n###############################################################\n";

}

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
    }
sub mean {
    my $sum=0;
    foreach my $n(@_){
    	$sum+=$n;
    	}
    my $total=scalar @_;
    my $mean=$sum/$total;
    return $mean;  	
}
sub minmax{
    my @unsorted=@_;
	my @sorted = sort { $a <=> $b } @unsorted;
	my $min = $sorted[0];
	my $max = $sorted[scalar @sorted -1];
	return ($min, $max);
	}
sub uniq{
    my @input=@_;
    my @unique = do { my %seen; grep { !$seen{$_}++ } @input };
    }

#ASSIGN INPUTS TO VARIABLES
our ($cutoff,$fastaFile, $outfile,$help,$ref,$indir,$weight_ceiling);
GetOptions(
'r:s' => \$ref,
'f:s' => \$fastaFile,
'd:s'=>\$indir,
'c:i'=>\$cutoff,
'o:s' => \$outfile,
'h'=> \$help,
'w:i' => \$weight_ceiling,
);

# Set defaults
if (!$weight_ceiling){$weight_ceiling=999999;}
if (!$cutoff){$cutoff=15;}

# If help option is specified or required files are not specified:

if ($help) {
    print print_usage();
	print "\n";
	exit;
}
if (!$indir and (scalar @ARGV==0)){
	print  "\nERROR: Please correctly specify input files or directory\n";
    print print_usage();
	print  "\n";
	exit;
}
if (!$fastaFile or !$ref){
	print  "\nERROR: Please correctly specify reference genome fasta and genbank files\n";
	print "Most genomes (in fasta and gbk format) are available at NCBI\n";
    print print_usage();
	print "\n";
	exit;
}

# Redirect STDOUT to log.txt. Anything print OUTed to the terminal will go into the log file
if (! $outfile){
	$outfile="summary.txt";
}

open OUT, ">",$outfile;

#Not sure if I'll need this but sometimes funky data inputs have hidden characters
sub cleaner{
    my $line=$_[0];
    chomp($line);
    $line =~ s/\x0d{0,1}\x0a{0,1}\Z//s;
    return $line;
}


#Get the input files out of the input directory, or take off of command line
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

#print OUT "Gathering data overview for Tn-Seq experiment\n\n";
#print OUT "Begin time: ",get_time(),"\n\n";

#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S). 
#These are the "results" files from calc_fitness.pl. Insertion location, fitness, etc.
#Go through each file from the commandline (ARGV array) and read each line as an array
#into select array if values satisfy the cutoff


#Store ALL insertion locations in this array. Later, get unique insertions
my @insertions_all;
#Store all genes with valid insertions here
my @genes_insertions;
#all lines that satisfied cutoff
my @unsorted;
#array to hold all positions of insertions. Going to use this later to match up with TA sites
my @insertPos;

#Markers
my $rows=-1;
my $last=0;

print OUT "Library description\n\n";
my @header=("library","file_path","ins","ins.f","genes.ins");
print OUT join ("\t",@header),"\n";

for (my $i=0; $i<$num; $i++){
	#Temp arrays for library
	my(@insertions_all_lib,@genes_insertions_lib,@insertPos_lib);
    my $file=$files[$i];
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
        my $coord= $line[0];
        push (@insertions_all_lib,$coord);
         #Average counts must be greater than cutoff (minimum allowed)
        my $avg = ($c1+$c2)/2;
        if ($avg > $cutoff) {
        	my @select=($coord,$w,$avg,$locus);
            my $select=\@select;
            push(@unsorted,$select);
            push(@insertPos_lib,$line[0]);   #keep track of actual insertion site position
            push (@genes_insertions_lib,$locus);
            $last=$select[0];
            $rows++;
        }
        if ($avg >= $weight_ceiling) { $avg = $weight_ceiling } # Maximum weight
    }
    close DATA;
    push (@insertions_all,@insertions_all_lib);
    @genes_insertions_lib= uniq @genes_insertions_lib;
    push (@genes_insertions,@genes_insertions_lib);
    push (@insertPos,@insertPos_lib);
    my @stat=($i+1,$file,scalar @insertions_all_lib,scalar @insertPos_lib,scalar @genes_insertions_lib);
    print OUT join("\t",@stat),"\n";
}

@insertPos = sort { $a <=> $b } @insertPos;
@insertPos= uniq @insertPos;
@genes_insertions= uniq @genes_insertions;
@insertions_all=uniq @insertions_all;
my $totalAll=scalar @insertions_all;
my $total=scalar @insertPos;
my $temp="1-".$num;
my @all_stat=($temp,"NA",$totalAll,$total,scalar @genes_insertions);
print OUT join("\t",@all_stat),"\n";

#Genome description: #TA sites, distance between TA sites, #TA sites in ORFS
print OUT "\n-------------------------\n";
print OUT "\nGenome description\n\n";
print OUT "File for genome: ", $fastaFile,"\n";

my @sites;
#First read fasta file into a string
my $seqio = Bio::SeqIO->new(-file => $fastaFile, '-format' => 'Fasta');
my $fasta;
while(my $seq = $seqio->next_seq) {
	$fasta = $seq->seq;
}
#Just in case $fasta file is in lowercase, change it to uppercase
$fasta=uc $fasta;

#Get genomic coordinate for TA sites:
my $x="TA";
my $offset=0;
my @indices;
my $result=index($fasta,$x,$offset);
while ($result !=-1){
	push (@indices,$result);
	$offset=$result+1;
	$result=index($fasta,$x,$offset);
}
my $countTA=scalar @indices;

#Get longest stretch with no TA sites
my @tempta=@indices;
my $prev=shift @tempta;
my $current=shift @tempta;
my $lg_dist_ta=$current-$prev;
foreach my $site(@tempta){
	$prev=$current;
	$current=$site;
	my $d=$current-$prev;
	if ($d>$lg_dist_ta){
		$lg_dist_ta=$d;
	}
}

#Get longest stretch of with no insertions
my @tempins=@insertPos;
$prev=shift @tempins;
$current=shift @tempins;
my $lg_dist_ins=$current-$prev;
foreach my $site(@tempins){
	$prev=$current;
	$current=$site;
	my $d=$current-$prev;
	if ($d>$lg_dist_ins){
		$lg_dist_ins=$d;
	}
}


my $genSize=length $fasta;
print OUT "$genSize\tGenome size\n";
print OUT "$countTA\tTotal number of TA sites\n\n";

my $sat=sprintf("%.2f", ($total/$countTA)*100);
my $satAll=sprintf("%.2f", ($totalAll/$countTA)*100);
my $inscov=sprintf("%.2f", ($total/$genSize)*100);
my $tacov=sprintf("%.2f", ($countTA/$genSize)*100);

#Get GC content of genome

my $sequence = ' ';
my $Ccount = 0;
my $Gcount = 0;
my $identifier = ' ';

my @nucleotides = split('', $fasta);

foreach my $nuc (@nucleotides) {
	if ($nuc eq 'G') {$Gcount++} 
	elsif ($nuc eq 'C') {$Ccount++}
}
my $sequencelength=length $fasta;

my $GCcontent = sprintf("%.2f",((($Gcount + $Ccount) / $sequencelength) * 100));
my $ATcontent =100-$GCcontent;

print OUT "$GCcontent%\tGC content of this genome\n";
print OUT "$ATcontent%\tAT content of this genome\n";

print OUT "$satAll%\tSaturation of TA sites before cutoff filter (allInsertions/TAsites)\n";
print OUT "$sat%\tSaturation of TA sites after cutoff filter (validInsertions/TAsites)\n";
print OUT "$inscov%\tGenome coverage by insertions (validInsertions/genomeSize)\n";
print OUT "$tacov%\tGenome coverage by TA sites (TAsites/genomeSize)\n";
print OUT "$lg_dist_ta\tLargest distance between TA sites\n";
print OUT "$lg_dist_ins\tLargest distance between insertions\n";

#Store everything to be print OUTed in array
my @table;

print OUT "\nGenes using the genbank annotation file\n\n";
###Get genbank file. Find all start and stop for genes
#See how many insertions fall into genes vs intergenic regions
#Get array of coordinates for all insertions then remove insertion if it is
#within a gene region
my $gb = Bio::SeqIO->new(-file=>$ref);
my $refseq = $gb->next_seq;

#store number of insertions in a gene here
my @geneIns;
my @allLengths;
my $blankGene=0; #Number of genes that don't have any insertions in them
my @genomeSeq=split('',$fasta);


#keep a copy to remove insertions that are in genes
my @insertPosCopy=@insertPos;

my @features = $refseq->get_SeqFeatures(); # just top level
foreach my $feature ( @features ) {
	if ($feature->primary_tag eq "gene"){
		my $start=$feature->start;
		my $end=$feature->end;
		my $length=$end-$start;
		push (@allLengths,$length);
		#turn this into a for loop
		my $i=0;
		my $ins=0;
		my $current=$insertPos[$i];;
		while ($current<=$end && $i<scalar @insertPos){
			if ($current>=$start){
				splice(@insertPosCopy, $i, 1);
				$ins++;			
			}
			$current=$insertPos[$i++];
		}
		if ($ins==0){$blankGene++}
		push (@geneIns,$ins);
	}
}
my $avgLength=sprintf("%.2f",mean(@allLengths));

my ($minLength, $maxLength) = minmax @allLengths;
my $avgInsGene=sprintf("%.2f",mean(@geneIns));





my ($minInsGene, $maxInsGene) = minmax @geneIns;
my $nonGeneIns=scalar @insertPosCopy;
my $totalIns=scalar @insertPos;
my $percNon=sprintf("%.2f",($nonGeneIns/$totalIns)*100);

print OUT "Length of a gene\n";
print OUT "$avgLength\tAverage","\n$minLength\tMininum","\n$maxLength\tMaximum\n";
print OUT "\nFor insertions in a gene:\n";
print OUT "$avgInsGene\tAverage","\n$minInsGene\tMininum","\n$maxInsGene\tMaximum\n";
print OUT "Number of genes that do not have any insertions: ",$blankGene,"\n";
print OUT "\n$nonGeneIns\tInsertions that are not in genes","\n$percNon% of all insertions\n";
#How many insertions are in genes and how many are in non-gene regions?



