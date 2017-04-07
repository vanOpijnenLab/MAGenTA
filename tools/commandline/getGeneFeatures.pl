#!/usr/bin/perl -w

#Margaret Antonio 16.04.22

#Get genbank features (coordinates, gene name, and gene function/description) for genes

#1. Get Perl module for gbk files
#2. Get important features into hash
#3. Read comp file with option for gene_id column
#4. Output full file

use Bio::DB::GenBank;
use strict;
use Getopt::Std;
use Bio::SeqIO;
use Set::Scalar;

# Get global options
our($opt_r,$opt_i,$opt_o,$opt_l);
getopt('riol');

#This is the output file with the full gene information
my $outfile = $opt_o || "impFeatures.csv";

#This is the reference genome file (genbank) from which features will be extracted
my $ref = $opt_r;
#This is the input file that needs the gene features
my $infile = $opt_i;
#The index of the gene locus
my $index=$opt_l || 0;

#Read in the genbank file with the Bio::DB::GenBank perl module.
#Full manual for module: http://search.cpan.org/dist/BioPerl/Bio/DB/GenBank.pm
#Genbank(.gb) files can't normally be read as a txt and csv files are
my $gb = Bio::SeqIO->new(-file=>$ref);
my $refseq = $gb->next_seq;

# Hashify features for easy lookup

open OUT,">",$outfile;

my %feature_summary;
my @features = $refseq->get_SeqFeatures(); # just top level
foreach my $feature ( @features ) {
	if ($feature->primary_tag eq "gene"){
		my $start=$feature->start;
		my $end=$feature->end;
		my $dir=$feature->strand;
		
		my @locus = $feature->get_tagset_values('locus_tag');
		my $locus = $locus[0] || '';

		my @gene = $feature->get_tagset_values('gene');
		my $gene = $gene[0] || '';
		
		$feature_summary{$locus}{g} =  $gene; 
		$feature_summary{$locus}{s} =  $start;
		$feature_summary{$locus}{e} =  $end;
		$feature_summary{$locus}{d} =  $dir;
	}
	elsif($feature->primary_tag eq "CDS"){
		my $start=$feature->start;
		my $end=$feature->end;
		my $dir=$feature->strand;
		my @locus = $feature->get_tagset_values('locus_tag');
		my $locus = $locus[0] || '';
				my @gene = $feature->get_tagset_values('gene');
		my $gene = $gene[0] || '';
		if (!$feature_summary{$locus}){
			print $locus;
			$feature_summary{$locus}{g} =  $gene; 
			$feature_summary{$locus}{s} =  $start;
			$feature_summary{$locus}{e} =  $end;
			$feature_summary{$locus}{d} =  $dir;
			}
			
		my @product = $feature->get_tagset_values('product');
		my $product = $product[0] || '';
		#replace commas in section with ; because making a csv file
		$product =~ s/,/;/g;
		
		my @note = $feature->get_tagset_values('note');
		my $note = $note[0] || '';
		$note =~ s/,/;/g;
		
		$feature_summary{$locus}{p} = $product;
		$feature_summary{$locus}{n} = $note;

		
		#Send all of the relevant features to an output file
		#my @feat=($feature->primary_tag,$locus,$gene,$start,$end,$dir,$product,$note);
		#print OUT join(",",@feat);
		#print OUT "\n";
		
		#Store features in a multi hash so they can be retrieved
		#and appended to gene fitness comparison file

	}
	#Now go retrieve next feature in the genbank file
}
close OUT;

#Now that we have a multi hash of all the features of interest for a gene
#Append these features to a file that has the gene locus

open (IN,"<",$infile) or print "Cannot open input file!\n" and die;
my @full;
my $header=<IN>;
while(<IN>){
	chomp $_;
	my @fields=split(",",$_);
	my $key=$fields[$index];
	if ($feature_summary{$key}){
		my @new=($feature_summary{$key}{s},$feature_summary{$key}{e},
		$feature_summary{$key}{d},$feature_summary{$key}{p});
		push (@fields,@new);
	}
	push (@full,\@fields);
}
open OUTFILE,">","19F_glucdapto_fullFeat3.csv";

print OUTFILE $header;
foreach (@full){
	my @entry=@{$_};
	print OUTFILE join(",",@entry);
	print OUTFILE "\n";
}

close OUTFILE;
		
