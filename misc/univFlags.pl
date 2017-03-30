
calcFitness.pl
   
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
    
dataOverview.pl

    print "\nRequired:\n";
    print "In the command line (without a flag), input the name(s) of the file(s) containing fitness values for individual insertion mutants.\n";
    print "\n perl dataOverview.pl <OPTIONS> --indir path/to/directory/with/results/files \n\n";
    print "\n Example:\n";
    print "perl ../Blueberries/dataOverview.pl -i ../9-daptomycin/data/19FGluc/ -f ../0-genome/19F_012469.fasta -r ../0-genome/19F_012469.gbk -l\n";

    print "\nOPTIONS:\n";
	print "-r\tThe name of the reference genome file, in GenBank format. Needed for wig and txt file creation\n";
	print "-l\t Send all output to a log file instead of the terminal\n";
	print "-f\tFasta file for genome\n";
    print "-i\t Directory containing all input files (results files from calc fitness script\n";
    print "-h\t Print usage\n";
    print "-c\t Cutoff average(c1+c2)>c\n";
    print "-o\t Outfile\n";

aggregate_MLA.pl

    print "\n";
    print "Usage: ./aggregate.pl -o outfile -d indirectory -r genbankFile -f fastaFile\n\n";
    print "Description: Finds average fitness of insertion mutations within annotated genes\n\n";
    print "Option List:\n\n";
    print " -o\tOutput file for aggregated data. (Required)\n";
    print " -d\tDirectory containing input files. Make sure / is included after name\n";
    print " -r\tCheck for missing genes in the data set - provide a reference genome in\n";
    print " \tgenbank format. Missing genes will be sent to stdout.\n";
    print " -m\tPlace a mark in an extra column for this set of genes. Provide a file\n";
    print " \twith a list of genes seperated by newlines.\n";
    print " -x\tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    print " -b\tBlanks: Exclude -b % of blank fitness scores (scores where c2 = 0) (default: 0 = 0%)\n";
    print " -w\tUse weighted algorithm to calculate averages, variance, sd, se\n";
    print " -l\tWeight ceiling: maximum value to use as a weight (default: 999,999)\n";
    print "-g\tGene list with coordinates so essentiality (pvalue) can be calculated\n";
    print "-f\tFasta file\n";
    print "-r\tRound to this number of decimals\n";
    print "-n\tNull distribution";
    print "\n";
  
  
  slidingwindow.pl  
'wig:s' => \$wig,
'ref:s' => \$ref_genome,
'cutoff:i'=>\$cutoff,
'in:s' => \$infile,
'csv:s'  => \$csv,
'step:i' => \$step,
'size:i' => \$size,
'txtg:s' => \$txtg,
'txt:s' => \$txt,
'random:s' =>\$random,
'round:i' =>\$round,
'fasta:s' => \$fasta,
'outdir:s' =>\$outdir,
'log' => \$log,
'usage' => \$h,
'tan'=>\$tan,
'indir:s'=>\$indir,
'inc:i'=>\$inc,
'sort:i'=>\$sortby,
'wc:i'=>\$weight_ceiling,
'nw'=>\$noweight,
'c:s'=>\$custom,