# Data Files
#### Sample data files for a Tn-Seq experiment of TIGR4 and Taiwan-19f (S. pneumoniae strains) in glucose and daptomycin.

#### Data files shown as examples
Each directory contains the results files from the Calculate Fitness tool. There are 6 results files for each strain-condition pairing (i.e. 19F in daptomycin) representing each of the 6 libraries for this particular experiment. These Calculate Fitness results files (label e.g. L1_19F_Dapto.csv) are inputs to several analysis tools: SingleFitness, RegionFitness, and AggregateFitness. In each strain-condition data directory, there is also the output file from the SingleFitness tool, which is an aggregate of all the library files (L1-L6) for that strain-condition pairing. The SingleFitness output file is used primarily for the Gviz multi-track visualization in R. 

#### Where to find corresponding genome files
Corresponding genome files for the two strains shown in these data files are located in the genome directory. Genome files include: Fasta, Genbank, and Gene files for both Taiwan-19F and TIGR4 strains.

#### Download and decompress tar.gz files
Each strain-comparison directory is compressed. Download the tar.gz file and double click the file to uncompress or uncompress the file in the terminal with:

```
tar -zxvf archive_name.tar.gz
```

