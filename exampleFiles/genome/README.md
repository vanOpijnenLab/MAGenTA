# Genome Files

Example genome files corresponding the strains represented in the example data files are included here. Genbank and Fasta files can also be found on NCBI genome. The "gene" files included (tigr4_nc003028_genes.txt and 19F_012469_genes.txt) contain genome, start and end coordinates, length, strand, and gene id for each gene in the genome. These files can be generated from a Genbank file and the [getCoordsGBK.py script](https://github.com/antmarge/Bioinformatics-Tools/blob/master/getCoordsGBK.py) here:
```
https://raw.githubusercontent.com/antmarge/Bioinformatics-Tools/master/getCoordsGBK.py
```
The genes files are formatted specifically to create a [GRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) object which is necessary for Gviz multi-track visualization.
