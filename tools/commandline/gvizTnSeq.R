# Instruction on how to input Tn-Seq data (output files from MAGenTA) into R for visualization
# Two methods shown here: (1) with Gviz (2) with standard R plotting

#Step 1: Install packages

# During installation step, allow updates or installing of
# dependency packages if prompted
# Run these two lines to check if necessary packages are
# installed If not, they will be installed.
packages <- c("Biostrings", "seqinr", "Gviz")
if (length(setdiff(packages, rownames(installed.packages()))) >
0) {
install.packages(setdiff(packages, rownames(installed.packages()))) }
# Load libraries for installed packages
library(seqinr) 
library(Biostrings) 
library(Gviz)

#Step 2: Specify input files. Sample file paths shown here.

# Set working directory
# setwd('~/Documents/lab/tvolab/magenta/')
# Specify file paths in the beginning of a session
# Fasta file for the reference genome, available on NCBI genomes
genomeFile = "19F_012469.fasta"
# Gene Features formatted for Gviz track from a Genbank file using the
# GetGeneFeatures tool
# Format: tab-delimited text file with gene_id, gene_name, start, end , strand, # and function fields
geneFile = "19F_012469_genes.txt"
# Single fitness files are created by the SingleFitness tool (Galaxy or command-line) # using input library csv files from the CalculateFitness tool
# Here, a control (glucose) file and experimental (glucose + daptomycin) are used ctrlSingleFitFile = "19fglucSingleFit.csv"
exptSingleFitFile = "19fdaptoSingleFit.csv"

# Step 3: Genome Axis and Sequence Tracks

# so Gviz doesn't look for NC_003028 in the UCSC site
options(ucscChromosomeNames = FALSE)
# MAGenTA input: None Number coordinates for genome
gTrack <- GenomeAxisTrack(littleTicks = TRUE)
# MAGenTA input: FASTA file (genomeFile).  Track for sequence (ATGC) in color
fcol <- c(A = "darkorange", C = "yellow", T = "darkred", G = "darkgreen") sTrack <- SequenceTrack(genomeFile, name = "Sequence", fontcolor = fcol)

# Step 4: Genome Annotation Track

# MAGenTA input: A formatted file containing gene id, start, end coordinates, etc.
# from Genbank file using script getCoordsGBK.py (geneFile) adds gene ids
# Make sure coordinate columns 2 and 3 are numeric.  If not then reassign as numeric:
# column is.numeric(genes[,2])
geneDF <- as.data.frame(read.table(file = geneFile, header = TRUE))
anTrack <- AnnotationTrack(start = geneDF$start, end = geneDF$end, chromosome = "genome", strand = geneDF$strand,
id = geneDF$id, showFeatureId = TRUE, name = "Taiwan-19F\nGenes", fill = "gray", fontcolor.item = "black")

# Step 5: Data Track

# MAGenTA input: single aggregate fitness per insertion site, from SingleFitness tool
makeJointAggData <- function(file1, file2) {
    # Adjust file format
insertFit <- read.csv(file1)
insertFit$seqnames = "genome"
insertFit$start = insertFit$pos
insertFit$end = insertFit$pos + 1
keepCols <- c("seqnames", "start", "end", "fitness")
insertFit <- insertFit[keepCols]
colnames(insertFit) <- c("seqnames", "start", "end", "control") insertFit2 <- read.csv(file2)
insertFit2$seqnames = "genome"
insertFit2$start = insertFit2$pos
insertFit2$end = insertFit2$pos + 1
keepCols <- c("seqnames", "start", "end", "fitness") insertFit2 <- insertFit2[keepCols]
merged <- merge(insertFit, insertFit2, by = c("seqnames", "start", "end"), all = TRUE)
colnames(merged) <- c("seqnames", "start", "end", "glucose", "glucose.daptomycin")
insertData <- makeGRangesFromDataFrame(merged, keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL, seqnames.field = "seqnames", start.field = "start", starts.in.df.are.0based = FALSE)
return(insertData) }
jointInsertData <- makeJointAggData(ctrlSingleFitFile, exptSingleFitFile)
count = length(jointInsertData)
trackTitle = "Insertion Fitness Cost\nfor Taiwan-19F mutants"
# Make track of insertion fitness. Histogram is used here, but other options are available # See Gviz documentation.
colors = c("royalblue4", "firebrick1")
aggTrack <- DataTrack(jointInsertData, chromosome = "genome",
groups = c("glucose", "glucose.daptomycin"), ylim = c(0,
2), type = c("p"), col.baseline = "gray", baseline = c(0, 0.5, 1, 1.5, 2), col = c("royalblue4", "firebrick1"), legend = TRUE, name = trackTitle, fontcolor.legend = "black", fontsize.legend = 12, box.legend = TRUE)

# Step 6: Plot Tracks

# Specify viewing window start and end coordinates
fromCoord = 1370000
toCoord = 1371700
# Change surrounding flanking region x = 1000
y = 1000
# Highlight Track
ht <- HighlightTrack(trackList = list(aggTrack), fill = "lightyellow",
col = "transparent", start = fromCoord, end = toCoord, chromosome = "genome")
# Plot all Tracks
plotTracks(list(gTrack, anTrack, ht), chromosome = "genome", background.title = "lightgray", fontcolor = "black", fontsize = 10, from = fromCoord - x, to = toCoord +
y, type = c("p", "smooth"), fontcol = "black",
col.title = "black", axis.col = "black", fontface.main = 0.8, fontfamily = "sans")

# SINGLE TRACK VISUALIZATION WITHOUT GVIZ

mainTitle = "Taiwan-19F Mutant Fitness in Glucose & Daptomycin"
visSingleFit <- function(ctrlFile, exptFile, x1, x2, h1, h2) { ctrl <- read.csv(ctrlFile)
expt <- read.csv(exptFile)
plot(ctrl$pos, ctrl$fitness, xlab = "Genomic Coordinate",
        ylab = "Insertion Fitness Cost", main = mainTitle, type = "h",
xlim = c(x1, x2), col = "darkblue")
legend("topright", 1.35, c("Glucose", "Glucose + Daptomycin"),
lty = c(1, 1), lwd = c(1.5, 2.5), col = c("darkblue", "red"), bty = "n", cex = 0.6)
abline(h = 1, col = "black")
lines(expt$pos, expt$fitness, col = "red", type = "h") rect(h1, 0, h2, 1.3, col = "#F4FA5840", border = NA)
}

visSingleFit(ctrlSingleFitFile, exptSingleFitFile, 1370000, 1372514,
    1370347, 1371514)

