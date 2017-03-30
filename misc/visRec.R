getHom<-function(gene_id,x,y){
  
  #Load libraries
  library(seqinr)
  library(Biostrings)
  library(Gviz)
  
  genomeFile="NC_003028.fa"
  geneFile="tigr4_genes.txt"
  singleFitFile="singleFit_tigr4_SDDM.txt"
  
  #genomeFile="NC_012469.fa"
  #geneFile="19Ftest.txt"
  #singleFitFile="singleFit_19F_SDDM.txt"
  
  #gene id was inputted so need to get start and end coordinates for that gene
  geneDF1<-as.data.frame(read.table(file = geneFile,header=TRUE))
  geneDF1$seqnames="genome"
  search<-geneDF1[grep(gene_id,geneDF1$id),]
  fromCoord<-as.numeric(search$start)
  toCoord<-as.numeric(search$end)
  print (fromCoord)
  print (toCoord)
  
  
  #so Gviz doesn't look for NC_003028 in the UCSC site
  options(ucscChromosomeNames = FALSE)
  
  #GenomeAxisTrack: Number coordinates for genome (no file input)
  gTrack <- GenomeAxisTrack()
  
  #SequenceTrack: Just fasta file for NC_003028. Modify it so first line at > reads ">NC_003028"
  fcol <- c(A = "darkorange", C = "yellow", T = "darkred",G = "darkgreen")
  sTrack<-SequenceTrack(genomeFile,name="Sequence",fontcolor=fcol)
  
  #AnnotationTrack1: from NC_003028.gbk. Got genes and coordinates from genbank file using script getCoordsGBK.py
  
  # is.numeric(genes[,2])
  # Make sure coordinate columns 2 and 3 are numeric. If not then reassign as numerically converted column
  #geneFile="tigr4_new.txt"
  geneDF1<-as.data.frame(read.table(file = geneFile,header=TRUE))
  anTrack<-AnnotationTrack(start=geneDF1$start,end=geneDF1$end,chromosome="genome",ignore.strand=FALSE,strand=geneDF1$strand,id=geneDF1$id,showFeatureId=TRUE, name="Genes",fill="gray",fontcolor.item="black")
  
  #DataTrack: single aggregate fitness per insertion site
  #singleFitFile="singleFit.txt"
  insertFit<-as.data.frame(read.table(file=singleFitFile))
  insertFit$seqnames="genome"
  insertFit<-insertFit[c("seqnames","V1","V2","V3")]
  colnames(insertFit)<-c("seqnames","start","end","fit")
  insertFit$fit=insertFit$fit;
  insertFit$seqnames="genome"
  insertData<-makeGRangesFromDataFrame(insertFit,keep.extra.columns=TRUE,seqinfo=NULL,seqnames.field="seqnames",start.field="start",starts.in.df.are.0based=FALSE)
  aggTrack<-DataTrack(insertData,chromosome="genome",ylim=c(0,2),type="histogram",
                      col.baseline="black",baseline=1,fill="darkred",
                      name="Insertions")
  #background.panel="#FFFFE0",
  
  #AnnotationTrack2: Fitness values that corresponding to DataTrack1
  output<-Wavg(insertFit,fromCoord,toCoord)
  avg<-round(as.numeric(output[1]),digits=4)
  count<-output[2]
  avgText<-paste0(avg,"  | Actual Aggregate Fitness= ",round(as.numeric(avg),digits=4))
  avgTrack<-AnnotationTrack(start=fromCoord,end=toCoord,chromosome="genome",id=avgText,showFeatureId=TRUE,fill="darkred", name="Window average")
  
  #Track for individual fitness values when the window is small enough
  fitTrack<-AnnotationTrack(insertFit,chromosome="genome",
                            showFeatureId=TRUE,showId=FALSE,id=insertFit$fit,fill="transparent",
                            fontcolor.item="black",col.frame="transparent",
                            background.title="transparent",col="transparent",name="Fitness Value")
  
  
  #######################################################################################################################################
  #Homologous gene in other strain
  
  matchFile ="tigr4-19F-conv.txt"
  matchDF<-as.data.frame(read.table(matchFile))
  match2<-matchDF[grep(gene_id,matchDF$V1),]
  
  gene2<-as.character(match2$V2)
  neg<-"No homolog exists in conversion file!"
  if(is.na(gene2))return(neg)
  else print(gene2)
  genomeFile2="NC_012469.fa"
  geneFile2="19F_012469_genes.txt"
  singleFitFile2="singleFit_19F_SDDM.txt"
  
  #gene id was inputted so need to get start and end coordinates for that gene
  geneDF<-as.data.frame(read.table(file = geneFile2,header=TRUE))
  geneDF$seqnames="genome2"
  #geneDF$width=abs(geneDF$width)
  #geneDF$id<-as.character(geneDF$id)
  search<-geneDF[grep(gene2,geneDF$id),]
  #print (search$start)
  #print (search$end)
  from2<-as.numeric(search$start)
  to2<-as.numeric(as.character(search$end))
  print (from2)
  print (to2)
  
  
  #so Gviz doesn't look for NC_003028 in the UCSC site
  options(ucscChromosomeNames = FALSE)
  
  #GenomeAxisTrack: Number coordinates for genome (no file input)
  gTrack2 <- GenomeAxisTrack()
  
  #SequenceTrack: Just fasta file for NC_003028. Modify it so first line at > reads ">NC_003028"
  fcol <- c(A = "darkorange", C = "yellow", T = "darkred",G = "darkgreen")
  sTrack2<-SequenceTrack(genomeFile2,name="Sequence",fontcolor=fcol)
  
  #AnnotationTrack1: from NC_003028.gbk. Got genes and coordinates from genbank file using script getCoordsGBK.py
  
  # is.numeric(genes[,2])
  # Make sure coordinate columns 2 and 3 are numeric. If not then reassign as numerically converted column
  #geneFile="tigr4_new.txt"
  geneDF<-as.data.frame(read.table(file = geneFile2,header=TRUE))
  geneDF$seqnames="genome2"
  #geneDF$end<-as.numeric(geneDF$end)
  anTrack2<-AnnotationTrack(start=geneDF$start,width=geneDF$width,chromosome="genome2",strand=geneDF$strand,
                            id=geneDF$id,showFeatureId=TRUE, name="Genes",fill="gray",
                            fontcolor.item="black")
  
  #DataTrack: single aggregate fitness per insertion site
  #singleFitFile="singleFit.txt"
  insertFit2<-as.data.frame(read.table(file=singleFitFile2))
  insertFit2$seqnames="genome2"
  insertFit2<-insertFit2[c("seqnames","V1","V2","V3")]
  colnames(insertFit2)<-c("seqnames","start","end","fit")
  insertFit2$fit=insertFit2$fit-1;
  insertFit2$seqnames="genome2"
  insertData2<-makeGRangesFromDataFrame(insertFit2,keep.extra.columns=TRUE,seqinfo=NULL,strand.field="width",seqnames.field="seqnames",start.field="start",starts.in.df.are.0based=FALSE)
  aggTrack2<-DataTrack(insertData2,chromosome="genome2",ylim=c(-.5,.5),type="histogram",
                       col.baseline="black",baseline=0,fill="darkred",
                       name="Insertions")
  #background.panel="#FFFFE0",
  
  #AnnotationTrack2: Fitness values that corresponding to DataTrack1
  output<-Wavg(insertFit2,from2,to2)
  avg<-round(as.numeric(output[1]),digits=4)
  count<-output[2]
  avgText2<-paste0(avg,"  | Actual Aggregate Fitness= ", round(as.numeric(avg+1),digits=4))
  avgTrack2<-AnnotationTrack(start=from2,end=to2,chromosome="genome2",id=avgText2,showFeatureId=TRUE,fill="darkblue", name="Window average")
  
  #Track for individual fitness values when the window is small enough
  fitTrack2<-AnnotationTrack(insertFit,chromosome="genome2",
                             showFeatureId=TRUE,showId=FALSE,id=insertFit$fit,fill="transparent",
                             fontcolor.item="black",col.frame="transparent",
                             background.title="transparent",col="transparent",name="Fitness Value")
  
  
  
  
  
  #######################################################################################################################################
  
  #Keep track pollution low:
  print(paste0("# of insertions in the window: ", count))
  if (count>15 || (toCoord-fromCoord>10000)){
    ht1 <- HighlightTrack(trackList = list(aggTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "genome",fill="#FFFACC",col="transparent")
    ht2 <-HighlightTrack(trackList = list(aggTrack2,avgTrack2), start =from2,end= to2,chromosome = "genome2",fill="#FFFACC",col="transparent")
  }
  else{
    ht1 <- HighlightTrack(trackList = list(aggTrack,avgTrack,fitTrack), start =fromCoord,end= toCoord,chromosome = "genome",fill="#FFFACC",col="transparent")
    ht2 <-HighlightTrack(trackList = list(aggTrack2,avgTrack2), start =from2,end= to2,chromosome = "genome2",fill="#FFFACC",col="transparent")
  }
  #plotTracks(list(gTrack,sTrack,anTrack,aggTrack,fitTrack),
  #background.title="darkred",fontsize=17,chromosome="genome",
  #from=fromCoord,to=toCoord)
  
  #TO CHANGE SURROUNDING WINDOWS ALTER:
  if (missing(x))x=50
  if (missing(y))y=50
  
  gTrack
  sTrack
  anTrack
  insertData
  aggTrack
  
  #plotTracks(list(gTrack,sTrack,anTrack,ht1),chromosome="genome",background.title="darkred",fontsize=17,from=fromCoord-x,to=toCoord+y,main="19F from Galaxy")
  par(mfrow=c(2,1))
  plotTracks(list(gTrack,sTrack,anTrack,ht1),chromosome="genome",background.title="darkred",fontsize=14,from=fromCoord-x,to=toCoord+y, main="TIGR4 | SDDM experiment | 15.11.04 | MLA")
  plotTracks(list(gTrack2,sTrack2,anTrack2,ht2),chromosome="genome2",background.title="darkblue",fontsize=14,from=from2-x,to=to2+y, main="Taiwan 19F | SDDM experiment | 15.10.31 | MLA")
  
  confirm=paste0("Tracks plotted for genomic coordinates ",fromCoord," to ",toCoord," with a window (",fromCoord-x,",",toCoord+y,")")
  return(confirm)
  
  
}

Wavg<-function(fit,fromCoord,toCoord){
  sum=0; count=0;i=1;
  while (fit[i,2]<=fromCoord) i<-i+1;
  
  while (fit[i,2]<=toCoord){
    sum<-sum+fit[i,4]; count<-count+1;
    i<-i+1;
  }
  avg<-(sum/count);
  output<-list(avg,count)
  return(output);
}


getTracks2<-function(fromCoord,toCoord,x,y){
  #Load libraries
  library(seqinr)
  library(Biostrings)
  library(Gviz)
  
  genomeFile="NC_003028.fa"
  geneFile="tigr4_genes.txt"
  singleFitFile="singleFit_tigr4.txt"
  
  #genomeFile="NC_012469.fa"
  #geneFile="19F_genes.txt"
  #singleFitFile="singleFit_19F_SDDM.txt"
  
  
  #so Gviz doesn't look for NC_003028 in the UCSC site
  options(ucscChromosomeNames = FALSE)
  
  #GenomeAxisTrack: Number coordinates for genome (no file input)
  gTrack <- GenomeAxisTrack()
  
  #SequenceTrack: Just fasta file for NC_003028. Modify it so first line at > reads ">NC_003028"
  fcol <- c(A = "darkorange", C = "yellow", T = "darkred",G = "darkgreen")
  sTrack<-SequenceTrack(genomeFile,name="Sequence",fontcolor=fcol)
  
  #AnnotationTrack1: from NC_003028.gbk. Got genes and coordinates from genbank file using script getCoordsGBK.py
  
  # is.numeric(genes[,2])
  # Make sure coordinate columns 2 and 3 are numeric. If not then reassign as numerically converted column
  #geneFile="tigr4_new.txt"
  geneDF<-as.data.frame(read.table(file = geneFile,header=TRUE))
  anTrack<-AnnotationTrack(start=geneDF$start,end=geneDF$end,chromosome="genome",strand=geneDF$strand,id=geneDF$id,showFeatureId=TRUE, name="Genes",fill="gray",fontcolor.item="black")
  
  #DataTrack: single aggregate fitness per insertion site
  #singleFitFile="singleFit.txt"
  insertFit<-as.data.frame(read.table(file=singleFitFile))
  insertFit$seqnames="genome"
  insertFit<-insertFit[c("seqnames","V1","V2","V3")]
  colnames(insertFit)<-c("seqnames","start","end","fit")
  insertFit$fit=insertFit$fit-1;
  insertFit$seqnames="genome"
  insertData<-makeGRangesFromDataFrame(insertFit,keep.extra.columns=TRUE,seqinfo=NULL,seqnames.field="seqnames",start.field="start",starts.in.df.are.0based=FALSE)
  aggTrack<-DataTrack(insertData,chromosome="genome",ylim=c(-1,1),type="histogram",
                      col.baseline="black",baseline=0,fill="darkred",
                      name="Insertions")
  #background.panel="#FFFFE0",
  
  #AnnotationTrack2: Fitness values that corresponding to DataTrack1
  output<-Wavg(insertFit,fromCoord,toCoord)
  avg<-round(as.numeric(output[1]),digits=4)
  count<-output[2]
  avgText<-paste0("Window Average= ",avg,"  |   Actual Aggregate Fitness= ", round(as.numeric(avg+1),digits=4))
  avgTrack<-AnnotationTrack(start=fromCoord,end=toCoord,chromosome="genome",id=avgText,showFeatureId=TRUE,fill="darkred", name="Window average")
  
  #Track for individual fitness values when the window is small enough
  fitTrack<-AnnotationTrack(insertFit,chromosome="genome",
                            showFeatureId=TRUE,showId=FALSE,id=insertFit$fit,fill="transparent",
                            fontcolor.item="black",col.frame="transparent",
                            background.title="transparent",col="transparent",name="Fitness Value")
  
  
  
  
  #Keep track pollution low:
  print(paste0("# of insertions in the window: ", count))
  if (count>15 || (toCoord-fromCoord>10000)){
    ht1 <- HighlightTrack(trackList = list(aggTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "genome")
  }
  else{
    ht1 <- HighlightTrack(trackList = list(aggTrack,fitTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "genome")
    
  }
  #plotTracks(list(gTrack,sTrack,anTrack,aggTrack,fitTrack),
  #background.title="darkred",fontsize=17,chromosome="genome",
  #from=fromCoord,to=toCoord)
  
  #TO CHANGE SURROUNDING WINDOWS ALTER:
  if (missing(x))x=50
  if (missing(y))y=50
  
  plotTracks(list(gTrack,sTrack,anTrack,ht1),chromosome="genome",background.title="darkred",fontsize=17,from=fromCoord-x,to=toCoord+y,main="TIGR4 from Penn expt")
  
  confirm=paste0("Tracks plotted for genomic coordinates ",fromCoord," to ",toCoord," with a window (",fromCoord-x,",",toCoord+y,")")
  return(confirm)
  
  
}